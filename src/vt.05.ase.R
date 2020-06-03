source("functions.R")
dird = "~/projects/reseq/data"
f_cfg = '~/projects/master.xlsx'

#{{{
yid = 'vt01'
diro = file.path(dird, '31_vnt_list')
fi = sprintf("%s/%s.txt", diro, yid)
ti = read_tsv(fi, col_names=c('sid')) %>%
    separate(sid, c('yid','Genotype'), sep="#", remove=F) %>%
    mutate(Genotype=str_to_upper(Genotype)) %>% select(-yid)

yids = c('ca20a3','rn14f','rn17b','rn17c','rn18b','rn18g','rn20a','rn20b','rn20d3','rn20e')
diri = '~/projects/barn/data/15_read_list'
tg = tibble(yid=yids) %>% mutate(fi = sprintf("%s/%s.tsv", diri, yid)) %>%
    mutate(data=map(fi, read_tsv)) %>%
    select(yid, data) %>% unnest(data) %>%
    distinct(yid, Genotype) %>%
    add_row(yid = 'rn18c',Genotype = 'W22xTeosinte')

to = tg %>% mutate(inbred = !str_detect(Genotype, 'x'))
to1 = to %>% filter(inbred) %>% mutate(Genotype=str_to_upper(Genotype)) %>%
    left_join(ti, by='Genotype')
to2 = to %>% filter(!inbred) %>%
    separate(Genotype, c('pa1','pa2'), sep='x', remove=F) %>%
    mutate(pa1 = str_to_upper(pa1), pa2 = str_to_upper(pa2)) %>%
    mutate(Genotype = str_c(pa1, pa2, sep='x')) %>%
    left_join(ti, by=c('pa1'='Genotype')) %>% rename(sid1=sid) %>%
    left_join(ti, by=c('pa2'='Genotype')) %>% rename(sid2=sid)
to1 %>% filter(is.na(sid)) %>% print(n=40)
to2 %>% filter(is.na(sid1) | is.na(sid2)) %>% print(n=40)
tp1 = to1 %>% filter(!is.na(sid)) %>% distinct(Genotype, inbred, sid)
tp1 = ti %>% mutate(Genotype=str_to_upper(Genotype), inbred=T) %>% select(Genotype,inbred, sid)
tp2 = to2 %>% filter(!is.na(sid1) & !is.na(sid2)) %>% distinct(Genotype, inbred, sid1, sid2)
tp = tp1 %>% bind_rows(tp2)

fo = sprintf("%s/35_vnt_ase/%s.tsv", dird, yid)
write_tsv(tp, fo, na='')
#}}}

read.bcftools.stats <- function(filename, sections = c("ID", "SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", 'PSC', 'PSI', "label")) {
    #{{{
    con <- file(filename, open = "r")
    lines <- readLines(con)
    obj <- list()
    # Loop the labels and convert to data frames
    for (x in sections) {
        if (x == "label") next
        header <- gsub("\\[[0-9]+\\]", "", unlist(strsplit(gsub("# ", "", lines[min(grep(paste(x, "\\t", sep=""), lines))]), "\t")))
        obj[[x]] <- read.delim(tc <- textConnection(lines[grepl(paste(x, "\\t", sep = ""), lines)]),
                                                              header = FALSE, skip = 1, col.names = header) %>% as_tibble()
        close(tc)
        obj[[x]][x] <- NULL
    }
    close(con)
    obj
    #}}}
}

merge.bcftools.stats <- function(ti, sections = c("ID", "SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", 'PSC', 'PSI', "label")) {
    #{{{
    obj <- list()
    for (x in names(ti$data[[1]])) {
        obj[[x]] = ti %>% mutate(data2 = map(data, x)) %>%
            select(rid, data2) %>% unnest()
    }
    obj
    #}}}
}

#{{{ vcf stats
yid = 'vt01'
fi = sprintf("%s/35_vnt_ase/%s.tsv", dird, yid)
ti = read_tsv(fi)
ti2 = ti %>%
    mutate(fs = sprintf("~/projects/reseq/data/ase_vcf/%s.txt", Genotype)) %>%
    filter(file.exists(fs)) %>%
    mutate(res = map(fs, read.bcftools.stats))

keys = c('records', 'no_ALTs', 'SNPs', 'indels', 'MNPs', 'others', 'multiallelic', 'multiallelic_SNP')
to = ti2 %>%
    mutate(x = map(res, 'SN')) %>% select(-fs, -res) %>% unnest() %>%
    select(-id) %>%
    mutate(key = str_replace(key, "^number of ", '')) %>%
    mutate(key = str_replace(key, ": *", '')) %>%
    mutate(key = str_replace(key, " sites$", '')) %>%
    mutate(key = str_replace(key, "[ -]", '_')) %>%
    filter(! key %in% c("samples")) %>%
    mutate(key = factor(key, levels = keys)) %>%
    spread(key, value) %>%
    arrange(-inbred, Genotype)

fo = sprintf("%s/35_vnt_ase/%s.stats.tsv", dird, yid)
write_tsv(to, fo, na='')
#}}}

