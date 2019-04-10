source("functions.R")
t_cfg

get_read_list <- function(ti) {
#{{{
if(yid == 'hmp2') {
    #{{{ hapmap2
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = SampleName,
                  Treatment = '',
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
    #}}}
} else if(yid == 'hmp3') {
    #{{{ hapmap3
    th = ti %>%
        separate("SampleName", c('pre', 'Genotype'), sep = "_", fill = 'left') %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = Genotype,
                  Treatment = '',
                  Replicate = '',
                  paired = paired)
    #}}}
} else if(yid == 'm282') {
    #{{{ 282set
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = str_sub(SampleName,8),
                  Treatment = '',
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
    #}}}
} else if(yid == 'gem31') {
    #{{{ german31
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = str_sub(SampleName,4),
                  Treatment = '',
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
    #}}}
} else if(yid == 'cs845') {
    #{{{
    th = ti %>%
        mutate(SampleName = str_replace_all(SampleName, "^(CRT|PGT)-(\\d+)-(\\d+)_", "\\1_\\2_\\3-")) %>%
        mutate(SampleName = str_replace_all(SampleName, "^(CRT|PGT)-(\\d+)_", "\\1_\\2-")) %>%
        separate(SampleName, c('idx','SampleName'), sep='-', extra='merge') %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = idx,
                  Treatment = SampleName,
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
    #}}}
} else {
    stop(sprintf("unknown yid: %s", yid))
}
th = sra_fill_replicate(th)
th
#}}}
}

yid = 'hmp3'
yid = 'm282'
yid = 'gem31'
yid = 'cs845'
fi = sprintf("%s/03_sra_list/%s.csv", dird, yid)
fi2 = sprintf("%s/03_sra_list/%s_exp.csv", dird, yid)
ti = read_sra_run(fi, fi2)

th = get_read_list(ti)
th %>% count(Tissue); th %>% count(Genotype); th %>% count(Replicate)
fo = sprintf("%s/05_read_list/%s.tsv", dird, yid)
write_tsv(th, fo)

#{{{ biomap
yid = 'biomap'
fi = file.path(dird, '03.raw.xlsx')
ti = read_xlsx(fi, sheet=yid)

# create sym-links, write read list
diro1 = sprintf("%s/%s/09_fq_interleaved", dirc, yid)
#if(!dir.exists(diro1)) system(sprintf("mkdir -p %s", diro1))
#
ndig = floor(log10(nrow(ti))) + 1
fmt_yid = sprintf("s%%0%dd", ndig)
ti = ti %>% fill(directory) %>%
    mutate(SampleID = sprintf(fmt_yid, 1:nrow(ti))) %>%
    mutate(f1 = sprintf("%s/%s.fastq.gz", directory, file)) %>%
    mutate(nf1 = sprintf("%s/%s.fq.gz", diro1, SampleID)) %>%
    mutate(f1 = normalizePath(f1), nf1 = normalizePath(nf1)) %>%
    mutate(cmd = sprintf("ln -sf %s %s", f1, nf1)) %>%
    mutate(tag = file.exists(f1))
sum(!ti$tag)

#map_int(ti$cmd, system)

th = ti %>% select(SampleID, Genotype, r0=f1) %>%
    mutate(Tissue='', Treatment='', Replicate = '', paired = T) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, r0)
th = sra_fill_replicate(th)
th %>% count(Genotype) %>% print(n=20)
fo = sprintf("%s/05_read_list/%s.tsv", dird, yid)
write_tsv(th, fo)
#}}}

#{{{ pc
yid = 'pc'
fi = sprintf("%s/05_read_list/%s.raw.tsv", dird, yid)
ti = read_tsv(fi)
diri = '/home/springer/shared/pcrisp/mutants_for_SNPs'
th1 = ti %>%
    transmute(SampleID = yid, Treatment = note,
              fi1 = file.path(diri, fq1), fi2 = file.path(diri, fq2),
              fv = file.exists(fi1) && file.exists(fi2)) %>%
    mutate(Genotype = Treatment, Tissue = '', Replicate = 1, paired = T) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, fi1, fi2, fv)
th1 %>% count(fv)
th = th1

th = th %>% mutate(r1=fi1, r2=fi2) %>% select(-fi1, -fi2, -fv)
fo = sprintf("%s/05_read_list/%s.tsv", dird, yid)
write_tsv(th, fo)
#}}}

#{{{ j01
yids = c('biomap', 'hmp2', 'hmp3')
fis = sprintf("%s/05_read_list/%s.tsv", dird, yids)
ti = tibble(yid=yids, fi=fis) %>%
    mutate(data = map(fi, read_tsv)) %>%
    select(yid, data) %>%
    unnest() %>%
    distinct(yid, Genotype) %>%
    mutate(gt=Genotype, sid=sprintf("%s_%s", yid, gt)) %>%
    select(sid,yid,gt)

ptn_gt = '(b73)|(mo17)|(w22$)|(ph207)|(b84)|(a682)'

ti1 = ti %>% filter(yid == yids[1])
ti2 = ti %>% filter(yid == yids[2], str_detect(str_to_lower(gt), ptn_gt))
ti3 = ti %>% filter(yid == yids[3], str_detect(str_to_lower(gt), ptn_gt) | gt %in% ti1$gt)

to = rbind(ti1,ti2,ti3) %>%
    select(SampleID = sid, yid, Genotype = gt)
yid = 'j01'
fo = sprintf("%s/05_read_list/%s.tsv", dird, yid)
write_tsv(to, fo)
#}}}

