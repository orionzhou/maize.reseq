source("functions.R")
dird = "~/projects/reseq/data"
f_cfg = '~/projects/master.xlsx'
t_ase = read_xlsx(f_cfg, sheet='barn', col_names=T) %>%
    replace_na(list(ase=F)) %>% filter(ase=='T')

#{{{
yid = 'vt01'
diro = file.path(dird, '31_vnt_list')
fi = sprintf("%s/%s.txt", diro, yid)
ti = read_tsv(fi, col_names=c('sid')) %>%
    separate(sid, c('yid','Genotype'), sep="#", remove=F) %>%
    mutate(Genotype=str_to_upper(Genotype)) %>% select(-yid)

diri = '~/projects/barn/data/15_read_list'
tg = t_ase %>% mutate(fi = sprintf("%s/%s.tsv", diri, yid)) %>%
    mutate(data=map(fi, read_tsv)) %>%
    select(yid, data) %>% unnest() %>%
    distinct(yid, Genotype)

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
tp2 = to2 %>% filter(!is.na(sid1) & !is.na(sid2)) %>% distinct(Genotype, inbred, sid1, sid2)
tp = tp1 %>% bind_rows(tp2)

fo = sprintf("%s/35_vnt_ase/%s.tsv", dird, yid)
write_tsv(tp, fo, na='')
#}}}


