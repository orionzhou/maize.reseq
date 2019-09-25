source("functions.R")
dird = "~/projects/reseq/phylo"
diro = file.path(dird, '05_sample_list')
fi = file.path(dird, '../data/11_geno_list/j01.tsv')
ti = read_tsv(fi)

#{{{ ph01
yid = 'ph01'
to = ti %>% select(sid)
fo = sprintf("%s/%s.txt", diro, yid)
write(to$sid, fo)
#}}}

#{{{ ph05
yid = 'ph05'
gts_x = c('B73','Mo17','W22','Ph207','PHB47','A682','B84')
ft = file.path(dird, '../data/08.vcfstats.tsv')
tt = read_tsv(ft)

gts = c(gts_nam25, gts_x, gts_biomap)
to  = ti %>%
    mutate(gt = str_to_upper(gt)) %>%
    filter(gt %in% str_to_upper(gts) | yid %in% c('dn14a','dn18a')) %>%
    inner_join(tt[,c('sid','avgDepth')], by='sid') %>%
    group_by(gt) %>%
    arrange(-avgDepth) %>%
    summarise(sid=sid[1], yid=yid[1], avgDepth=avgDepth[1])

fo = sprintf("%s/%s.txt", diro, yid)
write(to$sid, fo)
#}}}



