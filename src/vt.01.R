source("functions.R")
dird = "~/projects/reseq/data"
diro = file.path(dird, '31_vnt_list')
fi = file.path(dird, '../data/11_geno_list/j07.tsv')
ti = read_tsv(fi)

#{{{ vt01
yid = 'vt01'
gts_x = c('B73','Mo17','W22','Ph207','PHB47','A682','B84','Teosinte',
          'W64A','H84','H99','Oh43','B37')
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
to %>% filter(yid %in% c("dn14a",'dn15a','dn19a')) %>% print(n=50)

fo = sprintf("%s/%s.txt", diro, yid)
write(to$sid, fo)
#}}}



