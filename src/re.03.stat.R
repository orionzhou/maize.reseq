source("functions.R")
t_cfg %>% filter(libtype=='dnaseq')

#{{{ collect all vcf stats
yids = c('dn12a','dn14a','dn15a','dn17a','dn17b','dn18a','dn19a')
th = t_cfg %>% filter(yid %in% yids) %>%
    select(yid,author,study,n)

vcfstat_cols = c('fpath','id','sample','nRefHom','nNonRefHom',
      'nHets','nTransitions','nTransversions','nIndels',
      'avgDepth','nSingletons','nHapRef','nHapAlt','nMissing')

tp = th %>% mutate(fp = sprintf("%s/raw/%s/vcfstats.tsv", dird, yid)) %>%
    mutate(stat = map(fp, read_tsv, col_names = vcfstat_cols)) %>%
    unnest() %>%
    select(yid,author,study,n,sid=sample,nRefHom,nNonRefHom,nHets,avgDepth,nMissing)
tp %>% count(yid,author,study,n)
tp %>% group_by(yid,study,n) %>%
    summarise(d5 = sum(avgDepth>=5)) %>%
    ungroup()

fo = file.path(dird, '08.vcfstats.tsv')
write_tsv(tp, fo)
#}}}

fi = file.path(dird, '08.vcfstats.tsv')
ti = read_tsv(fi)
ti %>% count(yid,author,study,n)
ti %>% group_by(yid,study,n) %>%
    summarise(d5 = sum(avgDepth>=5)) %>%
    ungroup()

#{{{ j01: 5 studies
yid = 'j01'
yids = c('dn12a','dn14a','dn17a','dn17b','dn18a')
fo = sprintf("%s/11_geno_list/%s.tsv", dird, yid)
to = ti %>% filter(yid %in% yids, avgDepth>=5) %>%
    separate(sid, c("yid",'gt'), sep='#', remove=F) %>%
    select(sid, yid, gt)
write_tsv(to, fo)
#}}}

#{{{ j07: 7 studies
yid = 'j07'
yids = c('dn12a','dn14a','dn15a','dn17a','dn17b','dn18a','dn19a')
fo = sprintf("%s/11_geno_list/%s.tsv", dird, yid)
to = ti %>% filter(yid %in% yids, avgDepth>=5) %>%
    separate(sid, c("yid",'gt'), sep='#', remove=F) %>%
    select(sid, yid, gt)
to %>% count(yid)
write_tsv(to, fo)
#}}}


