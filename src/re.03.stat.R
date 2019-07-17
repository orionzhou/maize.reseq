source("functions.R")
t_cfg

#{{{ collect vcf stats
yids = c('dn12a','dn14a','dn17a','dn17b','dn18a')
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


yid = 'j01'
fo = sprintf("%s/11_geno_list/%s.tsv", dird, yid)
to = ti %>% filter(avgDepth>=5) %>%
    separate(sid, c("yid",'gt'), sep='#', remove=F) %>%
    select(sid, yid, gt)
write_tsv(to, fo)

