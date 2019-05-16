source("functions.R")
t_cfg

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

fo = file.path(dird, '10.vcfstats.tsv')
write_tsv(tp, fo)
