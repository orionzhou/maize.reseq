source("functions.R")
dirw = file.path(dird, '21_qc')
t_cfg

yid = 'j07'
fi = file.path(dird, 'raw', yid, 'bcftools_stats.rds')
ti = readRDS(fi)

tp = ti$PSC %>% mutate(nSites=nRefHom+nNonRefHom+nHets+nMissing) %>%
    group_by(sample) %>%
    summarise(nRefHom=sum(nRefHom), nNonRefHom=sum(nNonRefHom),
        nHets=sum(nHets), nMissing=sum(nMissing),
        avgDepth=sum(nSites*average.depth)/sum(nSites),
        nSites=sum(nSites),
        nSingletons=sum(nSingletons)) %>%
    ungroup() %>%
    separate(sample, c("yid", "genotype"), sep="#") %>%
    select(yid,genotype,nRefHom,nNonRefHom,nHets,nMissing,nSingletons,nSites,avgDepth) %>%
    inner_join(t_cfg[,c('yid','alias','author')], by='yid') %>% select(-yid) %>%
    select(study=alias,author,everything())

fo = file.path(dirw, yid, '01.bcfstats.tsv')
write_tsv(tp, fo)
