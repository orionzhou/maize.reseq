source("functions.R")
require(treeio)
require(ggtree)
require(tidytree)
dird = "~/projects/reseq/phylo"

#{{{ hmp3 tree
dirw = file.path(dird, '20_hmp')
fi = file.path(dirw, '25.treefile')
tree = read.newick(fi)

gts = c('all','b73','mo17','w22','ph207','x')
gt_col = c('black',pal_locuszoom()(length(gts)-1))
names(gt_col) = gts
tp = as_tibble(tree) %>%
    filter(!is.na(label)) %>%
    mutate(gt = 'all') %>%
    mutate(gt = ifelse(str_detect(str_to_lower(label), 'b73'), 'b73', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(label), 'mo17'), 'mo17', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(label), 'w22$'), 'w22', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(label), 'ph207'), 'ph207', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(label), 'b84'), 'x', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(label), 'a682'), 'x', gt)) %>%
    mutate(col = gt_col[gt])

p1 = ggtree(tree) +
    geom_tiplab(node=tp$node, size=2, color=tp$col) +
    scale_y_continuous(expand=c(0,2))
fo = file.path(dirw, '25.pdf')
ggsave(p1, file=fo, width=15, height=70, limitsize=F)
#}}}

#{{{
yid = 'ph01'
dirw = file.path(dird, '11_qc', yid)
diri = file.path(dird, 'raw', yid)
fi = file.path(dirw, '35.treefile')
tree = read.newick(fi)

th = t_cfg %>% select(study=alias, yid) %>% filter(!is.na(study))
ft = '~/projects/reseq/data/21_qc/j01/01.bcfstats.tsv'
tt = read_tsv(ft) %>% inner_join(th, by='study') %>%
    mutate(taxa = str_c(yid, genotype, sep='_'))

gts = c('all','h')
gt_col = c('black',pal_nejm()(1))
names(gt_col) = gts
tp = tibble(taxa = tree$tip.label) %>%
    left_join(tt, by='taxa') %>%
    mutate(lab = sprintf("%s  %4.01fx", genotype, avgDepth)) %>%
    mutate(gt = 'all') %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'b73'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'mo17'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'w22$'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'ph207'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'oh43'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'phb47'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'a682'), 'h', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'b84'), 'h', gt)) %>%
    mutate(col = gt_col[gt])
#
studies = unique(tp$study)
st_col = pal_aaas()(5)
names(st_col) = studies
p1 = ggtree(tree) %<+%
    tp +
    geom_tiplab(aes(label=study,col=study), size=2.5, hjust=0, align=T, linesize=.5) +
    geom_tiplab(aes(label=lab,col=gt), size=2.5, hjust=0, align=T, linesize=.5, offset=.012, linetype='blank') +
    #geom_text(aes(label=avgDepth), size=2.5, hjust=0, align=T, linesize=.5) +
    scale_color_manual(values=c(gt_col, st_col)) +
    scale_x_continuous(expand=expand_scale(mult=c(0,.3))) +
    scale_y_continuous(expand=c(0,1))
fo = file.path(dirw, '05.pdf')
ggsave(p1, file=fo, width=8, height=35, limitsize = F)
#}}}


