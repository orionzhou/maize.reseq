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
dirw = file.path(dird, 'raw', yid)
fi = file.path(dirw, '35.treefile')
tree = read.newick(fi)

fi = '~/projects/reseq/data/raw/j01/38.stat.tsv'
tt = read_tsv(fi, skip=1) %>%
    select(sid=`[3]sample`, avgDepth=`[10]average depth`)

gts = c('all','b73','mo17','w22','ph207','x')
gt_col = c('black',pal_locuszoom()(length(gts)-1))
names(gt_col) = gts
tp = tibble(taxa = tree$tip.label) %>%
    mutate(id = taxa) %>%
    left_join(tt, by=c('taxa'='sid')) %>%
    mutate(lab = sprintf("%s [%4.01fx]", taxa, avgDepth)) %>%
    mutate(gt = 'all') %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'b73'), 'b73', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'mo17'), 'mo17', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'w22$'), 'w22', gt)) %>%
    mutate(gt = ifelse(str_detect(str_to_lower(taxa), 'ph207'), 'ph207', gt)) %>%
    mutate(col = gt_col[gt])
#
p1 = ggtree(tree) %<+%
    tp +
    geom_tiplab(aes(label=lab,col=gt), size=2.5, hjust=0, align=T, linesize=.5) +
    #geom_text(aes(label=avgDepth), size=2.5, hjust=0, align=T, linesize=.5) +
    scale_color_manual(values=gt_col) +
    scale_x_continuous(expand=expand_scale(mult=c(0,.2))) +
    scale_y_continuous(expand=c(0,1))
fo = file.path(dirw, '36.pdf')
ggsave(p1, file=fo, width=8, height=10, limitsize=F)
#}}}
