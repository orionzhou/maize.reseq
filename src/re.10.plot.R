source("functions.R")
require(treeio)
require(ggtree)
require(tidytree)
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
