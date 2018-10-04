#{{{ load required libraries, define common variables
require(grid)
require(tidyverse)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(cluster)
require(Hmisc)
require(ggsignif)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(scales)
require(pheatmap)
require(yaml)
options(stringsAsFactors = FALSE)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))
#
dird = '~/projects/maize.reseq/data'
dirp = dird
dirc = '/scratch.global/zhoux379/maize.reseq'
f_cfg = file.path(dird, '01.cfg.tsv')
t_cfg = read_tsv(f_cfg)
f_yml = file.path(dird, '11.cfg.yaml')
Sys.setenv("R_CONFIG_FILE" = f_yml)
#}}}




