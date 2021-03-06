require(devtools)
load_all("~/git/rmaize")
dird = '~/projects/reseq/data'
dirp = dird
dirc = '/scratch.global/zhoux379/reseq'
t_cfg = read_gspread(lib='dnaseq')
f_biomap = '~/projects/barn/data/15_read_list/dn18a.tsv'
gts_biomap = read_tsv(f_biomap) %>% distinct(Genotype) %>% pull(Genotype)
gts_nam25 = c('B97','CML52','CML69','CML103','CML228','CML247','CML277',
    'CML322','CML333','HP301','Il14H','Ki3','Ki11','Ky21','M37W','M162W',
    'Mo18W','Ms71','NC350','NC358','Oh7B','Oh43','P39','Tx303','Tzi8')
gts_main = c('B73','Mo17','W22','PH207','PHB47','A682','B84')



