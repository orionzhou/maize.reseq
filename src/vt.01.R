source("functions.R")
dird = "~/projects/reseq/data"
diro = file.path(dird, '31_vnt_list')

#{{{ vt01
yid = 'vt01'
fi = file.path(dird, '../data/11_geno_list/j08.tsv')
ti = read_tsv(fi)
ti %>% filter(yid %in% c("dn14a",'dn15a','dn19a')) %>% print(n=50)

fo = sprintf("%s/%s.txt", diro, yid)
write(ti$sid, fo)
#}}}



