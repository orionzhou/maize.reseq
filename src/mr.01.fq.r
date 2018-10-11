#{{{
source("mr.fun.r")
source("sra.R")
t_cfg
#}}}

get_read_list <- function(ti) {
#{{{
if(sid == 'hmp2') {
    #{{{ hapmap2
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = '',
                  Genotype = SampleName,
                  Treatment = '',
                  Replicate = '',
                  paired = paired) %>% 
        arrange(SampleID)
    #}}}
} else if(sid == 'hmp3') {
    #{{{ hapmap3
    th = ti %>% 
        separate("SampleName", c('pre', 'Genotype'), sep = "_", fill = 'left') %>%
        transmute(SampleID = Run, 
                  Tissue = '', 
                  Genotype = Genotype, 
                  Treatment = '',
                  Replicate = '',
                  paired = paired)
    #}}}
} else {
    stop(sprintf("unknown sid: %s", sid))
}
th = sra_fill_replicate(th)
th
#}}}
}

sid = 'hmp3'
fi = sprintf("%s/03_sra_list/%s.csv", dird, sid)
fi2 = sprintf("%s/03_sra_list/%s_exp.csv", dird, sid)
ti = read_sra_run(fi, fi2)

th = get_read_list(ti)
#th = th %>% filter(Genotype %in% c("Mo17", "PH207"))
th %>% count(Tissue); th %>% count(Genotype); th %>% count(Replicate)
fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)

#{{{ widiv
sid = 'widiv'
fi = sprintf("%s/05_read_list/%s.raw.txt", dird, sid)
ti = read_tsv(fi)
diri = '/home/hirschc1/shared/reads_for_SRA_submission/wi_reseq_set2'
th1 = ti %>%
    transmute(SampleID = `2=SampleId`, Genotype = `4=sampleName`, 
              fi = file.path(diri, `8=fileName`),
              fv = file.exists(fi)) %>%
    mutate(SampleID = sprintf("s%s", str_sub(as.character(SampleID),4))) %>%
    mutate(Tissue = '', Treatment = '', Replicate = 1, paired = T) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, fi, fv)
#
diri = '/home/springer/zhoux379/data/widiv_mo17'
fi = file.path(diri, 'filelist.txt')
ti = read_tsv(fi, col_names = F)
th2 = ti %>% 
    mutate(SampleID = sprintf("s%02d", 1:nrow(ti)),
           Tissue = '', Genotype = "Mo17", 
           Treatment = '', Replicate = 1:nrow(ti), paired = T, 
           fi = file.path(diri, X1), fv = file.exists(fi)) %>% select(-X1)
th = rbind(th1, th2)
table(th$fv)

diro = sprintf("%s/cache/%s/09_fq_interleaved", dird, sid)
if(!file.exists(diro)) system(sprintf("mkdir -p %s", diro))
for (i in 1:nrow(th)) {
    fl = sprintf("%s/%s.fq.gz", diro, th$SampleID[i])
    cmd = sprintf("ln -sf %s %s", th$fi[i], fl)
    system(cmd)
}

th = th %>% select(-fi, -fv)
fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
#}}}


