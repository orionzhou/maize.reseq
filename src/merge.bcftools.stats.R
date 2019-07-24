#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'merge bcftools stats files')
parser$add_argument("fi", nargs='+', help = "stats file(s)")
parser$add_argument("fo", nargs=1, help = "output file")
parser$add_argument("-s", "--simple", action="store_true",
                    help = "simple output [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo

require(tidyverse)

read.bcftools.stats <- function(filename, sections = c("ID", "SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", 'PSC', 'PSI', "label")) {
    #{{{
    con <- file(filename, open = "r")
    lines <- readLines(con)
    obj <- list()
    # Loop the labels and convert to data frames
    for (x in sections) {
        if (x == "label") next
        header <- gsub("\\[[0-9]+\\]", "", unlist(strsplit(gsub("# ", "", lines[min(grep(paste(x, "\\t", sep=""), lines))]), "\t")))
        obj[[x]] <- read.delim(tc <- textConnection(lines[grepl(paste(x, "\\t", sep = ""), lines)]),
                                                              header = FALSE, skip = 1, col.names = header) %>% as_tibble()
        close(tc)
        obj[[x]][x] <- NULL
    }
    close(con)
    obj
    #}}}
}

merge.bcftools.stats <- function(ti, sections = c("ID", "SN", "TSTV", "SiS", "AF", "QUAL", "IDD", "ST", "DP", 'PSC', 'PSI', "label")) {
    #{{{
    obj <- list()
    for (x in names(ti$data[[1]])) {
        obj[[x]] = ti %>% mutate(data2 = map(data, x)) %>%
            select(rid, data2) %>% unnest()
    }
    obj
    #}}}
}

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(rid = str_replace(fname, '[\\.]\\w+$', ''))

to = ti %>% mutate(data = map(fi, read.bcftools.stats)) %>% select(rid, data)

x = merge.bcftools.stats(to)

if(args$simple) {
    to = x$SN %>% select(-id) %>%
        group_by(key) %>%
        summarise(value = sum(value)) %>% ungroup() %>%
        mutate(key = str_replace(as.character(key), ':$', ''))
    write_tsv(to, fo)
} else {
    saveRDS(x, file=fo)
}

