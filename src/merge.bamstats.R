#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Merge bam_stat output')
parser$add_argument("fi", nargs='+', help = "bam_stat result file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="bamstats.tsv",
                    help = "output file [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo

require(dplyr)
require(readr)
require(stringr)
require(tidyr)

nfile = length(fis)
tn = tibble()
for (i in 1:nfile) {
    fi = fis[i]
    ti = read_tsv(fi, col_names = c("type","cnt"), col_types='cd')
    if(ncol(ti) != 2)
        stop(sprintf("not 2 cols: %s", fi))
    fname = basename(fi)
    sid = str_replace(fname, '[\\.]tsv$', '')
    colnames(ti)[2] = sid

    if(i == 1) {
        tn = ti
    } else {
        if(! identical(tn$type, ti$type))
           stop('not identical 1st column')
        tn = cbind(tn, ti[,-1])
    }
}

to = tn %>% 
    gather(sid, cnt, -type) %>%
    spread(type, cnt)
write_tsv(to, fo)
