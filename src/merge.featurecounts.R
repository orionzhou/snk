#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Merge featurecounts output')
parser$add_argument("fi", nargs='+', help = "featurecount result file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="featurecounts.tsv",
                    help = "output file [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo

require(dplyr)
require(readr)

nfile = length(fis)
tn = tibble()
for (i in 1:nfile) {
    fi = fis[i]
    ti = read_tsv(fi, skip = 1)
    if(ncol(ti) != 7)
        stop(sprintf("not 7 cols: %s", fi))
    cnames = colnames(ti)[-c(1:6)]
    res = strsplit(cnames, split = ":")
    sids = sapply(res, "[", 2)
    tc = ti[,-c(2:6)]
    colnames(tc) = c('gid', sids)

    if(i == 1) {
        tn = tc
    } else {
        if(! identical(tn$gid, tc$gid))
           stop('not identical gids')
        tn = cbind(tn, tc[,-1])
    }
}

write_tsv(tn, fo)
