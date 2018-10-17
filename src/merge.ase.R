#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Merge ASE output')
parser$add_argument("fi", nargs='+', help = "ASE tabular file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="ase.tsv",
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
    ti = read_tsv(fi)
    if(!identical(colnames(ti), c('gid','n0','n1','ncft')))
        stop(sprintf("not 4 cols: %s", fi))
    fname = basename(fi)
    sid = str_replace(fname, '[\\.]tsv$', '')
    ti = ti %>% mutate(SampleID = sid) %>% 
        select(SampleID, everything())
    tn = tn %>% bind_rows(ti)
}

to = tn
write_tsv(to, fo)
