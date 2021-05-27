#!/usr/bin/Rscript
library(PropCIs)

RESULTS_DIR <- "/home/fa_bula/PCAWG/results/"
freqFiles <- list.files(RESULTS_DIR, pattern = "*freq.csv$")

for (fileName in freqFiles) {
    df = read.csv(file.path(RESULTS_DIR, fileName), sep="\t", header=TRUE)
    targetsColumn <- which(colnames(df) == "targets")
    mutationsColumn <- which(colnames(df) == "mutations")
    df$percent <- 100 * df$mutations / df$targets
    ll = apply(df, 1, function(x) exactci(as.integer(x[mutationsColumn]), as.integer(x[targetsColumn]), 0.95)$conf.int[1]*100)
    df$percent_ll = ll
    ul = apply(df, 1, function(x) exactci(as.integer(x[mutationsColumn]), as.integer(x[targetsColumn]), 0.95)$conf.int[2]*100)
    df$percent_ul = ul
    nms <- colnames(df)
    df <- df[, c(nms[nms!='APOBEC_enrichment'],'APOBEC_enrichment')]
    write.table(df, file.path(RESULTS_DIR, gsub(".csv", "_ci.csv", fileName)), sep="\t", row.names=FALSE)
}
