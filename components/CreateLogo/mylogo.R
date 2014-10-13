# read file
args <- commandArgs(trailingOnly = TRUE)
wm_file <- args[1]
format = args[2] # svg or ps
outfile = args[3]
library(seqLogo)
text <- readLines(wm_file)
wm <- matrix(nrow=0, ncol=4)
started <- FALSE
for (line in text) {
    if (grepl("^[[:digit:]]+", line)) {
       vals <- as.numeric(unlist(strsplit(line, "\t"))[2:5])
       vals <- vals / sum(vals)
       wm <- rbind(wm, vals)
    }
    else if (grepl("^NA", line)) {
         title <- paste(unlist(strsplit(line, " "))[-1], sep=" ")
    }
    else if (grepl("^P0|^PO", line)) {
         letters <- unlist(strsplit(line, "[[:space:]]+"))[2:5]
         if (all(letters != c("A", "C", "G", "T"))) {
            stop("Only A C G T order is accepted!")
         }
    }
    else if (grepl("^//", line)) {
         if (started) {
            started <- FALSE
            if (format == "svg") {
                svg(paste(outfile, "svg", sep="."))
            }
            if (format == "ps") {
                postscript(paste(outfile, "ps", sep="."))
            }
            seqLogo(t(wm), xfontsize=20, yfontsize=20)
            grid.text(title, y=.95, gp=gpar(fontsize=30))
            dev.off()
         }
         else {
            started <- TRUE
         }
    }
}
