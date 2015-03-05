#!/usr/bin/env Rscript

# Nan Wang
# University of Southern California
# ===
# Andrew R Wood
# Genetics of Complex Traits
# University of Exeter Medical School
# Email: a.r.wood@ex.ac.uk
# Date:  25/02/2015

setwd('/Users/Nan/Desktop/')

results <- 'di_all.txt'
output  <- 'di_all.tiff'
    
library(ggplot2)
library(dplyr)
library(data.table)

association_data<-fread(results, header=TRUE)
setnames(association_data,names(association_data),tolower(names(association_data)))
setkey(association_data,chr)

bak = copy(association_data)
# relabel chr X and XY = 23 if necessary
if ("X" %in% association_data$chr) {
    #system.time(association_data[,chr:= revalue(association_data$chr, c("X" = "23"))])    #revalue from plyr
    association_data["X",chr:= "23"]    #save time from 6s to 0.003s
}
if ("XY" %in% association_data$chr) {
    #association_data[,chr:= revalue(association_data$chr, c("XY" = "23"))]
    association_data["XY",chr:= "23"]
}

tbl.assoc = tbl_dt(association_data)

max_pos <- tbl.assoc %>%
    group_by(chr) %>%
    summarise(max(pos)) %>%
    arrange(as.numeric(chr))
setnames(max_pos,names(max_pos),c('chr','max_pos'))
max_pos$max_pos = as.numeric(max_pos$max_pos)

max_pos <- max_pos %>%
    mutate(cs = cumsum(max_pos)) %>%
    mutate(lcs = lag(cs,1,default=0))

tbl.assoc <- tbl.assoc %>%
    left_join(max_pos) %>%
    mutate(x_pos = pos+lcs)


# generate the plot
tiff(filename=output, compression="lzw", width = 1200, height = 600, units = "px")
print(ggplot(tbl.assoc, aes(x=x_pos, y=-log10(p), color=chr)) 
      + geom_point(shape=20, size = 1) 
      + xlab("") + theme(axis.text.x = element_blank(), legend.position="bottom") 
      + scale_color_manual(values=c("blue4", "orange3")) 
      + geom_hline(aes(yintercept=-log10(5e-08))))
dev.off()

rm(list=ls())

