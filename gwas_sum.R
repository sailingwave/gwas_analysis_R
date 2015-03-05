#!/usr/bin/env Rscript

# Nan Wang
# University of Southern California
# Email: sailingwave@gmail.com


#=== for devlopment

test = T

if(test){
    setwd('/Users/Nan/Desktop/')    
    results <- 'di_all.txt'
    output  <- 'di_all.tiff'
}else{
    args <- commandArgs(trailingOnly = TRUE)
    
    results <- args[1]
    output  <- args[2]
}

#==================#

library(ggplot2)
library(dplyr)
library(data.table)

association_data<-fread(results, header=TRUE)
setnames(association_data,names(association_data),tolower(names(association_data)))
setkey(association_data,chr)

# relabel chr X and XY = 23 if necessary
if ("X" %in% association_data$chr) {
    association_data["X",chr:= "23"]
}
if ("XY" %in% association_data$chr) {
    association_data["XY",chr:= "23"]
}

tbl.assoc = tbl_dt(association_data)

if(test){
    bak <- copy(tbl.assoc)
    set.seed(123)
    tbl.assoc = sample_frac(tbl.assoc,0.0001, replace = TRUE)
}

max_pos <- tbl.assoc %>%
    group_by(chr) %>%
    summarise(max(pos)) %>%
    arrange(as.numeric(chr))
setnames(max_pos,names(max_pos),c('chr','max_pos'))
max_pos$max_pos = as.numeric(max_pos$max_pos)
x_lab = paste0('chr',max_pos$chr)

max_pos <- max_pos %>%
    mutate(cs = cumsum(max_pos),lcs = lag(cs,1,default=0),chr_mk = lcs+max_pos/2)
#cumulative sum, lag cs, tick mark positions for each chr

tbl.assoc <- tbl.assoc %>%
    left_join(max_pos) %>%
    mutate(x_pos = pos+lcs)

# Manhattan plot
n_chr = length(max_pos$chr)
pallete = rep(c("blue4", "orange3"),times = ceiling(n_chr/2))[1:n_chr]

man_plot <- ggplot(tbl.assoc, aes(x=x_pos, y=-log10(p), color=chr)) +
    geom_vline(xintercept=max_pos$lcs[-1],linetype="dotted", size=.8,color='white') +
    geom_hline(aes(yintercept=-log10(5e-08)),color="red") +
    geom_point(shape=20, size = 3)+ xlab("") +
    scale_x_continuous(breaks=max_pos$chr_mk,labels=x_lab) +
    scale_color_manual(values=pallete,guide=FALSE) +
    theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank()) +
    theme(axis.title.y=element_text(size=16)) +
    theme(axis.text.y=element_text(size=12))


tiff(filename=output, compression="lzw", width = 1200, height = 600, units = "px")
    print(man_plot)
dev.off()



# Calculate GC lambda

lambda = median(qnorm(data$p/2)^2,na.rm=T)/qchisq(0.5,1) 
print(paste("File: ", results, ": Lamda = ", lambda, sep=""))


rm(list=ls())

