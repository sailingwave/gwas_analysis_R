#!/usr/bin/env Rscript

# Nan Wang
# University of Southern California
# Email: sailingwave@gmail.com

# Quick summarizing large GWAS result data sets by taking advantage of data.table and dplyr
# Thanks to Andrew R Wood from University of Exeter Medical School for the original structure

# usage: Rscript gwas_sum.R arg1 arg2 arg3
# arg1: GWAS summary result file with columns chr, pos and p
# arg2 (optional): output file name of manhattan plot
# arg3 (optional): output file name of qq plot

#=== for devlopment

test = F

if(test){
    setwd('/Users/Nan/Desktop/')    
    results <- 'test.txt'
    output_mh  <- 'test.tiff'
    output_qq  <- 'test.tiff'
}else{
    args <- commandArgs(trailingOnly = TRUE)
    
    results <- args[1]
    output_mh  <- args[2]
    output_qq  <- args[3]
}


#==================#

library(ggplot2)
library(dplyr)
library(data.table)

fast_plot = T    #whether plot a small fraction of SNPs below a p value cutoff
p_cut = 1e-2
snp_frac = 0.05

#=== Data preparation

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

if(test){    #for developing
    bak <- copy(tbl.assoc)
    set.seed(123)
    tbl.assoc = sample_frac(bak,0.0001, replace = F)
}

#=== for fast plot

if(fast_plot){
    #bak <- copy(tbl.assoc)
    set.seed(888)
    sampl_snp <- bak %>%
        filter(p >= p_cut) %>%
        sample_frac(snp_frac, replace = F)
    
    tbl.assoc <- bak %>%
        filter(p < p_cut) %>%
        rbind(sampl_snp)
}

#===========#

chr_sep = 2e7    #separation among chromosomes, in base pairs; can be set to 0 if no separation is wanted

#max bp for each chr
max_pos <- tbl.assoc %>%
    group_by(chr) %>%
    summarise(max(pos)) %>%
    arrange(as.numeric(chr))
setnames(max_pos,names(max_pos),c('chr','max_pos'))
max_pos$max_pos = as.numeric(max_pos$max_pos)
x_lab = paste0('chr',max_pos$chr)    #for mh plotting

max_pos <- max_pos %>%
    mutate(max_p_sep = max_pos+chr_sep,cs = cumsum(max_p_sep),lcs = lag(cs,1,default=0),
           chr_mk = lcs+max_pos/2,chr_sep_line = lcs-chr_sep/2) %>%
    select(chr,lcs,chr_mk,chr_sep_line)
#max+chr_sep, cumulative sum, lag cs, tick mark positions for each chr, dotted-line position for separating chr 

#merge
tbl.assoc <- tbl.assoc %>%
    left_join(max_pos,by='chr') %>%
    mutate(x_pos = pos+lcs) %>%
    arrange(chr,pos)


#=== Manhattan plot

pallete <- as.numeric(unique(tbl.assoc$chr))%%2 %>% as.character()
pallete[pallete == '1'] <- "blue4"
pallete[pallete == '0'] <- "orange3"



if(!is.na(output_mh)){
    cat("Generating Manhattan plot ... ")
    
    runtime = system.time(
    {
        mh_plot <- ggplot(tbl.assoc, aes(x=x_pos, y=-log10(p), color=chr)) +
            geom_vline(xintercept=max_pos$chr_sep_line[-1],linetype="dotted", size=.8,color='white') +
            geom_hline(aes(yintercept=-log10(5e-08)),color="red") +
            geom_point(shape=20, size = 3)+ xlab("") +
            scale_x_continuous(breaks=max_pos$chr_mk,labels=x_lab) +
            scale_color_manual(values=pallete,guide=FALSE) +
            theme(panel.grid.major.x = element_blank(),panel.grid.minor = element_blank()) +
            theme(axis.title.y=element_text(size=16)) +
            theme(axis.text.y=element_text(size=12))
        
        
        tiff(filename=output_mh, compression="lzw", width = 1200, height = 600, units = "px")
            print(mh_plot)
        dev.off()
    })[3]
    
    cat(paste(" in",runtime," s\n"))
}


#=== QQ plot

if(!is.na(output_qq)){
    cat("Generating QQ plot ... ")
    
    runtime = system.time(
    {
        alpha = 0.05    #alpha level for confidence interval
        
        n_snp = nrow(tbl.assoc)
        exp = -log10(1:n_snp/n_snp)
        obs = -log10(sort(tbl.assoc$p))
        
        #expected (null) p values ~ unif(0,1), of which the order stat ~ beta(i,N-i+1)
        ci_up=-log10(qbeta(1-alpha/2, 1:n_snp, (n_snp-1:n_snp+1)))
        ci_low=-log10(qbeta(alpha/2, 1:n_snp, (n_snp-1:n_snp+1)))
        
        lab_exp = expression(paste("Expected -log"[10], plain(P)))
        lab_obs = expression(paste("Observed -log"[10], plain(P)))
        
        qq_plot <- ggplot(tbl.assoc) +
            theme_bw() +
            geom_ribbon(aes(x=exp, ymin=ci_up, ymax=ci_low), fill="grey", alpha=.5) +
            geom_point(aes(exp, obs), shape=1, size=3) +
            geom_abline(intercept=0, slope=1, alpha=0.8 ,color='red') +     
            xlab(lab_exp) + ylab(lab_obs) +
            theme(axis.title=element_text(size=16)) +
            theme(axis.text=element_text(size=12))
        
        
        tiff(output_qq, compression="lzw")
            print(qq_plot)
        dev.off()
    })[3]
    
    cat(paste(" in",runtime," s\n"))
}


#=== Calculate GC lambda

lambda = median(qchisq(tbl.assoc$p,df = 1,lower.tail = F),na.rm = T)/qchisq(0.5,df = 1) 
print(paste0("File: ", results, " => GC Lamda = ", lambda))


#=== Finish

rm(list=ls())

