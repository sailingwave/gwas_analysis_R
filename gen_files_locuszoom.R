# Two purposes:
# 1. generate locuszoom input files
# 2. find the SNPs of interest with lowest p values according to a file 'interest_list'
# 06/21/2015
# 
# Required input files:
# 1. 'interest_list' with 3 columns: filename, chr, number
# 2. 'snp_id' with 2 columns: cptid, MarkerName
# 3. input files for gwas_sum.r with 3 columns: chr, pos, p


library('data.table')
library('dplyr')

uniq_sig_window = 20000    #up/down bp window size for a 'unique' signal
pref_out = 'locuszoom_'    #prefix for output
log_file = 'locuszoom_fileprep.log'

sink(log_file,split = T)


interest_list <- read.table('snp_interest.txt',header=T,stringsAsFactors = F)
snp_id <- fread('snp_id.txt',header=T,showProgress = F)
setnames(snp_id,names(snp_id),c('cptid','MarkerName'))
setkey(snp_id,cptid)

files = unique(interest_list$filename)

for(f in files){
    cat(paste("==> Processing",f,"...\n"))
    
    sub_list = interest_list[interest_list$filename==f,]
    max_num = max(sub_list[,'number'])
    
    res <- fread(f,showProgress = F)
    tbl_res <- tbl_dt(res)
    chrs = sub_list[,'chr']
    sub_set <- tbl_res %>%
        filter(chr %in% chrs) %>%
        mutate(cptid=paste(chr,pos,sep=":"))    
    setkey(sub_set,cptid)
    
    sub_set <- merge(snp_id,sub_set)
    
    #output input files for locuszoom
    out_data <- sub_set%>%
        select(MarkerName,p)
    setnames(out_data,names(out_data),c('MarkerName','P-value'))
    
    out_file = paste0(pref_out,gsub('\\.[^.]*$','',f),'.txt')
    write.table(out_data,out_file,quote=F,row.names=F,sep="\t")
    
    #find index SNPs
    for(i in 1:max_num){
        cat(paste(" -> for SNP ",i,"of each chr ...\n"))
        
        if(i > 1){
            rm_snp = selected[,c('chr','pos'),with=F]
            setnames(rm_snp,names(rm_snp),c('chr','pos1'))
            #remove previously-selected signals
            sub_set <- merge(sub_set,rm_snp,by='chr') %>%
                filter(pos<pos1-uniq_sig_window | pos>pos1+uniq_sig_window) %>%
                select(-pos1)
        }
        
        selected <- sub_set %>%
            filter(chr %in% chrs) %>%
            group_by(chr) %>%
            filter(p==min(p)) %>%
            select(chr,MarkerName,p,pos) %>%
            arrange(chr)
        
        print(selected)        
        
        selected <- selected %>%
            group_by(chr) %>%
            filter(pos==min(pos))    #break ties if there is any
        
        sub_list$number = sub_list$number-1
        chrs = sub_list[sub_list$number>0,'chr']
    }
}

sink()
