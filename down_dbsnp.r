# Scrape essential information of queried SNP from dbSNP
# Nan Wang
# Univeristy of Southern California
# 3/10/2015

library('XML')
library('RCurl')

down_dbsnp <- function(snp){
    url = paste0('http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=',snp)
    dbsnp = getURL(url)
    
    doc = htmlParse(dbsnp)
    tabs = readHTMLTable(doc,header = F,stringsAsFactors = F)
    
    info <<- as.data.frame(tabs[[11]])[2:14,1:2]    #general info
    
    div <<- data.frame(tabs$Diversity)[-1,-12]   #diversity
    names(div) = div[1,]
    div = div[div$Pop != '',]
    div = div[-1,]
    #div[,c(2,10,11)]
}

down_dbsnp('rs8078000')
down_dbsnp('rs13342692')

a = sapply(tabs[[6]],function(x) grep('Chr Pos',unlist(x)))
which(sapply(a,length)!=0)


