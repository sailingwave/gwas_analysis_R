# Scrape essential information of queried SNP from dbSNP
# Nan Wang
# Univeristy of Southern California
# 3/10/2015

library('XML')
library('RCurl')

down_dbsnp <- function(snp){
    #from NCBI website
    url = paste0('http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=',snp)
    dbsnp = getURL(url)
    
    doc <- htmlParse(dbsnp)
    tabs <- readHTMLTable(doc,header = F,stringsAsFactors = F)
    
    info <- as.data.frame(tabs[[11]])[2:14,1:2]    #general info
    
    div <- data.frame(tabs$Diversity)[-1,-12]   #diversity
    names(div) = div[1,]
    div = div[div$Pop != '',]
    div = div[-1,]
    #div[,c(2,10,11)]
    
    #from NCBI eutils
    id = gsub('rs','',snp)
    url = paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&report=XML&id=',id)
    
    file = getURL(url)
    doc <- xmlParse(file)
    rt = xmlRoot(doc)
    ns = xmlNamespaces(doc, simplify = TRUE)
    names(ns)[names(ns) == ""] <- 'x'        #name x as the default namespace
    
    assembly = xpathSApply(rt,"//x:Assembly",xmlAttrs,namespaces = ns)
    chr = xpathSApply(rt,"//x:Component",xmlGetAttr,name="chromosome",namespaces = ns)
    pos = xpathSApply(rt,"//x:MapLoc",xmlAttrs,namespaces = ns)
    
    assem_info <- t(rbind(assembly,chr,pos))[,c(1:3,6,11,12,15)]
    
    result <- list(info=info,diversity=div,assembly=assem_info)
}



#usage

info = down_dbsnp('rs8078000')    #create info and div in the global env
info[[3]]


#dev
down_dbsnp('rs8078000')

#a = sapply(tabs,function(x) grep('7044260',unlist(x)))
#which(sapply(a,length)!=0)

#test email address



