# PubMed text data mining
# validate genes

genes0=read.csv("genes_identified.csv")
genes0=genes0[,2]
genes=as.matrix(genes0)

g=as.vector(genes)
mode(g)="character"
# install.packages("org.Hs.eg.db")
library(org.Hs.eg.db)
g_info=select(org.Hs.eg.db,keys = g,columns=c("ENTREZID","SYMBOL","GENENAME"),keytype="ENTREZID")
gnam=paste(g_info$SYMBOL,g_info$GENENAME,sep=" ")
g_matrix=as.matrix(gnam)
rownames(g_matrix)= g_info$ENTREZID

vdg <- matrix(NA, nrow = length(g_matrix), ncol = 3);
# vdg <- matrix(NA, nrow = length(g_matrix), ncol = 1);
rownames(vdg) <- rownames(g_matrix);
colnames(vdg) <- c("names","all cancers PMIDs","cancer PMIDs")
# colnames(vdg) <- c("names")

# install.packages("easyPubMed")
library(easyPubMed)

for ( i in 1 : nrow(g_matrix)){
    gname=g_matrix[i,1]
    
    # gene name
    # url_1="https://www.ncbi.nlm.nih.gov/gene"
    # g=(paste(url_1,gene,sep="/"))
    # gURL=browseURL(g)
    # r=readLines(g)
    # gn=r[11]
    # gname=substr(gn,12,nchar(gn)-22)
    
    # validate
    # url_2="https://www.ncbi.nlm.nih.gov/pubmed/?term="
    # pm=(paste(url_2,gname,sep=""))
    
    # all cancers: text data mining
    dami_query_string=(paste(gname,"cancer",sep=" AND "))
    dami_on_pubmed=get_pubmed_ids(dami_query_string)
    num_pmid=(dami_on_pubmed$Count)
    if(num_pmid>1){
        pmid1=(unlist(dami_on_pubmed$IdList))
        pmid=pmid1[1]
    }
    if(num_pmid==1){
        pmid=(unlist(dami_on_pubmed$IdList))
    }
    if(num_pmid==0){
        pmid=0
    }
    
    # specific cancer: text data mining
    dami_query_string2=(paste(gname,"breast cancer prognosis marker",sep=" AND "))
    # dami_query_string2=(paste(gname,"skin cancer prognosis marker",sep=" AND "))
    # dami_query_string2=(paste(gname,"lung cancer prognosis marker",sep=" AND "))
    dami_on_pubmed2=get_pubmed_ids(dami_query_string2)
    num_pmid2=(dami_on_pubmed2$Count)
    if(num_pmid2>1){
        pmid12=(unlist(dami_on_pubmed2$IdList))
        pmid2=pmid12[1]
    }
    if(num_pmid2==1){
        pmid2=(unlist(dami_on_pubmed2$IdList))
    }
    if(num_pmid2==0){
        pmid2=0
    }
        
    vdg[i,1] <- gname
    vdg[i,2] <- pmid
    vdg[i,3] <- pmid2
}

write.csv(vdg,"valid_gene.csv")
summary(vdg[,2]!=0)
summary(vdg[,3]!=0)
