# PubMed text data mining
# validate pathways

path0=read.csv("path_identified.csv")
path0=path0[,2]
path=as.matrix(path0)

vdp <- matrix(NA, nrow = length(path), ncol = 3);
rownames(vdp) <- path;
colnames(vdp) <- c("BC","basal PMIDs","luminal PMIDs")

# install.packages("easyPubMed")
library(easyPubMed)

for ( i in 1 : nrow(path)){
    pid=path[i,1]

    # pathway name
    url_1="https://www.genome.jp/kegg-bin/show_pathway?hsa0"
    p=(paste(url_1,pid,sep=""))
    # pURL=browseURL(p)
    r=readLines(p,n=4)
    pn=r[4]
    pname=substr(pn,12,nchar(pn)-22)

    # validate
    # url_2="https://www.ncbi.nlm.nih.gov/pubmed/?term="
    # pm=(paste(url_2,gname,sep=""))

    # all cancers: text data mining
    dami_query_string=(paste(pname,"basal-like breast cancer prognosis marker",sep=" AND "))
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
    dami_query_string2=(paste(pname,"luminal breast cancer prognosis marker",sep=" AND "))
    # dami_query_string2=(paste(pname,"skin cancer prognosis marker",sep=" AND "))
    # dami_query_string2=(paste(pname,"lung cancer prognosis marker",sep=" AND "))
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

    vdp[i,1] <- pname
    vdp[i,2] <- pmid
    vdp[i,3] <- pmid2
}

write.csv(vdp,"valid_path.csv")
