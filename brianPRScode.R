library(snpStats)
library(data.table)
library(stringr)

matchSNP<-function(snps,table){
  matches<-list()
  for (snp in snps) matches[[snp]]<-grep(snp,colnames(table))
  return(unlist(matches))}

setwd("/Volumes/Brian Mac Backup/REQUITE")

#already edited RA_OR_European file
snps_info<-read.csv(file="RA_OR_European.csv", row.names=1)
sample_info<-read.csv(file="excite_ids_associated_data.csv")

setwd("/Volumes/RAPPER_Oncoarray/REQUITE Probabilities")

#Change name for specific studies
list_files<-list.files(pattern="excite_ids.oncoarray*")[-16]
samples<-as.character(unlist(sample_info$Id2))


snp_select<-list()
allele_dose<-list()

#read.impute has a problem with chr23, remove chr23 from analysis
alleleDose<-function(list_files){
  
  lapply(list_files, function(file){
    
    chr_code<-strsplit(file,'_')[[1]][4]
    chr_number<-substr(chr_code,4,str_count(chr_code))
    
    print(chr_number)
    
    excite_chr<-read.impute(file= file, rownames=samples,nsnp=NULL,snpcol=2)
    
    #SNPs on chr12 given as positions not rs, use strsplit index 2 for this chromosome...
    if(chr_number=='12') {
      colnames(excite_chr)<-unlist(lapply(colnames(excite_chr),function(x) strsplit(x,':')[[1]][2]))	
    } else {
      colnames(excite_chr)<-unlist(lapply(colnames(excite_chr),function(x) strsplit(x,':')[[1]][1]))}
    
    snps_chr<-snps_info$Modified_SNP[snps_info$Chr==chr_number]
    matches<-which(colnames(excite_chr)%in%snps_chr)
    
    print(matches)
    
    
    #some problematic loci...
    #two matches to rs2236668, second is "rs2236668:45650009:T:CG" so remove...
    if(chr_number==21) {matches<-matches[-6]	}
    
    #two matches; one missing SNP rs9373594
    if (chr_number==6) {matches<-matches[-2]}
    
    #rs73081554 is coded as "chr3_58302935_C_T"...
    if (chr_number==3) {matches<-c(matches,453986)}
    #rs67250450 is coded as "chr7_28174986_C_T"...
    if (chr_number==7) {matches<-c(matches,264776)}
    #rs67250450 is coded as "chr7_28174986_C_T"...
    if (chr_number==10){matches<-c(matches,493233)} else {
      matches<-matches}
    
    allele_dose[[chr_number]]<-t(as(excite_chr[,matches],"numeric"))
    
  })
}



complete_allele_dose<-do.call(rbind.data.frame,allele_dose)

#need to rename some of the SNPs...
rownames(complete_allele_dose)[23]<-"rs71508903"#chr10
rownames(complete_allele_dose)[32]<-"rs773125"#chr12
rownames(complete_allele_dose)[33]<-"rs1633360"#chr12
rownames(complete_allele_dose)[34]<-"rs10774624"#chr12
rownames(complete_allele_dose)[35]<-"rs9603616"#chr13
rownames(complete_allele_dose)[60]<-"rs4239702"#chr20
rownames(complete_allele_dose)[72]<-"rs73081554"#chr3
rownames(complete_allele_dose)[93]<-"rs67250450"#chr7



prs<-apply(complete_allele_dose,2,sum,na.rm=T)
names(prs)<-colnames(complete_alle_dose)

present<-rownames(snps_info)%in%rownames(complete_allele_dose)

lor<-log(as(unlist(lapply(snps_info$OR.95.CI[present],function(x) strsplit(as.character(unlist(x)),'[()]')[[1]][1])),"numeric"))
names(lor)<-rownames(complete_allele_dose)
wprs<-apply(complete_allele_dose,2,function(x) sum(lor%*%x))
names(wprs)<-colnames(complete_allele_dose)


excite_allele_dose<-list(complete_allele_dose,prs,wprs)
names(excite_allele_dose)<-c('allele_dose','prs','wprs')

