suppressMessages(library(stringr))
suppressMessages(library(dplyr))

Readdat <- function(inpath,idpath){
  print('Start reading and merging')
  ExpDat <- read.table(inpath,header = T,sep = "\t")
  colnames(ExpDat) <- gsub(colnames(ExpDat), pattern = '[.]', replacement = '_')
  names(ExpDat)[1] <- 'Ensembl_ID'
  symbol_id <- read.table(idpath,header = T,sep = '\t')
  names(symbol_id)[1] <- 'Ensembl_ID'
  ExpDat <- inner_join(symbol_id,ExpDat,by = 'Ensembl_ID')
  ExpDat <- aggregate(x = ExpDat[,7:ncol(ExpDat)],
                      by = list(ExpDat$gene),
                      FUN = median)
  rownames(ExpDat) <- ExpDat[,1]
  ExpDat <- ExpDat[,-1]
  ExpDat <- as.data.frame(ExpDat)
  print('successful')
  
  return(ExpDat)
}


ReadClinc <- function(inpath){
  print('Import clinical data')
  InClinc <- read.table(inpath,header = T,sep = '\t',quote = '')
  names(InClinc)[1] <- 'submitter'
  InClinc$submitter <- gsub(InClinc$submitter, pattern = '[-]', replacement = '_')
  InClinc$type <- ifelse(as.numeric(substr(InClinc$submitter,14,15)) > 9,'Normal','Tumor')
  rownames(InClinc) <- InClinc$submitter
  print('successful')
  
  return(InClinc)
}