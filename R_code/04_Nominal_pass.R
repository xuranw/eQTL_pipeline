## Pseudobulk
bed_files_fastqtl <- list.files('/sc/arion/projects/psychAD/data/share/cmu/pseudobulk',pattern='.bed.gz$',full.names = T)

read_fastQTL_write_QTLtools <- function(i){
  output_file_name <- gsub('bed.gz','qtltools.bed',bed_files_fastqtl[i]) %>% 
    gsub('/sc/arion/projects/psychAD/data/share/cmu/pseudobulk','/sc/arion/projects/psychAD/data/share/cmu/pseudobulk/PC5_indep/',.)
  
  fasqtl <- read_tsv(bed_files_fastqtl[i]) %>% 
    mutate(pid='.',strand='+') %>% 
    dplyr::select(`#Chr`:ID,pid,strand,everything())
  
  dir.create('/sc/arion/projects/psychAD/data/share/cmu/pseudobulk/PC5_indep/',recursive = TRUE, showWarnings = FALSE)
  write_tsv(fasqtl, output_file_name) 
}
lapply(1:length(bed_files_fastqtl), read_fastQTL_write_QTLtools)

### Sorted 

bed_files_fastqtl <- list.files('/sc/arion/projects/psychAD/data/share/cmu/sorted',pattern='.bed.gz$',full.names = T)

read_fastQTL_write_QTLtools <- function(i){
  output_file_name <- gsub('bed.gz','qtltools.bed',bed_files_fastqtl[i]) %>% 
    gsub('/sc/arion/projects/psychAD/data/share/cmu/sorted','/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep/',.)
  
  fasqtl <- read_tsv(bed_files_fastqtl[i]) %>% 
    mutate(pid='.',strand='+') %>% 
    dplyr::select(`#Chr`:ID,pid,strand,everything())
  
  dir.create('/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep/',recursive = TRUE, showWarnings = FALSE)
  write_tsv(fasqtl, output_file_name) 
}
lapply(1:length(bed_files_fastqtl), read_fastQTL_write_QTLtools)


### bMIND_small
bed_files_fastqtl <- list.files('/sc/arion/projects/psychAD/data/share/cmu/bmind_small',pattern='.bed.gz$',full.names = T)

read_fastQTL_write_QTLtools <- function(i){
  output_file_name <- gsub('bed.gz','qtltools.bed',bed_files_fastqtl[i]) %>% 
    gsub('/sc/arion/projects/psychAD/data/share/cmu/bmind_small','/sc/arion/projects/psychAD/data/share/cmu/bmind_small/PC10_indep/',.)
  
  fasqtl <- read_tsv(bed_files_fastqtl[i]) %>% 
    mutate(pid='.',strand='+') %>% 
    dplyr::select(`#Chr`:ID,pid,strand,everything())
  
  dir.create('/sc/arion/projects/psychAD/data/share/cmu/bmind_small/PC10_indep/',recursive = TRUE, showWarnings = FALSE)
  write_tsv(fasqtl, output_file_name) 
}
lapply(1:length(bed_files_fastqtl), read_fastQTL_write_QTLtools)


### bMIND_full
bed_files_fastqtl <- list.files('/sc/arion/projects/psychAD/data/share/cmu/bulk_bmind',pattern='.bed.gz$',full.names = T)

read_fastQTL_write_QTLtools <- function(i){
  output_file_name <- gsub('bed.gz','qtltools.bed',bed_files_fastqtl[i]) %>% 
    gsub('/sc/arion/projects/psychAD/data/share/cmu/bulk_bmind/','/sc/arion/projects/psychAD/data/share/cmu/bulk_bmind/PC5_indep/',.)
  
  fasqtl <- read_tsv(bed_files_fastqtl[i]) %>% 
    mutate(pid='.',strand='+') %>% 
    dplyr::select(`#Chr`:ID,pid,strand,everything())
  
  dir.create('/sc/arion/projects/psychAD/data/share/cmu/bulk_bmind/PC5_indep/',recursive = TRUE, showWarnings = FALSE)
  write_tsv(fasqtl, output_file_name) 
}
lapply(1:length(bed_files_fastqtl), read_fastQTL_write_QTLtools)





file.dir = '/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep/'
perm.txt = list.files('/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep/', pattern = 'permute.txt', full.names = T)
perm.22.txt = list.files('/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep/', pattern = '22_22', full.names = T)
perm.14.txt = list.files('/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep/', pattern = '14_22', full.names = T)

for(i in 1:4){
  temp = read.table(perm.txt[i], header = T)
  temp1 = read.table(perm.22.txt[i], header = F)
  temp2 = read.table(perm.14.txt[i], header = F)
  colnames(temp1) <- colnames(temp2) <- colnames(temp)
  gid1 = setdiff(temp1$phe_id, temp$phe_id)
  gid2 = setdiff(temp2$phe_id, temp$phe_id)
  temp = rbind(temp, temp1[temp1$phe_id %in% gid1, ], temp2[temp2$phe_id %in% gid2, ])
  write.table(temp, file = perm.txt[i], row.names = F, quote = F)
}


# bed_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/pseudobulk/PC5_indep',pattern='.qtltools.bed.gz$',full.names = T)
bed_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/sorted/PC5_indep',pattern='.qtltools.bed.gz$',full.names = T)
bed_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/bulk_bmind/PC5_indep',pattern='.qtltools.bed.gz$',full.names = T)
bed_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/bmind_small/PC10_indep',pattern='.qtltools.bed.gz$',full.names = T)

read_fastQTL_write_QTLtools_conditional <- function(i){
  output_file_name <- gsub('qtltools.bed.gz','qtltools.conditional.bed',bed_files[i])
  
  significant_hits <- read_delim(gsub('qtltools.bed.gz','permutations_all.significant.txt',bed_files[i]),delim = ' ',col_names=FALSE)
  
  fasqtl <- read_tsv(bed_files[i]) %>% 
    filter(ID%in%significant_hits$X1)
  
  threshold <- read_delim(gsub('qtltools.bed.gz','permutations_all.thresholds.txt',bed_files[i]),
                          delim = ' ',col_names=FALSE) %>% 
    filter(X1%in%significant_hits$X1)
  
  write_tsv(fasqtl,output_file_name) 
  write_delim(threshold,gsub('qtltools.bed.gz','permutations_all.thresholds.txt',bed_files[i]),col_names=FALSE,delim = ' ') 
}
lapply(1:length(bed_files), read_fastQTL_write_QTLtools_conditional)

