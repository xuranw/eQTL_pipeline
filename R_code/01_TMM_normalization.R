# TMM normalization of pseudo-bulk data

BiocManager::install("rtracklayer")
library(rtracklayer)
gtf <- rtracklayer::import('./data_store/Homo_sapiens.GRCh38.96.gtf') %>% 
  as.data.frame() %>% 
  dplyr::filter(type=='gene') %>% 
  mutate(TSS_start=ifelse(strand=='+',start,end),
         TSS_end=ifelse(strand=='+',start+1,end+1)) %>% 
  mutate(gene=paste0(gene_name,'_',gene_id)) %>% 
  dplyr::select(seqnames,TSS_start,TSS_end,gene) %>% 
  filter(seqnames%in%c(1:22))
gtf$seqnames <- droplevels(gtf$seqnames)
gtf$ensembl = unlist(lapply(strsplit(gtf$gene, '_'), function(x){x[2]}))
gtf$hgnc = unlist(lapply(strsplit(gtf$gene, '_'), function(x){x[1]}))

saveRDS(gtf, file = './new_data_store/gtf.RDS')

load('new_data_store/match_subject_names.RData')
# df.subject

##################### Pseudobulk #########################
# Preprocess: filter genes and samples
library(tidyverse)
library(SingleCellExperiment)
library(parallel)
sc.seu = readRDS('./new_data_store/sc_seurat_object.rds')
sc.seu = sc.seu[, which(sc.seu$subject_cmc != '#N/A')]
# Let's skip the function and repeat the procedure ourselves
cell_type_id <- cbind(as.character(sc.seu$celltype), sc.seu$subject_cmc) %>% unique() %>% as.data.frame()
cell_type_id <- cell_type_id[cell_type_id$V1 %in% c('In', 'Ex', 'Oli', 'Mic', 'Ast') & !(cell_type_id$V2 == '#N/A'), ]

sum_expression = NULL
for(i in 1:nrow(cell_type_id)){
  cnt.mtx = sc.seu@assays$RNA@counts[, sc.seu$celltype == cell_type_id$V1[i] & sc.seu$subject_cmc == cell_type_id$V2[i]]
  if(is.null(dim(cnt.mtx))){
    tot_counts = cnt.mtx
    n_expressed = as.numeric(cnt.mtx > 0)
    n_cells = 1
  }else{
    tot_counts <- Matrix::rowSums(cnt.mtx)
    n_expressed <- Matrix::rowSums(cnt.mtx > 0)
    n_cells <- ncol(cnt.mtx)
  }
  cnt.df = data.frame(counts = tot_counts, n_expressed = n_expressed, n_cells = n_cells) %>%
    rownames_to_column('ensembl') %>%
    mutate(perc.exprs = round(n_expressed*100/n_cells, digits = 2), 
           cell.type = cell_type_id$V1[i], 
           subject.cmc = cell_type_id$V2[i]) %>%
    dplyr::select(cell.type, subject.cmc, ensembl, counts, n_expressed, n_cells, perc.exprs)
  sum_expression = rbind(sum_expression, cnt.df)
}

sum_expression <- sum_expression %>% 
  group_by(cell.type, subject.cmc) %>% 
  mutate(counts_scaled_1M=counts*10^6/sum(counts)) %>% 
  ungroup() %>% 
  as_tibble()

sumstats_gene <- sum_expression %>% group_by(cell.type, ensembl) %>% 
  summarise(counts_per_celltype = sum(counts),
            n_expressed_per_celltype = sum(n_expressed),
            n_cells_per_celltype = sum(n_cells),
            perc_expressed_per_celltype = round(n_expressed_per_celltype*100/n_cells_per_celltype, digits=2),
            n_samples_expressed = sum(counts>0),
            perc_samples_expressed = round(n_samples_expressed*100/length(counts),digits=2),
            mean_cpm = mean(counts_scaled_1M)) %>% 
  group_by(ensembl) %>% 
  mutate(prop_max_cpm = mean_cpm/max(mean_cpm)) %>% 
  ungroup()
# We do have ambient records

# In total we have 99 subjects
# 1. Expressed in at least 10 individuals
# 2. Expressed with a mean CPM of at least 1
gene_to_keep_per_cell_type <- sumstats_gene %>% 
  filter(n_samples_expressed > 9, mean_cpm > 1) %>% 
  dplyr::select(cell.type, ensembl) %>% 
  unique()

sum_expression <- sum_expression %>% 
  inner_join(.,gene_to_keep_per_cell_type,by=c('cell.type','ensembl'))

(dplyr::count(sum_expression,cell.type,subject.cmc) %>% 
    dplyr::select(-subject.cmc) %>% 
    unique() %>% 
    arrange(-n))
# A tibble: 5 × 2
#cell.type     n
#<chr>     <int>
#1 Ex        17937
#2 In        17863
#3 Ast       17691
#4 Oli       16858
#5 Mic        9302

(n_cells_ind_cell_type <- sum_expression %>% 
    dplyr::select(cell.type, subject.cmc, n_cells) %>% 
    unique() %>% arrange(n_cells))
# Remove individuals with less than 10 cells for a given cell type.
keep <- filter(n_cells_ind_cell_type,n_cells>=10)
dplyr::count(keep, cell.type) %>% arrange(-n)
# A tibble: 5 × 2
#cell.type     n
#<chr>     <int>
#  1 Ex           97
#  2 In           94
#  3 Oli          92
#  4 Ast          89
#  5 Mic          17
sum_expression <- inner_join(sum_expression, keep, by=c('cell.type','subject.cmc','n_cells'))


get_fastqtl_pheno <- function(cell){
  #keep all expression for the cell type and patient that have genotype
  bed <- sum_expression %>% dplyr::filter(cell.type == cell) %>% 
    dplyr::select(subject.cmc, ensembl, counts) %>% 
    spread(subject.cmc, counts) %>% 
    column_to_rownames('ensembl') %>% 
    edgeR::DGEList(counts = .) %>% #TMM normalize the data using edgeR
    edgeR::calcNormFactors(.) %>% 
    edgeR::cpm(.) %>%
    as.data.frame() %>% 
    rownames_to_column('ensembl')
  bed <- bed %>% inner_join(.,gtf,by = 'ensembl') %>% 
    dplyr::select(seqnames,TSS_start,TSS_end,gene,everything()) %>% 
    arrange(seqnames,TSS_start)
  
  colnames(bed)[c(1,2,3,4)] <- c('#Chr','start','end','ID')
  bed <- bed[, -5]
  write_tsv(bed, path = paste0('./minerva_upload/pseudobulk/',make.names(cell),'.bed'))
} 

cells <- unique(sum_expression$cell.type)
TMM.pseudo.bulk = lapply(cells, get_fastqtl_pheno)
names(TMM.pseudo.bulk) = cells

p.gt = readRDS('./new_data_store/pca_gt.rds')
gt.bm = readRDS('./data_store/gt_bm_mtx.rds')
cells = c('Ex', 'In', 'Oli', 'Ast', 'Mic')

for(cell in cells){
  bed = read_tsv(paste0('./minerva_upload/pseudobulk/', make.names(cell), '.bed'))
  col4 = colnames(bed)[1:4]
  col.sub = colnames(bed)[-c(1:4)]
  col.sub = intersect(col.sub, df.subject$Individual_ID)
  cat(cell, '\t', length(col.sub), '\n')
  bed = bed[, c(col4, col.sub)]
  col.sub = df.subject$Genotyping_Sample_ID[match(col.sub, df.subject$Individual_ID)]
  colnames(bed) = c(col4, col.sub)
  write_tsv(bed, path = paste0('./minerva_upload/pseudobulk/',make.names(cell),'.bed'))
}


############################  Sorted #####################################
# Sorted data with normalization
geneCounts = readRDS('./new_data_store/cts_geneCounts_merged.RDS')  # counts
info = readRDS('./new_data_store/cts_info_merged.RDS')   # metadata
geneInfo = readRDS('./new_data_store/cts_geneInfo.RDS')
table(info$CellType)

library(reshape2)
#colnames(geneCounts) = info$Individual.ID
sorted.expression = as_tibble(melt(geneCounts))
colnames(sorted.expression) = c('ensembl', 'subject.cmc', 'counts')
sorted.expression$cell.type = info$CellType[match(sorted.expression$subject.cmc, info$Name)]
sorted.expression$subject.cmc = info$Individual.ID[match(sorted.expression$subject.cmc, info$Name)]

sorted.expression <- sorted.expression %>% 
  group_by(cell.type, subject.cmc) %>% 
  mutate(counts_scaled_1M=counts*10^6/sum(counts)) %>% 
  ungroup() %>% 
  as_tibble()

sumstats_gene <- sorted.expression %>% group_by(cell.type, ensembl) %>% 
  summarise(counts_per_celltype = sum(counts),
            n_samples_expressed = sum(counts>0),
            perc_samples_expressed = round(n_samples_expressed*100/length(counts),digits=2),
            mean_cpm = mean(counts_scaled_1M)) %>% 
  group_by(ensembl) %>% 
  mutate(prop_max_cpm = mean_cpm/max(mean_cpm)) %>% 
  ungroup()

# In total we have 101 subjects
# 1. Expressed in at least 9 individuals
# 2. Expressed with a mean CPM of at least 1
gene_to_keep_per_cell_type <- sumstats_gene %>% 
  filter(n_samples_expressed > 9, mean_cpm > 1) %>% 
  dplyr::select(cell.type, ensembl) %>% 
  unique()

sorted.expression <- sorted.expression %>% 
  inner_join(.,gene_to_keep_per_cell_type,by=c('cell.type','ensembl'))

(dplyr::count(sorted.expression,cell.type,subject.cmc) %>% 
    dplyr::select(-subject.cmc) %>% 
    unique() %>% 
    arrange(-n))
# A tibble: 4 × 2
#cell.type     n
#<fct>     <int>
# 1 MgAs      21855
# 2 GLU       21317
# 3 GABA      20728
# 4 Olig      19389


get_fastqtl_pheno_sorted <- function(cell){
  #keep all expression for the cell type and patient that have genotype
  bed <- sorted.expression %>% dplyr::filter(cell.type==cell) %>% 
    dplyr::select(subject.cmc, ensembl, counts) %>% 
    spread(subject.cmc, counts) %>% 
    column_to_rownames('ensembl') %>% 
    edgeR::DGEList(counts = .) %>% #TMM normalize the data using edgeR
    edgeR::calcNormFactors(.) %>% 
    edgeR::cpm(.) %>%
    as.data.frame() %>% 
    rownames_to_column('ensembl')
  
  bed <- bed %>% inner_join(.,gtf,by = 'ensembl') %>% 
    dplyr::select(seqnames,TSS_start,TSS_end,gene,everything()) %>% 
    arrange(seqnames,TSS_start)
  
  colnames(bed)[c(1,2,3,4)] <- c('#Chr','start','end','ID')
  bed <- bed[, -5]
  write_tsv(bed, path = paste0('./minerva_upload/sorted/', make.names(cell), '.bed'))
}    

cells.sorted <- unique(sorted.expression$cell.type)
TMM.sorted = lapply(cells.sorted, get_fastqtl_pheno_sorted)
names(TMM.sorted) = cells.sorted

cells.sorted = c('GABA', 'GLU', 'Olig', 'MgAs')
for(cell in cells.sorted){
  bed = read_tsv(paste0('./minerva_upload/sorted/', make.names(cell), '.bed'))
  col4 = colnames(bed)[1:4]
  col.sub = colnames(bed)[-c(1:4)]
  col.sub = intersect(col.sub, df.subject$Individual_ID)
  cat(cell, '\t', length(col.sub), '\n')
  bed = bed[, c(col4, col.sub)]
  col.sub = df.subject$Genotyping_Sample_ID[match(col.sub, df.subject$Individual_ID)]
  colnames(bed) = c(col4, col.sub)
  write_tsv(bed, path = paste0('./minerva_upload/sorted/',make.names(cell),'.bed'))
}


####################### Bulk #######################################
# Get TMM from bulk counts
load('data_store/bulk_match_snp_small_ready.RData')
bulk.expression = as_tibble(melt(bulk.count.mtx))
colnames(bulk.expression) = c('hgnc', 'subject.cmc', 'counts')

bulk.expression <- bulk.expression %>% 
  group_by(subject.cmc) %>% 
  mutate(counts_scaled_1M = counts*10^6/sum(counts)) %>% 
  ungroup() %>% 
  as_tibble()

sumstats_gene <- bulk.expression %>% group_by(hgnc) %>% 
  summarise(counts_per_celltype = sum(counts),
            n_samples_expressed = sum(counts>0),
            perc_samples_expressed = round(n_samples_expressed*100/length(counts),digits=2),
            mean_cpm = mean(counts_scaled_1M)) %>% 
  group_by(hgnc) %>% 
  mutate(prop_max_cpm = mean_cpm/max(mean_cpm)) %>% 
  ungroup()

gene_to_keep_per_cell_type <- sumstats_gene %>% 
  filter(mean_cpm > 1) %>% 
  dplyr::select(hgnc) %>% 
  unique()

bulk.expression <- bulk.expression %>% 
  inner_join(. , gene_to_keep_per_cell_type,by=c('hgnc'))

(dplyr::count(bulk.expression, subject.cmc) %>% 
    dplyr::select(-subject.cmc) %>% 
    unique() %>% 
    arrange(-n))

bed <- bulk.expression %>%
  dplyr::select(subject.cmc, hgnc, counts) %>% 
  spread(subject.cmc, counts) %>% 
  column_to_rownames('hgnc') %>% 
  edgeR::DGEList(counts = .) %>% #TMM normalize the data using edgeR
  edgeR::calcNormFactors(.) %>% 
  edgeR::cpm(.) %>%
  as.data.frame() %>% 
  rownames_to_column('hgnc')

bed <- bed %>% inner_join(.,gtf,by = 'hgnc') %>% 
  dplyr::select(seqnames,TSS_start,TSS_end,gene,everything()) %>% 
  arrange(seqnames,TSS_start)

colnames(bed)[c(1,2,3,4)] <- c('#Chr','start','end','ID')
bed <- bed[, c(1:4, 6:782)]

write_tsv(bed, path = paste0('./minerva_upload/bulk_bMIND/bulk.bed'))

bed = read_tsv('./minerva_upload/bulk_bMIND/bulk.bed')
col4 = colnames(bed)[1:4]
col.sub = colnames(bed)[-c(1:4)]
col.sub = intersect(col.sub, df.subject$Individual_ID)
cat(cell, '\t', length(col.sub), '\n')
bed = bed[, c(col4, col.sub)]
col.sub = df.subject$Genotyping_Sample_ID[match(col.sub, df.subject$Individual_ID)]
colnames(bed) = c(col4, col.sub)
write_tsv(bed, path = './minerva_upload/bulk_bMIND/bulk.bed')

############################## bMIND ##########################################
A = readRDS('new_data_store/A_tmm_bisque.rds')
cells.bmind = c('In', 'Ex', 'Oli', 'Ast', 'Mic')
cell = cells.bmind[1]
ensembl.names = rownames(A)

load('new_data_store/match_subject_names.RData')
load('./new_data_store/ind_idv.RData')

for(cell in cells.bmind){
  # Get TMM for each cell type
  # cell = cells.bmind[i]
  A.temp <- as_tibble(A[, cell, ])
  id.sub = which(rowMeans(A.temp) > 1) # Average CPM > 1
  A.temp = A.temp[id.sub, ]
  # colnames(A.temp) = df.subject$Genotyping_Sample_ID[match(colnames(A.temp), df.subject$Individual_ID)]
  A.temp$ensembl = ensembl.names[id.sub]
  
  bed <- A.temp %>% inner_join(., gtf, by = 'ensembl') %>%
    dplyr::select(seqnames,TSS_start,TSS_end,gene,everything()) %>% 
    arrange(seqnames,TSS_start)
  colnames(bed)[c(1,2,3,4)] <- c('#Chr','start','end','ID')
  
  bed <- bed[, c(1:4, which(colnames(bed) %in% df.subject$Individual_ID))]
  cat(cell, 'has ', nrow(bed), 'genes \n')
  write_tsv(bed, path = paste0('./minerva_upload/bulk_bMIND/', make.names(cell), '.bed'))
  
  # Get the independent samples
  # col4 = colnames(bed)[1:4]
  # col.sub = colnames(bed)[-c(1:4)]
  # id = which(df.subject$Individual_ID[match(col.sub, df.subject$Genotyping_Sample_ID)] %in% ind.idv)
  
  # bed = bed[, c(col4, col.sub[id])]
  # bed = bed[rowMeans(bed[, -c(1:4)]) > 1, ]
  
  # write_tsv(bed, path = paste0('./minerva_upload/bmind_small/', make.names(cell), '.bed'))
}


