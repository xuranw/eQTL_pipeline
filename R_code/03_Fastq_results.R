## Get FastQTL results
library(qvalue)
library(tidyverse)
# Sorted only need PC80
combine_fastqtl_out <- function(fastqtl_out_files){
  d <- tibble(files=fastqtl_out_files,
              cell_type = basename(files) %>% 
                gsub('.quantile.txt.gz..+','',.) %>% 
                gsub('\\.',' ',.) %>% 
                gsub('OPCs   COPs','OPCs / COPs',.)) %>% 
    mutate(PCs=dirname(files) %>% basename(.)) %>% 
    mutate(file_content = map(files, read.table, header=F)) %>% 
    unnest(file_content) %>% 
    filter(!is.na(V11)) %>% 
    dplyr::select(-files) %>% 
    setNames(c("cell_type","PCs","pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope","ppval", "bpval")) %>% 
    group_by(cell_type,PCs) %>% 
    mutate(adj_p=qvalue(bpval)$qvalues) %>% 
    ungroup() %>% 
    arrange(cell_type,bpval)
  return(d)
}

fastqtl_out_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/sorted',pattern='.gz.[0-9]',full.names = T,recursive = TRUE)
d <- combine_fastqtl_out(fastqtl_out_files)
# dir.create('/sc/arion/projects/psychAD/data/cmu/fastqtl_output/sorted/',showWarnings = FALSE,recursive = TRUE)
write_tsv(d, file = '/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/sorted_eqtl.allPCs.txt')
#d0 = read_tsv('/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/sorted_eqtl.allPCs.txt')
selected_pcs <- d %>% 
  filter(adj_p < 0.05) %>% 
  dplyr::count(PCs) %>% filter(n==max(n))

# selected_pcs 5

d %>% 
  filter(PCs==selected_pcs$PCs) %>% 
  dplyr::select(-PCs) %>% 
  write_tsv(.,paste0('/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/sorted_eqtl.',selected_pcs$PCs,'.txt'))

###### Pseudobulk  #####
fastqtl_out_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/pseudobulk',pattern='.gz.[0-9]',full.names = T,recursive = TRUE)
d <- combine_fastqtl_out(fastqtl_out_files)
write_tsv(d, file = '/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/pseudobulk_eqtl.allPCs.txt')

selected_pcs <- d %>% 
  filter(adj_p < 0.05) %>% 
  dplyr::count(PCs) %>% filter(n==max(n))
# Select PC5

d %>% 
  filter(PCs==selected_pcs$PCs) %>% 
  dplyr::select(-PCs) %>% 
  write_tsv(.,paste0('/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/pseudobulk_eqtl.',selected_pcs$PCs,'.txt'))


##### bMIND small #####
# We do not include Mic as well 
fastqtl_out_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/bmind_small',pattern='.gz.[0-9]',full.names = T,recursive = TRUE)
fastqtl_out_files <- fastqtl_out_files[-grep('Mic', fastqtl_out_files)]

d <- combine_fastqtl_out(fastqtl_out_files)
write_tsv(d, file = '/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/bmind_small_eqtl.allPCs.txt')

selected_pcs <- d %>% 
  filter(adj_p < 0.05) %>% 
  dplyr::count(PCs) %>% filter(n==max(n))

# PC 10

d %>% 
  filter(PCs==selected_pcs$PCs) %>% 
  dplyr::select(-PCs) %>% 
  write_tsv(.,paste0('/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/bmind_small_eqtl.',selected_pcs$PCs,'.txt'))

#### bMIND full #####
fastqtl_out_files <- list.files('/sc/arion/projects/psychAD/data/share/cmu/bulk_bmind/',pattern='.gz.[0-9]',full.names = T,recursive = TRUE)
d <- combine_fastqtl_out(fastqtl_out_files)

write_tsv(d, file = '/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/bulk_bmind_eqtl.allPCs.txt')

selected_pcs <- d %>% 
  filter(adj_p < 0.05) %>% 
  dplyr::count(PCs) %>% filter(n==max(n))

# selected_pcs 5

d %>% 
  filter(PCs==selected_pcs$PCs) %>% 
  dplyr::select(-PCs) %>% 
  write_tsv(.,paste0('/sc/arion/projects/psychAD/data/share/cmu/fastqtl_output/bulk_bmind_eqtl.',selected_pcs$PCs,'.txt'))


