set.seed(1000)
setwd("C:/Users/timot/Dropbox/chpt_2/rer_converge")
#BiocManager::install("fgsea")
library(fgsea)
library(dplyr)
library(tidyr)
library(tibble)

input_names <- c("polyandry","anti_repeat",
                 "anti_totipotent","col_size",
                 "crossover_per_chromosome",
                 "queen_age","anti_polygyny")
#read in rer converge results for each trait
input_names <- sapply(input_names,
                      function(x)paste("rer_results_",x,sep = ""),
                      USE.NAMES = F)

#the rer converge results for each predictive trait and gwrr
rer_trait_list <- rep(list(NULL),7)
for (i in 1:7) {
  rer_trait_list[[i]] <- read.table(paste("rer_data/",
                                          input_names[i],
                                          sep = ""))
}

#list of orthogroups which are rho positive
# and ranked in descending order of -log10(permulation p-val)
#for each predictive trait and gwrr
ranked_list <- rep(list(NULL),7)
for(i in 1:7){
  #list of genes which are evolving more rapidly with GWRR
  rer_res_pos <- rer_trait_list[[i]][rer_trait_list[[i]]$Rho>0,] %>% 
    filter(!is.na(Rho))
  #make the ranked named -log10(p) values for the rho pos genes
  rer_res_pos_p <- rer_res_pos$permpval + .00001
  names(rer_res_pos_p) <- rownames(rer_res_pos)
  ranked_list[[i]] <- -log10(sort(rer_res_pos_p,decreasing=F))
}



#make list where GO term is element name 
#and an element is each gene tied to that GO term
genes_go_terms <- read.table("ids_go_list.tsv",fill = T)

#each go term in a separate row; if it corresponds to one or more genes
#it is represented in that many rows
go_surjective <- genes_go_terms %>%  
  separate_rows(V2, sep=",") %>% 
  filter(V2 != "")
#make list entries containing genes corresponding to a single go term
go_list <- go_surjective %>% 
  group_by(V2) %>%  
  summarize(genes = list(V1)) %>% 
  deframe()

fgsea_list <- rep(list(NULL),7)
for (i in 1:7) {
#gene set enrichment analysis
fgsea_list[[i]] <- fgsea(
  pathways = go_list,  # Named list of GO terms to genes
  stats = ranked_list[[i]],                 # Named vector of transformed p-values or scores
  minSize = 10,                  # Minimum gene set size
  maxSize = 500                  # Maximum gene set size
)
}

for (i in 1:7) {
  print(fgsea_list[[i]]$pathway[which(fgsea_list[[i]]$padj < .05)])
}
