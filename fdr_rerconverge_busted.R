set.seed(1000)
setwd("C:/Users/timot/Dropbox/chpt_2/rer_converge")
#the busted p-vals for each single copy ortholog which had all 13 species present in dna alignments
busted_pvals <- read.table("orth_ids_pval_all_single_copy.txt")

input_names <- c("polyandry","anti_repeat",
                 "anti_totipotent","col_size",
                 "crossover_per_chromosome",
                 "queen_age","anti_polygyny")
input_names <- sapply(input_names,
                      function(x)paste("rer_results_",x,sep = ""),
                      USE.NAMES = F)
rer_trait_list <- rep(list(NULL),7)

for (i in 1:7) {
  rer_trait_list[[i]] <- read.table(paste("rer_data/",
                                          input_names[i],
                                          sep = ""))
}

#list of single copy orthogroups that are rho+ significant for each trait
sig_rho_pos_list <- rep(list(NULL),7)
busted_fdr_rho_pos_list <- rep(list(NULL),7)
index_val <- 0
#ids for significant rho+ gwrr genes
for (i in rer_trait_list) {
  trait_rho_pos <- which(i$Rho>0)
  trait_sig <- which(i$permpval<.05)
  index_val <- index_val +1
  sig_rho_pos_list[[index_val]] <- rownames(i[intersect(trait_sig,
                                                        trait_rho_pos),])
  #get the busted p-values of each significant rho+ ortholog for a given trait
  p_val_locs <- na.omit(match(sig_rho_pos_list[[index_val]],busted_pvals$V1))
  #p-values for rho+ significant  each trait with orthid
  busted_rho_pos_subset <- busted_pvals[p_val_locs,]
  #location of fdr corrected busted p-vals < .05
  fdr_locs <- which(p.adjust(busted_rho_pos_subset$V2,method = "BH")<.05)
  busted_fdr_rho_pos_list[[index_val]] <- busted_rho_pos_subset$V1[fdr_locs]
  }
#so these are the orthogroups for each trait that are rer rho+ significant
#before fdr correction, and also busted significant after FDR correction


#list of orthids that have a positive rho value with a significant permpval
#and a significant busted p val for GWRR 
#and each of the four pedictive traits of strongest model

#which genes are rho positive fdr significant for each trait
fdr_rho_sig <- rep(list(NULL),7)
index_val <- 0
for (i in rer_trait_list) {
  trait_rho_pos <- which(i$Rho>0)
  #corrected p-val below
  trait_sig_fdr <- which(i$permpvaladj<.05)
  index_val <- index_val+1
  fdr_rho_sig[[index_val]] <- rownames(i[intersect(trait_rho_pos,trait_sig_fdr),])
}


rho_pos_anti_polygyny <-
  rownames(rer_trait_list[[7]])[intersect(which(rer_trait_list[[7]]$Rho<0),which(rer_trait_list[[7]]$permpval <.05))]
#25 rho pos genes significant for busted prior to fdr correction
intersect(busted_ids_sig,rho_pos_anti_polygyny)
