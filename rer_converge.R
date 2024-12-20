set.seed(1000)
library(devtools)
library(RColorBrewer)
library(gplots)
library(phytools)
library(geiger)
library(knitr)
library(RcppArmadillo)
library(weights)
library(phangorn) 
library(impute)
library(RERconverge)
setwd("C:/Users/timot/Dropbox/chpt_2/rer_converge")
#make the trees for each orthogroup have the same structure...
#estimatePhangornTreeAll(alndir ="single_copy_all_aln_header", pattern = ".fa",
#                        treefile = "species_tree_rooted_sans_aas.txt", 
#                        output.file = "input_rerconverge.txt", submodel = "LG",
#                        type = "AA", format = "fasta", k = 4)
#read in gene trees of orthologues
ort_trees <- readTrees("input_rerconverge.txt", reestimateBranches = F, 
                       minSpecs = 10)
#rename from p rugosus to p barbatus to match ort_trees
#(p barbosus was a proxy genome/bc rugosus lacked annotation file)
trait_dat_rer <- read.csv("trait_dat.csv")
trait_dat_rer$Species <- trait_dat_rer$Species <- sub(" ","_",trait_dat_rer$Species)
trait_dat_rer$Species[10] <- "Pogonomyrmex_barbosus"
#run rer converge for GWRR and all predictors
#is there increase in evolutionary rate with increase in polyandry?
polyandry_dat <- trait_dat_rer$Polyandry
#decrease in repeat content
anti_repeat <- 100 - trait_dat_rer$repeat_content
#decrease in worker totipotency (3 is full totipotency 1 is sterile workers)
anti_totipotent <- 4 - trait_dat_rer$Totipotency
#increase in colony size
col_size <- as.numeric(trait_dat_rer$Colony.Size)
#then there is also GWRR (crossovers per chromosome)
cross_per_genome <- (trait_dat_rer$Recombination.Rate*trait_dat_rer$genome_size)/100
#divided by the number of chromosomes is the number of crossovers per chromosome
crossover_per_chromosome <- (cross_per_genome)/trait_dat_rer$chromosome_number
#increase in queen age
queen_age <- trait_dat_rer$mean_queen_age
#decrease in polygyny
anti_polygyny <- 3 - as.numeric(trait_dat_rer$Polygyny)

trait_list <- list(polyandry_dat,anti_repeat,
                   anti_totipotent,col_size, 
                   crossover_per_chromosome,queen_age,anti_polygyny)

charpaths <- rep(list(NULL),7)
RER_list <- rep(list(NULL),7)
res_list <- rep(list(NULL),7)
for (i in 7) {
names(trait_list[[i]]) <- trait_dat_rer$Species
charpaths[[i]] =char2Paths(trait_list[[i]], ort_trees)
RER_list[[i]] = getAllResiduals(ort_trees,
                                useSpecies=names(trait_list[[i]]), 
                                transform = "sqrt", 
                                weighted = T, scale = T)


res_list[[i]]=correlateWithContinuousPhenotype(RER_list[[i]], 
                                               charpaths[[i]], 
                                               min.sp = 10, 
                                               winsorizeRER = 3, 
                                               winsorizetrait = 3)
}
#permulation option of rerconverge
perms_list <- rep(list(NULL),7)
mt <- ort_trees$masterTree
mt=root.phylo(mt, outgroup="Vespula_vulgaris", resolve.root=T)
for (i in 7) {
perms_list[[i]]=RERconverge::getPermsContinuous(numperms =  1000, 
                                                      traitvec =  trait_list[[i]], 
                                                      RERmat =  RER_list[[i]],
                                                      trees = ort_trees,
                                                      mastertree =  mt,
                                                      calculateenrich = F,
                                                      )

}

#get the permulation p-vals
corpermpvals_list <- rep(list(NULL),7)
for (i in 7) {
corpermpvals_list[[i]]=RERconverge::permpvalcor(res_list[[i]], 
                                                      perms_list[[i]])
res_list[[i]]$permpval=corpermpvals_list[[i]][match(rownames(res_list[[i]]), 
                                                                names(corpermpvals_list[[i]]))]
res_list[[i]]$permpvaladj=p.adjust(res_list[[i]]$permpval, 
                                         method="BH")
}

plotRers(res_list[[5]],"OG0003699")
plotTreeHighlightBranches(ort_trees$masterTree,hlspecies = c(NULL))
plotTreeHighlightBranches(ort_trees$trees$OG0003699,
                          hlspecies = c("Pogonomyrmex_barbosus",
                                        "Cardiocondyla_obscurior",
                                        "Solenopsis_invicta",
                                        "Acromyrmex_echinatior"),
                          hlcols=c("blue"))
plotTreeHighlightBranches(ort_trees$trees$OG0006386,
                          hlspecies = c(NULL),
                          hlcols=c("blue"))
#write out the RERconverge results
output_names <- c("polyandry","anti_repeat",
                   "anti_totipotent","col_size",
                   "crossover_per_chromosome","queen_age","anti_polygyny")
output_names <- sapply(output_names,
                       function(x)paste("rer_results_",x,sep = ""),
                       USE.NAMES = F)


for (i in 7) {
  write.table(res_list[[i]],paste("rer_data/",
                                  output_names[i],
                                  sep = ""))
}


