set.seed(1000)
library(caper)
library(dplyr)
library(car)
library(stringr)
library(phytools)

hym_tree <- read.newick("species_tree_rooted_sans_aas.txt")
#remove _CDS_aas suffix from the tip labels of the tree
hym_tree$node.label <- NULL
trait_dat <- read.csv("trait_dat.csv")
trait_dat
trait_dat$Species <- trait_dat$Species <- sub(" ","_",trait_dat$Species)
#now we will calculate the number of crossovers per chromosome
#centimorgans per genome divided by one hundred is the number of crossovers per genome
cross_per_genome <- (trait_dat$Recombination.Rate*trait_dat$genome_size)/100
#divided by the number of chromosomes is the number of crossovers per chromosome
crossover_per_chromosome <- (cross_per_genome)/trait_dat$chromosome_number
#remove recombination rate in cM/Mb and add in new response variable
rec_data_set <- cbind(trait_dat[,-6],crossover_per_chromosome)

#scaling the predictors to have mean 0 sd 1
scaled_preds <- scale(rec_data_set[,2:7],center = T, scale = T)
scaled_rec_dat <- data.frame(rec_data_set$Species, 
                             scaled_preds,
                             rec_data_set$crossover_per_chromosome)
colnames(scaled_rec_dat)[8] <- "crossover_per_chromosome"
colnames(scaled_rec_dat)[1] <- "Species"
#all 13 species have the spellings the same in newick tree and trait data set:
#intersect(hym_tree$tip.label,scaled_rec_dat$Species)

#the phylogenetically aware data that has been scaled
hym_dat_scaled <- comparative.data(phy = hym_tree, 
                                   data = scaled_rec_dat, 
                                   names.col = Species, 
                                   vcv = TRUE, na.omit = FALSE, 
                                   warn.dropped = TRUE)

#the names of the predictive variables of interest
preds_rec <- colnames(rec_data_set)[2:7]
pgls_preds <- c() 
for (i in 1:6) {
  pred_combos <- combn(preds_rec,i)
#the single variable models don't need the "+" sign
if(i > 1){
pred_formulas <- apply(pred_combos,MARGIN=2,
      FUN=function(x) paste(x,collapse=" + "))
}else{
  pred_formulas <- pred_combos
}
  pgls_preds <- c(pgls_preds, pred_formulas)
}
#the equations to be put into the caper pgls regressions
pgls_eqns <- paste("crossover_per_chromosome~ ",
                   pgls_preds,sep = "")
#pgls_eqns

#get AIC vals for each pgls model
aic_vals <- c()
index_val <- 0
for (i in pgls_eqns) {
  index_val <- 1 + index_val
  skip_to_next <- FALSE
  pgls_model <- tryCatch({pgls(as.formula(i), 
                                 data = hym_dat_scaled, lambda = "ML")},
                           error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next==TRUE){
    print(index_val)
    next}
  aic_iter <- pgls_model$aic
  aic_vals <- c(aic_vals,aic_iter)
  #print(pgls_eqns)
}
#one of the 63 (eqn 22) models couldn't converge
#"crossover_per_chromosome~ Polyandry + Polygyny + Colony.Size"
model_summaries <- data.frame(pgls_preds[-22],pgls_eqns[-22],aic_vals)
colnames(model_summaries)[1:2] <- c("pgls_preds","pgls_eqns")
#sort by aic
model_summaries_order <- model_summaries[order(model_summaries$aic_vals,
                                               decreasing = FALSE),]

#leave one out cross validation of r-sq values
rsq_top_models <- c()
#for each of the top 12 formulas run loocv
for (i in 1:length(model_summaries_order$pgls_eqns)) {
  #leave one species out for each GWRR data point
  predicted_rates <- c()
  for (j in 1:13) {
    #dataset missing datapoint
    one_missing <- scaled_rec_dat[-j,]
    #make data set missing one data point
    cv_hym_dat <- comparative.data(phy = hym_tree, 
                                   data = one_missing, 
                                   names.col = Species, 
                                   vcv = TRUE, na.omit = FALSE, 
                                   warn.dropped = TRUE)
    #some iterations in leave-one-out have trouble convergenging with ML approach...
    skip_to_next <- FALSE
    cv_pgls_best <- tryCatch({pgls(as.formula(model_summaries_order$pgls_eqns[i]), 
                                   data = cv_hym_dat, lambda = "ML")},
                             error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next==TRUE){next}
    #to predict a rate, of an observation, we need to know which/how many predictors we have
    pred_traits <- str_split(model_summaries_order$pgls_preds[i],
                             " \\+ " )
    #predicted intercept plus each number of predictors coeffs*observed value for indv 
    pred_int <- cv_pgls_best[[1]]$coef[1]
    pred_effect <- c()
    for (k in 1:length(pred_traits[[1]])) {
      #k+1 because first coefficient of model is intercept
      trait_coeff <- cv_pgls_best[[1]]$coef[k+1]
      trait_val <- eval(parse(text=(paste("scaled_rec_dat[j,]$",
                                          pred_traits[[1]][k],sep = ""))))
      #regression coeffecient for predictive trait times the observed value of the left out individual for that trait
      pred_effect <-  c(pred_effect,trait_coeff*trait_val)
      #print(pred_effect)
    }
    #predicted recombination rate for individual left out of regression model construction
    pred_rate_indv <- pred_int+sum(pred_effect)
    #predicted rates of each left out individual
    predicted_rates <- c(predicted_rates,pred_rate_indv)
  }
  #calculate r-sq value
  y_i <- rec_data_set$crossover_per_chromosome
  rss <- sum((y_i-predicted_rates)^2)
  y_bar <- mean(rec_data_set$crossover_per_chromosome)
  sst <- sum((y_i-y_bar)^2)
  r_sq <- 1 - (rss/sst)
  rsq_top_models <- c(rsq_top_models,r_sq)
}
rsq_top_models
rsq_top_models_sort <- order(rsq_top_models,decreasing = T)

top_model <- pgls(crossover_per_chromosome~Polyandry+Polygyny+Colony.Size+Totipotency, 
                  data = hym_dat_scaled, lambda = "ML")
summary(top_model)
#top 10 models are within 2 aic of the lowest AIC model
model_summaries_order$pgls_preds[1:10]
#top 10 rsq models
top_10_rsq <- rsq_top_models_sort[1:10]
model_summaries_order$pgls_preds[top_10_rsq]
#how many of the top 10 overlap for both?                                 
intersect(model_summaries_order$pgls_preds[1:10],
          model_summaries_order$pgls_preds[top_10_rsq])

library(ggplot2)
library(dplyr)
# heatmap for predictor presence across 10 strongest models by AIC and loocv-rsq
data <- data.frame(
  Trait = c("Polyandry", "Polygyny", "Totipotency", "Colony size", "Repeat content", "Queen age"),
  Lowest_AIC = c(100, 60, 100, 40, 50, 30),
  Highest_LOOCV_rsq = c(100, 50, 90, 70, 30, 50)
)



# Convert data to long format
data_long <- tidyr::pivot_longer(data, cols = c(Lowest_AIC, Highest_LOOCV_rsq),
                                 names_to = "Metric", values_to = "Value")

# Calculate order of Traits by 'Highest_LOOCV_rsq'
trait_order <- data %>%
  arrange(Highest_LOOCV_rsq) %>%
  pull(Trait)

# Reorder Traits in data_long
data_long$Trait <- factor(data_long$Trait, levels = trait_order)


# Order traits by ascending LOOCV and ascending AIC
trait_order <- data %>%
  arrange(Highest_LOOCV_rsq, Lowest_AIC) %>%
  pull(Trait)

# Reorder Traits in data_long
data_long$Trait <- factor(data_long$Trait, levels = trait_order)

custom_labels <- c(
  "Lowest_AIC" = "Low AIC Models",
  "Highest_LOOCV_rsq" = "High LOOCV r-sq Models"
)

# Create the heatmap
heatmap_plot <- ggplot(data_long, aes(x = Metric, y = Trait, fill = Value)) +
  geom_tile(color = "white") +  # Add gridlines between tiles
  scale_fill_gradient(low = "white", high = "blue") +  # Customize color gradient
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 0, hjust = .5),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +  # Rotate x-axis labels
  labs(fill = "% of models \n containing predictor", title = "Trait Prevelance Across Strongest Models",
       x = "Performance Criterea", y = "Predictor of GWRR")+
  scale_fill_gradient(low = "#E0F7FA", high = "#0077BE",  # Adjust low to light blue
                      limits = c(30, 100),
                      breaks = c(30,65, 100))+
  scale_x_discrete(labels = custom_labels)+
  geom_text(aes(label=Value))

# Display the heatmap
print(heatmap_plot)









