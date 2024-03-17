################################################################################
### Empirical Validation
### Predictive ability between predicted and observed values
###
### This script will calculate the predictive ability (correlation) of predicted and observed parameters 
### from SoyNAM project/package. This will include:
### 1. Observed values obtained from ASRemlR (varG and rG)
### 2. Predicted values from PopVar using pheno data as BLUEs
### 3. Data visualization
### Author: Cleiton Wartha
### > session_info()
### version  R version 4.3.3 (2024-02-29 ucrt)
### os       Windows 10 x64 (build 19045)  system   x86_64, mingw32
### ui       RStudio language (EN)  tz       America/Chicago
### collate  English_United States.utf8  ctype    English_United States.utf8
### rstudio  2023.06.0+421 Mountain Hydrangea (desktop)
### package      * version    date (UTC) lib source
### boot         * 1.3-29     2024-02-19 [2] CRAN (R 4.3.3)
### cowplot      * 1.1.3      2024-01-22 [1] CRAN (R 4.3.3)
### devtools     * 2.4.5      2022-10-11 [1] CRAN (R 4.3.2)
### dplyr        * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
### forcats      * 1.0.0      2023-01-29 [1] CRAN (R 4.3.2)
### fs             1.6.3      2023-07-20 [1] CRAN (R 4.3.2)
### ggplot2      * 3.4.4      2023-10-12 [1] CRAN (R 4.3.2)
### ggpmisc      * 0.5.5      2023-11-15 [1] CRAN (R 4.3.3)
### ggridges     * 0.5.6      2024-01-23 [1] CRAN (R 4.3.2)
### gridExtra    * 2.3        2017-09-09 [1] CRAN (R 4.3.3)
### gtools         3.9.5      2023-11-20 [1] CRAN (R 4.3.3)
### Matrix         1.6-5      2024-01-11 [1] CRAN (R 4.3.3)
### MatrixModels   0.5-3      2023-11-06 [1] CRAN (R 4.3.3)
### moments      * 0.14.1     2022-05-02 [1] CRAN (R 4.3.1)
### RColorBrewer * 1.1-3      2022-04-03 [1] CRAN (R 4.3.1)
### Rcpp           1.0.11     2023-07-06 [1] CRAN (R 4.3.2)
### remotes        2.4.2.1    2023-07-18 [1] CRAN (R 4.3.2)
### rlang          1.1.2      2023-11-04 [1] CRAN (R 4.3.2)
### tidyr        * 1.3.0      2023-01-24 [1] CRAN (R 4.3.2)
### tidyverse    * 2.0.0      2023-02-22 [1] CRAN (R 4.3.2)
################################################################################
#Install packages and load libraries
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('moments')) install.packages('moments'); library(moments)
if (!require('boot')) install.packages('boot'); library(boot)
if (!require('cowplot')) install.packages('cowplot'); library(cowplot)
if (!require('ggpmisc')) install.packages('ggpmisc'); library(ggpmisc) 
if (!require('devtools')) install.packages('devtools'); library(devtools)
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library(RColorBrewer)
if (!require('ggridges')) install.packages('ggridges'); library(ggridges)
if (!require('grid')) install.packages('grid'); library(grid)
if (!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
options(scipen = 999) #prevents use of scientific notation

##Load object with outpuf files from scripts #1 and #2
load("data/obj3.RData")
# pheno: phenotypic BLUEs obtained from ASRemlR 4.2 #script1
# predVal: predicted values obtained from PopVar::pop_predict2 #script1
# varG.uv: Observed genetic variances for each family obtained from univariate models fitted with ASRemlR 4.2 #script2
# rG.bv: Observed genetic correlations for each family obtained from bivariate models fitted with ASRemlR 4.2 #script2

pheno <- pheno %>%  #convert phenotypic BLUEs from wide to long format for summarize
         pivot_longer(cols = height:yield, names_to = "trait")
DisTab <- pheno %>% group_by(trait, family) %>%
          summarise(Skew = skewness(value), Kurtosis = kurtosis(value)) #obtain estimates of Skewness and Kurtosis for each family-trait combination
#write.csv(DisTab, file = "Skewness_Kurtosis_Progeny_Phenotypic_Value_Distribution_Traits_by_Family.csv", row.names = F) #save file
# Relevant traits
traits <- c( "height","lodging","maturity","oil","protein","size","yield")
# Create a vector to replace trait names for data visualization
trait.replace <- setNames(c("Plant Height", "Lodging", "Days to Maturity", "Oil", "Protein", "Seed Size", "Seed Yield"), traits)
#Data.frame with parental lines and group information
ParID <- data.frame(parent = c('TN05-3027','4J105-3-4','5M20-2-5-2','CL0J095-4-6','CL0J173-6-8','HS6-3976','Prohio','LD00-3309','LD01-5907','LD02-4485','LD02-9050','Magellan',
                               'Maverick','S06-13640','NE3001','Skylla','U03-100612','LG03-2979','LG03-3191','LG04-4717','LG05-4292','LG05-4317','LG05-4464','LG05-4832','LG90-2550','LG92-1255',
                               'LG94-1128','LG94-1906','LG97-7012','LG98-1605','LG00-3372','LG04-6000','PI398881','PI427136','PI437169B','PI518751','PI561370','PI404188A','PI574486'),
                    family = c(2,3,4,5,6,8,9,10,11,12,13,14,15,17,18,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,41,42,48,50,54,64),
                    Group= c('EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','BX','BX','BX','BX','BX','BX',
                             'BX','BX','BX','BX','BX','BX','BX','BX','BX','PI','PI','PI','PI','PI','PI','PI'))
#Combine information for density plot
DisTab <- left_join(pheno, DisTab, by = join_by(family, trait)) %>%
  left_join(., ParID, by= "family") %>%
  mutate(Group = as.factor(Group)) %>%
  mutate(Trait = str_replace_all(trait, trait.replace)) %>%
  mutate(annotation = str_c("Skew = ", round(Skew, 4), "\n", "Kurtosis = ", round(Kurtosis, 4))) %>%
  mutate(Flag = case_when(Skew > 1 ~ "HighSkew",
                          Skew < -1 ~ "LowSkew",
                          Kurtosis > 5 ~ "HighKurt",
                          TRUE ~ NA))
#Create file to Flag departures from normality
Flagdat <- DisTab %>% select(family, trait, Flag) %>%
  distinct(., .keep_all = T)

## Plot per trait
for (t in 1:length(traits)){
       print(traits[t])
      DisTabSub <- DisTab %>% filter(trait == trait[t]) #subset data for single trait
      #Create density plot and add skewness and kurtosis values
      density.plot <- ggplot(data= DisTabSub, aes(x = value)) +
      geom_density(alpha=0.15, aes(fill= Group)) + #fill="#999999",
        scale_fill_manual(values = c('#999999','#E69F00','#56B4E9')) +
      facet_wrap(~ family, ncol = 4, scale = "free") +
      theme_classic() +
        labs(title=trait.replace[t], x = "Progeny Phenotypic Value", y ="Density") +
        theme(plot.title = element_text(face= "bold", hjust= 0.5), 
              axis.title.x= element_text(face= "bold"),
              axis.title.y= element_text(face= "bold"),
              legend.position = "none") +
      geom_text(data = distinct(DisTabSub, family, Skew, Kurtosis, annotation), aes(x = Inf, y = Inf, label = annotation), 
                 size = 3, vjust = "inward", hjust = "inward")
  ggsave(filename = paste0("Density_Plot_Progeny_Phenotypic_Value_",trait.replace[t],".png"), plot = density.plot, height = 15, width = 12, dpi = 600)  
}
#Summary table of number of RILs per family
RILSum <- pheno %>% group_by(family) %>%
  summarise(N = n_distinct(strain))

################################################################################
###Prepare files to estimate predictive ability, i.e., correlation between predicted and observed values
################################################################################
#Calculate mu per family using BLUEs 
tail.p= 0.1 #proportion selected 10% ksp= 1.75
Obsmu <- pheno %>% group_by(family, trait) %>%
  summarise(mu = mean(value),
            musp_low = mean(value[value <= quantile(value, probs= tail.p)]),
            musp_high = mean(value[value >= quantile(value, probs= 1 - tail.p)])) 
ObsmuspLow <- Obsmu %>% filter(trait == c("height", "lodging", "maturity")) %>% #select traits with desired lower values
            select(family, trait, musp_low) %>% #select only desired columns
            rename(musp = musp_low) #rename to be able to rbind
ObsmuspHigh <- Obsmu %>% filter(trait != c("height", "lodging", "maturity")) %>% #select traits with desired higher values
  select(family, trait, musp_high) %>% #select only desired columns
  rename(musp = musp_high) #rename to be able to rbind
Obsmusp <- rbind(ObsmuspLow, ObsmuspHigh) #combine musp files

Obsmu <- left_join(Obsmu, Obsmusp) %>% #join files
          select(-musp_low, -musp_high) #remove redundant columns

################################################################################
####Summary tables across all families for observed mu, varG and rG###############
#Obtain file for ObsmuvarG
ObsmuvarG <- varG.uv %>% select(family, trait, varG, std.error) %>%
  left_join(., Obsmu) %>% #combine obs values for Summary table
        pivot_longer(cols = varG:musp, names_to = "parameter")
#Summary table for mu and varG
ObsmuvarGSmry <- ObsmuvarG %>% group_by(trait, parameter) %>%
  summarize_at(vars(value), list(mean = mean, min = min, max = max)) %>% #funs(mean, min, max)
  arrange(parameter) #alphabetical order
#write.csv(ObsmuvarGSmry, file = "Descriptive_Summary_Observed_muvarG.csv", row.names = F) #save file
#Summary table for rG
ObsrGSmry <- rG.bv %>% select(family, trait1, trait2, rG, std.error) %>%
  group_by(trait1, trait2) %>%
  summarize_at(vars(rG), list(mean = mean, min = min, max = max))
#write.csv(ObsrGSmry, file = "Descriptive_Summary_Observed_rG.csv", row.names = F) #save file
################################################################################
### Predictions from deterministic equations PopVar              ###############
### Predicted values were obtained from MSI supercomputer        ###############
################################################################################
#Load file with results from SoyNAM package BLUEs
predVal <- predVal %>% select(parent1:cor_w_size) %>% #disregard columns with correlated progeny pred low and high superior progeny
          left_join(., ParID, by = c("parent2" = "parent")) #add family and group information
#Work on pred_musp
predmuspLow <-predVal %>% filter(trait %in% c("height", "lodging", "maturity")) %>% #select trait with desired lower values
  select(family, trait, pred_musp_low) %>% #select only desired columns
  rename(pred_musp = pred_musp_low) #rename to be able to rbind
predmuspHigh <-predVal %>% filter(trait %in% c("oil", "protein", "size", "yield")) %>% #select trait with desired lower values
  select(family, trait, pred_musp_high) %>% #select only desired columns
  rename(pred_musp = pred_musp_high) #rename to be able to rbind
predmusp <- rbind(predmuspLow, predmuspHigh) #combine musp files
#Bring results back into predVal file
predVal <- left_join(predVal, predmusp) %>% #recombine files
  select(-pred_musp_low, -pred_musp_high)
#More housekeeping
predValLong <- predVal %>% select(-parent1, -parent2) %>% #remove columns with parents 
  relocate(family, Group, .before = trait) %>% #bring Fam column to the front
  relocate(pred_musp, .after = pred_varG) %>% #bring pred_musp after pred_varG
  arrange(family) %>% # sort Fam column
  pivot_longer(cols = pred_mu:cor_w_size, names_to = "pred", values_to = "prediction") %>%
  filter(complete.cases(prediction)) #remove cases of cor(trait by same trait) 
 
#Subset the predicted mean, variance, and correlation
predmu <- predValLong %>% filter(grepl("^pred_mu$", pred)) %>% #filter rows with pred_mu - needs to not select pred_musp
    arrange(family, trait) #to have in same order as observed values by Fam and then trait in alphabetical order
predmusp <- predValLong %>% filter(grepl("pred_musp", pred)) %>% #filter rows with pred_mu
  arrange(family, trait) #to have in same order as observed values by family and then trait in alphabetical order
predvarG <- predValLong %>% filter(grepl("pred_varG", pred)) %>% #filter rows with pred_varG
    arrange(family, trait) #to have in same order as observed values by family and then trait in alphabetical order
predrG <- predValLong %>% filter(grepl("cor_w_", pred)) %>% #filter rows with predictions
    rename(trait2 = pred) %>% #rename column
    mutate(trait2 = gsub('cor_w_', '', trait2)) %>%#remove correlation portion from trait name
  rename(trait1 = trait) %>%
  arrange(family, trait1, trait2) %>% #arrange column order to match obsrG before removing duplicates
  filter(!duplicated(prediction))#filter duplicated rG values #WARNING: if two predicted values from different trait comb are =, they would be deleted
#please check file dimensions with rG.bv to have equal number of observations
  
#Combine Obs and Pred sets for predictive ability calculation
#mu
predObsmu <- left_join(Obsmu, predmu) %>%
  rename(., estimate= mu) %>% #rename column
  select(., -pred) #remove pred column
#musp
predObsmusp <- left_join(Obsmu, predmusp) %>%
  rename(., estimate= musp) %>% #rename column
  select(., -pred, -mu) #remove pred and mu column
#mu and musp to make joint figure
ObsmuLong <-  Obsmu %>% pivot_longer(cols = mu:musp, names_to = "parameter", values_to = "estimate")
predmuLong <- rbind(predmu, predmusp) %>% rename(parameter = pred) %>%
  mutate(parameter = gsub('pred_', '', parameter))
predObsmu_musp <- left_join(ObsmuLong, predmuLong, by = c("family", "trait", "parameter")) 
#varG
predObsvarG <- left_join(varG.uv, predvarG) %>%
  rename(., estimate= varG) %>% #rename column
  select(., -pred, -z.ratio, -bound, -X.ch) #remove columns
#write.csv(predObsvarG, "SoyNAM_Predicted_Observed_varG_39fam_7_traits.csv", row.names = F)
#rG
predObsrG <- left_join(rG.bv, predrG) %>%
  rename(., estimate= rG) %>% #rename column
    select(., -z.ratio, -bound, -X.ch)
#write.csv(predObsrG , "SoyNAM_Predicted_Observed_rG_39_fam_21_trait_combinations.csv", row.names = F)
### Measure prediction accuracy
# Calculate the correlation between predictions and observations
#Load bootstrapping function
source("src/bootstrap.R") #bootstrap by J Neyhart
boot_reps <- 10000
alpha <- 0.05
set.seed(2024)
#mu
predAccmu <- predObsmu %>% 
  group_by(trait) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup() %>%
  mutate(base_rounded = round(base, 2),
         base_signif = paste(base_rounded, annotation, sep= "")) %>%
  mutate(parameter = "mu") %>%
  relocate(parameter, .before = statistic)
  
muBias <- predObsmu %>%
        group_by(trait)%>%
  summarise(pred = mean(prediction),
            obs = mean(estimate))%>%
  mutate(mean_bias = ((pred - obs)/obs)*100) %>% #calculate the mean bias across all families
  select(trait, mean_bias)
predAccBiasmu <- left_join(predAccmu, muBias)
#write.csv(predAccBiasmu, "Prediction_Accuracy_Bias_mu_summary_table.csv", row.names = F)
#write.csv(predObsmu, "Obs_pred_mu_values_by_fam.csv", row.names = F)
#musp
predAccmusp <- predObsmusp %>% 
  group_by(trait) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup() %>%
  mutate(base_rounded = round(base, 2),
         base_signif = paste(base_rounded, annotation, sep= ""))%>%
  mutate(parameter = "musp") %>%
  relocate(parameter, .before = statistic)

muspBias <- predObsmusp %>%
  group_by(trait)%>%
  summarise(pred = mean(prediction),
            obs = mean(estimate))%>%
  mutate(mean_bias = ((pred - obs)/obs)*100) %>% #calculate the mean bias across all families
  select(trait, mean_bias)
predAccBiasmusp <- left_join(predAccmusp, muspBias)

#write.csv(predAccBiasmusp, "results/Prediction_Accuracy_Bias_musp_summary_table.csv", row.names = F)
# write.csv(predObsmu, "results/Obs_pred_mu_values_by_fam.csv", row.names = F)
predAccBiasmu_musp <- rbind(predAccBiasmu, predAccBiasmusp) #to use in figure

#varG
predAccvarG <- predObsvarG %>% 
    group_by(trait) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup() %>%
  mutate(base_rounded = round(base, 2),
         base_signif = paste(base_rounded, annotation, sep= ""))%>%
  mutate(parameter = "varG") %>%
  relocate(parameter, .before = statistic)

varGBiasValFam <- predObsvarG %>%
  mutate(Bias.calc = ((prediction - estimate)/estimate)*100) 

varGBias <- predObsvarG %>%
  group_by(trait) %>%
  summarise(pred = mean(prediction),
            obs = mean(estimate))%>%
  mutate(mean_bias = ((pred - obs)/obs)*100) %>% #calculate the mean bias across all families
  select(trait, mean_bias)
predAccBiasvarG <- left_join(predAccvarG, varGBias)

write.csv(predAccBiasvarG, "Prediction_Accuracy_Bias_VarG_summary_table.csv", row.names = F)
 #write.csv(predObsvarG, "Obs_pred_varG_values_by_fam.csv", row.names = F)

#Genetic correlations
### Measure prediction accuracy
# Calculate the correlation between predictions and observations
predAccrG <- predObsrG %>% 
  group_by(trait1, trait2) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", ""),
         trait_pair = str_c(trait1, " / ", trait2)) %>%
  ungroup() %>%
  distinct(base, .keep_all= TRUE) %>% #remove duplicated estimates of correlation for each pair of traits
  mutate(base_rounded = round(base, 2),
         base_signif = paste(base_rounded, annotation, sep= ""))%>%
  mutate(parameter = "rG") %>%
  relocate(parameter, .before = statistic)

rGBias <- predObsrG %>%
  group_by(trait1, trait2) %>%
  summarise(pred = mean(prediction),
            obs = mean(estimate))%>%
  mutate(mean_bias = ((pred - obs)/obs)*100) %>% #calculate the mean bias across all families
  select(trait1, trait2, mean_bias)
predAccBiasrG <- left_join(predAccrG, rGBias)

write.csv(predAccBiasrG, "Prediction_Accuracy_Bias_rG_summary_table.csv", row.names = F)
# write.csv(predObsrG, "Obs_pred_rG_values_by_fam.csv", row.names = F)

################################################################################
##### Data Visualization                                          ##############
################################################################################
#Obtain file for ObsmuvarG
ObsmuvarG <- varG.uv %>% select(family, trait, varG, std.error) %>%
  left_join(., Obsmu) %>%
  left_join(., ParID, by= "family") %>%
  mutate(Group = as.factor(Group)) %>%
  mutate(Trait = str_replace_all(trait, trait.replace))
ObsmuvarGCov <- ObsmuvarG %>%
  group_by(trait) %>%
  summarise(Covariance = cov(varG, mu, method = "pearson"))
#write.csv(ObsmuvarGCov, "results/Covariance_Observed_mean_genetic_variance_all_traits.csv", row.names = F)

#Supplementary figure
ObsmuvarG.cov.plot <- ggplot(ObsmuvarG, aes(x= mu, y = varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  #geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ Trait, ncol = 2, scales = "free") +
  theme_classic() + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")), coef.digits = 4, p.digits=2) + #add regression line values
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Observed" ~italic(mu)), y=expression("Observed" ~ italic(sigma[g]^2))) +
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
##Version 2 - regression line only for significant correlations
#Subset to remove traits that were not significant
ObsmuvarG.covSub <- ObsmuvarG %>% filter(trait %in% c("lodging", "maturity", "oil", "yield")) #significant traits 
ObsmuvarG.cov.plot <- ObsmuvarG.cov.plot + #new plot with regression line only for the significant cases
  geom_smooth(data= ObsmuvarG.covSub, se=F, method = "lm", color="darkblue")
ggsave(filename = "Covariance_Observed_mean_genetic_variance_all_traits.png", plot = ObsmuvarG.cov.plot, height = 12, width = 9, dpi = 1000)

#Figure 1. Plot the distributions of each trait and for each family
family.BLUE <- pheno %>% group_by(family, trait) %>%
  mutate(family.mean.est = mean(value)) %>% mutate(family= as.character(family)) %>% #convert family to character for bind_rows downstream
  ungroup()
# Combine 
family.BLUE.toplot <- bind_rows(
  mutate(family.BLUE, type = "By\nfamily"), 
  mutate(family.BLUE, family = "All", type = "All"))
# Relevant traits
traits <- c( "height","lodging","maturity","oil","protein","size","yield")
# Create a vector to replace trait names
trait.replace <- setNames(c("Plant Height", "Lodging", "Days to Maturity", "Oil", "Protein", "Seed Size", "Seed Yield"), traits)
# Create a vector of colors for each family
all.colors <- colorRampPalette(rev(brewer.pal(12, "Paired")))
# Plot per trait
g.family.density <- family.BLUE.toplot %>%
  mutate(trait = str_replace_all(trait, trait.replace)) %>%
  split(.$trait) %>%
  map(~{
    temp <- .
    # Order the family based on the family mean
    family.order <- distinct(temp, family, family.mean.est) %>% 
      filter(family != "All") %>% 
      arrange(family.mean.est) %>% 
      pull(family)
    # Convert family to factor
    temp$family <- factor(temp$family, levels = c(family.order, "All"))
    # Create a color scheme
    family.color <- setNames(c(all.colors(n = nlevels(temp$family) - 1), "grey"), levels(temp$family))
    # Plot
    temp %>%
      ggplot(aes(x = value, y = type, fill = family)) +
      geom_density_ridges(alpha = 0.2, size = 0.25) +
      facet_wrap(~ trait, ncol = 3, scale = "free") +
      scale_fill_manual( values = family.color, name = "Family", guide = "none") +
      theme_classic() +
      theme(axis.title = element_blank(), strip.text = element_text(size = 17), #increase text size of facet_wrap
            axis.text = element_text(size = 13)) #increase text size of axis
  })
# Cowplot
g.density.plot <- plot_grid(plotlist = g.family.density, ncol = 2, nrow = 4, align = "hv")  #create a grid with all charts for each trait
ggsave(filename = "Figure1_Validation_Fam_Trait_Density_Distribution.png", plot = g.density.plot,
       height = 15, width = 12, dpi = 1000)
#####################################
#Figure 2.1 Observed by predicted mean
predObsmu <- predObsmu %>%
  mutate(parameter = "mu") %>%
  mutate_at(c('trait', 'parameter'), as.factor) %>%
  left_join(., predAccBiasmu, by = c("trait", "parameter")) %>%
  mutate(annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))%>%
  mutate(trait = str_replace_all(trait, trait.replace)) #replace trait values by actual name

#Make plot
Mu.scatter.plot <- ggplot(predObsmu, aes(x= prediction, y = estimate)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5) +
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Predicted" ~italic(hat(mu))), y=expression("Observed" ~ italic(mu))) +
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=15), #change legend title font size
         legend.text = element_text(size=11),
         axis.title = element_text(size = 19),
         legend.position="bottom") #legend position #change legend text font size
ggsave(filename = "Figure2_1_ScatterPlot_MEAN_Obs_Predicted.png", plot = Mu.scatter.plot, height = 12, width = 9, dpi = 1000)

#Figure 2.2 Observed vs predicted superior progeny mean
predObsmusp <- predObsmusp %>%
  mutate(parameter = "musp")%>%
    mutate_at(c('trait', 'parameter'), as.factor) %>%
  left_join(., predAccBiasmusp, by = c("trait", "parameter")) %>%
  mutate(annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))%>%
  mutate(trait = str_replace_all(trait, trait.replace)) #replace trait values by actual name

Musp.scatter.plot <- ggplot(predObsmusp, aes(x= prediction, y = estimate)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(data = distinct(predObsmusp, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Predicted" ~italic(hat(mu)[sp])), y=expression("Observed" ~ italic(mu)[SP])) +
  theme( axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 14), #change font size of facet grid labels
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=11),
        axis.title = element_text(size = 19),
        legend.position="bottom") #legend position #change legend text font size
ggsave(filename = "Figure2_2_ScatterPlot_MUSP_Obs_Predicted.png", plot = Musp.scatter.plot, height = 12, width = 9, dpi = 1000)

#Figure 2.3 Observed by predicted mean and musp
# Create a vector to replace the parameters in graphs
 parameter<- c("mu",  "mu_sp")
 param_replace <- c("mu" = "mu", "mu_sp" = "musp")
predAccBiasmu_musp <-  predAccBiasmu_musp %>%
  mutate(trait = str_replace_all(trait, trait.replace),
         parameter = str_replace_all(parameter, param_replace))
######################
## Create figure with mu and musp combined
df1<- predObsmu_musp %>%
  mutate(trait = str_replace_all(trait, trait.replace)) %>% #replace trait values by actual name
  mutate_at(c('trait', 'parameter'), as.factor) %>%
      left_join(., predAccBiasmu_musp, by = c("trait", "parameter")) %>%
      mutate(parameter = str_replace_all(parameter, param_replace),
             parameter = factor(parameter, levels = param_replace),
             annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))
    # # Split by trait and parameter
     plot_list <- df1 %>%
       split(list(.$trait, .$parameter)) %>%
       map(~{
         df2 <- .
         ggplot(df2, aes(x = prediction, y = estimate)) +
           geom_smooth(method = "lm", se = FALSE, color="darkblue") + 
           geom_point(aes(colour= Group), size= 2, shape = "circle") + 
           scale_colour_manual(values= c('#999999','#E69F00','#56B4E9')) +
           geom_text(data = distinct(df2, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
                     parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5) + 
           ylab("Observed") +
           xlab("Predicted") + 
           facet_grid(trait ~ parameter, scales = "free", labeller = labeller(parameter = label_parsed), switch = "y") + 
           scale_y_continuous(breaks = scales::pretty_breaks(), labels = function(x) str_pad(x, width = 2, pad = "0")) + 
           theme_classic(base_size = 6) +
           theme(strip.placement = "outside", axis.title = element_blank(), strip.text.y = element_text(size = 11),
                 legend.position = "none", strip.text.x = element_text(size=15),
                 axis.text = element_text(size = 8))
       })
    # ## Re-order
     plot_list <- plot_list[c(1,8,2,9,3,10,4,11, 5,12, 6,13,7,14)]
     ## Edit parts of the grid
    ## Remove strips
     plot_list[c(3:14)] <- plot_list[c(3:14)] %>%
       map(~ . + theme(strip.text.x = element_blank(), strip.background.x = element_blank()))
     plot_list[c(2,4,6,8,10,12,14)] <- plot_list[c(2,4,6,8,10,12,14)] %>%
       map(~ . + theme(strip.text.y = element_blank(), strip.background.y = element_blank()))
    # # Create the grid by rows
     top <- plot_grid(plotlist = plot_list[1:2], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     second <- plot_grid(plotlist = plot_list[3:4], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     third <- plot_grid(plotlist = plot_list[5:6], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     fourth <- plot_grid(plotlist = plot_list[7:8], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     fifth <- plot_grid(plotlist = plot_list[9:10], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     sixth <- plot_grid(plotlist = plot_list[11:12], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     seventh <- plot_grid(plotlist = plot_list[13:14], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
    
    p_grid <- plot_grid(top, second, third, fourth, fifth, sixth, seventh, ncol = 1, rel_heights = c(1, 0.9, 0.9))
     # ## Add axis
     y_axis <- grid::textGrob(label = "Observed", gp = grid::gpar(fontsize = 20), rot = 90)
     x_axis <- grid::textGrob(label = "Predicted", gp = grid::gpar(fontsize = 20))
    # 
    # # Plot again
    mu.musp.arrange <- grid.arrange(arrangeGrob(p_grid, left = y_axis, bottom = x_axis))
    ggsave(filename = "Figure_2_3_ScatterPlot_MEAN_MUSP_Obs_Predicted.png", plot = mu.musp.arrange , height = 13, width = 9, dpi = 1000)
    
####################################
# Figure 3. Observed by predicted varG
predObsvarG <- predObsvarG %>%
      left_join(., Flagdat) %>%
 mutate(parameter = "varG") %>%
  mutate_at(c('trait', 'parameter'), as.factor) %>%
  left_join(., predAccBiasvarG, by = c("trait", "parameter")) %>%
  mutate(trait = str_replace_all(trait, trait.replace),
         annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))
##Needs to create individual plots due to drastic changes in variance scales impacting the x and y-axis
p.height <- predObsvarG %>% filter(trait == "Plant Height") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  theme_classic() +
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+ 
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=1.0, position=position_dodge(0.05))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Plant Height") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.57*")

p.lodging <- predObsvarG %>% filter(trait == "Lodging") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=.001, position=position_dodge(0.05))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Lodging") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.26 ")

p.maturity <- predObsvarG %>% filter(trait == "Days to Maturity") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=.5, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Days to Maturity") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.88*")

p.oil <- predObsvarG %>% filter(trait == "Oil") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=.003, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Oil") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.34 ")

p.protein <- predObsvarG %>% filter(trait == "Protein") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=.005, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Protein") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.17 ")

p.size <- predObsvarG %>% filter(trait == "Seed Size") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=.02, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Seed Size") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.64*")

p.yield <- predObsvarG %>% filter(trait == "Seed Yield") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=400, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Seed Yield") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.45*")
# Cowplot
varG.scatter.plot <- plot_grid(p.maturity, p.lodging, p.oil, p.height,
                              p.protein, p.size, p.yield, align='vh', ncol = 2, nrow = 4)
#create common x and y labels
y.grob <- textGrob(expression("Observed" ~ italic(sigma[G]^2)), gp=gpar(col="black", fontsize=22), rot=90)
x.grob <- textGrob(expression("Predicted" ~ italic(hat(sigma)[g]^2)), gp=gpar( col="black", fontsize=22))

#add to plot
varG.scatter.arrange <- grid.arrange(arrangeGrob(varG.scatter.plot, left = y.grob, bottom = x.grob))
ggsave(filename = "Figure_4_ScatterPlot_VAR_G_Obs_Predicted.png", plot = varG.scatter.arrange , height = 15, width = 12, dpi = 1000)

##########################
##Figure 4. Observed vs Predicted Correlation
#Create vector with trait names to replace the colnames
#Housekeeping to replace trait names
predObsrG <- predObsrG %>%
  mutate(trait_comb = paste(trait1, trait2, sep = "_"))%>%
  mutate(parameter = "rG")%>%
  mutate_at(c('trait1', "trait2", 'parameter'), as.factor) %>%
  left_join(., predAccBiasrG, by = c("trait1", "trait2", "parameter")) %>%
  mutate(annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))

trait_comb <- unique(predObsrG$trait_comb) #extract chr vector with 21 trait combinations
trait.comb.replace <- setNames(c("Plant Height / Lodging","Plant Height / Days to Maturity", "Plant Height / Oil", "Plant Height / Protein", "Plant Height / Seed Size", "Plant Height / Seed Yield",
                                 "Lodging / Days to Maturity", "Lodging / Oil", "Lodging / Protein", "Lodging / Seed Size", "Lodging / Seed Yield", "Days to Maturity / Oil",
                                 "Days to Maturity / Protein", "Days to Maturity / Seed Size", "Days to Maturity / Seed Yield", "Oil / Protein", "Oil / Seed Size",
                                 "Oil / Seed Yield", "Protein / Seed Size", "Protein / Seed Yield", "Seed size / Seed Yield"), trait_comb)
predObsrG <- predObsrG %>%
  mutate(trait_comb = str_replace_all(trait_comb, trait.comb.replace))
#Create scatterplot
Correl.scatter.plot <- ggplot(predObsrG, aes(x= prediction, y = estimate)) +
  geom_point(aes(colour= Group), size= 2.2, shape = "circle") +
  geom_errorbar(aes(ymin= estimate-std.error, ymax= estimate+std.error),
                width=.01, position=position_dodge(0.04))+
    facet_wrap(~ trait_comb, ncol = 3, scales = "free") +
  theme_classic() + 
  geom_text(data = distinct(predObsrG, trait_comb, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5) +
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) + #define colors to each Group
  labs(x= expression("Predicted"  ~italic(hat(r)[g])), y=expression("Observed"  ~ italic(r[G]))) +
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 11.5),
        strip.text.x = element_text(size = 15), #change font size of facet grid labels
        legend.position="bottom", #legend position
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=18),
        legend.key.height=unit(0, "cm"),      
        plot.margin = unit(c(1,0.5,0,0.5), "lines")) 

##Version 2 - regression line only for significant correlations
#Subset to remove traits that were not significant
predObsrGSub <- predObsrG %>% filter(!trait_comb %in% c("Days to Maturity / Protein", "Lodging / Protein", "Protein / Seed Yield")) #remove nonsignificant traits 
Correl.scatter.plot <- Correl.scatter.plot + #new plot with regression line only for the significant cases
  geom_smooth(data= predObsrGSub, se=F, method = "lm", color="darkblue")
ggsave(filename = "Figure5_ScatterPlot_CORREL_Obs_Predicted.png", plot = Correl.scatter.plot, height = 15, width = 12, dpi = 1000)

#######################################################################
####Relationship between predicted mean and variance  #################
#Make plot
predVal <- predVal %>%
  mutate(trait = str_replace_all(trait, trait.replace)) #replace trait values by actual name
options(scipen = 999)
MuvarG.scatter.plot <- ggplot(predVal, aes(x= pred_mu, y = pred_varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep = "*\", \"*")), coef.digits = 4, p.digits = 2) + #add regression line values
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Predicted" ~italic(hat(mu))), y=expression("Predicted" ~ italic(hat(sigma)[g]^2))) +
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size

##Version 2 - regression line only for significant correlations
#Subset to remove traits that were not significant
predValSub <- predVal %>% filter(!trait %in% c("Oil")) #remove nonsignificant traits 
MuvarG.scatter.plot <- MuvarG.scatter.plot + #new plot with regression line only for the significant cases
  geom_smooth(data= predValSub, se=F, method = "lm", color="darkblue")
ggsave(filename = "Figure3_ScatterPlot_MEAN_varG_Predicted.png", plot = MuvarG.scatter.plot, height = 12, width = 9, dpi = 1000)
