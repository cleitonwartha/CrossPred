## Empirical Validation
## Predictive ability between predicted and observed values
##
## This script will calculate the predictive ability (correlation) of predicted and observed parameters 
## from SoyNAM project/package. This will include:
## 1. Observed values obtained from ASRemlR (vG and rG)
## 2. Predicted values from PopVar using pheno data as BLUEs
## 3. Data visualization
## Author: Cleiton Wartha
## Last modified: February, 2024
##
################################################################################
#Install packages and load libraries
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('moments')) install.packages('moments'); library(moments)
if (!require('boot')) install.packages('boot'); library(boot)
if (!require('cowplot')) install.packages('cowplot'); library(cowplot)
if (!require('ggpmisc')) install.packages('ggpmisc'); library(ggpmisc) 
if (!require('devtools')) install.packages('devtools'); library(devtools)
devtools::install_github("karthik/wesanderson")
options(scipen = 999) #prevents use of scientific notation

##Load object with outpuf files from script #1 and #2
load("data/obj3.RData")
# pheno: phenotypic BLUEs obtained from ASRemlR 4.2 #script1
# predVal: predicted values obtained from PopVar::pop_predict2 #script1
# varG.uv: Observed genetic variances for each family obtained from univariate models fitted with ASRemlR 4.2 #script2
# rG.bv: Observed genetic correlations for each family obtained from bivariate models fitted with ASRemlR 4.2 #script2

pheno <- pheno %>%  #convert phenotypic BLUEs from wide to long format for summarize
         pivot_longer(cols = height:yield, names_to = "trait")
DisTab <- pheno %>% group_by(trait, family) %>%
          summarise(Skew = skewness(value), Kurtosis = kurtosis(value))
write.csv(DisTab, file = "results/Skewness_Kurtosis_Progeny_Phenotypic_Value_Distribution_Traits_by_Family.csv", row.names = F) #save file
# Relevant traits
traits <- c( "height","lodging","maturity","oil","protein","size","yield")
# Create a vector to replace trait names
trait.replace <- setNames(c("Plant Height", "Lodging", "Days to Maturity", "Oil", "Protein", "Seed Size", "Seed Yield"), traits)
#Data.frame with parental lines and group information
ParID <- data.frame(parent = c('TN05-3027','4J105-3-4','5M20-2-5-2','CL0J095-4-6','CL0J173-6-8','HS6-3976','Prohio','LD00-3309','LD01-5907','LD02-4485','LD02-9050','Magellan',
                               'Maverick','S06-13640','NE3001','Skylla','U03-100612','LG03-2979','LG03-3191','LG04-4717','LG05-4292','LG05-4317','LG05-4464','LG05-4832','LG90-2550','LG92-1255',
                               'LG94-1128','LG94-1906','LG97-7012','LG98-1605','LG00-3372','LG04-6000','PI398881','PI427136','PI437169B','PI518751','PI561370','PI404188A','PI574486'),
                    family = c(2,3,4,5,6,8,9,10,11,12,13,14,15,17,18,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,41,42,48,50,54,64),
                    group= c('EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','BX','BX','BX','BX','BX','BX',
                             'BX','BX','BX','BX','BX','BX','BX','BX','BX','PI','PI','PI','PI','PI','PI','PI'))
#Combine information for density plot
DisTab <- left_join(pheno, DisTab, by = join_by(family, trait)) %>%
  left_join(., ParID, by= "family") %>%
  mutate(group = as.factor(group)) %>%
  mutate(Trait = str_replace_all(trait, trait.replace)) %>%
  mutate(annotation = str_c("Skew = ", round(Skew, 4), "\n", "Kurtosis = ", round(Kurtosis, 4)))

## Plot per trait
for (t in 1:length(traits)){
       print(traits[t])
      DisTabSub <- DisTab %>% filter(trait == trait[t]) #subset data for single trait
      #Create density plot and add skewness and kurtosis values
      density.plot <- ggplot(data= DisTabSub, aes(x = value)) +
      geom_density(alpha=0.15, aes(fill= group)) + #fill="#999999",
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
  ggsave(filename = paste0("results/Density_Plot_Progeny_Phenotypic_Value_",trait.replace[t],"_fill_by_group.png"), plot = density.plot, height = 15, width = 12, dpi = 600)  
}
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
####Summary tables across all families for observed mu, vG and rG###############
#Obtain file for ObsmuVg
ObsmuvarG <- varG.uv %>% select(family, trait, varG, std.error) %>%
  left_join(., Obsmu) %>% #combine obs values for Summary table
        pivot_longer(cols = varG:musp, names_to = "parameter")
#Summary table for mu and Vg
ObsmuvarGSmry <- ObsmuvarG %>% group_by(trait, parameter) %>%
  summarize_at(vars(value), list(mean = mean, min = min, max = max)) %>% #funs(mean, min, max)
  arrange(parameter) #alphabetical order
#Summary table for rG
ObsrGSmry <- rG.bv %>% select(family, trait1, trait2, rG, std.error) %>%
  group_by(trait1, trait2) %>%
  summarize_at(vars(rG), list(mean = mean, min = min, max = max))

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
  relocate(family, group, .before = trait) %>% #bring Fam column to the front
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
#rG
predObsrG <- left_join(rG.bv, predrG) %>%
  rename(., estimate= rG) %>% #rename column
    select(., -z.ratio, -bound, -X.ch)
### Measure prediction accuracy
# Calculate the correlation between predictions and observations
#Load bootstrapping function
source("src/bootstrap.R") #bootstrap by J Neyhart
boot_reps <- 50000
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
#write.csv(predAccBiasmu, "results/Prediction_Accuracy_Bias_mu_summary_table.csv", row.names = F)
# write.csv(predObsmu, "results/Obs_pred_mu_values_by_fam.csv", row.names = F)
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

write.csv(predAccBiasmusp, "results/Prediction_Accuracy_Bias_musp_summary_table.csv", row.names = F)
# write.csv(predObsmu, "results/Obs_pred_mu_values_by_fam.csv", row.names = F)
predAccBiasmu_musp <- rbind(predAccBiasmu, predAccBiasmusp) #to use in figure

#Vg
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

varGBias <- predObsvarG %>%
  group_by(trait)%>%
  summarise(pred = mean(prediction),
            obs = mean(estimate))%>%
  mutate(mean_bias = ((pred - obs)/obs)*100) %>% #calculate the mean bias across all families
  select(trait, mean_bias)
predAccBiasvarG <- left_join(predAccvarG, varGBias)

write.csv(predAccBiasvarG, "results/Prediction_Accuracy_Bias_VarG_summary_table.csv", row.names = F)
 #write.csv(predObsvarG, "results/Obs_pred_Vg_values_by_fam.csv", row.names = F)

#for correlations
### Measure prediction accuracy
# Calculate the correlation between predictions and observations
# Boot reps
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
  group_by(trait1, trait2)%>%
  summarise(pred = mean(prediction),
            obs = mean(estimate))%>%
  mutate(mean_bias = ((pred - obs)/obs)*100) %>% #calculate the mean bias across all families
  select(trait1, trait2, mean_bias)
predAccBiasrG <- left_join(predAccrG, rGBias)

write.csv(predAccBiasrG, "results/Prediction_Accuracy_Bias_rG_summary_table.csv", row.names = F)
# write.csv(predObsrG, "results/Obs_pred_rG_values_by_fam.csv", row.names = F)


################################################################################
##### Data Visualization                                          ##############
################################################################################
#Obtain file for ObsmuVg
ObsmuvarG <- varG.uv %>% select(family, trait, varG, std.error) %>%
  left_join(., Obsmu) %>%
  left_join(., ParID, by= "family") %>%
  mutate(Group = as.factor(group)) %>%
  mutate(Trait = str_replace_all(trait, trait.replace))
ObsmuvarGCov <- ObsmuvarG %>%
  group_by(trait) %>%
  summarise(Covariance = cov(varG, mu, method = "pearson"))
#write.csv(ObsmuvarGCov, "results/Covariance_Observed_mean_genetic_variance_all_traits.csv", row.names = F)

ObsmuvarG.cov.plot <- ggplot(ObsmuvarG, aes(x= mu, y = varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
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
ggsave(filename = "results/Covariance_Observed_mean_genetic_variance_all_traits.png", plot = ObsmuvarG.cov.plot, height = 12, width = 9, dpi = 1000)
