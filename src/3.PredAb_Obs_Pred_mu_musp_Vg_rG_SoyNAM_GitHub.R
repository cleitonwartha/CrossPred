## Empirical Validation
## Predictive ability between predicted and observed values
##
## This script will calculate the predictive ability (correlation) of predicted and observed parameters 
## from SoyNAM project/package. This will include:
## 1. Observed values obtained from ASRemlR (vG and rG)
## 2. Predicted values from PopVar using pheno data as BLUEs
## 3. Data visualization
## Author: Cleiton Wartha
## Last modified: July, 2022
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




#Figure 1. Plot the distributions of each trait and for each family
family.BLUE <- pheno %>% group_by(family, trait) %>%
  mutate(family.mean.est = mean(value)) %>% mutate(family= as.character(family)) %>% #convert family to character for bind_rows downstream
  ungroup()

## Calculate the mean, min, and max for all traits for the VP
family.BLUE %>% 
  group_by(trait) %>% 
  summarise(mean = mean(value), min = min(value), max = max(value))

# # A tibble: 8 × 4
# trait       mean      min     max
# <chr>      <dbl>    <dbl>   <dbl>
#   1 fiber       4.83    4.45     5.11
# 2 height     85.5    46.3    126.  
# 3 lodging     1.17    0.317    2.90
# 4 maturity  120.    103.     131.  
# 5 oil        20.5    15.6     22.4 
# 6 protein    32.6    29.6     40.9 
# 7 size       14.1     9.67    20.5 
# 8 yield    2828.   1326.    3818.

# Combine
family.BLUE.toplot <- bind_rows(
  mutate(family.BLUE, type = "By\nfamily"), 
  mutate(family.BLUE, family = "All", type = "All"))

# Relevant traits
traits <- c( "height","lodging","maturity","oil","protein","size","yield")
# Create a vector to replance trait names
trait.replace <- setNames(c("Plant Height", "Lodging", "Days to Maturity", "Oil", "Protein", "Seed Size", "Seed Yield"), traits)
# Create a vector of colors for each family
all.colors <- colorRampPalette(rev(brewer.pal(12, "Paired")))
heat.colors <-  wesanderson::wes_palette("Zissou1")
darj.colors <- wesanderson::wes_palette("Darjeeling1")
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
g.density.plot <- plot_grid(plotlist = g.family.density, ncol = 4, nrow = 2, align = "hv")  #create a grid with all charts for each trait
ggsave(filename = "results/Figure1/Validation_Fam_Trait_Density_Distribution_Paired_2Col_12_08_Paired_alphabetical.png", plot = g.density.plot,
       height = 7, width = 14, dpi = 1000)
#size values for figure in ggsave above, also change ncol arg in plot_grid
#4 height = 8, width = 24,
#2 height = 15, width = 12,
#1 height = 20, width = 6,

#####################################
#Figure 2. Observed by predicted mean
predObsmu <- predObsmu %>%
  mutate(parameter = "mu")%>%
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
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  #  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*"))) + #add regression line values
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Predicted" ~italic(hat(mu))), y=expression("Observed" ~ italic(mu))) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=15), #change legend title font size
         legend.text = element_text(size=11),
         axis.title = element_text(size = 19),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure2/ScatterPlot_MEAN_Obs_Predicted_1Col_color_point_07_26_printed_r.png", plot =Mu.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure2/ScatterPlot_MEAN_Obs_Predicted_2Col_color_point_11_28_printed_r.png", plot = Mu.scatter.plot, height = 12, width = 9, dpi = 1000)

#Figure 2.1 Observed vs predicted superior progeny mean
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
  #labs(x= expression(~italic(mu)[SP]), y = "")+
  theme( axis.text = element_text(size = 11),
        strip.text.x = element_text(size = 14), #change font size of facet grid labels
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=11),
        axis.title = element_text(size = 19),
        legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure2/ScatterPlot_MUSP_Obs_Predicted_1Col_color_point_07_26_printed_r.png", plot =Musp.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure2/ScatterPlot_MUSP_Obs_Predicted_2Col_color_point_11_28_printed_r.png", plot = Musp.scatter.plot, height = 12, width = 9, dpi = 1000)

#Figure 2.2 Observed by predicted mean and musp
# Create a vector to replace the parameters in graphs
 parameter<- c("mu",  "mu_sp")
# param.replace <- setNames(c("mu","mu[SP]"),parameter)
 param_replace <- c("mu" = "mu", "mu_sp" = "musp")
# "varG", "rG"
#"sigma[G]^2", "r[G]"

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
           ylab("Observation") +
           xlab("Prediction") + 
           facet_grid(trait ~ parameter, scales = "free", labeller = labeller(parameter = label_parsed), switch = "y") + 
           scale_y_continuous(breaks = scales::pretty_breaks(), labels = function(x) str_pad(x, width = 2, pad = "0")) + 
           theme_classic(base_size = 6) +
           theme(strip.placement = "outside", axis.title = element_blank(), strip.text.y = element_text(size = 11),
                 legend.position = "none", strip.text.x = element_text(size=15),
                 axis.text = element_text(size = 8))
         
       })
     
    # ## Re-order
     plot_list <- plot_list[c(1,8,2,9,3,10,4,11, 5,12, 6,13,7,14)]
     #plot_list <- plot_list[c(1,9,2,10,3,11,4,12, 5,13, 6,14,7,15,8,16)]
     
     # # Extract the y axis with strips
     # y_axis_plot <- plot_list[c(1,4,7)] %>% 
    # #   map(~ . + theme(axis.title = element_blank())) %>%
    # #   plot_grid(plotlist = ., ncol = 1)
    # 
    # 
    # ## Edit parts of the grid
    # # Remove strips
     plot_list[c(3:14)] <- plot_list[c(3:14)] %>%
       map(~ . + theme(strip.text.x = element_blank(), strip.background.x = element_blank()))
     plot_list[c(2,4,6,8,10,12,14)] <- plot_list[c(2,4,6,8,10,12,14)] %>%
       map(~ . + theme(strip.text.y = element_blank(), strip.background.y = element_blank()))
    # 
    # # Create the grid by rows
     top <- plot_grid(plotlist = plot_list[1:2], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     second <- plot_grid(plotlist = plot_list[3:4], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     third <- plot_grid(plotlist = plot_list[5:6], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     fourth <- plot_grid(plotlist = plot_list[7:8], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     fifth <- plot_grid(plotlist = plot_list[9:10], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     sixth <- plot_grid(plotlist = plot_list[11:12], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
     seventh <- plot_grid(plotlist = plot_list[13:14], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
    # eighth <- plot_grid(plotlist = plot_list[15:16], ncol = 2, align = "h", rel_widths = c(1, 0.9, 0.9))
    # 
    # # First plot
     #grobs <- ggplotGrob(plot_list[[1]])$grobs
     #legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
    p_grid <- plot_grid(top, second, third, fourth, fifth, sixth, seventh, ncol = 1, rel_heights = c(1, 0.9, 0.9))
    
    # 
    # ## Add axis
     y_axis <- grid::textGrob(label = "Observation", gp = grid::gpar(fontsize = 20), rot = 90)
     x_axis <- grid::textGrob(label = "Prediction", gp = grid::gpar(fontsize = 20))
    # 
    # # Plot again
    mu.musp.arrange <- grid.arrange(arrangeGrob(p_grid, left = y_axis, bottom = x_axis))
    ggsave(filename = "results/Figure2/ScatterPlot_MEAN_MUSP_Obs_Predicted_color_point_2Col_11_28_printed_r.png", plot = mu.musp.arrange , height = 13, width = 9, dpi = 1000)
    

####################################
# Figure 3. Observed by predicted varG
predObsvarG <- left_join(predObsvarG, ObsParam) %>% #join to bring CI information to plot
  select(Fam, trait, Group, estimate, prediction, CI_varG) %>% #keep CI information
 mutate(parameter = "varG")%>%
  mutate_at(c('trait', 'parameter'), as.factor) %>%
  left_join(., predAccBiasvarG, by = c("trait", "parameter")) %>%
  mutate(trait = str_replace_all(trait, trait.replace),
         annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))

##Needs to create individual plots due to drastic changes in variance scales

# p.fiber <- predObsvarG %>% filter(trait == "Fiber Content") %>% #filter trait
#   ggplot(aes(x=prediction , y = estimate)) +
#   geom_point(aes(colour= Group), size= 3, shape = "circle") +
#   geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
#   #facet_wrap(~ trait, ncol = 2, scales = "free") +
#   theme_classic() + 
#   geom_text(aes(x = Inf, y = -Inf, label = annotation), 
#             parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
#   scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
#   geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG), width=.00005)+  #add 95% CI bars adjust ymin and ymax with calculated values
#   labs(title = "Fiber Content") +
#   theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
#         plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.63*")

p.height <- predObsvarG %>% filter(trait == "Plant Height") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() +
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+ 
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=1.0, position=position_dodge(0.05))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Plant Height") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.57*")

p.lodging <- predObsvarG %>% filter(trait == "Lodging") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=.001, position=position_dodge(0.05))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Lodging") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.26 ")

p.maturity <- predObsvarG %>% filter(trait == "Days to Maturity") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=.5, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Days to Maturity") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.88*")

p.oil <- predObsvarG %>% filter(trait == "Oil") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=.003, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Oil") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.34 ")

p.protein <- predObsvarG %>% filter(trait == "Protein") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=.005, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Protein") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.17 ")

p.size <- predObsvarG %>% filter(trait == "Seed Size") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=.02, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Seed Size") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.64*")

p.yield <- predObsvarG %>% filter(trait == "Seed Yield") %>% #filter trait
  ggplot( aes(x=prediction , y = estimate)) +
  geom_point(aes(colour= Group), size= 3, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") +  #Remove SE bar from regression line
  #facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  geom_text(aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'))+
  geom_errorbar(aes(ymin= estimate-CI_varG, ymax= estimate+CI_varG),
                width=400, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
  labs(title = "Seed Yield") +
  theme(legend.position="none",axis.title = element_blank(), axis.text = element_text(size = 13),
        plot.title = element_text(size =18, hjust=0.5)) # + annotate("text", x = Inf, y = -Inf, hjust=1, vjust=-0.8, label = "r = 0.45*")

# Cowplot

varG.scatter.plot <- plot_grid(p.maturity, p.lodging, p.oil, p.height,
                              p.protein, p.size, p.yield, align='vh', ncol = 4, nrow = 2)
#create common x and y labels
y.grob <- textGrob(expression("Observed" ~ italic(sigma[G]^2)), gp=gpar(col="black", fontsize=22), rot=90)
x.grob <- textGrob(expression("Predicted" ~ italic(hat(sigma)[g]^2)), gp=gpar( col="black", fontsize=22))
# grobs <- ggplotGrob(p.yield)$grobs
# legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

#add to plot
varG.scatter.arrange <- grid.arrange(arrangeGrob(varG.scatter.plot, left = y.grob, bottom = x.grob))
ggsave(filename = "results/Figure3/ScatterPlot_VAR_G_Obs_Predicted_color_point_2Col_12_08_printed_r.png", plot = varG.scatter.arrange , height = 7, width = 15, dpi = 1000)


##########################
##Figure 4. Observed vs Predicted Correlation
#Create vector with trait names to replace the colnames
#Housekeeping to replace trait names
predObsrG <- predObsrG %>%
  mutate(trait_comb = paste(trait, trait2, sep = "_"))%>%
  mutate(parameter = "rG")%>%
  mutate_at(c('trait', "trait2", 'parameter'), as.factor) %>%
  left_join(., predAccBiasrG, by = c("trait", "trait2", "parameter")) %>%
  mutate(annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))


trait_comb <- unique(predObsrG$trait_comb) #extract chr vector with 28 trait combinations
trait.comb.replace <- setNames(c("Plant Height / Lodging","Plant Height / Days to Maturity", "Plant Height / Oil", "Plant Height / Protein", "Plant Height / Seed Size", "Plant Height / Seed Yield",
                                 "Lodging / Days to Maturity", "Lodging / Oil", "Lodging / Protein", "Lodging / Seed Size", "Lodging / Seed Yield", "Days to Maturity / Oil",
                                 "Days to Maturity / Protein", "Days to Maturity / Seed Size", "Days to Maturity / Seed Yield", "Oil / Protein", "Oil / Seed Size",
                                 "Oil / Seed Yield", "Protein / Seed Size", "Protein / Seed Yield", "Seed size / Seed Yield"), trait_comb)

predObsrG <- predObsrG %>%
  mutate(trait_comb = str_replace_all(trait_comb, trait.comb.replace)) %>%
  filter(trait_comb %in% c("Oil / Seed Yield", "Plant Height / Seed Yield", "Protein / Seed Yield"))
#Create scatterplot
Correl.scatter.plot <- ggplot(predObsrG, aes(x= prediction, y = estimate)) +
  geom_point(aes(colour= Group), size= 2.2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  geom_errorbar(aes(ymin= ci_rG_lower, ymax= ci_rG_upper),
                width=.01, position=position_dodge(0.04))+  #add 95% CI bars adjust ymin and ymax with calculated values
    facet_wrap(~ trait_comb, ncol = 4, scales = "free") +
  theme_classic() + 
  geom_text(data = distinct(predObsrG, trait_comb, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
            parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5) +
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) + #define colors to each group
  labs(x= expression("Predicted"  ~italic(hat(r)[g])), y=expression("Observed"  ~ italic(r[G]))) +
    theme(axis.title = element_text(size = 23), axis.text = element_text(size = 11.5),
        strip.text.x = element_text(size = 15), #change font size of facet grid labels
        legend.position="bottom", #legend position
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=18),
        legend.key.height=unit(0, "cm"),      
        plot.margin = unit(c(1,0.5,0,0.5), "lines")) 

#ggsave(filename = "results/Figure4/ScatterPlot_CORREL_Obs_Predicted_4Col_color_point_07_26.png", plot = Correl.scatter.plot, height = 15, width = 12, dpi = 1000)
ggsave(filename = "results/Figure4/ScatterPlot_CORREL_Obs_Predicted_3Col_color_point_01_11_printed_r.png", plot = Correl.scatter.plot, height = 13, width = 16, dpi = 1000)


#######################################################################
####Relationship between predicted mean and variance  #################
#Make plot
predVal <- predVal %>%
  mutate(trait = str_replace_all(trait, trait.replace)) #replace trait values by actual name
options(scipen = 999)
MuvarG.scatter.plot <- ggplot(predVal, aes(x= pred_mu, y = pred_varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  # geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
  #           parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  #stat_cor()+
  #stat_regline_equation() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")), coef.digits = 5) + #add regression line values
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Predicted" ~italic(hat(mu))), y=expression("Predicted" ~ italic(hat(sigma)[g]^2))) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure5/ScatterPlot_MEAN_varG_Predicted_1Col_color_point_11_28_printed_r.png", plot =MuvarG.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure5/ScatterPlot_MEAN_varG_Predicted_2Col_color_point_11_28_reg_line.png", plot = MuvarG.scatter.plot, height = 12, width = 9, dpi = 1000)

#Musp by varG
MuspvarG.scatter.plot <- ggplot(predVal, aes(y= pred_musp, x = pred_varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  # geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
  #           parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(y= expression("Predicted" ~italic(hat(mu)[sp])), x=expression("Predicted" ~ italic(hat(sigma)[g]^2))) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure5/ScatterPlot_MEAN_varG_Predicted_1Col_color_point_11_28_printed_r.png", plot =MuspvarG.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure5/ScatterPlot_MUSP_varG_Predicted_2Col_color_point_11_28_reg_line.png", plot = MuspvarG.scatter.plot, height = 12, width = 9, dpi = 1000)

#Musp by varG
MuspBymu.scatter.plot <- ggplot(predVal, aes(y= pred_musp, x = pred_mu)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  # geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
  #           parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(y= expression("Predicted" ~italic(hat(mu)[sp])), x=expression("Predicted" ~ italic(hat(mu)))) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure5/ScatterPlot_MEAN_MU_Predicted_1Col_color_point_11_28_printed_r.png", plot =MuspBymu.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure5/ScatterPlot_MUSP_MU_Predicted_2Col_color_point_11_28_reg_line.png", plot = MuspBymu.scatter.plot, height = 12, width = 9, dpi = 1000)

#################Observed Plots  ########################
ObsmuvarGWide <- ObsmuvarG %>% pivot_wider(names_from = parameter, values_from = value) %>%
  merge(., ParID) %>% #to bring group info
  mutate(trait = str_replace_all(trait, trait.replace)) #replace trait values by actual name
#remove scientific notation
options(scipen=999)

MUVARG.scatter.plot <- ggplot(ObsmuvarGWide, aes(x= mu, y = varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  stat_poly_eq(aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "*\", \"*")), coef.digits = 5) + #add regression line values
  # geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
  #           parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(x= expression("Observed" ~italic(mu)), y=expression("Observed" ~ italic(sigma)[G]^2)) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure6/ScatterPlot_MEAN_varG_Predicted_1Col_color_point_11_28_printed_r.png", plot =MUVARG.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure6/ScatterPlot_MEAN_varG_Observed_2Col_color_point_07_01_reg_line.png", plot = MUVARG.scatter.plot, height = 12, width = 9, dpi = 1000)

#MUSP and varG
MUSPVARG.scatter.plot <- ggplot(ObsmuvarGWide, aes(y= musp, x = varG)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  # geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
  #           parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
  labs(y= expression("Observed" ~italic(mu)[SP]), x=expression("Observed" ~ italic(sigma)[G]^2)) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure6/ScatterPlot_MEAN_varG_Predicted_1Col_color_point_11_28_printed_r.png", plot =MUSPVARG.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure6/ScatterPlot_MUSP_varG_Observed_2Col_color_point_11_28_reg_line.png", plot = MUSPVARG.scatter.plot, height = 12, width = 9, dpi = 1000)

#MUSP and MU
MUSPMU.scatter.plot <- ggplot(ObsmuvarGWide, aes(y= musp, x = mu)) +
  geom_point(aes(colour= Group), size= 2, shape = "circle") +
  geom_smooth(se=F, method = "lm", color="darkblue") + #Remove SE bar from regression line
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  theme_classic() + 
  # geom_text(data = distinct(predObsmu, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
  #           parse = TRUE, size = 4, hjust = 1.1, vjust = -0.5)+
  scale_colour_manual(values= c('#999999','#E69F00','#56B4E9'), labels= c("Exotic Ancestry", "Elite Line", "Plant Introduction")) +
   labs(y= expression("Observed" ~italic(mu)[SP]), x=expression("Observed" ~ italic(mu))) +
  #labs(x= expression(~italic(mu)), y = "")+
  theme( axis.text = element_text(size = 11),
         strip.text.x = element_text(size = 14), #change font size of facet grid labels
         legend.title = element_text(size=18), #change legend title font size
         legend.text = element_text(size=14),
         axis.title = element_text(size = 20),
         legend.position="bottom") #legend position #change legend text font size
#first option change lab row and facet_wrap(ncol=1)
#ggsave(filename = "results/Figure6/ScatterPlot_MEAN_MU_Predicted_1Col_color_point_11_28_printed_r.png", plot =MUSPMU.scatter.plot, height = 20, width = 6, dpi = 1000)
ggsave(filename = "results/Figure6/ScatterPlot_MUSP_MU_Observed_2Col_color_point_11_28_reg_line.png", plot = MUSPMU.scatter.plot, height = 12, width = 9, dpi = 1000)

#Correlation of predictive abilities
predAb <- read.csv("results/Pred_abilities_mu_varG_musp.csv")
cor.test(predAb$Mean, predAb$Var)
# Pearson's product-moment correlation
# 
# data:  predAb$Mean and predAb$Var
# t = -0.46379, df = 5, p-value = 0.6623
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.8293143  0.6492654
# sample estimates:
#        cor 
# -0.2030909
cor.test(predAb$Mean, predAb$Superior)
# Pearson's product-moment correlation
# 
# data:  predAb$Mean and predAb$Superior
# t = 2.1358, df = 5, p-value = 0.08577
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1299223  0.9497579
# sample estimates:
#       cor 
# 0.6907147 
cor.test(predAb$Var, predAb$Superior)
# Pearson's product-moment correlation
# 
# data:  predAb$Var and predAb$Superior
# t = 1.3709, df = 5, p-value = 0.2287
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3799204  0.9154200
# sample estimates:
#       cor 
# 0.5226766