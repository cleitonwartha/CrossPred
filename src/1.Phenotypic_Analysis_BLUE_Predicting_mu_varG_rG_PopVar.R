################################################################################
### Empirical Validation - Phenotypic analysis
###
### This script will run statistical analyses of the phenotypic data from the SoyNAM project/package
### Phenotypic analysis of validation family data: obtain BLUEs
### Author: Cleiton Wartha
### 
### > sessionInfo()
### R version 4.3.3 (2024-02-29 ucrt)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 19045) 
### Matrix products: default 
### locale:
### [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8
### LC_MONETARY=English_United States.utf8 LC_NUMERIC=C LC_TIME=English_United States.utf8
### time zone: America/Chicago
### tzcode source: internal
### attached base packages:
### [1] parallel  stats     graphics  grDevices utils     datasets  methods   base
### other attached packages:
### [1] PopVar_1.3.0  data.table_1.14.8 doParallel_1.0.17 iterators_1.0.14  foreach_1.5.2     lme4_1.1-35.1     asreml_4.2.0.257  Matrix_1.6-5      SoyNAM_1.6.2     
### [9] lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2       readr_2.1.4       tidyr_1.3.0       tibble_3.2.1     
### [17] ggplot2_3.4.4     tidyverse_2.0.0 
################################################################################

#Install packages and load libraries
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('SoyNAM')) install.packages('SoyNAM'); library(SoyNAM)
if (!require('asreml')) install.packages('asreml'); library(asreml)
if (!require('lme4')) install.packages('lme4'); library(lme4)
if (!require('parallel')) install.packages('parallel'); library(parallel)
if (!require('doParallel')) install.packages('doParallel'); library(doParallel)
if (!require('data.table')) install.packages('data.table'); library(data.table)
if (!require('PopVar')) install.packages('PopVar'); library(PopVar)
asreml.options(ai.sing = TRUE, fail= "soft")
options(scipen = 999) #prevent the triggering of scientific notation for integers
#setwd(fs::path_home()) #set the working directory 

#Define function to estimate check effect 
CHECK <-  function(trait){ test=reshape2::dcast(data.check.qa.raw,environ+spot~strain,value.var=trait,mean)
rownames(test)=test[,2];E=test[,1];test=test[,-c(1,2)];test=data.matrix(test);test[is.nan(test)]=NA;
X=function(X) unlist(tapply(X,E,FUN=function(x){m=mean(x,na.rm=T);SD=sd(x,na.rm=T);return((x-m)/SD)}))
MEAN=apply(test,2,X);C=rowMeans(MEAN,na.rm=T);names(C)=rownames(test);C[is.nan(C)]=0;return(C)}

#Load raw pheno data stored at SoyBase.org (SoyNAM populations ONLY - not the checks data)
pheno.raw <- as.data.frame(data.table::fread("https://www.soybase.org/SoyNAM/data/phenotypes/SoyNAM%20all%20familes%20phenotype%20data.txt", fill=T))

##Load SoyNAM progeny and parental marker data for cross prediction - imputed and coded -1, 0, 1 format
load("data/geno.RData")
dim(geno)

#Housekeeping of the Soybase pheno data
pheno.raw <- pheno.raw %>%
  rename(environ= Env, strain = "Corrected Strain", family= FamNo, set = FamSet,
         yield= "Yld (kg/ha)", maturity= "Days to Mat", height= "Ht (cm)", 
         lodging = Lod,protein= Protein, oil= Oil, size = "100 sdwt (g)") %>%
  mutate(spot = paste0(set, "_", environ)) %>% #create spot column for check adjustment
  mutate(ge = paste0(environ, "_", strain)) %>% #create ge column for concatenated values for the strain x environ interaction
  mutate_at(c("yield", "maturity", "height", "lodging", "protein", "oil", "size"), ~as.numeric(na_if(., "."))) #convert missing values from dots to NA

#Load phenotypic check data from SoyNAM package (cleaned but without set information)
data(soybase, envir=environment(), package="SoyNAM") #soybase is the quality assured data set
gen.qa <- NULL #remove geno file to clear memory

#Select all the checks and lines from SoyNAM package to create separate files for analysis
checks.chr <- as.character(unique(data.check.qa$strain))
line.chr <- as.character(unique(data.line.qa$strain))
#Check file
data.check.qa.raw <- pheno.raw %>% filter(strain %in% checks.chr)%>% #select only relevant columns
  dplyr::select(strain, environ, ge, family, spot, set, height, lodging, maturity, oil,
                protein, size, yield)
#RILs file  
data.line.qa.raw <-  pheno.raw %>% filter(!strain %in% checks.chr)%>% #select only relevant columns
  dplyr::select(strain, environ, ge, family, spot, set, height, lodging, maturity, oil,
                protein, size, yield) 
#Loop through traits to extract check correction for each set
data.line.qa.rawC <- data.line.qa.raw #Create file to store check info
for (i in c(7:13)){
  trait <- colnames(data.line.qa.raw)[i]
  print(trait)
  check = CHECK(trait);set = as.character(data.line.qa.raw[ , "spot"]) #calculate average correction from checks in each set
  C = check[set] #bring estimates back to the raw file format
  data.line.qa.rawC <- cbind(data.line.qa.rawC, C) #combine values into the raw file
  colnames(data.line.qa.rawC)[colnames(data.line.qa.rawC) == 'C'] <- paste0("C.", trait) #rename column
}

#Housekeeping to remove not needed columns
data.line <- data.line.qa.rawC %>%
  select(-spot, -set) %>%
  mutate_at(c('strain', 'environ'), as.factor) #convert model terms as factor
#Clean up RAM deleting unused files
data.check.qa <-NULL; data.check.qa.raw <- NULL
data.line.qa <- NULL; data.line.qa.raw <- NULL
data.line.qa.rawC <- NULL

################################################################################
############ MTMV models in AsremlR
################################################################################
#Univariate model to estimate genetic variances
traits <- c( "height", "lodging", "maturity", "oil", "protein", "size", "yield")

#Define formula for fixed effect of each individual trait - check correction
height.uv <- height ~ C.height + strain + environ
lodging.uv <- lodging ~ C.lodging + strain + environ
maturity.uv <- maturity ~ C.maturity + strain + environ
oil.uv <- oil ~ C.oil + strain + environ
protein.uv <- protein ~ C.protein + strain + environ
size.uv <- size ~ C.size + strain + environ
yield.uv <- yield ~ C.yield + strain + environ
#Create a list with the formulas for each trait
model.uv <-  list(height.uv, lodging.uv, maturity.uv, oil.uv,
                  protein.uv, size.uv, yield.uv)
blue.long <-c() #Create empty object to store results from loop
startTime <- Sys.time() #calculate start time
  #Loop over all traits
  for (t in 1:length(traits)){ #sequence of seven traits from the model
    #Define formula for specific trait combination
    form <- model.uv[[t]]
    print(traits[t])
    #Fit the univariate model
    uv.asr <- asreml(fixed = form, #adding the interaction to the model does not change the BLUEs and results in not enough workspace in regular personal computers
                     ai.sing= T,
                     workspace= "6gb",
                     data = data.line)
    #Extract the BLUEs from the model
    blue.asr <- as.data.frame(summary(uv.asr, coef = TRUE)$coef.fixed) # BLUEs
    int  <- blue.asr[rownames(blue.asr)%in% c("(Intercept)"), "solution"]  #retrieve intercept value
    #Create tidy long format file to store the results from all traits
    blue.asr <- blue.asr %>% rownames_to_column()  %>% 
      filter(str_detect(rowname, "strain_")) %>%  #filter rownames with strains
      mutate(strain = str_replace(rowname, "strain_", "")) %>% 
      mutate(value = solution + int) %>%
      mutate(family =gsub('DS1.-','', strain)) %>% #create Fam column
      mutate(family =gsub('...$','', family, perl = T)) %>%
      mutate(family =as.numeric(family)) %>%
      rename(std.error= 'std error') %>%
      select(strain, family, value, std.error) %>%
      mutate(trait = traits[t]) %>%
      relocate(c("family","trait"), .before = value)
    
    blue.long <- rbind(blue.long, blue.asr) #store results from each loop
  }
endTime <- Sys.time() #calculate end time of BLUE calculation
print(endTime - startTime) # print recorded time 7.70 min

#Convert BLUE file from long to wide format
pheno <- blue.long %>% select(-std.error) %>% #remove standard error not needed for downstream analysis
  pivot_wider(names_from = trait, values_from = value)

################################################################################
#Load SoyNAM composite map with genetic distances cM  - composite linkage map based on the RILs from all 39 families
#Supporting Information File 8 - Supplementary Table S5 from the article Genetic Characterization of the Soybean Nested Association Mapping Population
#https://doi.org/10.3835/plantgenome2016.10.0109
map <- as.data.frame(readxl::read_excel(path = "tpg2plantgenome2016100109-sup-0008.xlsx", skip = 1)) %>% #skip row 1 with concatenated header title to have a data frame
       rename(ssID= "NCBI ssID", SNPID = "SNP ID", Chr = "Linkage Group", Pos = "Linkage map position (cM)") %>%
       filter(ssID %in% colnames(geno))  %>% # remove three markers not present in the QCed progeny marker file Gm12_15420968_G_T, Gm12_17864392_T_C, and Gm13_20431930_G_A
       select(-SNPID) #remove SNPID column 

ParID <- data.frame(strain = c('TN05-3027','4J105-3-4','5M20-2-5-2','CL0J095-4-6','CL0J173-6-8','HS6-3976','Prohio','LD00-3309','LD01-5907','LD02-4485','LD02-9050','Magellan',
        'Maverick','S06-13640','NE3001','Skylla','U03-100612','LG03-2979','LG03-3191','LG04-4717','LG05-4292','LG05-4317','LG05-4464','LG05-4832','LG90-2550','LG92-1255',
        'LG94-1128','LG94-1906','LG97-7012','LG98-1605','LG00-3372','LG04-6000','PI398881','PI427136','PI437169B','PI518751','PI561370','PI404188A','PI574486'),
        family = c(2,3,4,5,6,8,9,10,11,12,13,14,15,17,18,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,41,42,48,50,54,64),
        group= c('EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','EL','BX','BX','BX','BX','BX','BX',
        'BX','BX','BX','BX','BX','BX','BX','BX','BX','PI','PI','PI','PI','PI','PI','PI'))
#Housekeeping on the pheno file to only have progenies that have marker data
pheno <- pheno %>% subset(strain %in% row.names(geno))
#Include parental lines in the pheno file because all genotypes must be part of both files
#To attend mixed.solve(y = y, Z = M, method = "REML") : nrow(Z) == n
ParPheno <- data.frame(matrix(NA, nrow = 40, ncol = length(pheno))) #create empty data frame with NAs
names(ParPheno) <- names(pheno)
ParPheno$strain <- c('Parent_IA3023','TN05-3027','4J105-3-4','5M20-2-5-2','CL0J095-4-6','CL0J173-6-8','HS6-3976','Prohio','LD00-3309','LD01-5907','LD02-4485','LD02-9050','Magellan',
                     'Maverick','S06-13640','NE3001','Skylla','U03-100612','LG03-2979','LG03-3191','LG04-4717','LG05-4292','LG05-4317','LG05-4464','LG05-4832','LG90-2550','LG92-1255',
                     'LG94-1128','LG94-1906','LG97-7012','LG98-1605','LG00-3372','LG04-6000','PI398881','PI427136','PI437169B','PI518751','PI561370','PI404188A','PI574486')
pheno <- rbind(pheno, ParPheno) #combine pheno file with progeny and parents
fam <- unique(ParID$family) #Create a vector with the 39 families


mrkPar <- subset(geno, row.names(geno) %in% ParPheno$strain)
################################################################################
#####            Run PopVar::pop_predict2 function
#### Generates predictions of the genetic variance and genetic correlation in
#### bi-parental populations using a set of deterministic equations instead of simulations.
################################################################################
## With the large training set N ~ 5000 - use parallel computing for optimization
# Number of cores
nCores<- detectCores() - 1
# Starting parallel
registerDoParallel(nCores)
startTime <- Sys.time() #calculate start time
################################################################################
### WARNING: foreach loop can take several hours due to the training pop > 5000 ind
### Recommendation is to run in a supercomputer cluster to have more cores available
### Alternatively, computational time reduces dramatically if run with smaller training pop size
### The output object predVal with the predicted values is available in obj3.RData
### for downstream analyzes in the following script #3
################################################################################

# Foreach loop over each family c(2:6,8:15,17,18,22:34,36:42,48,50,54,64)
predVal <-foreach(i = c(2:6,8:15,17,18,22:34,36:42,48,50,54,64), .packages = c("tidyverse","PopVar"), #family numbers hard coded in the numeric vector
                  .combine = rbind) %dopar% {
          PhenoSub <- subset(pheno, family != i | is.na(family)) %>% #remove family i and keep NAs
                      select(-family) %>% #remove Fam column to only have strains and phenotypes
                      select(strain, maturity, yield)
          DropFam <- pheno %>%filter(family == i) #select strains from family to drop
          DropFam <- DropFam$strain #strain IDs as character
          #Subset geno file
          mrkSub <- subset(geno, !row.names(geno) %in% DropFam) #remove lines from the family
          #Select parental line for the given family being predicted
          ParIDPred <- ParID %>% filter(family == i) #Select family being predicted
          ParIDPred <- ParIDPred$strain #convert to character vector
          parents <- c("Parent_IA3023", ParIDPred) #character vector of parents
          #Run PopVar::pop_predict2
          out <- PopVar::pop_predict2(M = mrkSub, y.in= PhenoSub, map.in = map, parents = parents,
                                      tail.p = 0.1, models = c("rrBLUP"))
          }
endTime <- Sys.time() #calculate end time of BLUE calculation
print(endTime - startTime) # print recorded time - 1.46h for only maturity and yield in seven families 2,3,4,5,6,8 (one fam per core)
#write.csv(predVal, "Pred_PopVar_SoyNAM_mu_Vg_rG_ALL_fam.csv", row.names = F) #save output file
##################    End of the run   #########################################
