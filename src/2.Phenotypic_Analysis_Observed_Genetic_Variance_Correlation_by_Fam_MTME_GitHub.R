################################################################################
### Empirical Validation - Phenotypic analysis
###
### This script will run statistical analyses of the phenotypic data from the SoyNAM project/package
### Phenotypic analysis of validation family data: obtain Vg and rG
### Author: Cleiton Wartha
### > sessionInfo()
### R version 4.3.3 (2024-02-29 ucrt)
### Platform: x86_64-w64-mingw32/x64 (64-bit)
### Running under: Windows 10 x64 (build 19045)
### 
### Matrix products: default
### locale:
### [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
### [5] LC_TIME=English_United States.utf8    
### time zone: America/Chicago
### tzcode source: internal
### attached base packages:
###   [1] stats     graphics  grDevices utils     datasets  methods   base     
### other attached packages:
### [1] gtools_3.9.5      data.table_1.14.8 lme4_1.1-35.1     asreml_4.2.0.257  Matrix_1.6-5      lubridate_1.9.3   forcats_1.0.0     stringr_1.5.1    
### [9] dplyr_1.1.4       purrr_1.0.2       readr_2.1.4       tidyr_1.3.0       tibble_3.2.1      ggplot2_3.4.4     tidyverse_2.0.0   SoyNAM_1.6.2 
################################################################################

#Install packages and load libraries
if (!require('SoyNAM')) install.packages('SoyNAM'); library(SoyNAM)
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('asreml')) install.packages('asreml'); library(asreml)
if (!require('lme4')) install.packages('lme4'); library(lme4)
if (!require('data.table')) install.packages('data.table'); library(data.table)
if (!require('gtools')) install.packages('gtools'); library(gtools)
if (!require('ASRgenomics')) install.packages('ASRgenomics'); library(ASRgenomics)
options(scipen = 999) #prevent the triggering of scientific notation for integers
#setwd(fs::path_home()) #set the working directory 
#Define function to estimate check effect 
CHECK <-  function(trait){ test=reshape2::dcast(data.check.qa.raw,environ+spot~strain,value.var=trait,mean)
rownames(test)=test[,2];E=test[,1];test=test[,-c(1,2)];test=data.matrix(test);test[is.nan(test)]=NA;
X=function(X) unlist(tapply(X,E,FUN=function(x){m=mean(x,na.rm=T);SD=sd(x,na.rm=T);return((x-m)/SD)}))
MEAN=apply(test,2,X);C=rowMeans(MEAN,na.rm=T);names(C)=rownames(test);C[is.nan(C)]=0;return(C)}

#Load raw pheno data stored at SoyBase.org (SoyNAM populations ONLY - not the checks data)
pheno.raw <- as.data.frame(data.table::fread("https://www.soybase.org/SoyNAM/data/phenotypes/SoyNAM%20all%20familes%20phenotype%20data.txt", fill=T))

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
  select(-spot, -set) 
#Clean up RAM deleting unused files
data.check.qa <-NULL; data.check.qa.raw <- NULL
data.line.qa <- NULL; data.line.qa.raw <- NULL
data.line.qa.rawC <- NULL

################################################################################
############ MTMV models in AsremlR
################################################################################
#Univariate model to estimate genetic variances
traits <- c( "height", "lodging", "maturity", "oil", "protein", "size", "yield")
fam <- unique(data.line$family)
#Define formula for fixed effect of each individual trait - check correction
height.uv <- height ~ C.height 
lodging.uv <- lodging ~ C.lodging
maturity.uv <- maturity ~ C.maturity
oil.uv <- oil ~ C.oil
protein.uv <- protein ~ C.protein
size.uv <- size ~ C.size
yield.uv <- yield ~ C.yield
#Create a list with the formulas for each trait
model.uv <-  list(height.uv, lodging.uv, maturity.uv, oil.uv,
                  protein.uv, size.uv, yield.uv)
varComp.uv <-c() #Create empty object to store results from loop

startTime <- Sys.time() #calculate start time
###Loop through the families and traits
for (f in fam){ #numerical vector of the families
  data.sub <- data.line %>% filter(family == f) %>% #create file for ASRemlR function
    mutate_at(c('strain', 'environ'), as.factor) #convert random terms as factor
  
  #Loop over all trait combinations
  for (t in 1:7){ #sequence of seven traits from the model
    #Define formula for specific trait combination
    form <- model.uv[[t]]
    #Fit the univariate model
    uv.asr <- asreml(fixed = form, 
                     random = ~ strain + environ, #+ environ:strain - interaction argument removed because of singularities in the Average Information matrix
                     maxit = 1000,
                     data = data.sub)
    #Retrieve the variance components from the model
    varComp <- summary(uv.asr)$varcomp %>% rownames_to_column(.) %>%
      data.frame(family= f,  trait=traits[t], .)
    varComp.uv <-rbind(varComp.uv, varComp) #store results from each loop
  }
}
endTime <- Sys.time() #calculate end time
print(endTime - startTime) # print recorded time 32.83392 secs

#Retrieve only genetic variances from univariate model
varG.uv <- varComp.uv %>% filter(rowname == "strain") %>% # filter only the genetic variances
  rename(varG = component) %>% select(-rowname)

####### Bivariate model to obtain the genetic correlations
trait.comb <- gtools::combinations(length(traits), 2) #create trait combinations
trait.comb.chr <- data.frame(trait1= traits[trait.comb[, 1]], trait2= traits[trait.comb[, 2]])
varComp.bv <-c() #Create empty object to store results from loop
asreml.options(ai.sing = T, fail= "soft") #arguments needed because some of the variance components were near zero in the rG model below
startTime <- Sys.time() #calculate start time
###Loop through the families and trait combinations
for (f in fam){ #numerical vector of the families
  data.sub <- data.line %>% filter(family == f) %>% #create file for ASRemlR function
    mutate_at(c('strain', 'environ'), as.factor) #convert random terms as factor
  for (t in 1:nrow(trait.comb.chr)){   #Loop over all trait combinations
    #Define formula for specific trait combination accounting for the respective trait covariates
    form <- as.formula(paste0("cbind(",paste0(trait.comb.chr[t, 1]),",", paste0(trait.comb.chr[t, 2]), ") ~ trait + at(trait, 1):C.", paste0(trait.comb.chr[t, 1]), "+ at(trait, 2):C.", paste0(trait.comb.chr[t, 2])))
    #Fit the bivariate model
    bv.asr <- asreml(fixed = form, #fixed component takes the checks correction for each trait
                     random = ~ corgh(trait):strain + trait:environ, # + trait:environ:strain, #interaction argument removed because of singularities in the Average Information matrix
                     residual = ~ id(units):us(trait),
                     na.action = na.method(x = "include", y = "include"),
                     data = data.sub)
    #Retrieve the variance components from the model
    varComp <- summary(bv.asr)$varcomp %>% rownames_to_column(.) %>%
      data.frame(family= f,  trait1=trait.comb.chr[t, 1], trait2=trait.comb.chr[t, 2], .)
    varComp.bv <-rbind(varComp.bv, varComp) #store results from each loop
  }
}
endTime <- Sys.time() #calculate end time
print(endTime - startTime) # print recorded time 2.86 min #wGxE 7.274104 mins
#Retrieve only genetic correlations from bivariate model
rG.bv <- varComp.bv %>% filter(grepl('trait:strain', rowname),
                               grepl('.cor', rowname)) %>% #filter the genetic components
  rename(rG = component) %>% select(-rowname)
####################################
#Full multi-trait multi-environment model fitted jointly 
# at(family) included to estimate parameters within a family
#The model below can become very complex to fit because of the large number of parameters to
#be estimated jointly with large number of traits and families - not reaching model convergence
## It is an elegant model that might work with fewer indiv/families and traits but did not work with the SoyNAM data set
##returned variances and correlations have the same value for all families instead of individual estimates for each family
data.mv <- data.line %>% mutate_at(c('strain', 'family', 'environ', "ge"), as.factor) #convert random terms as factor
startTime <- Sys.time() #calculate start time
mv.asr <- asreml(fixed = cbind(height, lodging, maturity, oil, protein, size, yield) ~ trait + at(trait, 1):C.height + at(trait, 2):C.lodging +
                   at(trait, 3):C.maturity + at(trait, 4):C.oil + at(trait, 5):C.protein + at(trait, 6):C.size + at(trait, 7):C.yield, 
                 random = ~ corgh(trait):at(family):strain + trait:at(family):environ, # + trait:at(family):ge, ##interaction argument removed because of singularities in the Average Information matrix
                 #residual = ~ id(units):us(trait), #generic residual structure
                 residual = ~ dsum(~units:us(trait)|family),
                 workspace.mem = "8gb", na.action = na.method(x = "include", y = "include"),
                 maxit = 1000,
                 data = data.mv)
endTime <- Sys.time() #calculate end time
print(endTime - startTime) #3.16 mins
#Retrieve the variance components from the .asr object
#Genetic correlation - multivariate multi-trait model
rG.mv <- summary(mv.asr)$varcomp %>% rownames_to_column(., var = "Component") %>%
  filter(grepl('trait:strain', Component),
         grepl('.cor', Component)) %>% #filter the genetic components
  separate(Component, into= c("Component","trait1", "trait2"), sep= "trait!") %>%
  mutate_all(~gsub(":!|.cor", "", .)) %>%
  mutate(family = readr::parse_number(Component)) %>%
  rename(rG = component) %>% mutate_at(c("rG"), as.numeric) %>%
  select(-Component) %>%
  relocate(family, .before = trait1)
#Genetic variance - multivariate multi-trait model
varG.mv <- summary(mv.asr)$varcomp %>% rownames_to_column(., var = "Component") %>%
  filter(grepl('trait:strain!trait_', Component)) %>%
  separate(Component, into= c("Component","trait"), sep= "_") %>%
  mutate(family = readr::parse_number(Component)) %>%
  rename(varG = component) %>% mutate_at(c("varG"), as.numeric) %>%
  select(-Component) %>%
  relocate(family, .before = trait)
#################################################################################
###  Phenotypic reliability using the average variance of a difference between a pair
###  of BLUPs (vBLUP) or the predicted error variance (PEV) - these methods are more suitable for unbalanced data
################################################################################
Reliability.uv <-c() #Create empty object to store results from loop
startTime <- Sys.time() #calculate start time
###Loop through the families and traits
for (f in fam){ #numerical vector of the families
  data.sub <- data.line %>% filter(family == f) %>% #create file for ASRemlR function
    mutate_at(c('strain', 'environ', "ge"), as.factor) #convert random terms as factor
  for (t in 1:7){ #sequence of seven traits from the model
    #Define formula for specific trait combination
    form <- model.uv[[t]]
    #Fit the univariate model
    uv.asr <- asreml(fixed = form, 
                     random = ~ strain + environ,  #+ environ:strain, #interaction argument removed because of singularities in the Average Information matrix
                     data = data.sub)
    varG <- summary(uv.asr)$varcomp %>% rownames_to_column(.) %>% filter(rowname == "strain") %>% select(component) %>% as.numeric() #genetic variance
    #Retrieve the predicted error variance
    PEV <- summary(uv.asr, coef=TRUE)$coef.random %>% data.frame() %>% rownames_to_column(.) %>%
      filter(grepl('strain', rowname)) %>%
      mutate(PEV = std.error^2) %>%cbind(., varG) %>%
      mutate(Reliability.PEV = mean(1 - (PEV/varG))) %>% #average across all BLUPs
      distinct(Reliability.PEV) #retrieve single mean value
    #Retrieve the mean variance of a difference of two BLUPs
    vBLUP <- predict(uv.asr, "strain", vcov=TRUE)$avsed[2]^2 #mean avsed (average standard error of the difference between a pair of BLUPs) squared to obtain the mean variance
    Reliability <- summary(uv.asr)$varcomp %>% rownames_to_column(.) %>%
      filter(rowname == "strain") %>% cbind(., vBLUP) %>% #add the variance of difference of two BLUPs
    mutate(Reliability.vBLUP = 1 - (vBLUP/(2*component))) %>%
      cbind(., PEV) %>% #add the reliability estimate coming from PEV
      select(-rowname) %>% rename(varG = component) %>% #housekeeping
      data.frame(family= f,  trait=traits[t], .)
    Reliability.uv <-rbind(Reliability.uv, Reliability) #store results from each loop
  }
}
endTime <- Sys.time() #calculate end time 1.671527 min
print(endTime - startTime) # print recorded time
#Summary of the phenotypic reliabilities obtained by two different methods - based on vBLUP or PEV
Reliab.Smry <- Reliability.uv %>% group_by(trait) %>% #group results by trait
  summarise(Mean.Reliability.vBLUP = mean(Reliability.vBLUP),
            Min.Reliability.vBLUP = min(Reliability.vBLUP),
            Max.Reliability.vBLUP = max(Reliability.vBLUP),
            Mean.Reliability.PEV = mean(Reliability.PEV),
            Min.Reliability.PEV = min(Reliability.PEV),
            Max.Reliability.PEV = max(Reliability.PEV))
#write.csv(Reliab.Smry, "Phenotypic_Reliability_Summary_vBLUP_PEV.csv", row.names = F)
################################################################################
####### Alternative method to obtain the observed varG and rG by including the Genomic Relationship Matrix (GRM)
####### Method recommended for when all the phenotyped individuals have marker data available
####### Not all SoyNAM RILs have marker data that passed quality control and therefore the GRM was not included to avoid undersampling of indiv per family
################################################################################
##Load SoyNAM progeny and parental marker data
load("data/geno.RData"); dim(geno) #Marker scores in the -1,0,1 format
varComp.uv <-c() #Create empty object to store results from loop
startTime <- Sys.time() #calculate start time
###Loop through the families and traits
for (f in fam){ #numerical vector of the families
  data.sub <- data.line %>% filter(family == f) %>% mutate_at(c('strain', 'environ', "ge"), as.factor) #create file for ASRemlR function
  strains <- unique(data.sub$strain)
  cat("Family:", f,'\n')
  cat("Number of phenotyped strains:", length(strains),'\n')
  M.sub <- geno[rownames(geno) %in% strains, ]  #filter only marker data of genotypes that have marker data
  cat("Number of genonotyped strains:", length(rownames(M.sub)),'\n')
  #Second filtering step to remove pheno data of individuals without marker data
  #This step is needed to avoid levels of the factor strain with pheno data and no information in the GRM - causing singularities in the AI matrix
  data.sub <- data.line %>% filter(strain %in% rownames(M.sub)) %>% droplevels() %>%  #create file for ASRemlR function
    mutate_at(c('strain', 'environ', "ge"), as.factor) #convert random terms as factor
  cat("Number of phenotyped strains used in the analysis:", length(unique(data.sub$strain)),'\n')
  #GRM <- rrBLUP::A.mat(M.sub, shrink= T) #alternative method to obtain the GRM using the rrBLUP package - could be used in lieu of the ASRgenomics G.matrix() and G.tuneup()
  GRM <- ASRgenomics::G.matrix(M.sub+1, method='VanRaden')$G  #Plus 1 operation done to convert from -1,0,1 to 0,1,2 marker score format required by ASRgenomics package
  GRM <- ASRgenomics::G.tuneup(G=GRM, bend=TRUE)$Gb #Bending step needed to ease the matrix inversion step - consists on adjusting the original G matrix to obtain a near positive definite matrix, which is done by making all negative or very small eigenvalues slightly positive.
  Ginv <- ASRgenomics::G.inverse(GRM, sparseform = T)$Ginv  #Advisable to input the inverse of the GRM in the sparse format - function uses chol2inv() - inverse from Cholesky Decomposition
  #ASReml-R does not handle the GRM well or with the same efficiency internally - it relies on solve() - incurring in additional computational time
  #Loop over all trait combinations
  for (t in 1:7){ #sequence of seven traits from the model
    #Define formula for specific trait combination
    form <- model.uv[[t]]
    #Fit the univariate model
    uv.asr <- asreml(fixed = form, 
                     random = ~ vm(strain, Ginv) + environ,  #+ environ:strain, #interaction argument removed because of singularities in the Average Information matrix
                     maxit = 100,
                     data = data.sub)
    #Retrieve the variance components from the model
    varComp <- summary(uv.asr)$varcomp %>% rownames_to_column(.) %>%
      data.frame(family= f,  trait=traits[t], .)
    varComp.uv <-rbind(varComp.uv, varComp) #store results from each loop
  }
}
endTime <- Sys.time() #calculate end time 1.671527 min
print(endTime - startTime) # print recorded time

#Retrieve only genetic variances from univariate model
varG.uvGRM <- varComp.uv %>%  #filter only the genetic variances
  filter(rowname == "vm(strain, Ginv)") %>%   
  rename(varG = component) %>% select(-rowname)
################################################################################
####### Bivariate model to obtain the genetic correlations
trait.comb <- gtools::combinations(length(traits), 2) #create trait combinations
trait.comb.chr <- data.frame(trait1= traits[trait.comb[, 1]], trait2= traits[trait.comb[, 2]])
varComp.bv <-c() #Create empty object to store results from loop
startTime <- Sys.time() #calculate start time
###Loop through the families and trait combinations
for (f in fam){ #numerical vector of the families
  data.sub <- data.line %>% filter(family == f) %>% mutate_at(c('strain', 'environ', "ge"), as.factor) #create file for ASRemlR function
  strains <- unique(data.sub$strain)
  cat("Family:", f,'\n')
  cat("Number of phenotyped strains:", length(strains),'\n')
  M.sub <- geno[rownames(geno) %in% strains, ]  #filter only marker data of genotypes that have marker data
  cat("Number of genonotyped strains:", length(rownames(M.sub)),'\n')
  #Second filtering step to remove pheno data of individuals without marker data
  #This step is needed to avoid levels of the factor strain with pheno data and no information in the GRM - causing singularities in the AI matrix
  data.sub <- data.line %>% filter(strain %in% rownames(M.sub)) %>% droplevels() %>%  #create file for ASRemlR function
    mutate_at(c('strain', 'environ', "ge"), as.factor) #convert random terms as factor
  cat("Number of phenotyped strains used in the analysis:", length(unique(data.sub$strain)),'\n')
  #GRM <- rrBLUP::A.mat(M.sub, shrink= T) #alternative method to obtain the GRM using the rrBLUP package - could be used in lieu of the ASRgenomics G.matrix() and G.tuneup()
  GRM <- ASRgenomics::G.matrix(M.sub+1, method='VanRaden')$G  #Plus 1 operation done to convert from -1,0,1 to 0,1,2 marker score format required by ASRgenomics package
  GRM <- ASRgenomics::G.tuneup(G=GRM, bend=TRUE)$Gb #Bending step needed to ease the matrix inversion step - consists on adjusting the original G matrix to obtain a near positive definite matrix, which is done by making all negative or very small eigenvalues slightly positive.
  Ginv <- ASRgenomics::G.inverse(GRM, sparseform = T)$Ginv  #Advisable to input the inverse of the GRM in the sparse format - function uses chol2inv() - inverse from Cholesky Decomposition
  #ASReml-R does not handle the GRM well or with the same efficiency internally - it relies on solve() - incurring in additional computational time
  for (t in 1:nrow(trait.comb.chr)){   #Loop over all trait combinations
    cat(t,'\n')
    #Define formula for specific trait combination accounting for the respective trait covariates
    form <- as.formula(paste0("cbind(",paste0(trait.comb.chr[t, 1]),",", paste0(trait.comb.chr[t, 2]), ") ~ trait + at(trait, 1):C.", paste0(trait.comb.chr[t, 1]), "+ at(trait, 2):C.", paste0(trait.comb.chr[t, 2])))
    #Fit the bivariate model
    bv.asr <- asreml(fixed = form, #fixed component takes the checks correction for each trait
                     random = ~ corgh(trait):vm(strain, Ginv) + trait:environ, # + trait:environ:strain, #interaction removed because the model did not converge in several instances with estimates near 1 or -1
                     residual = ~ id(units):us(trait),
                     maxit = 1000,
                     na.action = na.method(x = "include", y = "include"),
                     data = data.sub)
    #Retrieve the variance components from the model
    varComp <- summary(bv.asr)$varcomp %>% rownames_to_column(.) %>%
      data.frame(family= f,  trait1=trait.comb.chr[t, 1], trait2=trait.comb.chr[t, 2], .)
    varComp.bv <-rbind(varComp.bv, varComp) #store results from each loop
  }
}
endTime <- Sys.time() #calculate end time
print(endTime - startTime) # print recorded time  5.065258 min
#Retrieve only genetic correlations from bivariate model
rG.bvGRM <- varComp.bv %>% filter(grepl('trait:vm', rowname), 
                               grepl('.cor', rowname)) %>% #filter the genetic correlation components
  rename(rG = component) %>% select(-rowname)

######Phenotyped versus genotyped individuals per SoyNAM family
fam = rownames(geno)
fam = gsub("DS1.-", "", fam)
fam = gsub("...$", "", fam, perl = T)
fam = as.numeric(fam)
GenoN <- fam %>% data.frame() %>% rename(family = ".") %>% drop_na(family) %>%
  group_by(family) %>%  summarise(N.geno = n())
PhenoN <- data.line %>% group_by(family) %>%
  summarise(N.pheno = n_distinct(strain))
PhenoGenoN <- left_join(PhenoN, GenoN) %>%
  mutate(Diff = N.pheno - N.geno) %>% arrange(desc(Diff))
#write.csv(PhenoGenoN, "N_Count_Phenotyped_Genotyped_RIL_per_Family.csv", row.names = F)
