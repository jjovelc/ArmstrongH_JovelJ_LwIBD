library(ZIBR)
library(tidyverse)
library(nlme)

setwd('/Users/juanjovel/jj/data_analysis/heatherArmstrong/living_w_IBD/ZIBR')

data <- read.table("all_samples_mpa_tax_shortNames_37.tsv", header = T, 
                   row.names = 1, sep = '\t')

metadata <- read.table("metadata_37.tsv", header = T, 
                   sep = '\t', row.names = 1)

# Extract columns from the 6th onward with an average of 10 or more
data_counts <- data[, 6:ncol(data)][, sapply(data[, 6:ncol(data)], mean) >= 10]

# Combine the first 5 columns with the filtered columns
#data <- cbind(data[, 1:5], filtered_cols)
data_counts <- data_counts/colSums(data_counts)

cat('samples','taxa',dim(data_counts),'\n')
data_counts[1:5,1:3]

### Filter low abundant bacterial data
filter.index1 <- apply(data_counts,2,function(X){sum(X>0)>0.4*length(X)})
filter.index2 <- apply(data_counts,2,function(X){quantile(X,0.9)>1})
taxa.filter <- data_counts[,filter.index1 & filter.index2]
taxa.filter <- 100*sweep(taxa.filter, 1, rowSums(taxa.filter), FUN="/")
cat('after filter:','samples','taxa',dim(taxa.filter),'\n')
cat(colnames(taxa.filter),'\n')
head(rowSums(taxa.filter))

### 
taxa.data.lwIBD <- taxa.filter
dim(taxa.data.lwIBD)

#### create covariates, 
#### timePoint, remission
reg.cov.lwIBD <- data.frame(Sample = rownames(taxa.data.lwIBD), stringsAsFactors = FALSE) %>% 
  left_join(tibble::rownames_to_column(metadata, var = 'Sample'), by = 'Sample') %>%
  dplyr::select(Sample, timePoint, patient, severity, diagnosis, sex) %>%
  mutate(Treat = ifelse(severity == 'flare', 1, 0),
         patient = paste('S', patient, sep = ''),
         Time = ifelse(timePoint == 'wk0', 0, ifelse(timePoint == 'wk26', 1, ifelse(timePoint == 'wk52', 2, NA))),
         timePoint.X.Severity = Time * Treat) %>%
  as.data.frame


### take out first time point
reg.cov.lwIBD.t1   <-  subset(reg.cov.lwIBD,Time==0)
rownames(reg.cov.lwIBD.t1) <- reg.cov.lwIBD.t1$patient
reg.cov.lwIBD.t23 <-  subset(reg.cov.lwIBD,Time!=0)

reg.cov.lwIBD.t23 <- data.frame(
  baseline.sample = reg.cov.lwIBD.t1[reg.cov.lwIBD.t23$patient, 'Sample'],
  baseline.patient = reg.cov.lwIBD.t1[reg.cov.lwIBD.t23$patient, 'patient'],
  reg.cov.lwIBD.t23,
  stringsAsFactors = FALSE
)


#### Fit ZIBR model on the real data
spe.all <- colnames(taxa.data.lwIBD)
p.species.list.zibr <- list()
p.species.list.lme <- list()

for (spe in spe.all){
  #spe = 'g__Collinsella'
  ###### create covariates
  X <- data.frame(
    Baseline=taxa.data.lwIBD[reg.cov.lwIBD.t23$baseline.sample, spe]/100,
    #reg.cov.t234[,c('log.days','Delivery','Delivery.X.log.days')]
    reg.cov.lwIBD.t23[,c('Time','Treat')]
  )
  rownames(X) <- reg.cov.lwIBD.t23$Sample
  Z <- X
  patient.ind <- reg.cov.lwIBD.t23$patient
  time.ind   <- reg.cov.lwIBD.t23$Time
  ####
  cat(spe,'\n')
  Y <- taxa.data.lwIBD[reg.cov.lwIBD.t23$Sample, spe]/100
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  ####################
  ## fit linear random effect model with arcsin square transformation on Y
  tdata <- data.frame(Y.tran=asin(sqrt(Y)),X,SID=patient.ind)
  lme.fit <- lme(Y.tran ~ Baseline + Time + Treat,random=~1| SID, data = tdata)
  coef.mat <- summary(lme.fit)$tTable[-1,c(1,5)]
  p.species.list.lme[[spe]] <- coef.mat[,2]
  ###################
  if (sum(Y>0)<10 | sum(Y==0) <10 | max(Y)<0.01){
    print('skip')
    p.species.list.zibr[[spe]] <- p.species.list.lme[[spe]]
    next
  }else{
    est <- zibr(logistic.cov=X,beta.cov=Z,Y=Y,
                subject.ind=patient.ind,
                time.ind=time.ind,
                quad.n=30,verbose=TRUE)
    p.species.list.zibr[[spe]] <- est$joint.p
    
  }
  #break
}

### adjust p values
p.species.zibr <- t(as.data.frame(p.species.list.zibr))
p.species.zibr.adj <-
  tibble::rownames_to_column(as.data.frame(p.species.zibr),var = 'Species') %>% 
  mutate_each(funs(p.adjust(.,'fdr')),-Species)

write.table(p.species.zibr.adj,
            file=paste('4_Results/Real_Data_Estimation_Results_flare_ZIBR.tsv',
                       sep=''), sep = '\t', quote = F, row.names = F)


p.species.lme <- t(as.data.frame(p.species.list.lme))
p.species.lme.adj <-  
  add_rownames(as.data.frame(p.species.lme),var = 'Species') %>% mutate_each(funs(p.adjust(.,'fdr')),-Species)
write.table(p.species.lme.adj,
            file=paste('4_Results/Real_Data_Estimation_Results_flare_LME.tsv',
                       sep=''), sep = '\t', quote = F, row.names = F)



