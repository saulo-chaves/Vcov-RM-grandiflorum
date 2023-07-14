rm(list = ls())
# Analysis of multi-harvest data through mixed models: 
# an application in Theobroma grandiflorum breeding

# Crop Science

## Loading the packages
require(asreml)
require(tidyverse)
require(gghighlight)
require(ggrepel)
require(ComplexHeatmap)
require(nadiv)
require(patchwork)

## Import dataset
dataset = read.table("dataset.txt", header = TRUE)
dataset = dataset %>% group_by(Years, Plots, Replicates, Hybrids) %>% 
  summarise(pf = mean(pf,na.rm=T))
pedigree = read.table('pedigree.txt', header= T)

## Transform numeric variables into factors
dataset <- transform(dataset, Years = factor(Years), 
                     Plots = factor(Plots), 
                     Replicates = factor(Replicates), 
                     Hybrids = factor(Hybrids),
                     Hybrids2 = factor(Hybrids))

## Relationship matrices
Amat = as.matrix(makeA(pedigree = pedigree))

##### Modelling the additive genetic effects ------

#CS
mod1 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat)+ vm(Hybrids,Amat):idv(Years) + Plots,
              data = dataset, maxiter=100, 
              na.action = na.method(x="include", y = "include"))
sum1 = summary(mod1)$varcomp
aic1 = summary(mod1)$aic
logl1 = summary(mod1)$loglik

predm1_vcov = predict(mod1, classify = "Hybrids", vcov = T)
predm1_sed = predict(mod1, classify = "Hybrids", sed = T)

PEV = mean(diag(predm1_vcov$vcov)) 
MVdelta = mean((predm1_sed$sed^2)[upper.tri(predm1_sed$sed^2,diag = F)]) 

acc1 = sqrt(1-(PEV/sum1[1,1])) 
her1 = 1-(MVdelta/(2*summary(mod1)$varcomp[1,1])) 

#DIAG
mod2 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):idh(Years) + Plots,
              data = dataset, maxiter=100, 
              na.action = na.method(x="include", y = "include"))
mod2 = update(mod2)
sum2 = summary(mod2)$varcomp
aic2 = summary(mod2)$aic
logl2 = summary(mod2)$loglik

acc2 = NULL
her2 = NULL
for(i in 1:nlevels(dataset$Years)){
  predm2_vcov = predict(mod2, classify = "Hybrids:Years", 
                        level=list(Years = i), vcov = T)
  predm2_sed = predict(mod2, classify = "Hybrids:Years", 
                       level=list(Years = i), sed = T)
  
  predm2_sed$pvals    #EBLUPs
  
  PEV = mean(diag(predm2_vcov$vcov)) 
  MVdelta = mean((predm2_sed$sed^2)[upper.tri(predm2_sed$sed^2,diag = F)])
  
  acc2[i] = sqrt(1-(PEV/sum2[grep('Hybrids',rownames(sum2)),1][i]))
  her2[i] = 1-(MVdelta/(2*sum2[grep('Hybrids',rownames(sum2)),1][i]))
}
mean(acc2,na.rm = T)
mean(her2,na.rm = T)


#AR1 
mod3 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):ar1v(Years) + Plots,
              data = dataset, maxiter = 100, 
              na.action = na.method(x="include", y = "include"))
#mod3 = update(mod3)
sum3 = summary(mod3)$varcomp
aic3 = summary(mod3)$aic
logl3 = summary(mod3)$loglik

predm3_vcov = predict(mod3, classify = "Hybrids", vcov = T)
predm3_sed = predict(mod3, classify = "Hybrids", sed = T)

PEV = mean(diag(predm3_vcov$vcov)) 
MVdelta = mean((predm3_sed$sed^2)[upper.tri(predm3_sed$sed^2,diag = F)]) 

acc3 = sqrt(1-(PEV/sum3[grep('var',rownames(sum3)),1])) 
her3 = 1-(MVdelta/(2*sum3[grep('var',rownames(sum3)),1])) 


#AR1H
mod4 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):ar1h(Years) + Plots,
              data = dataset, maxiter = 100, 
              na.action = na.method(x="include", y = "include"))
mod4 = update(mod4)
sum4 = summary(mod4)$varcomp
aic4 = summary(mod4)$aic
logl4 = summary(mod4)$loglik

acc4 = NULL
her4 = NULL
for(i in 1:nlevels(dataset$Years)){
  predm4_vcov = predict(mod4, classify = "Hybrids:Years", 
                        level=list(Years = i), vcov = T)
  predm4_sed = predict(mod4, classify = "Hybrids:Years", 
                       level=list(Years = i), sed = T)
  
  PEV = mean(diag(predm4_vcov$vcov)) 
  MVdelta = mean((predm4_sed$sed^2)[upper.tri(predm4_sed$sed^2,diag = F)])
  
  acc4[i] = sqrt(1-(PEV/sum4[grep('Years_',rownames(sum4)),1][i]))
  her4[i] = 1-(MVdelta/(2*sum4[grep('Years_',rownames(sum4)),1][i]))
}
mean(acc4,na.rm = T)
mean(her4,na.rm = T)

#CORH
mod5 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):corh(Years) + Plots,
              data = dataset,maxiter = 100, 
              na.action = na.method(x="include", y = "include"))

sum5 = summary(mod5)$varcomp
aic5 = summary(mod5)$aic
logl5 = summary(mod5)$loglik

acc5 = NULL
her5 = NULL
for(i in 1:nlevels(dataset$Years)){
  predm5_vcov = predict(mod5, classify = "Hybrids:Years", 
                        level=list(Years = i), vcov = T)
  predm5_sed = predict(mod5, classify = "Hybrids:Years", 
                       level=list(Years = i), sed = T)
  
  PEV = mean(diag(predm5_vcov$vcov)) 
  MVdelta = mean((predm5_sed$sed^2)[upper.tri(predm5_sed$sed^2,diag = F)])
  
  acc5[i] = sqrt(1-(PEV/sum5[grep('Years_',rownames(sum5)),1][i]))
  her5[i] = 1-(MVdelta/(2*sum5[grep('Years_',rownames(sum5)),1][i]))
}
mean(acc5,na.rm = T)
mean(her5,na.rm = T)

#US - Did not converge
mod6 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):us(Years) + Plots,
              data = dataset, maxiter= 100, 
              na.action = na.method(x="include", y = "include"))
sum6 = summary(mod6)$varcomp
aic6 = summary(mod6)$aic
logl6 = summary(mod6)$loglik

#FA1
mod7 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):fa(Years,1) + Plots,
              data = dataset, maxiter= 100, 
              na.action = na.method(x="include", y = "include"))
sum7 = summary(mod7)$varcomp
aic7 = summary(mod7)$aic
logl7 = summary(mod7)$loglik

lambda = summary(mod7)$varcomp[grep('fa1',rownames(summary(mod7)$varcomp)),1]
psi = diag(summary(mod7)$varcomp[grep('var',rownames(summary(mod7)$varcomp)),1])
Gvcov = lambda %*% t(lambda) + psi
lambdacross = tcrossprod(lambda)
colnames(lambdacross) = rownames(lambdacross) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lambdacross[i,i] + lambdacross[j,j]) - lambdacross[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod7 = lambda_ASV/Gvcov_ASV

acc7 = NULL
her7 = NULL
for(i in 1:nlevels(dataset$Years)){
  predm7_vcov = predict(mod7, classify = "Hybrids:Years", 
                        level=list(Years = i), vcov = T)
  predm7_sed = predict(mod7, classify = "Hybrids:Years", 
                       level=list(Years = i), sed = T)
  
  predm7_sed$pvals    #EBLUPs
  
  PEV = mean(diag(predm7_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predm7_sed$sed^2)[upper.tri(predm7_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc7[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her7[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc7, na.rm = T)
mean(her7, na.rm = T)


#FA2
mod8 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):fa(Years,2)+ Plots,
              data = dataset, maxiter= 100)
mod8=update(mod8)
sum8 = summary(mod8)$varcomp
aic8 = summary(mod8)$aic
logl8 = summary(mod8)$loglik

fa1.loadings = summary(mod8)$varcomp[grep('fa1', rownames(summary(mod8)$varcomp)),1]
fa2.loadings = summary(mod8)$varcomp[grep('fa2', rownames(summary(mod8)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2')
psi = diag(summary(mod8)$varcomp[grep('var',rownames(summary(mod8)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)

Gvcov = lamblamb.star + psi

colnames(lamblamb.star) = rownames(lamblamb.star) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lamblamb.star[i,i] + lamblamb.star[j,j]) - lamblamb.star[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod8 = lambda_ASV/Gvcov_ASV

acc8 = NULL
her8 = NULL
for(i in 1:nlevels(dataset$Years)){
  predmod8_vcov = predict(mod8, classify = "Hybrids:Years", 
                          level=list(Years = i), vcov = T)
  predmod8_sed = predict(mod8, classify = "Hybrids:Years", 
                         level=list(Years = i), sed = T)
  
  PEV = mean(diag(predmod8_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod8_sed$sed^2)[upper.tri(predmod8_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc8[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her8[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc8, na.rm = T)
mean(her8, na.rm = T)

#FA3

mod9 = asreml(pf ~ Years + Replicates + Years:Replicates,
              random = ~vm(Hybrids,Amat):fa(Years,3)+ Plots,
              data = dataset, maxiter= 100, 
              na.action = na.method(x="include", y = "include"))
mod9=update(mod9)
sum9 = summary(mod9)$varcomp
aic9 = summary(mod9)$aic
logl9 = summary(mod9)$loglik

fa1.loadings = summary(mod9)$varcomp[grep('fa1', rownames(summary(mod9)$varcomp)),1]
fa2.loadings = summary(mod9)$varcomp[grep('fa2', rownames(summary(mod9)$varcomp)),1]
fa3.loadings = summary(mod9)$varcomp[grep('fa3', rownames(summary(mod9)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings, fa3.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2','fa3')
psi = diag(summary(mod9)$varcomp[grep('var',rownames(summary(mod9)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
Gvcov = lamblamb.star + psi
colnames(lamblamb.star) = rownames(lamblamb.star) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lamblamb.star[i,i] + lamblamb.star[j,j]) - lamblamb.star[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod9 = lambda_ASV/Gvcov_ASV

acc9 = NULL
her9 = NULL
for(i in 1:nlevels(dataset$Years)){
  predmod9_vcov = predict(mod9, classify = "Hybrids:Years", 
                          level=list(Years = i), vcov = T)
  predmod9_sed = predict(mod9, classify = "Hybrids:Years", 
                         level=list(Years = i), sed = T)
  
  PEV = mean(diag(predmod9_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod9_sed$sed^2)[upper.tri(predmod9_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc9[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her9[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc9, na.rm = T)
mean(her9, na.rm = T)

#FA4

mod10 = asreml(pf ~ Years + Replicates + Years:Replicates,
               random = ~vm(Hybrids,Amat):fa(Years,4)+ Plots,
               data = dataset, maxiter= 100)
mod10=update(mod10)
sum10 = summary(mod10)$varcomp
aic10 = summary(mod10)$aic
logl10 = summary(mod10)$loglik

fa1.loadings = summary(mod10)$varcomp[grep('fa1', rownames(summary(mod10)$varcomp)),1]
fa2.loadings = summary(mod10)$varcomp[grep('fa2', rownames(summary(mod10)$varcomp)),1]
fa3.loadings = summary(mod10)$varcomp[grep('fa3', rownames(summary(mod10)$varcomp)),1]
fa4.loadings = summary(mod10)$varcomp[grep('fa4', rownames(summary(mod10)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings, fa3.loadings, fa4.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2','fa3','fa4')
psi = diag(summary(mod10)$varcomp[grep('var',rownames(summary(mod10)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
Gvcov = lamblamb.star + psi

colnames(lamblamb.star) = rownames(lamblamb.star) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lamblamb.star[i,i] + lamblamb.star[j,j]) - lamblamb.star[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod10 = lambda_ASV/Gvcov_ASV

acc10 = NULL
her10 = NULL
for(i in 1:nlevels(dataset$Years)){
  predmod10_vcov = predict(mod10, classify = "Hybrids:Years", 
                           level=list(Years = i), vcov = T)
  predmod10_sed = predict(mod10, classify = "Hybrids:Years", 
                          level=list(Years = i), sed = T)
  
  PEV = mean(diag(predmod10_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod10_sed$sed^2)[upper.tri(predmod10_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc10[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her10[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc10, na.rm = T)
mean(her10, na.rm = T)

##### Modelling the residual effects ------

# IDH

mod11 = asreml(pf ~ Years + Replicates + Years:Replicates,
               random = ~vm(Hybrids,Amat):fa(Years,3)+ Plots,
               residual = ~dsum(~id(units)|Years),
               data = dataset, maxiter= 100, 
              na.action = na.method(x="include", y = "include"))
mod11=update(mod11)
sum11 = summary(mod11)$varcomp
aic11 = summary(mod11)$aic
logl11 = summary(mod11)$loglik

fa1.loadings = summary(mod11)$varcomp[grep('fa1', rownames(summary(mod11)$varcomp)),1]
fa2.loadings = summary(mod11)$varcomp[grep('fa2', rownames(summary(mod11)$varcomp)),1]
fa3.loadings = summary(mod11)$varcomp[grep('fa3', rownames(summary(mod11)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings, fa3.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2','fa3')
psi = diag(summary(mod11)$varcomp[grep('var',rownames(summary(mod11)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)

Gvcov = lamblamb.star + psi

colnames(lamblamb.star) = rownames(lamblamb.star) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lamblamb.star[i,i] + lamblamb.star[j,j]) - lamblamb.star[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod11 = lambda_ASV/Gvcov_ASV

acc11 = NULL
her11 = NULL
for(i in 1:nlevels(dataset$Years)){
  predmod11_vcov = predict(mod11, classify = "Hybrids:Years", 
                           level=list(Years = i), vcov = T)
  predmod11_sed = predict(mod11, classify = "Hybrids:Years", 
                          level=list(Years = i), sed = T)
  
  PEV = mean(diag(predmod11_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod11_sed$sed^2)[upper.tri(predmod11_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc11[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her11[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc11, na.rm = T)
mean(her11, na.rm = T)

# AR1

mod12 = asreml(pf ~ Years + Replicates + Years:Replicates,
               random = ~vm(Hybrids,Amat):fa(Years,3)+ Plots,
               residual = ~ar1(Years):Plots,
               data = dataset, maxiter= 120, 
              na.action = na.method(x="include", y = "include"))
mod12=update(mod12)
sum12 = summary(mod12)$varcomp
aic12 = summary(mod12)$aic
logl12 = summary(mod12)$loglik

fa1.loadings = summary(mod12)$varcomp[grep('fa1', rownames(summary(mod12)$varcomp)),1]
fa2.loadings = summary(mod12)$varcomp[grep('fa2', rownames(summary(mod12)$varcomp)),1]
fa3.loadings = summary(mod12)$varcomp[grep('fa3', rownames(summary(mod12)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings, fa3.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2','fa3')
psi = diag(summary(mod12)$varcomp[grep('var',rownames(summary(mod12)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)

Gvcov = lamblamb.star + psi

colnames(lamblamb.star) = rownames(lamblamb.star) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lamblamb.star[i,i] + lamblamb.star[j,j]) - lamblamb.star[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod12 = lambda_ASV/Gvcov_ASV

acc12 = NULL
her12 = NULL
for(i in 1:nlevels(dataset$Years)){
  predmod12_vcov = predict(mod12, classify = "Hybrids:Years", 
                           level=list(Years = i), vcov = T)
  predmod12_sed = predict(mod12, classify = "Hybrids:Years", 
                          level=list(Years = i), sed = T)
  
  PEV = mean(diag(predmod12_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod12_sed$sed^2)[upper.tri(predmod12_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc12[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her12[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc12, na.rm = T)
mean(her12, na.rm = T)

# AR1H

mod13 = asreml(pf ~ Years + Replicates + Years:Replicates,
               random = ~vm(Hybrids,Amat):fa(Years,3)+ Plots,
               residual = ~ar1h(Years):Plots,
               data = dataset, maxiter= 100)
mod13=update(mod13)
sum13 = summary(mod13)$varcomp
aic13 = summary(mod13)$aic
logl13 = summary(mod13)$loglik

fa1.loadings = summary(mod13)$varcomp[grep('fa1', rownames(summary(mod13)$varcomp)),1]
fa2.loadings = summary(mod13)$varcomp[grep('fa2', rownames(summary(mod13)$varcomp)),1]
fa3.loadings = summary(mod13)$varcomp[grep('fa3', rownames(summary(mod13)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings, fa3.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2','fa3')
psi = diag(summary(mod13)$varcomp[grep('var',rownames(summary(mod13)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)

Gvcov = lamblamb.star + psi
corr = cov2cor(Gvcov)
colnames(lamblamb.star) = rownames(lamblamb.star) = levels(dataset$Years)
colnames(Gvcov) = rownames(Gvcov) = levels(dataset$Years)

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(lamblamb.star[i,i] + lamblamb.star[j,j]) - lamblamb.star[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

fullsv = list()
for (i in levels(dataset$Years)) {
  semivar = matrix(nrow = nlevels(dataset$Years), ncol = 3,
                   dimnames = list(levels(dataset$Years), c('i','j','semivar')))
  for (j in levels(dataset$Years)) {
    semivar[j,] = c(i,j,
                    0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  }
  fullsv[[i]] = semivar
}

fullsv = as.data.frame(do.call(rbind,fullsv)) %>% 
  mutate(semivar = as.numeric(semivar)) %>% 
  pivot_wider(names_from = i, values_from = semivar) %>% 
  column_to_rownames(var = 'j')

Gvcov_ASV = 2/(nrow(fullsv) * (nrow(fullsv) - 1)) * sum(fullsv[upper.tri(fullsv)])

expvar_mod13 = lambda_ASV/Gvcov_ASV 
corr = cov2cor(Gvcov)

acc13 = NULL
her13 = NULL
for(i in 1:nlevels(dataset$Years)){
  predmod13_vcov = predict(mod13, classify = "Hybrids:Years", 
                           level=list(Years = i), vcov = T)
  predmod13_sed = predict(mod13, classify = "Hybrids:Years", 
                          level=list(Years = i), sed = T)
  
  PEV = mean(diag(predmod13_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod13_sed$sed^2)[upper.tri(predmod13_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc13[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her13[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc13, na.rm = T)
mean(her13, na.rm = T)

fa1.scores = mod13$coefficients$random[grep("Comp1",row.names(mod13$coefficients$random)),1];
names(fa1.scores) = sub("vm(Hybrids, Amat)_","",names(fa1.scores),fixed=T)
names(fa1.scores) = sub(":fa(Years, 3)_Comp1","",names(fa1.scores),fixed=T)
fa2.scores = mod13$coefficients$random[grep("Comp2",row.names(mod13$coefficients$random)),1];
names(fa2.scores) = sub("vm(Hybrids, Amat)_","",names(fa2.scores),fixed=T)
names(fa2.scores) = sub(":fa(Years, 2)_Comp1","",names(fa2.scores),fixed=T)
fa3.scores = mod13$coefficients$random[grep("Comp3",row.names(mod13$coefficients$random)),1];
names(fa3.scores) = sub("vm(Hybrids, Amat)_","",names(fa3.scores),fixed=T)
names(fa3.scores) = sub(":fa(Years, 3)_Comp3","",names(fa3.scores),fixed=T)

fa.scores = rbind(as.matrix(fa1.scores),as.matrix(fa2.scores),as.matrix(fa3.scores))

fa.scores.star = -kronecker(t(svd.mat.loadings$v), diag(nlevels(factor(pedigree$gen))))%*%fa.scores
rownames(fa.scores.star) = rownames(fa.scores)

fa1.scores.star = fa.scores.star[1:nlevels(factor(pedigree$gen)),1]
fa2.scores.star = fa.scores.star[(nlevels(factor(pedigree$gen))+1):(nlevels(factor(pedigree$gen))*2),1]
fa3.scores.star = fa.scores.star[((nlevels(factor(pedigree$gen))*2)+1):(nlevels(factor(pedigree$gen))*3),1]

EBLUPs_marg = (kronecker(mat.loadings.star,diag(nlevels(factor(pedigree$gen))))) %*% fa.scores.star
EBLUPs_marg = data.frame("Year" = rep(levels(dataset$Years),each = nlevels(factor(pedigree$gen))),
                         "Hyb" = rep(levels(factor(pedigree$gen)),nlevels(dataset$Years)),
                         "EBLUP_marg" = EBLUPs_marg)
EBLUPs_marg$FL_1 = rep(mat.loadings.star[,1], each = nlevels(factor(pedigree$gen)))

OP = mean(mat.loadings.star[,1]) * as.matrix(fa1.scores.star)
OP = OP[order(rownames(OP, vector)),,drop = FALSE]
rest = (EBLUPs_marg$EBLUP_marg -
          (kronecker(mat.loadings.star[,1],diag(nlevels(factor(pedigree$gen))))) %*% fa1.scores.star)^2
STA = data.frame("gen" = rep(levels(factor(pedigree$gen)), nlevels(dataset$Years)),
                 "ST" = rest)
STA = STA %>% group_by(gen) %>% summarise(ST = sqrt(mean(ST))) %>% arrange(gen)
plot1 = data.frame("gen"=levels(factor(pedigree$gen)),
                   "OP" = OP,
                   "ST" = STA$ST,
                   'ind' = ((OP - mean(OP))/sd(OP))-((STA$ST - mean(STA$ST))/sd(STA$ST)))

ggplot(data=plot1)+
  geom_point(aes(x = ST, y = OP, color=gen),size=2)+
  geom_vline(xintercept = 0,linetype="dashed", colour = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Overall Performance") + xlab("Root Mean Square Deviation")+
  gghighlight(max(ST)<=1.8,max(OP)>=.0,use_direct_label = F)+
  geom_label_repel(aes(x = ST, y = OP, label = gen), size = 4,
                   box.padding = 1)+
  theme(legend.position = 'none')
#ggsave('OPST.svg',dpi = 1000, width = 5, height = 5)


### Table and Figures ------

# Model selection

modsel = data.frame('mod' = paste0('M',seq(1:13)),
                    'Rstr' = c('IDV','IDV','IDV','IDV','IDV','IDV','IDV','IDV',
                               'IDV','IDV','DIAG', 'AR1', 'AR1H'),
                    'Gstr' = c('CS','DIAG','AR1','AR1H','CSH','US','FA1',
                               'FA2','FA3','FA4','FA3','FA3','FA3'),
                    'aic' = c(aic1,aic2,aic3,aic4,aic5,NA,aic7,aic8,aic9,
                              aic10,aic11,aic12,aic13),
                    'expvar' = c(NA,NA,NA,NA,NA,NA,expvar_mod7, expvar_mod8,
                                 expvar_mod9,expvar_mod10,expvar_mod11,expvar_mod12,
                                 expvar_mod13)*100,
                    'loglik' = c(logl1,logl2,logl3,logl4,logl5,NA,logl7,logl8,logl9,
                                 logl10,logl11,logl12,logl13),
                    'accuracy' = c(acc1,mean(acc2,na.rm=T),acc3,
                                   mean(acc4,na.rm=T),mean(acc5,na.rm=T),NA,
                                   mean(acc7,na.rm=T),mean(acc8),mean(acc9),
                                   mean(acc10),mean(acc11),mean(acc12),mean(acc13)))

# Genetic correlation heatmap

col_fun = colorRampPalette(c('red','green'))

rownames(corr) = colnames(corr) = levels(dataset$Years)

Heatmap(corr,col=col_fun(5), rect_gp = gpar(col = "white", lwd = 1, type = "none"),
        column_dend_height = unit(2, "cm"), 
        row_dend_width = unit(2, "cm"),
        heatmap_legend_param = list(title="Correlation",at=c(-1,0,1),labels=c("-1","0","1")),
        row_order = sort(rownames(corr)), show_row_names = F,
        column_order = sort(colnames(corr)),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j)  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          if(i >= j) grid.text(sprintf("%.2f", corr[i, j]), x, y, gp = gpar(fontsize = 15))})


# Kinship heatmap

Heatmap(
  Amat, col = col_fun(20), show_column_names = T, show_row_names = T,
  heatmap_legend_param = list(title = ""), 
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6)
)


# Ranking comparison

a = coef(mod1)$random[grep('vm', rownames(coef(mod1)$random)),][1:nlevels(factor(pedigree$gen))]
names(a) = sub("vm(Hybrids, Amat)_","", names(a), fixed= T)
a = rownames_to_column(as.data.frame(a))
colnames(a) = c("Hyb", "value"); a = arrange(a, Hyb)
b = plot1 %>% group_by(gen) %>% summarise(ind) %>% 
  arrange(gen) %>% rename(Hyb = gen, value = ind)
cor(a$value, b$value, method = 'spearman')  #Spearman correlation = 0.49

blups2 = rbind(a,b)
blups2$Model = rep(c('BLUP M1', 'OPST M13'), each = nrow(blups2)/2)

ggplot(data = blups2, aes(x = Model, y=value, group = Hyb, color=Hyb))+
  geom_line(alpha = 1/2, size = 2)+
  geom_point(alpha = 1/1.5, size = 2) + xlab(" ") + ylab('Breeding values')+
  theme(legend.title = element_blank(), legend.position = 'none',
        axis.text.y = element_text(size =12),axis.text.x = element_text(size =13),
        axis.title.y = element_text(size = 13))+
  labs(caption = paste('Ranking correlation =',
                       round(cor(a$value, b$value, method = 'spearman'),3)))


d = coef(mod11)$random[grep('_Y',rownames(coef(mod11)$random)),]
names(d) = sub("vm(Hybrids, Amat)_","",names(d),fixed=T)
names(d) = sub("fa(Years, 3)_","",names(d),fixed=T)

d = rownames_to_column(as.data.frame(d))
d = d %>% separate(rowname, c('Hyb', 'Yr'), sep = ':') %>% group_by(Hyb) %>% 
  summarise(blup = mean(d))


mean(d[d$Hyb %in% c('H140','H117','H143','H129','286x215','434','186x434',
                    'H163','513','215'),]$blup) / mean(dataset$pf,na.rm=T)  *100

mean(d[d$Hyb %in% c('H165','H162','H169','H140','1074','H123','H117',
                    '215','H143','H144'),]$blup) / mean(dataset$pf,na.rm=T) *100

# Relative gain

pm1 = predict(mod1, classify = 'Hybrids')$pvals

pm11 = predict(mod11, classify = 'Hybrids')$pvals

Gm1 = (mean((pm1[order(pm1$predicted.value, decreasing = T),][1:10,])$predicted.value) -
         (mean(pm1$predicted.value)))/(mean(pm1$predicted.value)) * 100

Gm11 = (mean((pm11[order(pm1$predicted.value, decreasing = T),][1:10,])$predicted.value) -
          (mean(pm11$predicted.value)))/(mean(pm11$predicted.value)) * 100

Gind = (mean((pm11[pm11$Hybrids %in% c('215', '286x215', '186x434', '513', 'H129',
                                       'H117', 'H125', 'H143', 'H120', 'H124'),])$predicted.value) -
          (mean(pm11$predicted.value)))/(mean(pm11$predicted.value)) * 100






