rm(list = ls())
# Analysis of multi-harvest data through mixed models: 
# an application in Theobroma grandiflorum breeding

# Crop Science

##>> Loading the packages
require(asreml)
require(tidyverse)
require(gghighlight)
require(ggrepel)
require(ComplexHeatmap)
require(nadiv)
require(patchwork)

##>> Import dataset
dataset = read.table("Data/dataset.txt", header = TRUE)
dataset = dataset %>% group_by(Harvests, Plots, Replicates, Hybrids) %>% 
  summarise(pf = mean(pf,na.rm=T))
pedigree = read.table('Data/pedigree.txt', header= T)

##>> Transform numeric variables into factors
dataset <- transform(dataset, Harvests = factor(Harvests), 
                     Plots = factor(Plots), 
                     Replicates = factor(Replicates), 
                     Hybrids = factor(Hybrids),
                     Hybrids2 = factor(Hybrids))

##>> Relationship matrices
Amat = as.matrix(makeA(pedigree = pedigree))


##### Modelling the additive effects ------

#CS
mod1 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat)+ vm(Hybrids,Amat):idv(Harvests) + Plots,
              data = dataset, maxiter=100)
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
mod2 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):idh(Harvests) + Plots,
              data = dataset, maxiter=100)
mod2 = update(mod2)
sum2 = summary(mod2)$varcomp
aic2 = summary(mod2)$aic
logl2 = summary(mod2)$loglik

acc2 = NULL
her2 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predm2_vcov = predict(mod2, classify = "Hybrids:Harvests", 
                         level=list(Harvests = i), vcov = T)
  predm2_sed = predict(mod2, classify = "Hybrids:Harvests", 
                        level=list(Harvests = i), sed = T)
  
  predm2_sed$pvals    #EBLUPs
  
  PEV = mean(diag(predm2_vcov$vcov)) 
  MVdelta = mean((predm2_sed$sed^2)[upper.tri(predm2_sed$sed^2,diag = F)])
  
  acc2[i] = sqrt(1-(PEV/sum2[grep('Hybrids',rownames(sum2)),1][i]))
  her2[i] = 1-(MVdelta/(2*sum2[grep('Hybrids',rownames(sum2)),1][i]))
}
mean(acc2,na.rm = T)
mean(her2,na.rm = T)


#AR1 
mod3 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):ar1v(Harvests) + Plots,
              data = dataset, maxiter = 100)
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
mod4 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):ar1h(Harvests) + Plots,
              data = dataset, maxiter = 100)
mod4 = update(mod4)
sum4 = summary(mod4)$varcomp
aic4 = summary(mod4)$aic
logl4 = summary(mod4)$loglik

acc4 = NULL
her4 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predm4_vcov = predict(mod4, classify = "Hybrids:Harvests", 
                        level=list(Harvests = i), vcov = T)
  predm4_sed = predict(mod4, classify = "Hybrids:Harvests", 
                       level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predm4_vcov$vcov)) 
  MVdelta = mean((predm4_sed$sed^2)[upper.tri(predm4_sed$sed^2,diag = F)])
  
  acc4[i] = sqrt(1-(PEV/sum4[grep('Harvests_',rownames(sum4)),1][i]))
  her4[i] = 1-(MVdelta/(2*sum4[grep('Harvests_',rownames(sum4)),1][i]))
}
mean(acc4,na.rm = T)
mean(her4,na.rm = T)

#CORH
mod5 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):corh(Harvests) + Plots,
              data = dataset,maxiter = 100)

sum5 = summary(mod5)$varcomp
aic5 = summary(mod5)$aic
logl5 = summary(mod5)$loglik

acc5 = NULL
her5 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predm5_vcov = predict(mod5, classify = "Hybrids:Harvests", 
                        level=list(Harvests = i), vcov = T)
  predm5_sed = predict(mod5, classify = "Hybrids:Harvests", 
                       level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predm5_vcov$vcov)) 
  MVdelta = mean((predm5_sed$sed^2)[upper.tri(predm5_sed$sed^2,diag = F)])
  
  acc5[i] = sqrt(1-(PEV/sum5[grep('Harvests_',rownames(sum5)),1][i]))
  her5[i] = 1-(MVdelta/(2*sum5[grep('Harvests_',rownames(sum5)),1][i]))
}
mean(acc5,na.rm = T)
mean(her5,na.rm = T)

#US - Did not converge
mod6 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):us(Harvests) + Plots,
              data = dataset, maxiter= 100)
sum6 = summary(mod6)$varcomp
aic6 = summary(mod6)$aic
logl6 = summary(mod6)$loglik

#FA1
mod7 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):fa(Harvests,1) + Plots,
              data = dataset, maxiter= 100)
sum7 = summary(mod7)$varcomp
aic7 = summary(mod7)$aic
logl7 = summary(mod7)$loglik

lambda = summary(mod7)$varcomp[grep('fa1',rownames(summary(mod7)$varcomp)),1]
psi = diag(summary(mod7)$varcomp[grep('var',rownames(summary(mod7)$varcomp)),1])

Gvcov = lambda %*% t(lambda) + psi
expvar1 = (sum(diag(lambda %*% t(lambda)))/sum(diag(lambda %*% t(lambda) + psi)))*100
#Variância explicada: 70.81%

acc7 = NULL
her7 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predm7_vcov = predict(mod7, classify = "Hybrids:Harvests", 
                         level=list(Harvests = i), vcov = T)
  predm7_sed = predict(mod7, classify = "Hybrids:Harvests", 
                        level=list(Harvests = i), sed = T)
  
  predm7_sed$pvals    #EBLUPs
  
  PEV = mean(diag(predm7_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predm7_sed$sed^2)[upper.tri(predm7_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc7[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her7[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc7, na.rm = T)
mean(her7, na.rm = T)


#FA2
mod8 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):fa(Harvests,2)+ Plots,
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
expvar2.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + psi)))*100
#Variância explicada: 83.14065

Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi

acc8 = NULL
her8 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predmod8_vcov = predict(mod8, classify = "Hybrids:Harvests", 
                         level=list(Harvests = i), vcov = T)
  predmod8_sed = predict(mod8, classify = "Hybrids:Harvests", 
                        level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predmod8_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod8_sed$sed^2)[upper.tri(predmod8_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc8[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her8[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc8, na.rm = T)
mean(her8, na.rm = T)

##### Modelling the residual effects ------

# IDH

mod9 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):fa(Harvests,2)+ Plots,
              residual = ~dsum(~id(units)|Harvests),
              data = dataset, maxiter= 100)
mod9=update(mod9)
sum9 = summary(mod9)$varcomp
aic9 = summary(mod9)$aic
logl9 = summary(mod9)$loglik

fa1.loadings = summary(mod9)$varcomp[grep('fa1', rownames(summary(mod9)$varcomp)),1]
fa2.loadings = summary(mod9)$varcomp[grep('fa2', rownames(summary(mod9)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2')
psi = diag(summary(mod9)$varcomp[grep('var',rownames(summary(mod9)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar2.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + psi)))*100

Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi

acc9 = NULL
her9 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predmod9_vcov = predict(mod9, classify = "Hybrids:Harvests", 
                          level=list(Harvests = i), vcov = T)
  predmod9_sed = predict(mod9, classify = "Hybrids:Harvests", 
                         level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predmod9_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod9_sed$sed^2)[upper.tri(predmod9_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc9[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her9[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc9, na.rm = T)
mean(her9, na.rm = T)


# AR1

mod10 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):fa(Harvests,2)+ Plots,
              residual = ~ar1(Harvests):Plots,
              data = dataset, maxiter= 100)
mod10=update(mod10)
sum10 = summary(mod10)$varcomp
aic10 = summary(mod10)$aic
logl10 = summary(mod10)$loglik

fa1.loadings = summary(mod10)$varcomp[grep('fa1', rownames(summary(mod10)$varcomp)),1]
fa2.loadings = summary(mod10)$varcomp[grep('fa2', rownames(summary(mod10)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2')
psi = diag(summary(mod10)$varcomp[grep('var',rownames(summary(mod10)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar2.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + psi)))*100

Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi

acc10 = NULL
her10 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predmod10_vcov = predict(mod10, classify = "Hybrids:Harvests", 
                          level=list(Harvests = i), vcov = T)
  predmod10_sed = predict(mod10, classify = "Hybrids:Harvests", 
                         level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predmod10_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod10_sed$sed^2)[upper.tri(predmod10_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc10[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her10[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc10, na.rm = T)
mean(her10, na.rm = T)


# AR1H

mod11 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
              random = ~vm(Hybrids,Amat):fa(Harvests,2)+ Plots,
              residual = ~ar1h(Harvests):Plots,
              data = dataset, maxiter= 100)
mod11=update(mod11)
sum11 = summary(mod11)$varcomp
aic11 = summary(mod11)$aic
logl11 = summary(mod11)$loglik

fa1.loadings = summary(mod11)$varcomp[grep('fa1', rownames(summary(mod11)$varcomp)),1]
fa2.loadings = summary(mod11)$varcomp[grep('fa2', rownames(summary(mod11)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2')
psi = diag(summary(mod11)$varcomp[grep('var',rownames(summary(mod11)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar2.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + psi)))*100

Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi
corr = cov2cor(Gvcov)

acc11 = NULL
her11 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predmod11_vcov = predict(mod11, classify = "Hybrids:Harvests", 
                          level=list(Harvests = i), vcov = T)
  predmod11_sed = predict(mod11, classify = "Hybrids:Harvests", 
                         level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predmod11_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod11_sed$sed^2)[upper.tri(predmod11_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc11[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her11[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc11, na.rm = T)
mean(her11, na.rm = T)

fa1.scores = mod11$coefficients$random[grep("Comp1",row.names(mod11$coefficients$random)),1];
names(fa1.scores) = sub("vm(Hybrids, Amat)_","",names(fa1.scores),fixed=T)
names(fa1.scores) = sub(":fa(Harvests, 2)_Comp1","",names(fa1.scores),fixed=T)
fa2.scores = mod11$coefficients$random[grep("Comp2",row.names(mod11$coefficients$random)),1];
names(fa2.scores) = sub("vm(Hybrids, Amat)_","",names(fa2.scores),fixed=T)
names(fa2.scores) = sub(":fa(Harvests, 2)_Comp1","",names(fa1.scores),fixed=T)

fa.scores = rbind(as.matrix(fa1.scores),as.matrix(fa2.scores))

fa.scores.star = -kronecker(t(svd.mat.loadings$v), diag(nlevels(factor(pedigree$gen))))%*%fa.scores
rownames(fa.scores.star) = rownames(fa.scores)

fa1.scores.star = fa.scores.star[1:nlevels(factor(pedigree$gen)),1]
fa2.scores.star = fa.scores.star[(nlevels(factor(pedigree$gen))+1):(nlevels(factor(pedigree$gen))*2),1]

EBLUPs_marg = (kronecker(mat.loadings.star,diag(nlevels(factor(pedigree$gen))))) %*% fa.scores.star
EBLUPs_marg = data.frame("Harvest" = rep(levels(dataset$Harvests),each = nlevels(factor(pedigree$gen))),
                         "Hyb" = rep(levels(factor(pedigree$gen)),nlevels(dataset$Harvests)),
                         "EBLUP_marg" = EBLUPs_marg)
EBLUPs_marg$FL_1 = rep(mat.loadings.star[,1], each = nlevels(factor(pedigree$gen)))

OP = mean(mat.loadings.star[,1]) * as.matrix(fa1.scores.star)
OP = OP[order(rownames(OP, vector)),,drop = FALSE]
rest = (EBLUPs_marg$EBLUP_marg -
          (kronecker(mat.loadings.star[,1],diag(nlevels(factor(pedigree$gen))))) %*% fa1.scores.star)^2
STA = data.frame("gen" = rep(levels(factor(pedigree$gen)), nlevels(dataset$Harvests)),
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
  ylab("Overall Performance") + xlab("Stability")+
  gghighlight(max(ST)<=1.5,max(OP)>=.3,use_direct_label = F)+
  geom_label_repel(aes(x = ST, y = OP, label = gen), size = 4,
                   box.padding = 1)+
  theme(legend.position = 'none')
#ggsave('OPST.svg',dpi = 1000, width = 5, height = 5)


# US

mod12 = asreml(pf ~ Harvests + Replicates + Harvests:Replicates,
               random = ~vm(Hybrids,Amat):fa(Harvests,2) + Plots,
               residual = ~us(Harvests):Plots,
               data = dataset, maxiter= 100)
mod12=update(mod12)
sum12 = summary(mod12)$varcomp
aic12 = summary(mod12)$aic
logl12 = summary(mod12)$loglik

fa1.loadings = summary(mod12)$varcomp[grep('fa1', rownames(summary(mod12)$varcomp)),1]
fa2.loadings = summary(mod12)$varcomp[grep('fa2', rownames(summary(mod12)$varcomp)),1]
mat.loadings = as.matrix(cbind(fa1.loadings,fa2.loadings))
svd.mat.loadings=svd(mat.loadings) 
mat.loadings.star = mat.loadings%*%svd.mat.loadings$v *-1
colnames(mat.loadings.star) = c("fa1",'fa2')
psi = diag(summary(mod12)$varcomp[grep('var',rownames(summary(mod12)$varcomp)),1])
lamblamb.star = mat.loadings.star %*% t(mat.loadings.star)
expvar2.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + psi)))*100

Gvcov = mat.loadings.star %*% t(mat.loadings.star) + psi

acc12 = NULL
her12 = NULL
for(i in 1:nlevels(dataset$Harvests)){
  predmod12_vcov = predict(mod12, classify = "Hybrids:Harvests", 
                           level=list(Harvests = i), vcov = T)
  predmod12_sed = predict(mod12, classify = "Hybrids:Harvests", 
                          level=list(Harvests = i), sed = T)
  
  PEV = mean(diag(predmod12_vcov$vcov)) # Prediction Error Variance
  MVdelta = mean((predmod12_sed$sed^2)[upper.tri(predmod11_sed$sed^2,diag = F)]) #Mean pairwise prediction error variance
  
  acc12[i] = sqrt(1-(PEV/Gvcov[i,i])) # Accuracy
  her12[i] = 1-(MVdelta/(2*Gvcov[i,i])) # Heritability
}

mean(acc12, na.rm = T)
mean(her12, na.rm = T)




### Tables and Figures ------

# Model selection

modsel = data.frame('mod' = paste('model',seq(1:12)),
                    'Rstr' = c('IDV','IDV','IDV','IDV','IDV','IDV','IDV','IDV',
                               'DIAG', 'AR1', 'AR1H', 'US'),
                    'Gstr' = c('CS','DIAG','AR1','AR1H','CSH','US','FA1',
                               'FA2','FA2','FA2','FA2','FA2'),
                    'aic' = c(aic1,aic2,aic3,aic4,aic5,NA,aic7,aic8,aic9,
                              aic10,aic11,NA),
                    'npar'= c(attr(aic1,"parameters"),attr(aic2,"parameters"),
                              attr(aic3,"parameters"),attr(aic4,"parameters"),
                              attr(aic5,"parameters"),NA,
                              attr(aic7,"parameters"),attr(aic8,"parameters"),
                              attr(aic9,"parameters"),attr(aic10,"parameters"),
                              attr(aic11,"parameters"),NA),
                    'loglik' = c(logl1,logl2,logl3,logl4,logl5,NA,logl7,logl8,logl9,
                                 logl10,logl11,NA),
                    'accuracy' = c(acc1,mean(acc2,na.rm=T),acc3,
                                   mean(acc4,na.rm=T),mean(acc5,na.rm=T),NA,
                                   mean(acc7,na.rm=T),mean(acc8),mean(acc9),
                                   mean(acc10),mean(acc11),NA))


# Heatmap - genetic correlation

col_fun = colorRampPalette(c('yellow','red'))

rownames(corr) = colnames(corr) = levels(dataset$Harvests)

Heatmap(corr,col=col_fun(20), rect_gp = gpar(col = "white", lwd = 1, type = "none"),
        column_dend_height = unit(2, "cm"), 
        row_dend_width = unit(2, "cm"),
        heatmap_legend_param = list(title="Correlation",at=c(-1,0,1),labels=c("-1","0","1")),
        row_order = sort(rownames(corr)), show_row_names = F,
        column_order = sort(colnames(corr)),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j)  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          if(i >= j) grid.text(sprintf("%.2f", corr[i, j]), x, y, gp = gpar(fontsize = 15))})


# Heatmap - kinship

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
cor(a$value, b$value, method = 'spearman')  #Spearman correlation = 0.5146

blups2 = rbind(a,b)
blups2$Model = rep(c('BLUP M1', 'OPST M11'), each = nrow(blups2)/2)

ggplot(data = blups2, aes(x = Model, y=value, group = Hyb, color=Hyb))+
  geom_line(alpha = 1/2, size = 2)+
  geom_point(alpha = 1/1.5, size = 2) + xlab(" ") + ylab('Breeding values')+
  theme(legend.title = element_blank(), legend.position = 'right',
        axis.text.y = element_text(size =12),axis.text.x = element_text(size =13),
        axis.title.y = element_text(size = 13))


d = coef(mod11)$random[grep('_Ha',rownames(coef(mod11)$random)),]
names(d ) = sub("vm(Hybrids, Amat)_","",names(d),fixed=T)
names(d) = sub("fa(Harvests, 2)_","",names(d),fixed=T)

d = rownames_to_column(as.data.frame(d))
d = d %>% separate(rowname, c('Hyb', 'Harv'), sep = ':') %>% group_by(Hyb) %>% 
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

Gind = (mean((pm11[pm11$Hybrids %in% c('H140','H117','H143','H129','286x215','434','186x434',
                         'H163','513','215'),])$predicted.value) -
    (mean(pm11$predicted.value)))/(mean(pm11$predicted.value)) * 100

1- Gm1 / Gm11

# Latent regression (optional) - Top 5 and Bottom 5

num.gen = nlevels(predmod11_sed$pvals$Hybrids)

EBLUPs_marg$EBLUPs_marg_f2 = EBLUPs_marg$EBLUP_marg - 
  (kronecker(mat.loadings.star[,1],diag(num.gen))) %*% as.matrix(fa1.scores.star)
EBLUPs_marg$FL_2 = rep(mat.loadings.star[,2], each = num.gen)


reg1 = ggplot(subset(EBLUPs_marg, Hyb %in% c('H140','215','H117','H129',"H143",
                                      'H123','174x186','H166','186x1074',"H152")))+
  geom_point(aes(x=FL_1,y=EBLUP_marg,color = Hyb))+
  geom_smooth(aes(x=FL_1,y=EBLUP_marg, color = Hyb),method=lm, fill = NA)+
  facet_wrap(~factor(Hyb, levels = c('H140','215','H117','H129',"H143",
                                     'H123','174x186','H166','186x1074',"H152")),
             ncol = 5, as.table = T, dir = 'h', scales = 'free_x')+
  xlab("Loadings for first factor")+ylab("EBLUP")+
  labs(color = "Genotypes", tag = '(A)')+
  theme(legend.position = "none")

reg2 = ggplot(subset(EBLUPs_marg, Hyb %in% c('H140','215','H117','H129',"H143",
                                      'H123','174x186','H166','186x1074',"H152")))+
  geom_point(aes(x=FL_2,y=EBLUPs_marg_f2,color = Hyb))+
  geom_smooth(aes(x=FL_2,y=EBLUPs_marg_f2, color = Hyb),method=lm, fill = NA)+
  facet_wrap(~factor(Hyb, levels = c('H140','215','H117','H129',"H143",
                                     'H123','174x186','H166','186x1074',"H152")),
             ncol = 5, as.table = T, dir = 'h', scales = 'free_x')+
  xlab("Loadings for second factor")+
  ylab(expression(EBLUP[m] - lambda[1]^'*' * f[1]^'*'))+
  labs(color = "Genotypes", tag = "(B)")+
  theme(legend.position = "none")


reg1/reg2






