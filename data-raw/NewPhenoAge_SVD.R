#Make sure you have your DNAm data in a matrix with columns as CpGs and rows as observations
library(glmnet)
CV= cv.glmnet(datMeth_InCHIANTI,Pheno_InCHIANTI$PredAge, nfolds=10,alpha=0.5, family="gaussian")
fit = glmnet(datMeth_InCHIANTI,Pheno_InCHIANTI$PredAge, family="gaussian", alpha=0.5, nlambda=100)
lambda.Min = CV$lambda.min

#This is plotting the error as a function of the lambda value. It will also show how many CpGs are selected as a function of lambda
plot(CV, main="CV Elastic Net")

#This Predicts the Age in the Training data
DNAmPhenoAge_Train=as.numeric(predict(fit,datMeth_InCHIANTI,type="response",s=lambda.Min))
DNAmPhenoAge_Test=as.numeric(predict(fit,datMeth_WHI,type="response",s=lambda.Min))

#Plot the predicted Age versus the DNAmPhenoAge
verboseScatterplot(Pheno_WHI$agewhi,DNAmPhenoAge_Test,xlab="Age",ylab="New DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")
verboseScatterplot(Pheno_WHI$agewhi,Pheno_WHI$DNAmPhenoAge,xlab="Age",ylab="Original DNAmPhenoAge",main="Validation Set")
abline(0,1,col="red")

#Plot New vs. Original DNAmPhenoAge
verboseScatterplot(DNAmPhenoAge_Test,Pheno_WHI$DNAmPhenoAge,xlab="New DNAmPhenoAge",ylab="Original DNAmPhenoAge",main="Training Set")
abline(0,1,col="red")

library(survival)
coxph(Surv(Data$ENDFOLLOWALLDY, Data$DEATHALL) ~ DNAmPhenoAge_Test+Pheno_WHI$agewhi)
coxph(Surv(Data$ENDFOLLOWALLDY, Data$DEATHALL) ~ Pheno_WHI$DNAmPhenoAge+Pheno_WHI$agewhi)

