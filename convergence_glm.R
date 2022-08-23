#' ---
#' title: "Performance of the JSDM model across different simulated scenarios"
#' author: "Giovanni Bianco"
#' date: "August 2022"
#' output: html_document
#' ---
#'
#'
#' To explore whether the value space characteristics of a dataset can affect the convergence of the JSDM, 
#' we ran the JSDM several times on simulated datasets representing species pools 
#' with different degrees of relatedness, species specific traits and responses to 
#' environmental variables 
#'
#'
# load packages
library(dplyr)
library(tidyr)
library(stringr)
library(knitr)
#'
# load dataset 
res<-read.csv("demo_looped_results.csv")
#'
# add a binary variable for convergence -> 1 converged; 0 did not converge 
# non convergence occurs if Rhat>1.1 & n_eff>100
#'
res<-mutate(res, conv=1)

head(res)

for (i in 1:nrow(res)){
  if (res$Rhat[i]>=1.1 & res$n_eff[i]<100){
    res$conv[i]<-0
  }
}
#'
#'
#' We assess whether the parameter space of species responses to the environment (betas),
#' of the species specific traits (zetas), and of the phylogenetic signal (rho) can 
#' determine whether the model converges or not. Convergence is treated as a binomial variable,
#' where 1 stands for successfully converged model and 0 for unsuccesful convergence. As indicator of 
#' the parameter space of the zetas and betas we used their mean and standard deviation.
#' 
#' 
# fix replicate ID
replicates<-rep(1:103,each=70)
res$replicate<-replicates
#'
# data frame for glm
df<-data.frame(matrix(ncol = 7,nrow = 103))
colnames(df)<-c("replicate","mean_beta","sd_beta","mean_z","sd_z","rho","convergence")
df$replicate<-seq(1:103)
#'
for (i in 1:nrow(df)){
  
  tmp<-subset(res,replicate==i)
  
  df$mean_beta[i]<-mean(tmp$TRUE.[grepl("beta",tmp$par)])
  df$sd_beta[i]<-sd(tmp$TRUE.[grepl("beta",tmp$par)])
  df$mean_z[i]<-mean(tmp$TRUE.[grepl("z",tmp$par)])
  df$sd_z[i]<-sd(tmp$TRUE.[grepl("z",tmp$par)])
  df$rho[i]<-tmp$TRUE.[grepl("rho",tmp$par)]
  df$convergence[i]<-tmp$conv[1]
  
}
#'
df$convergence<-as.logical(df$convergence)
#'
#'
# histograms of all parameter
hist(df$mean_beta)
hist(df$mean_z)
hist(df$sd_beta)
hist(df$sd_z)
hist(df$rho)
#'
#'Now we can fit a binomial glm that can estimate the role of different parameters 
#' in determining whether the JSDM converges or not. 
#'
#model selection
#'
# complex to simple 
model1<-glm(convergence~mean_beta+mean_z+sd_beta+sd_z+rho+rho:sd_z+rho:sd_beta,data=df,family = binomial)
summary(model1)
#'
model2<-glm(convergence~mean_beta+mean_z+sd_beta+sd_z+rho+rho:sd_beta,data=df,family = binomial)
summary(model2)
#'
model3<-glm(convergence~mean_beta+mean_z+sd_z+rho,data=df,family = binomial)
summary(model3)
#'
model4<-glm(convergence~mean_z+sd_z+rho,data=df,family = binomial)
summary(model4)
#'
model5<-glm(convergence~mean_z+sd_z,data=df,family = binomial)
summary(model5)
anova(model4,model5,test="Chisq")

#' Model 4 is the best model, including rho, mean zeta and standard deviation 
#' of zeta but no interactions

# summary
final_model<-glm(convergence~mean_z+sd_z+rho,data=df,family = binomial)
summary(final_model)
#'
#' We can use the best fitting model to produce figures
#' 
# plotting real vs fitted
estimated_convergence<-round(fitted(final_model))
plot(jitter(as.numeric(df$convergence))~jitter(estimated_convergence))
#'
# probability of convergence as a function of explanatory variables 
plot(fitted(final_model)~df$mean_z)
plot(fitted(final_model)~df$sd_z)
plot(fitted(final_model)~df$rho)
#'
#'
# predict convergence using model fit and produce a graph 
rho<-seq(0,1,0.1)
mean_z<-seq(-0.2,.2,0.05)
sd_z<-seq(0,1,0.1)
#'
newdata<-as.data.frame(expand.grid(rho,mean_z,sd_z))
colnames(newdata)<-c("rho","mean_z","sd_z")
str(newdata)
#'
ypredict<-predict(final_model,newdata = newdata,type = "response")
max(ypredict)
min(ypredict)
fitted.data<-cbind(newdata,ypredict)
#'
plot(fitted.data$ypredict~fitted.data$sd_z)
#'
plot(fitted.data$ypredict~fitted.data$mean_z)
#'
plot(c(0,1),c(0,1),type = "n",xlab="rho",ylab = "conv prob")
#'
plot(ypredict~rho)
#'
#' The only apparent pattern is in the graph with the std dev of z parameters 
#' because its effect in the model output was almost 10x greater than the effect of 
#' rho or the effect of mean z so we can conclude that it is the main driver of 
#' wheter the model converges or not. However to produce graphs for the other 2 
#' explanatory variables i can set the others to their mean values 
#' so that the graph is not really influenced by their effects.
#' 
mean_rho<-mean(df$rho)
mean_sd_z<-mean(df$sd_z)
high_z<-0.2
low_z<--0.2

# probability of convergence as a function of rho with high z values
newdata<-as.data.frame(cbind(rho,rep(high_z,length(rho)),rep(mean_sd_z,length(rho))))
colnames(newdata)<-c("rho","mean_z","sd_z")
str(newdata)

plot(c(0,1),c(0,1),type = "n",xlab="rho",ylab = "conv prob")
lines(rho,predict.glm(final_model,newdata = newdata,type = "response"))

# with low z values 

newdata<-as.data.frame(cbind(rho,rep(low_z,length(rho)),rep(mean_sd_z,length(rho))))
colnames(newdata)<-c("rho","mean_z","sd_z")
str(newdata)

plot(c(0,1),c(0,1),type = "n",xlab="rho",ylab = "conv prob")
lines(rho,predict.glm(final_model,newdata = newdata,type = "response"))

# probability as a function of sd_z with high mean value of z

newdata<-as.data.frame(cbind(rep(mean_rho,length(sd_z)),rep(high_z,length(sd_z)),sd_z))
colnames(newdata)<-c("rho","mean_z","sd_z")
str(newdata)

plot(c(0,1),c(0,1),type = "n",xlab="SD Z",ylab = "conv prob",main= "Convergence when Z= 0.2 and Rho=0.55")
lines(sd_z,predict.glm(final_model,newdata = newdata,type = "response"))

# with low z value 

newdata<-as.data.frame(cbind(rep(mean_rho,length(sd_z)),rep(low_z,length(sd_z)),sd_z))
colnames(newdata)<-c("rho","mean_z","sd_z")
str(newdata)

plot(c(0,1),c(0,1),type = "n",xlab="SD Z",ylab = "conv prob", main = "Convergence when Z=-0.2 and Rho=0.55")
lines(sd_z,predict.glm(final_model,newdata = newdata,type = "response"))

# probability as a function of mean_z

newdata<-as.data.frame(cbind(rep(mean_rho,length(mean_z)),mean_z,rep(mean_sd_z,length(mean_z))))
colnames(newdata)<-c("rho","mean_z","sd_z")
str(newdata)

plot(c(-0.2,0.2),c(0,1),type = "n",xlab="Z value",ylab = "conv prob", main = "Convergence ~ Mean Z when SD Z=-0.34 and Rho=0.55")
lines(mean_z,predict.glm(final_model,newdata = newdata,type = "response"))

