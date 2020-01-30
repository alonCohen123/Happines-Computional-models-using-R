
library(dfoptim)
library(tidyverse)
#get all data from csv file - change directory to the data file here
#dir <- "C:\Users\alon8\Downloads"
#setwd(dir)
alldata <- read.table("SubData.csv",header=F,sep=",");

#initialize the free weights and variables
trials<-60
weights <- c(0.5,0, 0, 0, 0)

names(weights) <- c("Gamma", "w0", "w1","w2","w3")

BIC = function(ll, trials, num_of_params) {
   2 * ll + log(trials) * num_of_params 
}

#############################################################
#############################################################
##########################MODEL A############################
#############################################################
#############################################################

#plot data and current predictions          
getregpred<- function(parms,data) {
  gama <- parms[1]
  getregpred <-c()
  for(trial in c(2:length(data[,1]))){
    #if we have an happiness value
    if(!is.nan(data[trial,8])){
      #sum of cr's by choise
      CR<-sum(data[1:trial,3]*(1-data[1:trial,6])*(gama^(trial-data[1:trial,2])))
      #sum of ev's by choise
      EV<- sum(((data[1:trial,4]+data[1:trial,5])/2)*data[1:trial,6]*(gama^(trial-data[1:trial,2])))
      #sum of rpe's by choise
      RPE<- EV-sum(data[1:trial,7]*data[1:trial,6]*(gama^(-data[1:trial,2])))
      result<- parms[2]+parms[3]*CR+parms[4]*EV+parms[5]*RPE
      getregpred <- c(getregpred,result)
    }
  }
  return(getregpred)
}                                           

#obtain current predictions and compute discrepancy
rmsd <-function(parms, data1) {
  gama=parms[1]
  #gama value supposed to be between 0 to 1
  if (gama>0 && gama <1){
    happiness <- alldata[2:length(data1[,1]),8]
    #get only the happiness value
    happiness<- na.omit(happiness)
    #calculate predictions
    preds<-getregpred(parms, data1)
    rmsd<-sqrt(sum((preds-happiness)^2)/length(preds))
  }else{
    rmsd<-100000
  }
}
#result_model_a <-c()
result_model_a <-matrix(rep(0,5*trials),ncol=5)
ll_model_a<-c()
bic_sub_a<-c()
for(participate in c(1:trials))
{
  optim_result_a <-optim(weights,rmsd, data1=alldata[((participate-1)*31+1):((participate-1)*31+30),])
  #bic score for each suject
  bic_sub_a<-c(bic_sub_a,BIC(unlist(optim_result_a[2]),trials,length(weights)))
  #get the ll value from optim result
  ll_model_a <- c(ll_model_a, optim_result_a[2])
  #add the new set of weights to the result list
  result_model_a[participate,] <-  optim_result_a$par 
}
ll_model_a <- unlist(ll_model_a)
ll_model_a <- mean(ll_model_a)
bic_model_a <- BIC(ll_model_a,trials,length(weights))

#result_model_a
mean_model_A <- c()
var_model_A <- c()
#calculate mean (avg) and variance for eah parameter
for(i in c(1:length(weights))){
  mean_model_A<-c(mean_model_A,mean(result_model_a[,i]))
  var_model_A<-c(var_model_A,var(result_model_a[,i]))
}

#############################################################
#############################################################
##########################MODEL B############################
#############################################################
#############################################################


#plot data and current predictions          
getregpred_b<- function(parms,data) {
  gama <- parms[1]
  getregpred <-c()
  for(trial in c(2:length(data[,1]))){
    #if we have an happiness value
    if(!is.nan(data[trial,8])){
      #sum of cr's by choise
      CR<-sum(data[1:trial,3]*(1-data[1:trial,6])*(gama^(trial-data[1:trial,2])))
      #sum of GR's by choise
      GR<- sum(data[1:trial,7])*data[1:trial,6]*(gama^(trial-data[1:trial,2]))
      result<- parms[2]+parms[3]*CR+parms[4]*GR
      #append(getregpred, result)
      getregpred <- c(getregpred,result)
    }
  }
  return(getregpred)
}                                           

#obtain current predictions and compute discrepancy
rmsd_b <-function(parms, data1) {
  gama=parms[1]
  #gama value supposed to be between 0 to 1
  if (gama>0 && gama <1){
    happiness <- alldata[2:length(data1[,1]),8]
    #get only the happiness value
    happiness<- na.omit(happiness)
    #calculate predictions
    preds<-getregpred_b(parms, data1)
    rmsd<-sqrt(sum((preds-happiness)^2)/length(preds))
  }else{
    rmsd<-100000;
  }
}
result_model_b <-matrix(rep(0,4*trials),ncol=4)
ll_model_b<-c()
bic_sub_b<-c()
for(participate in c(1:trials))
{
  optim_result_b <-optim(weights[1:4],rmsd_b, data1=alldata[((participate-1)*31+1):((participate-1)*31+30),])
  #bic score for each subject
  bic_sub_b<-c(bic_sub_b,BIC(unlist(optim_result_b[2]),trials,length(weights)))
  #get the ll value from optim result
  ll_model_b <- c(ll_model_b, optim_result_b[2])
  #add the new set of weights to the result list
  result_model_b[participate,] <-  optim_result_b$par
}
ll_model_b <- unlist(ll_model_b)
ll_model_b <- mean(ll_model_b)
bic_model_b <- BIC(ll_model_b,trials,length(weights)-1)

#result_model_b
mean_model_b <- c()
var_model_b <- c()
#calculate mean (avg) and variance for eah parameter
for(i in c(1:4)){
  mean_model_b<-c(mean_model_b,mean(result_model_b[,i]))
  var_model_b<-c(var_model_b,var(result_model_b[,i]))
}



#############################################################
#############################################################
##########################PLOT MODEL A#######################
#############################################################
#############################################################
dataA = data.frame(
  Data = c(mean_model_A),
  weights_model_A = factor(c('gamma','w0','w1','w2','w3')),
  SE = c(sd(mean_model_A)) / sqrt(trials)
)

p<-ggplot(dataA, aes(x = weights_model_A, y = Data , fill=weights_model_A))+
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = Data - SE, ymax = Data + SE),
                width = .2,
                # Width of the error bars
                position = position_dodge(.9))

p<- p + labs(title = "")
p + ylab('Mean Value')

#############################################################
#############################################################
####################PLOT VAR Weights A#######################
#############################################################
#############################################################
dataAV = data.frame(
  Data = c(var_model_A),
  weights_variance_model_A = factor(c('gamma','w0','w1','w2','w3')),
  SE = c(sd(var_model_A)) / sqrt(trials)
)

p<-ggplot(dataAV, aes(x = weights_variance_model_A, y = Data , fill=weights_variance_model_A))+
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = Data - SE, ymax = Data + SE),
                width = .2,
                # Width of the error bars
                position = position_dodge(.9))

p<- p + labs(title = "")
p + ylab('Variance Values')

#############################################################
#############################################################
##########################PLOT MODEL B#######################
#############################################################
#############################################################
dataB = data.frame(
  Data = c(mean_model_b),
  weights_model_B = factor(c('gamma','w0','w1','w2')),
  SE = c(sd(mean_model_b)) / sqrt(trials)
)

p<-ggplot(dataB, aes(x = weights_model_B, y = Data , fill=weights_model_B))+
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = Data - SE, ymax = Data + SE),
                width = .2,
                # Width of the error bars
                position = position_dodge(.9))

p<- p + labs(title = "")
p + ylab('Mean Values')


#############################################################
#############################################################
##########################VAR MODEL B########################
#############################################################
#############################################################
dataBV = data.frame(
  Data = c(var_model_b),
  weights_variance_model_B = factor(c('gamma','w0','w1','w2')),
  SE = c(sd(var_model_b)) / sqrt(trials)
)

p<-ggplot(dataBV, aes(x = weights_variance_model_B, y = Data , fill=weights_variance_model_B))+
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = Data - SE, ymax = Data + SE),
                width = .2,
                # Width of the error bars
                position = position_dodge(.9))

p<- p + labs(title = "")
p + ylab('Variance Values')

library(ggplot2)
#############################################################
#############################################################
####################BIC Score comparancy#####################
#############################################################
#############################################################
hist((bic_sub_a-bic_sub_b),xlab = "bic_sub_a-bic_sub_b" ,col = "grey",breaks = 30,xlim = c(-150,250))

#############################################################
#############################################################
####################BIC Model comparancy#####################
#############################################################
#############################################################

dataBIC = data.frame(
  Data = c(bic_model_a,bic_model_b),
  wssBIC = factor(c('BIC Model A','BIC Model B')),
  SE = c(sd(var_model_b)) / sqrt(trials)
)

p<-ggplot(dataBIC, aes(x = wssBIC, y = Data , fill=wssBIC))+
  geom_bar(position = position_dodge(), stat = "identity") +
  geom_errorbar(aes(ymin = Data - SE, ymax = Data + SE),
                width = .2,
                # Width of the error bars
                position = position_dodge(.9))

p<- p + labs(title = "")
p + ylab('Variance Values')


