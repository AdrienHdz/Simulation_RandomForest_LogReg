###########################################################################
######################## Initialisation ###################################
###########################################################################

library(caret)
library(MASS)
library(ROCR)

data =  read.csv("diabetes.csv")

# Creation of a list with 4 predictors
predictors = c("Glucose","BloodPressure","BMI","DiabetesPedigreeFunction")

# Applying a logistic regression on the dataset to generate the Betas (for the dependent variable)
model = glm(Outcome~
              +Glucose
            +BloodPressure
            +BMI
            +DiabetesPedigreeFunction, family = "binomial",data = data)

# Get the sigma m to generate independent variables X.

n= 10000
nsim = 100000
set.seed(2019)

s = (var(data$Glucose))/mean(data$Glucose)
a = (mean(data$Glucose)^2)/var(data$Glucose)
Glucose_baseline = rgamma(n,shape=a,scale=s)

s = (var(data$BloodPressure))/mean(data$BloodPressure)
a = (mean(data$BloodPressure)^2)/var(data$BloodPressure)
BloodPressure_baseline = rgamma(n,shape=a,scale=s)

s = (var(data$BMI))/mean(data$BMI)
a = (mean(data$BMI)^2)/var(data$BMI)
BMI_baseline = rgamma(n,shape=a,scale=s)
hist(BMI_baseline)
s = (var(data$DiabetesPedigreeFunction))/mean(data$DiabetesPedigreeFunction)
a = (mean(data$DiabetesPedigreeFunction)^2)/var(data$DiabetesPedigreeFunction)
DiabetesPedigreeFunction_baseline = rgamma(n,shape=a,scale=s)

mega_df <- data.frame("Glucose" = Glucose_baseline, 
                      "BloodPressure" = BloodPressure_baseline, 
                      "BMI" = BMI_baseline,
                      "DiabetesPedigreeFunction" = DiabetesPedigreeFunction_baseline, stringsAsFactors = FALSE)

###########################################################################
#################### Function to add the outliers #########################
###########################################################################

add.outlier.to.vector <- function(vector, amount) {
  set.seed(2019)
  cells.to.modify.upper <- sample(1:length(vector), amount, replace=F)
  
  mean.val <- mean(vector)
  sd.val <- sd(vector)
  max.val.upper <- mean.val + 7 * sd.val 
  min.val.upper <- mean.val + 5 * sd.val
  
  
  vector[cells.to.modify.upper] <- runif(amount, min=  min.val.upper, max= max.val.upper)
  
  return(vector)
}

###########################################################################
################## Function to simulate the data X ########################
##################          with outliers          ########################
###########################################################################

simulated_data <- function(df, predictors) {
  n= 1000
  nsim = 10000
  set.seed(2019)
  amount = round(nrow(df))
  predictors
  for (predictor in predictors) {
    switch(predictor, 
           Glucose={
             df$Glucose_1 = add.outlier.to.vector(df[,c(predictor)],amount*0.01)
             df$Glucose_5 = add.outlier.to.vector(df[,c(predictor)],amount*0.05)
           },
           BloodPressure={
             df$BloodPressure_1 = add.outlier.to.vector(df[,c(predictor)],amount*0.01)
             df$BloodPressure_5 = add.outlier.to.vector(df[,c(predictor)],amount*0.05)
           },
           BMI={
             df$BMI_1 = add.outlier.to.vector(df[,c(predictor)],amount*0.01)
             df$BMI_5 = add.outlier.to.vector(df[,c(predictor)],amount*0.05)
           },
           DiabetesPedigreeFunction={
             df$DiabetesPedigreeFunction_1 = add.outlier.to.vector(df[,c(predictor)],amount*0.01)
             df$DiabetesPedigreeFunction_5 = add.outlier.to.vector(df[,c(predictor)],amount*0.05)
           },
           {
             print('default')
           }
    )
  }
  return(df)
}

###########################################################################
############## Function to count the nbr of outliers ######################
###########################################################################

outliers_count <- function(df_vec, data_mean, data_sd) {
  mean_vector = data_mean
  sd_vector = data_sd
  max_outliers = sum(df_vec>(mean_vector + 3*sd_vector))
  min_outliers = sum(df_vec<(mean_vector - 3*sd_vector))
  
  return(c(min_outliers, max_outliers) )
  
}

###########################################################################
########################## Adding outliers ################################
###########################################################################

# Creation of the big data frame
mega_df = simulated_data(mega_df,predictors)


get_scenario <- function(i,model){
  if(i==1){df= mega_df[c("Glucose_1","BloodPressure","BMI","DiabetesPedigreeFunction")]
  }else if(i==2){df = mega_df[c("Glucose_5","BloodPressure","BMI","DiabetesPedigreeFunction")]
  }else if(i==3){df = mega_df[c("Glucose_1","BloodPressure_1","BMI","DiabetesPedigreeFunction")]
  }else if(i==4){df = mega_df[c("Glucose_5","BloodPressure_5","BMI","DiabetesPedigreeFunction")]
  }else if(i==5){df = mega_df[c("Glucose_1","BloodPressure_1","BMI_1","DiabetesPedigreeFunction")]
  }else if(i==6){df = mega_df[c("Glucose_5","BloodPressure_5","BMI_5","DiabetesPedigreeFunction")]
  }else if(i==7){df = mega_df[c("Glucose_1","BloodPressure_1","BMI_1","DiabetesPedigreeFunction_1")]
  }else if(i==8){df = mega_df[c("Glucose_5","BloodPressure_5","BMI_5","DiabetesPedigreeFunction_5")]
  }else{  df = mega_df[c("Glucose","BloodPressure","BMI","DiabetesPedigreeFunction")]}
  colnames(df) <- c("Glucose","BloodPressure","BMI","DiabetesPedigreeFunction")
  print(summary(df))
  set.seed(2019)
  y = rbinom(n,1,predict.glm(model,df, type="response"))
  df$Outcome = y
  return(df)
}




###########################################################################
############## Function to find an optimal cutoff      ####################
###########################################################################

opt.cut = function(perf, pred)
{
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(True_positive_rate = y[[ind]], False_positive_rate = x[[ind]],
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

###########################################################################
############## Function to divide the dataset  ############################
###########################################################################
separe_df= function(df,mode)
{
  
  # Split for dataset in 80% train and 20% test
  size.train=floor(n*0.80)
  size.test=n-size.train
  
  # Get the obersations ID assigned to the training and test set.
  set.seed(2019)
  id.train=sample(1:n,size.train,replace=FALSE)
  id.test=setdiff(1:n,id.train)
  
  df.train=df[id.train,]
  df.test=df[id.test,]
  # Dividing the sample
  if (mode == "train") {
    return(df.train)
  }
  else
  {
    return(df.test)
    
  }
  
}

###########################################################################
#################### Function to predict with a ###########################
#################### Random Forest (rf)         ###########################
###########################################################################


get_rf = function(df){
  
  # 4. The simulated dataset is ready for the logistic regression and Random Forest
  df$Outcome <- as.factor(df$Outcome)
  df.train = separe_df(df,"train")
  
  # Model 1 Random Forest
  library(randomForest)
  
  
  rf=randomForest(Outcome~.,data=df.train,ntree=600,mtry=2, importance=TRUE)
  return(rf)
}


###########################################################################
################### Function to predict with a  ###########################
################### logistic regression         ###########################
###########################################################################
get_log = function(df){
  
  # 4. The simulated dataset is ready for the logistic regression and Random Forest
  df$Outcome <- as.factor(df$Outcome)
  df.train = separe_df(df,"train")
  
  model_log =glm(Outcome~.,family="binomial",data=df.train)
  return(model_log)
}

###########################################################################
#################### Function to get the AUC    ###########################
###########################################################################
validation_model = function(model_mode, model, df.valid)
{
  if (model_mode == 'rf') {
    prob.predict=predict(model,newdata=df.valid, type="prob")[,2]
    
  }else if (model_mode == 'log') {
    prob.predict=predict.glm(model,newdata=df.valid ,type="response")
  }
  
  pred=prediction(prob.predict,df.valid$Outcome )
  perf=ROCR::performance(pred,"tpr","fpr")
  str(perf)
  par(mfrow=c(1,2))
  perf_lift=ROCR::performance(pred,"lift","rpp")
  plot(perf)
  abline(a=0,b=1)
  plot(perf_lift)
  cutoff = opt.cut(perf, pred)[3]
  test.pred = rep(0, nrow(df.valid))
  
  test.pred[prob.predict > cutoff] = 1
  
  M=table(test.pred, df.valid$Outcome,dnn=c("Prediction","Observation"))
  print(M)
  #print(cat("The cutoff is : ",cutoff, "\n"))
  auc.perf = ROCR::performance(pred,"auc")
  #print(cat("The cutoff is : ",cutoff, "\n"))
  return(data.frame(c(auc.perf@y.values,M)))
}

get_auc = function(df, rf_model, log_model){
  
  df$Outcome <- as.factor(df$Outcome)
  df.valid = separe_df(df,"valid")
  valid_rf = validation_model("rf",rf_model,df.valid)
  valid_log = validation_model("log",log_model, df.valid)
  
  return(data.frame(c(valid_rf,valid_log)))
}


for (i in 1:9) {
  assign(paste0("s_", i),get_scenario(i,model))
  assign(paste0("rf_S_", i),get_rf(get(paste0("s_", i))))
  assign(paste0("log_S_", i),get_log(get(paste0("s_", i))))
  
  s_i = get(paste0("s_", i))
  rf_i = get(paste0("rf_S_", i))
  log_i = get(paste0("log_S_", i))
  
  assign(paste0("auc_S_", i),get_auc(s_i,rf_i,log_i))
  
}


print(outliers_count(mega_df$Glucose,mean(mega_df$Glucose),sd(mega_df$Glucose)))
print(outliers_count(mega_df$Glucose_1,mean(mega_df$Glucose_1,),sd(mega_df$Glucose_1,)))
print(outliers_count(mega_df$Glucose_5,mean(mega_df$Glucose_5),sd(mega_df$Glucose_5)))

print(outliers_count(mega_df$BloodPressure,mean(mega_df$BloodPressure),sd(mega_df$BloodPressure)))
print(outliers_count(mega_df$BloodPressure_1,mean(mega_df$BloodPressure_1,),sd(mega_df$BloodPressure_1,)))
print(outliers_count(mega_df$BloodPressure_5,mean(mega_df$BloodPressure_5),sd(mega_df$BloodPressure_5)))

print(outliers_count(mega_df$DiabetesPedigreeFunction,mean(mega_df$DiabetesPedigreeFunction),sd(mega_df$DiabetesPedigreeFunction)))
print(outliers_count(mega_df$DiabetesPedigreeFunction_1,mean(mega_df$DiabetesPedigreeFunction_1,),sd(mega_df$DiabetesPedigreeFunction_1,)))
print(outliers_count(mega_df$DiabetesPedigreeFunction_5,mean(mega_df$DiabetesPedigreeFunction_5),sd(mega_df$DiabetesPedigreeFunction_5)))

print(outliers_count(mega_df$BMI,mean(mega_df$BMI),sd(mega_df$BMI)))
print(outliers_count(mega_df$BMI_1,mean(mega_df$BMI_1,),sd(mega_df$BMI_1,)))
print(outliers_count(mega_df$BMI_5,mean(mega_df$BMI_5),sd(mega_df$BMI_5)))

