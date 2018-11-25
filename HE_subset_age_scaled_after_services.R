# This real data analysis uses pre-HE serum creatinine levels as a covariate.

#####################################################
##         Real Data Code on Hypotension Data Set  ##
##                  Kevin Gunn                     ##
##                  5/27/2017                      ##
#####################################################

library("car")
library('ggplot2')
library('MASS')
library("randomForest")
#library("caret")
library("readr")
library("tidyr")
library("dplyr")
library("reshape2")
library("numDeriv")

library("stargazer")

set.seed(1234567)

# Load Dataset
HE_dataset_full_in <- read_csv("~/EMR-research/SQL/csv_files/HE_DS_ELIX.csv")
#View(HE_dataset_full_in)
problems(HE_dataset_full_in)
#The warnings are related to the way the data is read but not actually a problem. 

#########################################################################################

#KS_2D_Gauss performs kernel smoothing on Yt and outputs a vector for imputed data.  
KS_2D_Gauss = function(Yt.v, X_impute , X_label, const=c(2,2)){ 
  
  gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
  #Calculate bandwidth first as it is used in every iteration.
  n = length(Yt.v)
  m = dim(X_impute)[1]
  
  h1 = const[1]*n^(-1/5)
  h2 = const[2]*n^(-1/5)
  
  kde1.g = matrix(0 , ncol = m , nrow = n)
  kde2.g = matrix(0 , ncol = m , nrow = n)
  
  C_np_impute = rep(0, m)
  
  for(j in 1:m){
    kde1.g[,j] = (1/h1)*gauss((X_label[,c("X1")] - X_impute[j,c("X1")])/h1)
    kde2.g[,j] = (1/h2)*gauss((X_label[,c("X2")] - X_impute[j,c("X2")])/h2)
    C_np_impute[j] = sum(kde1.g[,j]*kde2.g[,j]*Yt.v) / sum(kde1.g[,j]*kde2.g[,j])
  }
  C_np_impute[is.nan(C_np_impute)] = 0
  return(C_np_impute)
}

#CV1 function - use directly with KS_2D_Gauss function.
cv.h <- function(Yt.in , x.in , seq.c){
  
  n = length(Yt.in)
  num.folds=5
  bandwidths <- expand.grid(seq.c,seq.c)
  bw_seq_length <- dim(bandwidths)[1]
  fold_MSEs <- matrix(0,nrow=num.folds,
                      ncol=bw_seq_length)
  
  colnames(fold_MSEs) <- 1:bw_seq_length
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    for (bw in 1:bw_seq_length) {
      mx.hat <- KS_2D_Gauss( y.train , x.test[,c("X1","X2")] , 
                             x.train[,c("X1","X2")] , const = c(bandwidths[bw,1] , bandwidths[bw,2]) ) 
      fold_MSEs[fold,paste(bw)] <- sum((y.test - mx.hat)^2)
    }
  }
  
  CV_MSEs = colMeans(fold_MSEs)
  best.bw = bandwidths[which.min(CV_MSEs),]
  return(as.numeric(best.bw))
}

hlscv.ks <- function(Yt.in , x.in , x.impute , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
  m = dim(x.impute)[1]
  
  fold_dfs = data.frame()
  imp_mat = matrix(0, ncol=num.folds , nrow = m)
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. needs to be readjusted with coefficients for k !=k'.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[,c("X1","X2")] , 
                          x.train[,c("X1","X2")] , const = const_in)
    
    mx.hat.imp <- KS_2D_Gauss(y.train , x.impute[,c("X1","X2")] , 
                              x.train[,c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    #print(head(x.test))
    #print(head(offLM))
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM)
    imp_mat[,fold] = mx.hat.imp 
    
  }
  
  #Need to change this part so it fits coefficients for the different folds.
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  offlm.model = lm(return_df$yt.test ~ X1 + X2 , data = return_df , offset = mx_k.hat,
                   weights = prop_score^-1)
  
  beta.off = offlm.model$coefficients
  
  #mu_hat = return_df$mx_k.hat + beta.off%*%t(return_df[c("int","X1","X2")])
  mu_hat = as.vector(rowMeans(imp_mat) + beta.off%*%t(x.impute))
  mu_all = list(mu_hat,beta.off)
  return(mu_all)
}

double_cv.ks1 <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=10
  #m = dim(x.impute)[1]
  
  fold_dfs = data.frame()
  #imp_mat = matrix(0, ncol=num.folds , nrow = m)
  #beta.off = rep(0,num.folds)
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. needs to be readjusted with coefficients for k !=k'.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    prop_score.train = prop_score[-test.rows]
    prop_score.test = prop_score[test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[,c("X1","X2")] , 
                          x.train[,c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM )
  }
  
  #Need to change this part so it fits coefficients for the different folds.
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  #print(return_df)
  mu_bc = c()
  for (k in 1:num.folds){
    # folds when k != k'
    fold_df = subset( return_df , fold==k, select = -fold )
    #print(fold_df)
    fold_lm.model = lm(fold_df$yt.test ~ X1 + X2 , data = fold_df , offset = mx_k.hat )
    beta.off = fold_lm.model$coefficients
    
    # k = k'
    kfold_df = subset( return_df , fold==k, select = -fold )
    mu_hat =  kfold_df$mx_k.hat + beta.off%*%t( kfold_df[,c("int","X1","X2")] )
    mu_bc = c(mu_bc,mu_hat)
  }
  mu_bc = mu_bc[order(as.numeric(rownames(fold_dfs)))]
  return(mu_bc)
}
double_cv.ks0 <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=10
  #m = dim(x.impute)[1]
  
  fold_dfs = data.frame()
  #imp_mat = matrix(0, ncol=num.folds , nrow = m)
  #beta.off = rep(0,num.folds)
  
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  #This for loop gets mx_k.hat. needs to be readjusted with coefficients for k !=k'.
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold,arr.ind=TRUE)
    x.test = x.in[test.rows,]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,]
    y.train = Yt.in[-test.rows]
    prop_score.train = prop_score[-test.rows]
    prop_score.test = prop_score[test.rows]
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[,c("X1","X2")] , 
                          x.train[,c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM )
  }
  
  #Need to change this part so it fits coefficients for the different folds.
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  #print(return_df)
  mu_bc = c()
  for (k in 1:num.folds){
    # folds when k != k'
    fold_df = subset( return_df , fold!=k, select = -fold )
    #print(fold_df)
    fold_lm.model = lm(fold_df$yt.test ~ X1 + X2 , data = fold_df , offset = mx_k.hat )
    beta.off = fold_lm.model$coefficients
    
    # k = k'
    kfold_df = subset( return_df , fold==k, select = -fold )
    mu_hat =  kfold_df$mx_k.hat + beta.off%*%t( kfold_df[,c("int","X1","X2")] )
    mu_bc = c(mu_bc,mu_hat)
  }
  mu_bc = mu_bc[order(as.numeric(rownames(fold_dfs)))]
  return(mu_bc)
}

###########################################################################################

# Influence Function
IF_sd <- function(XtX.inv,X,Y,A,mv,prop ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm1 = matrix(0,ncol=p,nrow=p)
  varm0 = matrix(0,ncol=p,nrow=p)
  n = length(Y)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # Variance  of IF
      #print((Y[i] - mv[i]))
      IF1 = X[i,]*(Y[i] - mv[i])*prop[i]^-1
      Var_IF1 = ( IF1%*%t(IF1) ) 
      varm1 = varm1 + Var_IF1
    } else{
      #print((Y[i] - mv[i]))
      IF0 = -X[i,]*(Y[i] - mv[i])*(1-prop[i])^-1 
      Var_IF0 = ( IF0%*%t(IF0) ) 
      varm0 = varm0 + Var_IF0
    } 
    
  }
  print(varm1+varm0)
  V1 = XtX.inv%*%varm1%*%t(XtX.inv) 
  V0 = XtX.inv%*%varm0%*%t(XtX.inv) 
  SE = sqrt(diag((V1+V0)/n^2 ) )
  #print(V1);print(V0)
  return(SE)
}

IF_ols_sd <- function(XtX.inv,X,Yt,mv ){
  
  p = dim(XtX.inv)[2]
  varm = matrix(0,ncol=p,nrow=p)
  n = length(Yt)
  
  for(i in 1:n){
    
    # Variance  of IF
    IF = X[i,]*(Yt[i] - mv[i])
    Var_IF = ( IF%*%t(IF) ) 
    varm = varm + Var_IF
    
  }
  print(varm)
  V = XtX.inv%*%varm%*%t(XtX.inv) 
  SE = sqrt(diag((V)/n^2 ) )
  return(SE)
}

# Influence Function with n/N -/> 0.
IF_sd_N <- function(XtX.inv,X,Y,A, mv1, mv0, prop, beta.in ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm1 = matrix(0,ncol=p,nrow=p)
  varm0 = matrix(0,ncol=p,nrow=p)
  n = length(Y)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # Variance  of IF
      #print((Y[i] - mv[i]))
      IF1 = X[i,]*(Y[i] - mv1[i])*prop[i]^-1 + X[i,]*(mv1[i] - mv0[i] - beta.in*X[i,]) 
      Var_IF1 = ( IF1%*%t(IF1) ) 
      varm1 = varm1 + Var_IF1
    } else{
      #print((Y[i] - mv[i]))
      IF0 = -X[i,]*(Y[i] - mv0[i])*(1-prop[i])^-1 + X[i,]*(mv1[i] - mv0[i] - beta.in*X[i,]) 
      Var_IF0 = ( IF0%*%t(IF0) ) 
      varm0 = varm0 + Var_IF0
    } 
    
  }
  print(varm1+varm0)
  V1 = XtX.inv%*%varm1%*%t(XtX.inv) 
  V0 = XtX.inv%*%varm0%*%t(XtX.inv) 
  SE = sqrt(diag((V1+V0)/n^2 ) )
  #print(V1);print(V0)
  return(SE)
}

# TR-OLS Influence Function with estimated propensity score.
IF_ols_sd_estprop <- function(XtX.inv,X,Yt,A,mv, Y, alpha.in, Hess_mat, prop){
  
  p = dim(XtX.inv)[2]
  varm = matrix(0,ncol=p,nrow=p)
  n = length(Yt)
  
  deriv_est = matrix(0 , ncol=p,nrow=n)
  for(i in 1:n){
    deriv_est[i,] =  X[i,]*( A[i]*Y[i]*exp(-sum(alpha.in*X[i,]))*(exp(sum(2*alpha.in*X[i,]))-1) - Y[i]*exp(sum(alpha.in*X[i,]) ))
  }
  IF_d = apply(deriv_est, 2 , mean)
  
  #print(IF_d)
  
  for(i in 1:n){
    
    # Variance  of IF
    IF = X[i,]*(Yt[i] - mv[i]) - IF_d*
      (Hess_mat %*% X[i,]*(A[i] - prop[i]))
    
    Var_IF = ( IF%*%t(IF) ) 
    
    varm = varm + Var_IF
    
  }
  
  V = XtX.inv%*%varm%*%t(XtX.inv) 
  SE = sqrt(diag((V)/n^2 ) )
  return(SE)
}

# SS Influence Function with estimated propensity score.
IF_sd_estprop_N <- function(XtX.inv,X,Y,A, mv1, mv0, prop, beta.in, alpha.in, IF_deriv, Hess_mat ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm1 = matrix(0,ncol=p,nrow=p)
  varm0 = matrix(0,ncol=p,nrow=p)
  n = length(Y)
  
  deriv_est1 = matrix(0 , ncol=p,nrow=n)
  deriv_est0 = matrix(0 , ncol=p,nrow=n)
  for(i in 1:n){
    
    if(A[i]==1){
      deriv_est1[i,] = IF_deriv[i]
    } else{
      deriv_est0[i,] = IF_deriv[i]
    }
  }
  deriv_est = rbind(deriv_est0,deriv_est1)
  IF_d = apply(deriv_est, 2, mean)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # Variance  of IF
      #print((Y[i] - mv[i]))
      IF1 = X[i,]*(Y[i] - mv1[i])*prop[i]^-1 + X[i,]*(mv1[i] - mv0[i] - beta.in*X[i,]) -
        IF_d*(Hess_mat %*% X[i,]*(A[i]- prop[i]) )
      Var_IF1 = ( IF1%*%t(IF1) ) 
      varm1 = varm1 + Var_IF1
    } else{
      #print((Y[i] - mv[i]))
      IF0 = -X[i,]*(Y[i] - mv0[i])*(1-prop[i])^-1 + X[i,]*(mv1[i] - mv0[i] - beta.in*X[i,]) -
        IF_d*(Hess_mat %*% X[i,]*(A[i]- prop[i]) )
      Var_IF0 = ( IF0%*%t(IF0) ) 
      varm0 = varm0 + Var_IF0
    } 
    
  }
  
  V1 = XtX.inv%*%varm1%*%t(XtX.inv) 
  V0 = XtX.inv%*%varm0%*%t(XtX.inv) 
  SE = sqrt(diag((V1+V0)/n^2 ) )
  #print(V1);print(V0)
  return(SE)
}

g <- function(x, theta) 1 / (1 + exp(-1 * x %*% theta))

logistic_loglik <- function(theta){
  sum(log(g(X.L, theta)) * Lin_ds$trt_ind) + sum((1 - Lin_ds$trt_ind) * log(1 - g(X.L, theta)))
}

# Influence Function with estimated propensity score.
IF_Curve_N <- function(XtX.inv,X,Y,A, mv1, mv0, prop, beta.in ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  n = length(Y)
  IF1 = matrix(0,ncol=p,nrow=n)
  IF0 = matrix(0,ncol=p,nrow=n)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # Variance  of IF
      #print((Y[i] - mv[i]))
      IF1[i,] = XtX.inv %*% X[i,]*(Y[i] - mv1[i])*prop[i]^-1 + X[i,]*(mv1[i] - mv0[i] - beta.in*X[i,]) 
      
    } else{
      #print((Y[i] - mv[i]))
      IF0[i,] = -XtX.inv %*% X[i,]*(Y[i] - mv0[i])*(1-prop[i])^-1 + X[i,]*(mv1[i] - mv0[i] - beta.in*X[i,]) 
    } 
    
  }
  
  IF = rbind(IF0,IF1)
  #print(V1);print(V0)
  return(IF)
}

###########################################################################################
# Estimates the propensity score in the bootstrap procedure.
bstrap_var_prop <- function(X.l1, Y1, X.l0, Y0 , X.U, hc1,hc0, numB ){
  
  n.l1 = length(Y1);n.l0 = length(Y0); n.u = dim(X.U)[1]; p = dim(X.U)[2]
  
  Beta_mat = matrix(0 , ncol=p , nrow = numB)
  Beta_mat2 = matrix(0 , ncol=p , nrow = numB)
  
  Ab = c( rep(0,n.l0), rep(1,n.l1) )
  #prop1 = prop[which(Ab==1)];prop0 = prop[which(Ab==0)]
  
  for( b in 1:numB){
    
    # Bootstrap sample for observed cases.
    samp.l1 = sample.int(n=n.l1,size=n.l1,replace=TRUE)
    samp.l0 = sample.int(n=n.l0,size=n.l0,replace=TRUE)
    X.l1b = X.l1[samp.l1,] ; X.l0b = X.l0[samp.l0,] ; Y1b = Y1[samp.l1] ; Y0b = Y0[samp.l0]
    
    X.lb = cbind(Ab, rbind(X.l0b, X.l1b) )
    
    #est prop
    #estimated propensity score.
    prop_modelb <- glm(Ab ~ X1 + X2 ,family=binomial(link='logit'),data=as.data.frame(X.lb))
    est.propb = prop_modelb$fitted.values
    propb = est.propb
    prop1b = propb[which(Ab==1)];prop0b = propb[which(Ab==0)]
    
    # Bootstrap sample for incomplete cases.
    samp.u = sample.int(n=n.u,size=n.u,replace=TRUE)
    X.Ub = X.U[samp.u,]
    
    #Predictions with bias reduction by least squares. Average over 10 cross-validated estimates.
    #CX_ls_mat = matrix(0,nrow=dim(X.Ub)[1],ncol=10)
    #for(k in 1:10){
    #  
    #  mx.hat_ls1 = hlscv.ks(Y1b ,x.in= X.l1b,
    #                        x.impute=X.Ub
    #                        ,const_in=hc1, prop_score = prop1b )[[1]]
    #  
    #  mx.hat_ls0 = hlscv.ks(Y0b , x.in = X.l0b,
    #                        x.impute = X.Ub
    #                        , const = hc0, prop_score = 1 - prop0b )[[1]]
    #  
    #  CX_ls_iter = mx.hat_ls1 - mx.hat_ls0
    #  CX_ls_mat[,k] = CX_ls_iter
    #  
    #} 
    #CX_ls = apply(CX_ls_mat,1,mean)
    
    mx.hat_ls1 = hlscv.ks(Y1b ,x.in= X.l1b,
                          x.impute=X.Ub
                          ,const_in=hc1, prop_score = prop1b )[[1]]
    
    mx.hat_ls0 = hlscv.ks(Y0b , x.in = X.l0b,
                          x.impute = X.Ub
                          , const_in = hc0, prop_score = 1 - prop0b )[[1]]
    
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    # Get regression coefficients for SS estimator.
    Reg.CX_ls = lm(CX_ls ~ X1 + X2  , data = as.data.frame(X.Ub))
    beta.CX_ls = Reg.CX_ls$coefficients
    print(beta.CX_ls)
    Beta_mat[b,] = beta.CX_ls
    
    # Get regression coefficients for complete cases.
    Yt_b = c(Y0b,Y1b); prop_b = c(prop0b,prop1b)
    Yt = Yt_b*(Ab - prop_b) / (prop_b*(1 - prop_b))
    Reg.Yt= lm(Yt ~ X1 + X2  , data = as.data.frame(X.lb))
    beta.Yt = Reg.Yt$coefficients
    print(beta.Yt)
    Beta_mat2[b,] = beta.Yt
    
  }
  
  # CX_ls beta
  beta.CX_ls.avg = apply(Beta_mat,2,mean)
  VB = var(Beta_mat)
  se = sqrt(diag(VB))
  
  # Yt Ols
  beta.ols.avg = apply(Beta_mat2,2,mean)
  VB_Yt = var(Beta_mat2)
  se_Yt = sqrt(diag(VB_Yt))
  
  se_list = list(CX_ls_beta_avg = beta.CX_ls.avg, CX_ls_se = se, beta.ols.avg = beta.ols.avg,
                 Yt_se = se_Yt)
  
  return(se_list)
}

# This version does not estimate the propensity score.
bstrap_var <- function(X.l1, Y1, X.l0, Y0 , X.U, hc1,hc0 , prop, numB ){
  
  n.l1 = length(Y1);n.l0 = length(Y0); n.u = dim(X.U)[1]; p = dim(X.U)[2]
  
  Beta_mat = matrix(0 , ncol=p , nrow = numB)
  Beta_mat2 = matrix(0 , ncol=p , nrow = numB)
  
  Ab = c( rep(0,n.l0), rep(1,n.l1) )
  prop1 = prop[which(Ab==1)];prop0 = prop[which(Ab==0)]
  
  for( b in 1:numB){
    
    # Bootstrap sample for observed cases.
    samp.l1 = sample.int(n=n.l1,size=n.l1,replace=TRUE)
    samp.l0 = sample.int(n=n.l0,size=n.l0,replace=TRUE)
    X.l1b = X.l1[samp.l1,] ; X.l0b = X.l0[samp.l0,] ; Y1b = Y1[samp.l1] ; Y0b = Y0[samp.l0]
    
    X.lb = cbind(Ab, rbind(X.l0b, X.l1b) )
    
    prop1b = prop1[samp.l1];prop0b = prop0[samp.l0]
    
    # Bootstrap sample for incomplete cases.
    samp.u = sample.int(n=n.u,size=n.u,replace=TRUE)
    X.Ub = X.U[samp.u,]
    
    #Predictions with bias reduction by least squares. Average over 10 cross-validated estimates.
    #CX_ls_mat = matrix(0,nrow=dim(X.Ub)[1],ncol=10)
    #for(k in 1:10){
    #  
    #  mx.hat_ls1 = hlscv.ks(Y1b ,x.in= X.l1b,
    #                        x.impute=X.Ub
    #                        ,const_in=hc1, prop_score = prop1b )[[1]]
    #  
    #  mx.hat_ls0 = hlscv.ks(Y0b , x.in = X.l0b,
    #                        x.impute = X.Ub
    #                        , const = hc0, prop_score = 1 - prop0b )[[1]]
    #  
    #  CX_ls_iter = mx.hat_ls1 - mx.hat_ls0
    #  CX_ls_mat[,k] = CX_ls_iter
    #  
    #} 
    #CX_ls = apply(CX_ls_mat,1,mean)
    
    mx.hat_ls1 = hlscv.ks(Y1b ,x.in= X.l1b,
                          x.impute=X.Ub
                          ,const_in=hc1, prop_score = prop1b )[[1]]
    
    mx.hat_ls0 = hlscv.ks(Y0b , x.in = X.l0b,
                          x.impute = X.Ub
                          , const_in = hc0, prop_score = 1 - prop0b )[[1]]
    
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    # Get regression coefficients for SS estimator.
    Reg.CX_ls = lm(CX_ls ~ X1 + X2  , data = as.data.frame(X.Ub))
    beta.CX_ls = Reg.CX_ls$coefficients
    print(beta.CX_ls)
    Beta_mat[b,] = beta.CX_ls
    
    # Get regression coefficients for complete cases.
    Yt_b = c(Y0b,Y1b); prop_b = c(prop0b,prop1b)
    Yt = Yt_b*(Ab - prop_b) / (prop_b*(1 - prop_b))
    
    #print(length(Yt_b));print(dim(X.lb));print(colnames(X.lb))
    
    Reg.Yt= lm(Yt ~ X1 + X2  , data = as.data.frame(X.lb))
    beta.Yt = Reg.Yt$coefficients
    print(beta.Yt)
    Beta_mat2[b,] = beta.Yt
    
  }
  
  # CX_ls beta
  beta.CX_ls.avg = apply(Beta_mat,2,mean)
  VB = var(Beta_mat)
  se = sqrt(diag(VB))
  
  # Yt Ols
  beta.ols.avg = apply(Beta_mat2,2,mean)
  VB_Yt = var(Beta_mat2)
  se_Yt = sqrt(diag(VB_Yt))
  
  se_list = list(CX_ls_beta_avg = beta.CX_ls.avg, CX_ls_se = se, beta.ols.avg = beta.ols.avg,
                 Yt_se = se_Yt)
  
  return(se_list)
}

# This version does not estimate the propensity score.
bstrap_value <- function(X.l1, Y1, X.l0, Y0 , X.U, hc1,hc0 , prop, Y_nosc1, Y_nosc0, numB ){
  
  n.l1 = length(Y1);n.l0 = length(Y0); n.u = dim(X.U)[1]; p = dim(X.U)[2]
  
  Beta_mat = matrix(0 , ncol=p , nrow = numB)
  Beta_mat2 = matrix(0 , ncol=p , nrow = numB)
  
  Ab = c( rep(0,n.l0), rep(1,n.l1) )
  prop1 = prop[which(Ab==1)];prop0 = prop[which(Ab==0)]
  
  Value_SS_vec = rep(0,numB)
  Value_OLS_vec = rep(0,numB)
  
  for( b in 1:numB){
    
    # Bootstrap sample for observed cases.
    samp.l1 = sample.int(n=n.l1,size=n.l1,replace=TRUE)
    samp.l0 = sample.int(n=n.l0,size=n.l0,replace=TRUE)
    X.l1b = X.l1[samp.l1,] ; X.l0b = X.l0[samp.l0,] ; Y1b = Y1[samp.l1] ; Y0b = Y0[samp.l0]
    Y1b_noscale = Y_nosc1[samp.l1] ; Y0b_noscale = Y_nosc0[samp.l0]
    Y_v = c(Y0b_noscale, Y1b_noscale)
    
    X.lb = cbind(Ab, rbind(X.l0b, X.l1b) )
    
    prop1b = prop1[samp.l1];prop0b = prop0[samp.l0]
    
    # Bootstrap sample for incomplete cases.
    samp.u = sample.int(n=n.u,size=n.u,replace=TRUE)
    X.Ub = X.U#[samp.u,]
    
    #Predictions with bias reduction by least squares. Average over 10 cross-validated estimates.
    #CX_ls_mat = matrix(0,nrow=dim(X.Ub)[1],ncol=10)
    #for(k in 1:10){
    #  
    #  mx.hat_ls1 = hlscv.ks(Y1b ,x.in= X.l1b,
    #                        x.impute=X.Ub
    #                        ,const_in=hc1, prop_score = prop1b )[[1]]
    #  
    #  mx.hat_ls0 = hlscv.ks(Y0b , x.in = X.l0b,
    #                        x.impute = X.Ub
    #                        , const = hc0, prop_score = 1 - prop0b )[[1]]
    #  
    #  CX_ls_iter = mx.hat_ls1 - mx.hat_ls0
    #  CX_ls_mat[,k] = CX_ls_iter
    #  
    #} 
    #CX_ls = apply(CX_ls_mat,1,mean)
    
    mx.hat_ls1 = hlscv.ks(Y1b ,x.in= X.l1b,
                          x.impute=X.Ub
                          ,const_in=hc1, prop_score = prop1b )[[1]]
    
    mx.hat_ls0 = hlscv.ks(Y0b , x.in = X.l0b,
                          x.impute = X.Ub
                          , const_in = hc0, prop_score = 1 - prop0b )[[1]]
    
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    # Get regression coefficients for SS estimator.
    Reg.CX_ls = lm(CX_ls ~ X1 + X2  , data = as.data.frame(X.Ub))
    beta.CX_ls = Reg.CX_ls$coefficients
    print(beta.CX_ls)
    Beta_mat[b,] = beta.CX_ls
    
    # Get regression coefficients for complete cases.
    Y_b = c(Y0b,Y1b); prop_b = c(prop0b,prop1b)
    Yt = Y_b*(Ab - prop_b) / (prop_b*(1 - prop_b))
    
    #print(length(Yt_b));print(dim(X.lb));print(colnames(X.lb))
    
    Reg.Yt= lm(Yt ~ X1 + X2  , data = as.data.frame(X.lb))
    beta.Yt = Reg.Yt$coefficients
    print(beta.Yt)
    Beta_mat2[b,] = beta.Yt
    
    # Treatment recommendation
    OLS_trt <- ifelse(beta.Yt%*%t(X.lb[,c("int","X1","X2")]) <= 0 ,1,0)
    SS_trt <- ifelse(beta.CX_ls%*%t(X.lb[,c("int","X1","X2")]) <= 0 ,1,0)
    
    # Value function Estimation
    value_ind_ss <- ifelse(SS_trt == Ab,1,0 )
    value_ind_ols <- ifelse(OLS_trt == Ab,1,0 )
    
    #print(sum(value_ind_ss));print(sum(value_ind_ols))
    
    Value_SS <- mean(value_ind_ss*Y_v*prop_b^-1)
    Value_OLS <- mean(value_ind_ols*Y_v*prop_b^-1)
    
    Value_SS_vec[b] <- Value_SS 
    Value_OLS_vec[b] <- Value_OLS 
  }
  
  # CX_ls beta
  beta.CX_ls.avg = apply(Beta_mat,2,mean)
  VB = var(Beta_mat)
  se = sqrt(diag(VB))
  
  # Yt Ols
  beta.ols.avg = apply(Beta_mat2,2,mean)
  VB_Yt = var(Beta_mat2)
  se_Yt = sqrt(diag(VB_Yt))
  
  se_val_list = list(CX_ls_beta_avg = beta.CX_ls.avg, CX_ls_se = se, beta.ols.avg = beta.ols.avg,
                     Yt_se = se_Yt, mean_value_ss = mean(Value_SS_vec), sd_value_ss = sd(Value_SS_vec),
                     mean_value_ols = mean(Value_OLS_vec),sd_value_ols = sd(Value_OLS_vec))
  
  return(se_val_list)
}

###########################################################################################

# Initial Clean up.
head(HE_dataset_full_in)

# Make elixhauser scores numeric.
HE_dataset_full_in[,2:4] <- lapply(HE_dataset_full_in[,2:4], as.numeric)
# Examine Patients that had Hypotensive episodes of at least 20 mins.
HE_dataset_full <- HE_dataset_full_in[which(HE_dataset_full_in$he_duration_mins >= 20),]

# Remove patients that do not have baseline creatinine values recorded.
HE_Cohort_Fluid_VP_0 <- subset(HE_dataset_full,!is.na(HE_dataset_full$pre_se_cr) )

# Remove patients that recieved both.
HE_Cohort_Fluid_VP_1 <- subset(HE_Cohort_Fluid_VP_0,!( ( iv_fluid_ind ==1) & ( vp_ind ==1 ) ) )

# Create missing indicator for each patient.
missing_ind_resp <- ifelse( is.na(HE_Cohort_Fluid_VP_1$chg_ser) , 1 , 0)

missing_ind_trt <- ifelse(( ( HE_Cohort_Fluid_VP_1$iv_fluid_ind ==0) & ( HE_Cohort_Fluid_VP_1$vp_ind ==0 ) ), 1 , 0)

missing_ind_count <- missing_ind_resp + missing_ind_trt 

missing_ind <- ifelse(missing_ind_count>0, 1, 0)

###########################
## Dataset with missing indicator ##
HE_Cohort_Fluid_VP2 <- cbind(HE_Cohort_Fluid_VP_1, missing_ind)
###########################
# remove 4 patients with missing total_urine_vol output.
HE_Cohort_Fluid_VP3 <- HE_Cohort_Fluid_VP2[which(!is.na(HE_Cohort_Fluid_VP2$tot_urine_vol)),]
HE_Cohort_Fluid_VP4 <- HE_Cohort_Fluid_VP3[which(!is.na(HE_Cohort_Fluid_VP3$mean_spo2)),]
HE_Cohort_Fluid_VP_elix <- HE_Cohort_Fluid_VP4[which(!is.na(HE_Cohort_Fluid_VP4$elixhauser_vanwalraven)),]
HE_Cohort_Fluid_VP_med_surg <- HE_Cohort_Fluid_VP_elix[which(!(HE_Cohort_Fluid_VP_elix$curr_service
                                                               %in% c("ENT","GU","GYN","OBS","TRAUM"))),]
HE_Cohort_Fluid_VP5 <- HE_Cohort_Fluid_VP_med_surg[-which(HE_Cohort_Fluid_VP_med_surg$pre_se_cr %in% 
                                                            sort(HE_Cohort_Fluid_VP_med_surg$pre_se_cr,decreasing = TRUE)[1:3]),]
############################################################################################################
# Normalize covariates.

norm_tuv <- scale(HE_Cohort_Fluid_VP5$tot_urine_vol)

norm_pre_cr <-  scale(HE_Cohort_Fluid_VP5$pre_se_cr)

norm_mean_map <-  scale(HE_Cohort_Fluid_VP5$mean_map)

norm_mean_spo2 <-  scale(HE_Cohort_Fluid_VP5$mean_spo2)

norm_mean_heart_rate <-  scale(HE_Cohort_Fluid_VP5$mean_heart_rate)

norm_age <-  scale(HE_Cohort_Fluid_VP5$age)

norm_chg_ser <- scale(HE_Cohort_Fluid_VP5$chg_ser)

log_elix <- log(-min(HE_Cohort_Fluid_VP5$elixhauser_vanwalraven) + 1 + HE_Cohort_Fluid_VP5$elixhauser_vanwalraven )
###########################################################################################################

## Dataset with scaled covariates.
HE_Cohort_Fluid_VP_services <- cbind(HE_Cohort_Fluid_VP5[,-29], norm_pre_cr, norm_tuv, norm_mean_map,
                                     norm_mean_spo2, norm_mean_heart_rate, norm_age, norm_chg_ser, log_elix )

# service indicator -> did patient recieve surgical or medical services?
service_ind <- ifelse(HE_Cohort_Fluid_VP_services$curr_service %in% c("MED","CMED","OMED","NMED"),1,0)

HE_Cohort_Fluid_VP <- cbind(HE_Cohort_Fluid_VP_services,service_ind)

# Seperate two treatment groups.
IV_Cohort <- subset( HE_Cohort_Fluid_VP, iv_fluid_ind ==1 )
VP_Cohort <- subset( HE_Cohort_Fluid_VP, vp_ind ==1 )


IV_Cohort_No_NA <- na.omit(IV_Cohort)
# which.min(IV_Cohort_No_NA$chg_ser) - Large baseline serum creatinine measure
#IV_Cohort_No_NA <- IV_Cohort_No_NA[which( IV_Cohort_No_NA$chg_ser != min(IV_Cohort_No_NA$chg_ser)),]

#norm_chg_ser <- scale(IV_Cohort_No_NA$chg_ser)

VP_Cohort_No_NA <- na.omit(VP_Cohort) 
# min value way off so most likely input into database incorrectly.
#VP_Cohort_No_NA <- VP_Cohort_No_NA[which( VP_Cohort_No_NA$chg_ser != min(VP_Cohort_No_NA$chg_ser)),]

trt_set1 <- rbind(IV_Cohort_No_NA, VP_Cohort_No_NA)

#trt_ind = 1 if IV fluid, trt_ind = 0 if vassopressor.
trt_ind <- c( rep(0,dim(IV_Cohort_No_NA)[1]) , rep(1,dim(VP_Cohort_No_NA)[1]) ) 

trt_set <- cbind(trt_set1, trt_ind)

###########################################################################################################
#############################################################################################################
###
# Missing Data set

### NOTE: HAVE TO DROP PATIENTS W/O BASELINE CREATININE MEASUREMENTS.
unlabeled_set <- subset(HE_Cohort_Fluid_VP, missing_ind == 1 )

##################################################################################################################
# subset 300 from labeled set.
n=300

# VP trt.
VP_count = length(which(trt_set$trt_ind==1))
# IV fluid trt
IV_count = length(which(trt_set$trt_ind==0))
#proportions
IV_prop = IV_count / (IV_count + VP_count)
VP_prop = VP_count / (VP_count + IV_count)

# subsample
IV_samp_num <- sample(dim(IV_Cohort_No_NA)[1], size = floor(n*IV_prop))
VP_samp_num <- sample(dim(VP_Cohort_No_NA)[1], size = ceiling(n*VP_prop))

# Simulated complete data.
IV_subset <- IV_Cohort_No_NA[IV_samp_num, ]
VP_subset <- VP_Cohort_No_NA[VP_samp_num, ]

# left out data.
IV_out <- IV_Cohort_No_NA[-IV_samp_num, ]
VP_out <- VP_Cohort_No_NA[-VP_samp_num, ]
hold_set <- cbind( rep(1,dim(IV_out)[1]+dim(VP_out)[1]), c( rep(0,dim(IV_out)[1]) , rep(1,dim(VP_out)[1]) )  , 
                   rbind(IV_out,VP_out) )

# Updated unlabeled set.
unlabeled_ds <- rbind(unlabeled_set, IV_out, VP_out )
int_unlabeled <- rep( 1, dim(unlabeled_ds)[1] )
unlabeled_ds_int <- cbind( unlabeled_ds, int = int_unlabeled )

#unlabeled set
unlabeled_SS <- as.matrix(unlabeled_ds_int[,c(37,40,31,36)])

colnames(unlabeled_SS) <- c("chg_ser","int","X1","X2")
rownames(unlabeled_SS) <- c()

## Supervised Regression with complete cases.
trt0_num = length(rep(0,dim(IV_subset)[1]))
trt_ind <- c( rep(0,dim(IV_subset)[1]) , rep(1,dim(VP_subset)[1]) ) 
trt_set1 <- rbind(IV_subset,VP_subset)
trt_set_in <- cbind(trt_set1, trt_ind)

trt_propensity <- glm(trt_ind ~ norm_pre_cr+norm_age
                      ,data=trt_set_in,family = binomial)
prop_score <-  trt_propensity$fitted.values

IV_subset_int <- cbind(IV_subset, rep(1,dim(IV_subset)[1]))
VP_subset_int <- cbind(VP_subset, rep(1,dim(VP_subset)[1]))

# Make subsets of dataset to work in semisupervised model.
IV_SS <- as.matrix(cbind(IV_subset_int[,c(14,40,31,36)],prop_score_in = prop_score[1:trt0_num] ) )
VP_SS <- as.matrix(cbind( VP_subset_int[,c(14,40,31,36)],prop_score_in = prop_score[-(1:trt0_num)] ) )
colnames(IV_SS) <- colnames(VP_SS) <- c("chg_ser","int","X1","X2","prop_score")
rownames(IV_SS) <- rownames(VP_SS) <- c()

trt_set_in_2covs <- as.data.frame(cbind(rbind(IV_SS,VP_SS), trt_ind))

hc1 = cv.h(VP_SS[,c("chg_ser")], x.in = VP_SS[,c("X1","X2")],seq.c = seq(0.1,10,0.5) )
hc0 = cv.h(IV_SS[,c("chg_ser")], x.in= IV_SS[,c("X1","X2")],seq.c = seq(0.1,10,0.5) )
print(hc1);print(hc0)

# Same approach as simulations
mx.hat_ls1 = hlscv.ks(VP_SS[,c("chg_ser")],x.in= VP_SS[,c("int","X1","X2")],x.impute=unlabeled_SS[,c("int","X1","X2")]
                      ,const_in=hc1, prop_score = VP_SS[,c("prop_score")])[[1]]

mx.hat_ls0 = hlscv.ks(IV_SS[,c("chg_ser")] , x.in= IV_SS[,c("int","X1","X2")],x.impute=unlabeled_SS[,c("int","X1","X2")]
                      , const_in = hc0, prop_score = 1 - IV_SS[,c("prop_score")])[[1]]

CX_ls = mx.hat_ls1 - mx.hat_ls0

X.imp = as.data.frame(cbind(CX_ls, unlabeled_SS))
# After examining leverage plots. Row 979 is a huge outlier. 
#X.imp = X.imp[-979,]
colnames(X.imp) <- c("CX_ls","chg_ser","int","norm_pre_cr","norm_age")

Reg.CX_ls = lm(-CX_ls ~ norm_pre_cr + norm_age , data = X.imp)
beta.CX_ls = Reg.CX_ls$coefficients
summary(Reg.CX_ls)
leveragePlots(Reg.CX_ls)
outlierTest(Reg.CX_ls)

# Decrease in serum creatinine = good kidney function.
#scale_Y = c(IV_SS[,1], VP_SS[,1])
Yt <- (trt_set_in_2covs$chg_ser*(trt_set_in_2covs$trt_ind - trt_set_in_2covs$prop_score)) / 
  (trt_set_in_2covs$prop_score*(1-trt_set_in_2covs$prop_score)) 

Lin_ds_subset <- as.data.frame(cbind(Yt, trt_set_in_2covs))

Lin_Reg <- lm(-Yt ~ X1 + X2 ,data=Lin_ds_subset )
summary(Lin_Reg)

beta.Yt = Lin_Reg$coefficients

###################################################################################################################
# Variance Estimation

Y_A1 = Lin_ds_subset[which(Lin_ds_subset$trt_ind == 1),2]
X.L1 = cbind( rep(1, length(Y_A1)) , as.matrix(Lin_ds_subset[which(Lin_ds_subset$trt_ind == 1),c(4,5)]) )
colnames(X.L1) <- c("int","X1","X2")

Y_A0 = Lin_ds_subset[which(Lin_ds_subset$trt_ind == 0),2]
X.L0 =  cbind( rep(1, length(Y_A0)) , as.matrix(Lin_ds_subset[which(Lin_ds_subset$trt_ind == 0),c(4,5)]) )
colnames(X.L0) <- c("int","X1","X2")

prop1 = prop_score[trt_ind==1]; prop0 = prop_score[trt_ind==0]

# Inverse of XtX
X.L = rbind(X.L0, X.L1)
Y.if = c(Y_A0, Y_A1)
prop = prop_score
X.m = as.matrix(X.imp[,-(1:2)])

XL_inv = t(as.matrix(X.L)) %*% as.matrix(X.L) / dim(X.L)[1]
Lambda_n_inv = solve(XL_inv)


Xm_inv = t(as.matrix(X.m)) %*% as.matrix(X.m) / dim(X.m)[1]
Lambda_N_inv = solve(Xm_inv)

X.all = rbind(X.m, X.L )
full_xtx = t(as.matrix(X.all)) %*% as.matrix(X.all) / dim(X.all)[1]
Lambda_all_inv = solve(full_xtx)

# IF estimation
# double estimates for IF.

mv1 = double_cv.ks1(-Y_A1, X.L1, hc1, prop1)
mv0 = double_cv.ks0(-Y_A0, X.L0, hc0, 1-prop0)

IF_sd_CXls = IF_sd(Lambda_all_inv, X.L, -Y.if, trt_ind, c(mv0, mv1) , prop )
print(IF_sd_CXls)

IF_sd_OLS = IF_ols_sd(Lambda_all_inv, X.L, -Yt, Lin_Reg$fitted.values)
print(IF_sd_OLS)

# p-value calculation
#pvalue2sided=2*pnorm(-abs(z))

beta.CX_ls / IF_sd_CXls
pvalue_ss = 2*pnorm(-abs(beta.CX_ls / IF_sd_CXls))

Lin_Reg$coefficients / IF_sd_OLS
pvalue_tr = 2*pnorm(-abs(Lin_Reg$coefficients / IF_sd_OLS))

# RE
IF_sd_OLS^2 / IF_sd_CXls^2

#CI
beta.CX_ls + 1.96*IF_sd_CXls
beta.CX_ls - 1.96*IF_sd_CXls

beta.Yt + 1.96*IF_sd_OLS
beta.Yt - 1.96*IF_sd_OLS
###################################################################################################################

# Treatment recommendation
OLS_trt <- ifelse(beta.Yt%*%t(X.all[,c("int","norm_pre_cr","norm_age")]) > 0 ,1,0)
SS_trt <- ifelse(beta.CX_ls%*%t(X.all[,c("int","norm_pre_cr","norm_age")]) > 0 ,1,0)

# Distribution Tests.
ks.test(X.L[,2],X.imp[,4])
ks.test(X.L[,3],X.imp[,5])

wilcox.test(trt_set_in$pre_se_cr,unlabeled_ds_int$pre_se_cr)
wilcox.test(trt_set_in$age,unlabeled_ds_int$age)

ks.test(trt_set_in$pre_se_cr,unlabeled_ds_int$pre_se_cr)
ks.test(trt_set_in$age,unlabeled_ds_int$age)

mean(trt_set_in$pre_se_cr);sd(trt_set_in$pre_se_cr)
mean(unlabeled_ds_int$pre_se_cr);sd(unlabeled_ds_int$pre_se_cr)

mean(trt_set_in$age);sd(trt_set_in$age)
mean(unlabeled_ds_int$age);sd(unlabeled_ds_int$age)

##################################################################################################

sum((OLS_trt)*(SS_trt))

sum((1-OLS_trt)*(1-SS_trt))

sum((1-OLS_trt)*SS_trt)

sum(OLS_trt*(1-SS_trt))

table_22 = cbind(c(sum((OLS_trt)*(SS_trt)),sum(OLS_trt*(1-SS_trt))), 
                 c(sum((1-OLS_trt)*(1-SS_trt)), sum((1-OLS_trt)*SS_trt)))

######################################################################################################
qplot(prop_score, binwidth = 0.03, geom = "histogram", fill=I("red"),xlab="Propensity Score")
