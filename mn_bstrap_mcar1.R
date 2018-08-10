#####################################################
##          Missing Data OTR Simulation            ##
##                  Kevin Gunn                     ##
##                  5/27/2017                     ##
#####################################################

#library('ggplot2')
library('MASS')
library("randomForest")
#library("caret")

set.seed(12345)

#########################################################################################
#########################################################################################
#KS
#KS_2D_Gauss performs kernel smoothing on Yt and outputs a vector for imputed data.  
#KS
#KS_2D_Gauss performs kernel smoothing on Yt and outputs a vector for imputed data.  
KS_2D_Gauss = function(Yt.v, X_impute , X_label, const=2){ 
  
  gauss <- function(x) 1/sqrt(2*pi) * exp(-(x^2)/2)
  #Calculate bandwidth first as it is used in every iteration.
  n = length(Yt.v)
  m = dim(X_impute)[1]
  
  h1 = const*n^(-1/3)
  h2 = const*n^(-1/3)
  
  kde1.g = matrix(0 , ncol = m , nrow = n)
  kde2.g = matrix(0 , ncol = m , nrow = n)
  
  C_np_impute = rep(0, m)
  
  for(j in 1:m){
    kde1.g[,j] = (1/h1)*gauss((X_label$X1 - X_impute$X1[j])/h1)
    kde2.g[,j] = (1/h2)*gauss((X_label$X2 - X_impute$X2[j])/h2)
    C_np_impute[j] = sum(kde1.g[,j]*kde2.g[,j]*Yt.v) / sum(kde1.g[,j]*kde2.g[,j])
  }
  C_np_impute[is.nan(C_np_impute)] = 0
  return(C_np_impute)
}

#CV1 function - use directly with KS_2D_Gauss function.
cv.h <- function(Yt.in , x.in , seq.c){
  
  n = length(Yt.in)
  num.folds=10
  bandwidths <- seq.c
  fold_MSEs <- matrix(0,nrow=num.folds,
                      ncol=length(bandwidths))
  
  colnames(fold_MSEs) = bandwidths
  case.folds <- rep(1:num.folds,length.out=n)
  case.folds <- sample(case.folds)
  
  for (fold in 1:num.folds) {
    test.rows = which(case.folds==fold)
    x.test = x.in[test.rows,-1]
    y.test = Yt.in[test.rows]
    x.train = x.in[-test.rows,-1]
    y.train = Yt.in[-test.rows]
    for (bw in bandwidths) {
      mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                            x.train[c("X1","X2")] , const = bw)
      fold_MSEs[fold,paste(bw)] <- sum((y.test - mx.hat)^2)
    }
  }
  
  CV_MSEs = colMeans(fold_MSEs)
  best.bw = bandwidths[which.min(CV_MSEs)]
  return(best.bw)
}

hlscv.ks <- function(Yt.in , x.in , x.impute , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=10
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
    mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                          x.train[c("X1","X2")] , const = const_in)
    
    mx.hat.imp <- KS_2D_Gauss(y.train , x.impute[c("X1","X2")] , 
                              x.train[c("X1","X2")] , const = const_in)
    
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
  offlm.model = lm(return_df$yt.test ~ X1 + X2 , data = return_df , offset = mx_k.hat,weights = prop_score^-1)
  
  beta.off = offlm.model$coefficients
  
  #mu_hat = return_df$mx_k.hat + beta.off%*%t(return_df[c("int","X1","X2")])
  mu_hat = as.vector(rowMeans(imp_mat) + beta.off%*%t(x.impute))
  mu_all = list(mu_hat,beta.off)
  return(mu_all)
}

beta_cv.ks <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
  
  fold_dfs = data.frame()
  
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
    mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                          x.train[c("X1","X2")] , const = const_in)
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    
    colnames(offLM) <- c("yt.test" , "mx_k.hat","X1","X2","fold")
    
    fold_dfs = rbind(fold_dfs , offLM)
    
  }
  
  #Need to change this part so it fits coefficients for the different folds.
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  offlm.model = lm(return_df$yt.test ~ X1 + X2 , data = return_df , offset = mx_k.hat, weights = prop_score^-1)
  
  beta.off = offlm.model$coefficients
  
  #mu_hat = return_df$mx_k.hat + beta.off%*%t(return_df[c("int","X1","X2")])
  return(beta.off)
}

double_cv.ks <- function(Yt.in , x.in , const_in, prop_score){
  
  ####
  ## Get mu_k.hat for use in the IF functions
  ####
  n = length(Yt.in)
  num.folds=5
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
    
    #kernel smoothing
    mx.hat <- KS_2D_Gauss(y.train , x.test[c("X1","X2")] , 
                          x.train[c("X1","X2")] , const = const_in)
    
    beta_k = beta_cv.ks(y.train,x.train[c("X1","X2")],const_in = const_in, prop_score = prop_score.train )
    
    #Add to data set to perform regression on residuals
    offLM = cbind(y.test , mx.hat , x.test , 
                  rep(fold , length(y.test)) )
    #print(head(x.test))
    #print(head(offLM))
    colnames(offLM) <- c("yt.test" , "mx_k.hat","int","X1","X2","fold")
    
    #bias corrected regression coefficients
    #offlm.model = lm(offLM$yt.test ~ X1 + X2 , data = offLM , offset = mx_k.hat )
    #beta.off = offlm.model$coefficients
    mu_hat = as.vector(offLM$mx_k.hat + beta_k%*%t(x.test))
    bcLM = cbind(offLM,mu_hat)
    
    fold_dfs = rbind(fold_dfs , bcLM )
    #imp_mat[,fold] = mx.hat.imp 
    
  }
  
  #Need to change this part so it fits coefficients for the different folds.
  return_df = fold_dfs[order(as.numeric(rownames(fold_dfs))),]
  
  #mu_hat = return_df$mx_k.hat + beta.off%*%t(return_df[c("int","X1","X2")]
  
  mu_all = return_df[,ncol(return_df)]
  return(mu_all)
}

###########################################################################################

# Influence Function
IF_sd <- function(XtX.inv,X,Y,A,mv,prop ){
  
  #prop1 vector of propensity scores for patients assigned to trt 1
  #prop0 vector of propensity scores for patients assigned to trt 0
  
  p = dim(XtX.inv)[2]
  varm = matrix(0,ncol=p,nrow=p)
  n = length(Y)
  
  for(i in 1:n){
    
    if(A[i]==1){
      # Variance  of IF
      #print((Y[i] - mv[i]))
      IF1 = X[i,]*(Y[i] - mv[i]) + X[i,]*( A[i] - prop[i] )
      Var_IF1 = ( IF1%*%t(IF1) )*prop[i]^(-2)
      varm = Var_IF1 + varm
    } else{
      #print((Y[i] - mv[i]))
      IF0 = X[i,]*(Y[i] - mv[i]) + X[i,]*( A[i] - prop[i] )
      Var_IF0 = ( IF0%*%t(IF0) )*(1-prop[i])^(-2)
      varm = Var_IF0 + varm
    }
    
  }
  
  var = ( XtX.inv%*%varm%*%t(XtX.inv) ) / n^2
  return(sqrt(diag(var)))
}

##########################################################################################
mn_bstrap_sample <- function(X.l, Y, A, prop_score, X.U, numB, hc1,hc0, m){
  
  n.l = length(Y); n.u = dim(X.U)[1]; p = dim(X.U)[2]
  
  Beta_mat = matrix(0 , ncol=p+1 , nrow = numB)
  
  for( b in 1:numB){
    
    # Bootstrap sample for observed cases.
    samp.l = sample.int(n=n.l,size=m,replace=TRUE)
    X.lb = X.l[samp.l,] ; Yb = Y[samp.l] ; Ab = A[samp.l]
    
    #est prop
    #estimated propensity score.
    prop_modelb <- glm(Ab ~ . ,family=binomial(link='logit'),data=as.data.frame(X.lb))
    est.propb = prop_modelb$fitted.values
    propb = est.propb
    
    # Bootstrap sample for incomplete cases.
    samp.u = sample.int(n=n.u,size=n.u,replace=TRUE)
    X.Ub = X.U[samp.u,]
    
    # Seperate subjects given trt 1 and trt 0.
    Y_A1b = Yb[which(Ab==1)];Y_A0b = Yb[which(Ab==0)]
    X.labelA1b = X.lb[Ab==1,];X.labelA0b = X.lb[Ab==0,]
    prop1b = propb[which(Ab==1)];prop0b = propb[which(Ab==0)]
    
    #Kernel Smoothing.
    X.l_1b = cbind(rep(1,dim(X.labelA1b)[1]),X.labelA1b)
    X.l_0b = cbind(rep(1,dim(X.labelA0b)[1]),X.labelA0b)
    X.i_KSb = cbind(rep(1,dim(X.Ub)[1]),X.Ub)
    
    #Predictions with bias reduction by least squares
    mx.hat_ls1 = hlscv.ks(Y_A1b,x.in= X.l_1b,x.impute=X.i_KSb,const_in=hc1, prop_score = prop1b )[[1]]
    mx.hat_ls0 = hlscv.ks(Y_A0b , X.l_0b, X.i_KSb, const = hc0, prop_score = 1 - prop0b )[[1]]
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    # Get regression coefficients
    formula.CX_ls  = as.formula(paste("CX_ls ~", paste(colnames(X.U), collapse="+")))
    Reg.CX_ls = lm(formula.CX_ls  , data = X.Ub)
    beta.CX_ls = Reg.CX_ls$coefficients
    #print(beta.CX_ls)
    Beta_mat[b,] = beta.CX_ls
    
  }
  
  return(Beta_mat)
}


get_m_MCAR <- function(X.l, Y, A, X.U, numB, hc1,hc0, beta_est){
  
  n = length(Y)
  N = dim(X.U)[1]
  #Let p be 0.85
  p = 0.80
  J = c(1,2,3,4,5)
  #Empty matrices
  ecdf_b0_mat = matrix(0,ncol=length(J),nrow=numB)
  ecdf_b1_mat = matrix(0,ncol=length(J),nrow=numB)
  ecdf_b2_mat = matrix(0,ncol=length(J),nrow=numB)
  
  # loop through to get bet best value of m and M.
  for (j in J) {
    mj = ceiling(p^j*n)
    
    # Get m-out-of-n bootstrap samples
    beta_samps = mn_bstrap_sample(X.l, Y, A, prop_score, X.U, 
                                  numB, hc1,hc0, mj)
    ecdf_beta_mat = sqrt(mj)*(beta_samps - beta_est) / numB
    
    ecdf_b0_mat[,j] = sort(ecdf_beta_mat[,1])
    ecdf_b1_mat[,j] = sort(ecdf_beta_mat[,2])
    ecdf_b2_mat[,j] = sort(ecdf_beta_mat[,3])
    
  }
  
  #Find mj that minimizes distance of b0 ecdf.
  Distance_b0 <- sapply(1:(length(J)-1), 
                        function(x) sqrt(sum((ecdf_b0_mat[(x+1),]-ecdf_b0_mat[x,])^2)))
  b0_min_J = which.min(Distance_b0)
  
  #Find mj that minimizes distance of b1 ecdf.
  Distance_b1 <- sapply(1:(length(J)-1), 
                        function(x) sqrt(sum((ecdf_b1_mat[(x+1),]-ecdf_b1_mat[x,])^2)))
  b1_min_J = which.min(Distance_b1)
  
  #Find mj that minimizes distance of b2 ecdf.
  Distance_b2 <- sapply(1:(length(J)-1), 
                        function(x) sqrt(sum((ecdf_b1_mat[(x+1),]-ecdf_b1_mat[x,])^2)))
  b2_min_J = which.min(Distance_b2)
  
  min_J = max( b0_min_J, b1_min_J, b2_min_J )
  
  m_j_return = ceiling(p^min_J * n)
  
  return(m_j_return)
}


mn_bstrap_var <- function(X.l, Y, A, X.U, numB, hc1,hc0, m ){
  
  n.l = length(Y); n.u = dim(X.U)[1]; p = dim(X.U)[2]
  
  Beta_mat = matrix(0 , ncol=p+1 , nrow = numB)
  
  for( b in 1:numB){
    
    # Bootstrap sample for observed cases.
    samp.l = sample.int(n=n.l,size=m,replace=FALSE)
    X.lb = X.l[samp.l,] ; Yb = Y[samp.l] ; Ab = A[samp.l]
    
    #est prop
    #estimated propensity score.
    prop_modelb <- glm(Ab ~ . ,family=binomial(link='logit'),data=as.data.frame(X.lb))
    est.propb = prop_modelb$fitted.values
    propb = est.propb
    
    # Bootstrap sample for incomplete cases.
    samp.u = sample.int(n=n.u,size=n.u,replace=TRUE)
    X.Ub = X.U[samp.u,]
    
    # Seperate subjects given trt 1 and trt 0.
    Y_A1b = Yb[which(Ab==1)];Y_A0b = Yb[which(Ab==0)]
    X.labelA1b = X.lb[Ab==1,];X.labelA0b = X.lb[Ab==0,]
    prop1b = propb[which(Ab==1)];prop0b = propb[which(Ab==0)]
    
    #Kernel Smoothing.
    X.l_1b = cbind(rep(1,dim(X.labelA1b)[1]),X.labelA1b)
    X.l_0b = cbind(rep(1,dim(X.labelA0b)[1]),X.labelA0b)
    X.i_KSb = cbind(rep(1,dim(X.Ub)[1]),X.Ub)
    
    #Predictions with bias reduction by least squares
    mx.hat_ls1 = hlscv.ks(Y_A1b,x.in= X.l_1b,x.impute=X.i_KSb,const_in=hc1, prop_score = prop1b )[[1]]
    mx.hat_ls0 = hlscv.ks(Y_A0b , X.l_0b, X.i_KSb, const = hc0, prop_score = 1 - prop0b )[[1]]
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    # Get regression coefficients
    formula.CX_ls  = as.formula(paste("CX_ls ~", paste(colnames(X.U), collapse="+")))
    Reg.CX_ls = lm(formula.CX_ls  , data = X.Ub)
    beta.CX_ls = Reg.CX_ls$coefficients
    print(beta.CX_ls)
    Beta_mat[b,] = beta.CX_ls
    
  }
  
  quantiles = apply(Beta_mat, 2,function(x){quantile(x,probs= c(0.025,0.975))})
  VB = var(Beta_mat)
  se = sqrt(diag(VB))
  se_quantile = list(se,quantiles)
  return(se)
}
#########################################################################################
bstrap_var <- function(X.l, Y.b, A.b, X.U, numB, hc1,hc0 ){
  
  n.l1 = length(Y.b[which(A.b==1)]);n.l0 = length(Y.b[which(A.b==0)])
  n.u = dim(X.U)[1]; p = dim(X.U)[2]
  
  # Estimate Y_A1 and Y_A0 separately. 
  Y_A1 = Y.b[which(A.b==1)];Y_A0 = Y.b[which(A.b==0)]
  X.labelA1 = X.l[which(A.b==1),];X.labelA0 = X.l[which(A.b==0),]
  
  Z1 = cbind(Y_A1,X.labelA1)
  rownames(Z1) <- c()
  Z0 = cbind(Y_A0,X.labelA0)
  rownames(Z0) <- c()
  
  Beta_mat = matrix(0 , ncol=p+1 , nrow = numB)
  
  for( b in 1:numB){
    
    samp.l1 = sample(x=n.l1,size=n.l1,replace=TRUE)
    samp.l0 = sample(x=n.l0,size=n.l0,replace=TRUE)
    samp.u = sample(x=n.u,size=n.u,replace=TRUE)
    
    # Estimate Y_A1 and Y_A0 separately. 
    Z1b = Z1[samp.l1,]
    Z1b = Z1b[order(as.numeric(rownames(Z1b))),]
    rownames(Z1b) <- c()
    Z0b = Z0[samp.l0,]
    Z0b = Z0b[order(as.numeric(rownames(Z0b))),]
    rownames(Z0b) <- c()
    Ab = c(rep(1,n.l1),rep(0,n.l0))
    X.lb = rbind(Z1b[,-1],Z0b[,-1])
      
    #estimated propensity score.
    prop_model <- glm(Ab ~ X1+X2 ,family=binomial(link='logit'),
                      data=as.data.frame(X.lb))
    est.prop = prop_model$fitted.values
    prop1b = est.prop[1:n.l1]
    prop0b = est.prop[-(1:n.l1)]
    
    X.Ub = X.U[samp.u,]
    ##################################################################################################
    
    #Kernel Smoothing approaches.
    Z.l_1b = cbind(rep(1,n.l1),Z1b)
    Z.l_0b = cbind(rep(1,n.l0),Z0b)
    X.i_KSb = cbind(rep(1,n.u),X.Ub)
    #X.l_KSb = cbind(rep(1,n),X.label)
    #Get h by grid search
    
    #Predictions with bias reduction by least squares
    mx.hat_ls1_b = hlscv.ks(Z.l_1b$Y_A1,x.in= Z.l_1b[,-2],x.impute=X.i_KSb,
                            const_in=hc1, prop_score = prop1b )[[1]]
    mx.hat_ls0_b = hlscv.ks(Z.l_0b$Y_A0 , x.in= Z.l_0b[,-2], x.impute=X.i_KSb, 
                            const = hc0, prop_score = 1 - prop0b )[[1]]
    CX_ls_b = mx.hat_ls1_b - mx.hat_ls0_b
    
    # Get regression coefficients
    formula.CX_ls_b  = as.formula(paste("CX_ls_b ~", paste(colnames(X.Ub), collapse="+")))
    Reg.CX_ls_b = lm(formula.CX_ls_b  , data = X.Ub)
    beta.CX_ls_b = Reg.CX_ls_b$coefficients
    #print(beta.CX_ls_b)
    Beta_mat[b,] = beta.CX_ls_b
    
  }
  quantiles = apply(Beta_mat, 2,function(x){quantile(x,probs= c(0.025,0.975))})
  VB = var(Beta_mat)
  se = sqrt(diag(VB))
  se_quantile = list(se,quantiles)
  return(se_quantile)
}
#########################################################################################

# Monte Carlo Simulation to get true values of regression coefficients
MC_coeff <-function(model = "SL", n.samp = 300000,p,g=NULL,theta,lambda,beta){
  
  
  X = mvrnorm(n=n.samp , mu = rep(0,p), Sigma = diag(p))
  
  colnames(X) <- paste("X",1:p,sep="")
  
  logit_coefs <- rep(1,ncol(X))
  true.prop = plogis( 0.5*X[,1] - 0.5*X[,2] )
  #true.prop = runif(n.samp,min = 0.4,max=0.6)
  #summary(true.prop)
  A = sapply(true.prop, function(x){ rbinom(1,1,x) })
  
  #Generate Y based on the given true model.
  if(model == "QL"){
    Y = as.vector((X%*%theta)^2 + A*lambda*(X%*%beta) + rnorm(n.samp,0,1))
    mu_X = (X%*%theta)^2; CX = X%*%beta
  }else if(model == "L2L"){
    Y = as.vector( X%*%beta + A*lambda*(X%*%beta) + rnorm(n.samp,0,1))
    mu_X = X%*%beta; CX = X%*%beta
  }else if(model == "L2NL3"){
    Y = as.vector(X%*%beta + A*lambda*((X%*%g)^3) + rnorm(n.samp,0,1))
    mu_X = X%*%beta; CX = (X%*%g)^3
  }else if(model == "QNL3"){
    Y = as.vector((X%*%theta)^2 + A*lambda*((X%*%g)^3) + rnorm(n.samp,0,1))
    mu_X = (X%*%theta)^2; CX = (X%*%g)^3
  }else if(model == "L2S"){
    Y = as.vector(X%*%beta + A*lambda*(sin(X%*%beta) ) + rnorm(n.samp,0,1))
    mu_X = X%*%beta; CX = sin(X%*%beta)
  }else if(model == "QS"){
    Y = as.vector((X%*%theta)^2 + A*lambda*( sin(X%*%beta) ) + rnorm(n.samp,0,1))
    mu_X = (X%*%theta)^2; CX = sin(X%*%beta)
  }
  
  
  Yt.trueP = Y*(A - true.prop) / (true.prop*(1 - true.prop))
  Yt = Yt.trueP
  
  
  formula.linReg = as.formula(paste("CX~", paste(colnames(X), collapse="+")))
  regM = lm(formula.linReg, data=as.data.frame(X))
  
  
  V_0 = mu_X + ifelse(CX>0,1,0)*CX
  V_b0 = mu_X + ifelse(cbind(rep(1,dim(X)[1]),X)%*%regM$coefficients>0,1,0)*CX
  
  ret_list = list(coefs = regM$coefficients , ux = mu_X , X = X, CX=CX,V_0=V_0,V_b0=V_b0 )
  return(ret_list)
}

#########################################################################################

#Simulation function

#########################################################################################


#n = 500;N=5000;p=2;mu_0 = 1;g=1;g2=1;lambda=1;model = "QL";beta = c(1,1);num.sims=50

Sim_OTR <- function(n , N , p=2, model = "SL" ,lambda=1,beta=NULL,g=NULL,num.sims=500){
  
  theta = rep(c(1,-1),floor(p/2))
  #g = c(1,rep(1,(p/2)),rep(0,(p/2)))
  #if(sum(beta)>0){
  #  g = rep(-0.5,p)
  #}else{
  #  g = rep(0.5,p)
  #}
  
  MSE.ols = rep(0,num.sims)
  MSE.np1 = rep(0,num.sims)
  MSE.np2 = rep(0,num.sims)
  
  
  #Store beta coefficients.
  Beta_mat_OLS = matrix(0 , ncol=p+1 , nrow = num.sims)
  
  Beta_mat1 = matrix(0 , ncol=p+1 , nrow = num.sims)
  imput_est_mat1.CX = matrix(0 , ncol = 4 , nrow= num.sims)
  
  Beta_mat2 = matrix(0 , ncol=p+1 , nrow = num.sims)
  imput_est_mat2.CX = matrix(0 , ncol = 4 , nrow= num.sims)
  
  
  #Store bootstrap/IF sd
  IF.CXls_mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  
  IF.mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  boot.mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  b_CI_list = list()
  
  #CP
  CP_CX.mat = matrix(0 , ncol=p+1 , nrow = num.sims)
  CP_CX.mat2 = matrix(0 , ncol=p+1 , nrow = num.sims)
  CP_CX.mat3 = matrix(0 , ncol=p+1 , nrow = num.sims)
  
  #pcd
  pcd_sim.RF = rep(0,num.sims)
  pcd_sim.CXLS = rep(0,num.sims)
  pcd_sim.OLS = rep(0,num.sims)
  # Value functions
  VF.RF = rep(0,num.sims)
  VF.CXLS = rep(0,num.sims)
  VF.OLS = rep(0,num.sims)
  VF.true = rep(0,num.sims)
  
  
  #Get beta.t
  MC_list = MC_coeff(model = model, n.samp = 500000,p=p,
                     g=g,theta=theta,lambda=lambda,beta=beta)
  beta.t = MC_list$coefs
  print(beta.t)
  
  for(sim in 1:num.sims){
    
    X = mvrnorm(n=n+N , mu = rep(0,p), Sigma = diag(p))
    
    colnames(X) <- paste("X",1:p,sep="")
    
    #Impute incomplete cases with Kernel Smoothing and Random Forests.
    X.imp = as.data.frame(X[(n+1):(N+n),]); X.label = as.data.frame(X[1:n,])
    
    logit_coefs <- rep(1,ncol(X))
    true.prop = plogis( 0.5*X[,1] - 0.5*X[,2] )
    #true.prop = runif(n.samp,min = 0.4,max=0.6)
    #summary(true.prop)
    A = sapply(true.prop, function(x){ rbinom(1,1,x) })
    
    #Generate Y based on the given true model.
    if(model == "QL"){
      Y = as.vector((X%*%theta)^2 + A*lambda*(X%*%beta) + rnorm(n+N,0,1))
      CX = X%*%beta
    }else if(model == "L2L"){
      Y = as.vector( (X%*%beta)*(1+X%*%theta) + A*lambda*(X%*%beta) + rnorm(n+N,0,1))
      mu_X = X%*%beta ; CX = X%*%beta
    }else if(model == "L2NL3"){
      Y = as.vector((X%*%beta)*(1+X%*%theta) + A*lambda*((X%*%g)^3) + rnorm(n+N,0,1))
      CX = (X%*%g)^3
    }else if(model == "QNL3"){
      Y = as.vector((X%*%theta)^2 + A*lambda*((X%*%g)^3) + rnorm(n+N,0,1))
      CX = (X%*%g)^3
    }else if(model == "L2S"){
      Y = as.vector((X%*%beta)*(1+X%*%theta) + A*lambda*(sin(X%*%beta) ) + rnorm(n+N,0,1))
      CX = sin(X%*%beta)
    }else if(model == "QS"){
      Y = as.vector((X%*%theta)^2 + A*lambda*( sin(X%*%beta) ) + rnorm(n+N,0,1))
      CX = sin(X%*%beta)
    }
    
    #need these to get propensity score.
    X.known = X[1:n,]
    A.known = A[1:n]
    Y.known = Y[1:n]
    #estimated propensity score.
    prop_model <- glm(A.known ~ . ,family=binomial(link='logit'),data=as.data.frame(X.known))
    est.prop = prop_model$fitted.values
    prop = est.prop[1:n]
    
    #values for response variable in transformed outcome linear regression.
    Yt = Y.known*(A.known - prop) / (prop*(1 - prop))
    #Yt = Yt.true[1:n] #prop = true.prop[1:n];Nprop = true.prop[-(1:n)]
    
    # Estimate Y_A1 and Y_A0 separately. 
    Y_A1 = Y.known[which(A.known==1)];Y_A0 = Y.known[which(A.known==0)];prop0 = prop[which(A.known==0)]
    X.labelA1 = X.label[A.known==1,];X.labelA0 = X.label[A.known==0,];prop1 = prop[which(A.known==1)]
    
    ###################################################################################################
    # Random Forest approach
    
    #Predict missing Y.
    formula.YA1 = as.formula(paste("Y_A1~", paste(colnames(X.label), collapse="+")))
    formula.YA0 = as.formula(paste("Y_A0~", paste(colnames(X.label), collapse="+")))
    
    #Random Forest Model
    RF.model_A1 = randomForest(formula.YA1 , data = X.labelA1)
    RF.model_A0 = randomForest(formula.YA0 , data = X.labelA0)
    mx.hat_A1 = predict(RF.model_A1,newdata = X.imp)
    mx.hat_A0 = predict(RF.model_A0,newdata = X.imp)
    CX_RF = mx.hat_A1 - mx.hat_A0
    
    ##################################################################################################
    
    #Kernel Smoothing approaches.
    X.l_1 = cbind(rep(1,dim(X.labelA1)[1]),X.labelA1);X.l_0 = cbind(rep(1,dim(X.labelA0)[1]),X.labelA0)
    X.i_KS = cbind(rep(1,N),X.imp)
    X.l_KS = cbind(rep(1,n),X.label)
    #Get h by grid search
    if(sim==1){
      hc1 = cv.h(Y_A1, X.l_1,seq.c = seq(0.1,5,0.2)) #0.75
      hc0 = cv.h(Y_A0, X.l_0,seq.c = seq(0.1,5,0.2))
      #hc = cv.h(Yt, X.l_KS,seq.c = seq(0.1,5,0.2),prop_score = prop)
      #best.rf = cv.rf(Yt, X.label, num.trees=seq(500,1000,250),num.ns = seq(5,15,5) )
      print(hc1);print(hc0)
    }
    
    #Predictions with bias reduction by least squares
    mx.hat_ls1 = hlscv.ks(Y_A1,x.in= X.l_1,x.impute=X.i_KS,const_in=hc1, prop_score = prop1 )[[1]]
    mx.hat_ls0 = hlscv.ks(Y_A0 , X.l_0, X.i_KS, const = hc0, prop_score = 1 - prop0 )[[1]]
    CX_ls = mx.hat_ls1 - mx.hat_ls0
    
    #Get mx.hat for complete case for IF estimation.
    comp_mx.hat_ls1 = hlscv.ks(Y_A1,x.in= X.l_1,x.impute=X.l_KS,const_in=hc1 , prop_score = prop1)[[1]]
    comp_mx.hat_ls0 = hlscv.ks(Y_A0 , X.l_0, X.l_KS, const = hc0, prop_score = 1 - prop0 )[[1]]
    
    ####################################################################################################
    #Regression predicitions.
    formula.Yt = as.formula(paste("Yt~", paste(colnames(X.label), collapse="+")))
    Reg.OLS = lm(formula.Yt , data = X.label)
    beta.OLS = Reg.OLS$coefficients
    
    formula.CX_RF  = as.formula(paste("CX_RF ~", paste(colnames(X.label), collapse="+")))
    Reg.CX_RF = lm(formula.CX_RF  , data = X.imp)
    beta.CX_RF = Reg.CX_RF$coefficients
    
    formula.CX_ls  = as.formula(paste("CX_ls ~", paste(colnames(X.label), collapse="+")))
    Reg.CX_ls = lm(formula.CX_ls  , data = X.imp)
    beta.CX_ls = Reg.CX_ls$coefficients
    
    # Relative Efficiency
    MSE.ols[sim] = sum((beta.OLS - beta.t)^2)
    MSE.np1[sim] = sum((beta.CX_RF - beta.t)^2)
    MSE.np2[sim] = sum((beta.CX_ls - beta.t)^2)
    
    # Record regression coefficients
    Beta_mat_OLS[sim,] = beta.OLS
    Beta_mat1[sim,] = beta.CX_RF
    Beta_mat2[sim,] = beta.CX_ls
    #Beta_mat5[sim,] = beta.CX_preds
    
    #Record imputation value
    imput_est_mat1.CX[sim,] = c(Q1_bias=quantile(CX_RF)[2]-quantile(CX[-(1:n)])[2] ,mean_bias=mean(CX_RF) - mean(CX[-(1:n)]) , 
                                Q3_bias=quantile(CX_RF)[4]-quantile(CX[-(1:n)])[4],MSE=mean((CX_RF - CX[-(1:n)])^2))
    imput_est_mat2.CX[sim,] = c(Q1_bias=quantile(CX_ls)[2]-quantile(CX[-(1:n)])[2] ,mean_bias=mean(CX_ls) - mean(CX[-(1:n)]) , 
                                Q3_bias=quantile(CX_ls)[4]-quantile(CX[-(1:n)])[4],MSE=mean((CX_ls - CX[-(1:n)])^2))
    
    ##########################################################################
    # Variance Estimation
    # First get double cv mu_hat estimates
    dcv.mx1 = double_cv.ks(Y_A1,x.in= X.l_1,const_in=hc1, prop_score = prop1 )
    dcv.mx0 = double_cv.ks(Y_A0 , X.l_0, const = hc0, prop_score = 1 - prop0 )
    #u1.est = hlscv.ks(Y_A1,x.in= X.l_1,x.impute=X.l_0,prop_score = prop1, const_in=hc1 )[[1]]
    #u0.est = hlscv.ks(Y_A0 , X.l_0, x.impute=X.l_1,prop_score = 1 - prop0, const = hc0 )[[1]]
    #
    
    full_X.imp = as.matrix(cbind(rep(1,N),X.imp))
    names(X.l_1) <- names(X.l_0) 
    full_X.label = as.matrix(rbind(X.l_1,X.l_0))
    
    GammaInv = solve(t(full_X.imp)%*%full_X.imp/N)
    
    ### IF sd 3
    #YA = c(Y_A1,Y_A0);mu = c(dcv.mx1,dcv.mx0);A_all = c(rep(1,length(Y_A1)),rep(0,length(Y_A0)))
    
    X.in = rbind( as.matrix(X.l_1) , as.matrix(X.l_0) );Y.in = c(Y_A1,Y_A0);prop.in = c(prop1,prop0)
    mv.input = c(dcv.mx1,dcv.mx0); A.in = c( rep(1,length(Y_A1)) , rep(0,length(Y_A0)) ) 
    ####################################################################################################
    # Combine these estimates to put in influence function.
    #u1 = c(dcv.mx1,u1.est);u0 = c(u0.est,dcv.mx0)
    #dcx = u1 - u0
    # Get double cv beta est.
    #formula.dcv_b  = as.formula(paste("dcx ~", paste(colnames(X.label), collapse="+")))
    #Reg.dcvb = lm(formula.dcv_b  , data = as.data.frame(X.in[,-1]))
    #beta.dcvb = Reg.dcvb$coefficients
    
    ####################################################################################################
    
    IF.CXls_sd = IF_sd(GammaInv, X.in, Y.in, A.in, mv.input,prop.in )
    IF.CXls_mat[sim,] = IF.CXls_sd
    
    ##############################################################
    #Coverage probability
    #OLS
    z = qnorm(0.975)
    #CX_ls RAL estimator
    lower.cx = beta.CX_ls - z*IF.CXls_sd
    upper.cx = beta.CX_ls + z*IF.CXls_sd
    CP_CX.mat[sim,] = 1*((lower.cx <= beta.t) & (beta.t <= upper.cx))
    ###############################################################
    #m-out-of-n Bootstrap
    #if(sim==1){
    #  m.val = get_m_MCAR(X.l=X.label,Y=Y.known,A=A.known,X.U=X.imp,
    #                            numB=500, hc1=hc1, hc0=hc0, beta_est = beta.CX_ls)
    #}
    #bVar = mn_bstrap_var(X.l=X.label,Y=Y.known,A=A.known,X.U=X.imp,
    #                     numB=500, hc1=hc1, hc0=hc0, m = 350)
    #n out of n bstrap
    bVar = bstrap_var( X.l=X.label,Y.b=Y.known,A.b=A.known,X.U=X.imp,
                       numB=250, hc1=hc1, hc0=hc0 )
    
    #print(bVar)
    boot.mat[sim,] = bVar[[1]]
    #CX_ls RAL estimator
    lower.cx = beta.CX_ls - z*bVar[[1]]
    upper.cx = beta.CX_ls + z*bVar[[1]]
    CP_CX.mat2[sim,] = 1*((lower.cx <= beta.t) & (beta.t <= upper.cx))
    
    b_CI_list[[sim]] = bVar[[2]] 
    CP_CX.mat3[sim,] = 1*((bVar[[2]][1,] <= beta.t) & (beta.t <= bVar[[2]][2,]))
    ###############################################################
    #pcd for the current simulation.
    int_X = cbind(rep(1,dim(X)[1]),X)
    BX_RF = beta.CX_RF%*%t(int_X)
    BX_CXLS = beta.CX_ls%*%t(int_X)
    BX_OLS = beta.OLS%*%t(int_X)
    BX_true = beta.t%*%t(int_X)
    
    Ind_RF = sapply(BX_RF,function(x)ifelse((x>0),1,0))
    Ind_CXLS = sapply(BX_CXLS,function(x)ifelse((x>0),1,0))
    Ind_OLS = sapply(BX_OLS,function(x)ifelse((x>0),1,0))
    Ind_true = sapply(BX_true,function(x)ifelse((x>0),1,0))
    
    pcd_sim.RF[sim] = 1 - mean(abs(Ind_RF - Ind_true))
    pcd_sim.CXLS[sim] = 1 - mean(abs(Ind_CXLS - Ind_true))
    pcd_sim.OLS[sim] = 1 - mean(abs(Ind_OLS - Ind_true))
    
    #VF 
    X_MC = MC_list$X
    int_MCX = cbind(rep(1,dim(X_MC)[1]),X_MC)
    MC_RF = beta.CX_RF%*%t(int_MCX)
    MC_CXLS = beta.CX_ls%*%t(int_MCX)
    MC_OLS = beta.OLS%*%t(int_MCX)
    
    
    Ind_MC_RF = sapply(MC_RF,function(x)ifelse((x>0),1,0))
    Ind_MC_CXLS = sapply(MC_CXLS,function(x)ifelse((x>0),1,0))
    Ind_MC_OLS = sapply(MC_OLS,function(x)ifelse((x>0),1,0))
    
    VF.RF[sim] = mean(as.vector(MC_list$ux + Ind_MC_RF*lambda*MC_list$CX))
    VF.CXLS[sim] = mean(as.vector(MC_list$ux + Ind_MC_CXLS*lambda*MC_list$CX))
    VF.OLS[sim] = mean(as.vector(MC_list$ux + Ind_MC_OLS*lambda*MC_list$CX))
  }
  
  #pcd
  pcd.RF.mean = mean(pcd_sim.RF);pcd.RF.sd = sd(pcd_sim.RF)
  pcd.CXLS.mean = mean(pcd_sim.CXLS);pcd.CXLS.sd = sd(pcd_sim.CXLS)
  pcd.OLS.mean = mean(pcd_sim.OLS);pcd.OLS.sd = sd(pcd_sim.OLS)
  #VF
  VF.RF.mean = mean(VF.RF);VF.RF.sd = sd(VF.RF)
  VF.CXLS.mean = mean(VF.CXLS);VF.CXLS.sd = sd(VF.CXLS)
  VF.OLS.mean = mean(VF.OLS);VF.OLS.sd = sd(VF.OLS )
  
  #mean
  betaOLS.mean = apply(Beta_mat_OLS,2,mean)
  beta.CX_RF.mean = apply(Beta_mat1,2,mean)
  beta.CX_ls.mean = apply(Beta_mat2,2,mean)
  
  #Empirical sd
  betaOLS.sd = apply(Beta_mat_OLS,2,sd)
  beta.CX_RF.sd = apply(Beta_mat1,2,sd)
  beta.CX_ls.sd = apply(Beta_mat2,2,sd)
  
  #Average IF/boot sd
  IF.CXls_avg = apply(IF.CXls_mat,2,mean)
  boot.CXls_avg = apply(boot.mat,2,mean)
  
  #CP
  CP_CX = apply(CP_CX.mat,2,mean)
  CP_CX_boot = apply(CP_CX.mat2,2,mean)
  CP_CX_boot_CI = apply(CP_CX.mat3,2,mean)
  #Relative Efficiency
  RE.RF = mean(MSE.ols)/mean(MSE.np1)
  
  RE.CX_ls = mean(MSE.ols) / mean(MSE.np2)
  
  ##########################################################################################
  
  RF_2_imput = as.matrix(apply(imput_est_mat1.CX,2,mean),ncol=4)
  rownames(RF_2_imput)<-c("Q1_bias","mean_bias","Q3_bias","MSE")
  CXls_imput = as.matrix(apply(imput_est_mat2.CX,2,mean),ncol=4)
  rownames(CXls_imput)<-c("Q1_bias","mean_bias","Q3_bias","MSE")
  
  RE = list( RE_RF = RE.RF,
             RE_CXls = RE.CX_ls,
             RF_2_imput = RF_2_imput,
             CXls_imput = CXls_imput,
             bias_sd_OLS = cbind(beta.t,bias=betaOLS.mean-beta.t,Esd=betaOLS.sd),
             bias_sd_RF = cbind(beta.t,bias=beta.CX_RF.mean-beta.t ,sd=beta.CX_RF.sd),
             bias_sd_CXls = cbind(beta.t,bias=beta.CX_ls.mean-beta.t ,Esd=beta.CX_ls.sd,
                                  Asd_IF = IF.CXls_avg,CP=CP_CX,
                                  boot_se = boot.CXls_avg,CP_boot = CP_CX_boot, CP_CI_boot = CP_CX_boot_CI ),
             pcd_OLS = c(pcd.OLS.mean,pcd.OLS.sd),
             VF_OLS = c(VF.OLS.mean,VF.OLS.sd),
             pcd_RF = c(pcd.RF.mean, pcd.RF.sd),
             VF_RF = c(VF.RF.mean,VF.RF.sd),
             pcd_CXls = c(pcd.CXLS.mean, pcd.CXLS.sd),
             VF_Cxls = c(VF.CXLS.mean,VF.CXLS.sd),
             VF_0 = c(mean_V0 = mean(MC_list$V_0),sd_V0 = sd(MC_list$V_0)),
             VF_B0 = c(mean_VB0 = mean(MC_list$V_b0),sd_VB0 = sd(MC_list$V_b0)),
             PCD_vals_RF = pcd_sim.RF,
             PCD_vals_CXLS = pcd_sim.CXLS,
             PCD_vals_oLS = pcd_sim.OLS,
             V_RF = VF.RF,
             V_CXLS = VF.CXLS,
             V_OLS = VF.OLS,
             Beta_CX = Beta_mat1,
             b_CI_list
             )
  
  return(RE)
}

#p=10
beta.in = rep(1,2);g.in = c(0.3,0.6); g_large = c(0.7,0.7)
###################################################################################
#Sim_L2L = Sim_OTR(n=500,N=5000, p=2, model = "L2L",beta = beta.in,num.sims = 500)
#print(Sim_LL)

Sim_QL = Sim_OTR(n=500,N=5000, p=2, model = "QL",beta = beta.in,num.sims=500)

#Sim_L2NL3 = Sim_OTR(n=500,N=5000, p=2, model = "L2NL3",beta=beta.in,g = g.in,num.sims=500)

#Sim_QNL3 = Sim_OTR(n=500,N=5000, p=2, model = "QNL3",beta=NULL,g = g.in,num.sims= 500)

#Sim_L2S = Sim_OTR(n=500,N=5000, p=2, model = "L2S",beta=beta.in, num.sims=500)

#Sim_QS = Sim_OTR(n=500,N=5000, p=2, model = "QS",beta=beta.in, num.sims=500)

#capture.output(Sim_LL,file="VAR_MCAR_SL.txt")
capture.output(Sim_QL,file="bstrap_estprop_MCAR_QL.txt")
#capture.output(Sim_LNL3,file="VAR_MCAR_SNL3.txt")
#capture.output(Sim_QNL3,file="mnbstrap_estprop_MCAR_QNL3.txt")
#capture.output(Sim_LS,file="VAR_MCAR_SS.txt")
#capture.output(Sim_QS,file="VAR_MCAR_QS.txt")

