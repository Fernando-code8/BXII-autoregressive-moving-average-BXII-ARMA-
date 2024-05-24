predict_bxii<-function(fit_bxii){
  alpha<-fit_bxii$coeff[1]
  phi<-fit_bxii$phi
  theta<-fit_bxii$theta
  #### Forecasting
  h1<-length(datatest)
  ar<-1:length(phi)
  ma<-1:length(theta)
  m<-max(length(ma),length(ar)) ## parei aqui
  y_t<-fit_bxii$stats$linkfun(c(yy[(n-m):(n-2)],dadosh))
  errorhat<- c(fit_bxii$errorhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <-y_prev <- c()
  
  for(i in 1:h1){
    ynew_prev[i] <- alpha  + (phi%*%y_t[m+i-ar])+(theta%*%errorhat[m+i-ma])
    errorhat[m+i]<- y_t[m+i]-ynew_prev[i] # predictor scale
    y_prev[i] <- fit_bxii$stats$linkinv(ynew_prev[i])
  }
  return(y_prev)
}

predict_bxii_wind<-function(fit_bxii){
  alpha<-fit_bxii$coeff[1]
  phi<-fit_bxii$phi
  theta<-fit_bxii$theta
  #### Forecasting
  h1<-length(datatest)
  ar<-1:length(phi)
  ma<-1:length(theta)
  m<-max(length(ma)) ## parei aqui
  y_t<-fit_bxii$stats$linkfun(c(yy[(n-m):(n-2)],dadosh))
  errorhat<- c(fit_bxii$errorhat[(n-m):(n-1)],rep(NA,h1))
  ynew_prev <-y_prev <- c()
  
  for(i in 1:h1){
    ynew_prev[i] <- alpha +(phi%*%y_t[m+i-ar])+(theta%*%errorhat[m+i-ma])
    errorhat[m+i]<- y_t[m+i]-ynew_prev[i] # predictor scale
    y_prev[i] <- fit_bxii$stats$linkinv(ynew_prev[i])
  }
  return(y_prev)
}