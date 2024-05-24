rm(list = ls())
source("bxiiarma.fit.r")
source("simu_bxiiarma.R")

library(beepr)

R <- 1000
vn <- c(
  #
 500
  # ,150,
  # 300
  # ,
  # 500
        )
# n=500

# alpha <- .05
# phi<-c(.04,.01)
# theta <- c(.02,.05)
# c<-2


## CENÁRIO CERTO
# alpha <- -1
# # phi<-c(.1,.4)
# theta <- c(0.3,0.2)
# c<-3

# alpha <- .5
# phi<-c(.27,-.6)
# theta <- c(0.5,.27)
# c<-2

alpha <- 1
# phi<-c(.5)
theta <- c(0.2)
c<-.5
model <- c(alpha,
             # phi,
            theta,
           c)

# p1<- 1:length(phi)
q1<- 1:length(theta)

set.seed(2024)
tempo.inicio = Sys.time()
for(n in vn){
# n=70  
# For the confidence intervals
alpha_erro <-0.05
quantil <-1-alpha_erro/2
z<-qnorm(quantil)
t<-qt(quantil, n-1)
  
# To save the results
estim <-ICi<-ICs<- err <- matrix(NA, nrow = R, ncol = length(model))
calpha<-cphi1<-cphi2<-ctheta1<-ctheta2<-cc<-0
i<-0


### simulation
while(i<R) {
  y <- simu.bxiiarma(n,
                      # phi=phi,
                     theta=theta,
                     alpha=alpha,c=c,tau = 0.5)
  # print(y)# hist(y)
  result <- try(bxiiarma.fit(y,
                              # ar=p1,
                              ma=q1,
                             mod=2,tau = 0.5,diag1=1), silent = T)
  if(class(result) == "try-error" )
  {
    result$conv<-1
    convergencia<-1
  }
  
  if(result$conv == 0)
  {
    convergencia<-  result$conv
  }else{
    
    warning("FUNCTION DID NOT CONVERGE!")
  }
  
  if(convergencia==0 && sum(is.na(result$stderror))==0)
  {
    i<-i+1
    print(c("i=",i))
    estim[i,] <- result$model[,1]
    err[i,] <- result$model[,2]
    ICi[i,]<-estim[i,]-(z*err[i,])
    ICs[i,]<-estim[i,]+ (z*err[i,])
    
    if (ICi[i,1]<=alpha && ICs[i,1]>=alpha)
    {
      calpha<-calpha+1
    }
    
    # if (ICi[i,2]<= phi[1] && ICs[i,2]>=phi[1])
    # {
    #   cphi1<-cphi1+1
    # }
    # if (ICi[i,3]<= phi[2] && ICs[i,3]>=phi[2])
    # {
    #   cphi2<-cphi2+1
    # }
    if (ICi[i,2]<= theta[1] && ICs[i,2]>=theta[1])
    {
      ctheta1<-ctheta1+1
    }
    # if (ICi[i,3]<= theta[2] && ICs[i,3]>=theta[2])
    # {
    #   ctheta2<-ctheta2+1
    # }
    if (ICi[i,3]<= c && ICs[i,3]>=c)
    {
      cc<-cc+1
    }
  } # fim convergencia MC
} #fim loop MC

### mean
m <- apply(estim, 2, mean)
### bias
bias <- (model-m)
### relative percentage bias
biasP <- bias/model *100
### SD
erro <- apply(estim, 2, sd)
### MSE
MSE <- apply(estim, 2, var)+bias^2
### 
TC<-c(calpha,
      # cphi1,cphi2,
      ctheta1,
      # ctheta2,
      cc)/R
## final results
results <- rbind(m, bias, biasP, erro, MSE,TC)
rownames(results) <- c("Mean", "Bias","RB%", "SE", "MSE","TC")
print(c("Tamanho da Amostra:",n))
print(round(results,4))
}
mnjtempo.fim = Sys.time()
tempo.exec = mnjtempo.fim-tempo.inicio
print(tempo.exec)
beep(8)



## TC1 gráfico 
postscript(file = "ACR",width = 6, height = 4,family = "Times")
# par(mar=c(5,6,4,1)+.1)
par(mfrow=c(1,1))
par(mar=c(2.8, 2.5, 0.5, 0.5))
par(mgp=c(1.3, 0.45, 0))

ARMA11<-c(0.9356,0.9417,0.9420,0.9430,0.9420)
ARMA10<-c(0.9478,0.9491,0.9481,0.9480,0.9489)
ARMA01<-c(0.9392,0.9455,0.9488,0.9482,0.9488)

plot(ARMA11,xaxs="r",type="l",
     ylab=expression("ACR"),
     xlab=expression(plain(italic(n))),
     axes=FALSE,
     ylim=c(0.92,0.97))

lines(ARMA11, lty=1,col=1)
lines(ARMA10, lty=1,col=7)
lines(ARMA01, lty=1,col=4)

points(ARMA11,pch=1,col=1)
points(ARMA10,pch=2,col=7)
points(ARMA01,pch=3,col=4)

abline(h=0.95,col=2,lty=3)

legend("top",c(expression(plain(BXII-ARMA(1,1))),
                    expression(plain(BXII-ARMA(1,0))), 
                    # expression(plain(theta[1])),
                    expression(BXII-ARMA(0,1))),
       lty=c(1,1,1), pch=c(1,2,4),col=c(1,2,4),bty="n",cex=1)
axis(1,1:5,c(70,150,300,500,1000))
box();axis(2)
dev.off()  


## TC2 gráfico 
postscript(file = "TC2",width = 6, height = 4,family = "Times")
# par(mar=c(5,6,4,1)+.1)
par(mfrow=c(1,1))
par(mar=c(2.8, 2.5, 0.5, 0.5))
par(mgp=c(1.3, 0.45, 0))

TC1<-c(0.9462,0.9494,0.9475,0.9480,0.9477)
TC2<-c(0.9389,0.9457,0.9446,0.9464,0.9468)
TC4<-c(0.9584,0.9521,0.9522,0.9496,0.9522)


plot(TC1,xaxs="r",type="l",
     ylab=expression("TC95%"),
     xlab=expression(plain(italic(n))),
     axes=FALSE,
     ylim=c(min(c(TC1,TC2,#TC3,
                  TC4)),max(c(TC1,TC2,#TC3,
                              TC4))))

lines(TC1, lty=1,col=1)
lines(TC2, lty=1,col=7)
# lines(TC3, lty=1,col=3)
lines(TC4, lty=1,col=4)

points(TC1,pch=1,col=1)
points(TC2,pch=2,col=7)
# points(TC3,pch=3,col=3)
points(TC4,pch=4,col=4)
abline(h=0.95,col=2,lty=3)


legend("bottomright",c(expression(plain(alpha)),
                       expression(plain(phi[1])), 
                       # expression(plain(theta[1])),
                       expression(plain(c))),
       lty=c(1,1,1), pch=c(1,2,4),col=c(1,2,4),bty="n",cex=1.25)
axis(1,1:5,c(70,150,300,500,1000))
box();axis(2)
dev.off()  

