# density function
dbxii<-function(y,mu,c,tau=.5)
{
  s=1
  d<-log(1/(1-tau))*(c*y^(c-1))/(s^(c)*log(1+(mu/s)^c))*
    (1+(y/s)^c)^(log(1-tau)/(log(1+(mu/s)^c))-1)
  d
}

# cumulative distribution function
pbxii<-function(y,mu,c,tau=.5)
{
  s=1
  p<- 1-(1+(y/s)^c)^(log(1-tau)/(log(1+(mu/s)^c)))
  p
}

# quantile function
qbxii<-function(u,mu,c,tau=.5)
{
  s=1
  q<- s*((1-u)^(log(1+(mu/s)^c)/log(1-tau))-1)^(1/c)
  q
}

# inversion method for randon generation
rbxii<-function(n,mu,c,tau=.5)
{
  s=1
  u<- runif(n)
  y<- s*((1-u)^(log(1+(mu/s)^c)/log(1-tau))-1)^(1/c)
  y
}

