library(MCMCpack)
library(ggplot2)
library(patchwork)
 

# True values
th1 = 6
th2 =-6
sigma = 1

set.seed(123)
n_data <- 500       

z <- rbinom(n_data, size = 1, prob = 0.5)
x <- rnorm(n_data, mean = ifelse(z == 1, th1, th2),sd=1)

hist(x,breaks=70)

# Uncoupled MCMC
th1_init =  8   # prior: N(0,100)
th2_init = 8    # prior: N(0,100)

tune = 2
iters = 1e5


chain=function(x,tune,th1_init,th2_init,sigma,iters)
{
  x_th1 = th1_init
  x_th2 = th2_init

  xvec=matrix(NA,nrow=iters,ncol=2)
  for (i in 1:iters) {
    print(i)
    can_th1= x_th1 + rnorm(1,0,tune)
    can_th2= x_th2 + rnorm(1,0,tune)
    
  
    log_numA = sum(log(dnorm(x,can_th1,sigma)/2 + dnorm(x,can_th2,sigma)/2)) + dnorm(can_th1,0,100,log=TRUE) +
               dnorm(can_th2,0,100,log=TRUE)   
    
    
    log_denomA = sum(log(dnorm(x,x_th1,sigma)/2 + dnorm(x,x_th2,sigma)/2)) + dnorm(x_th1,0,100,log=TRUE) +
      dnorm(x_th2,0,100,log=TRUE) 
    
    logA = log_numA - log_denomA
    
    if (log(runif(1))<logA){
      x_th1 = can_th1
      x_th2 = can_th2    
      }
    
    xvec[i, 1] = x_th1; xvec[i,2] = x_th2;
  }
  xvec
}


samps <- chain(x,tune,th1_init,th2_init,sigma,iters)

samps_df <- data.frame(samps)
samps_true <- data.frame(samps[(iters*0.95):iters,])

p1 <- ggplot(data=samps_true,aes(seq((iters*0.95):iters),samps_true[,1]))+geom_line()+labs(title="th1")+xlab("iterations")+ylab("value")
p2 <- ggplot(data=samps_true,aes(seq((iters*0.95):iters),samps_true[,2]))+geom_line()+labs(title="th2")+xlab("iterations")+ylab("value")
p1+p2


z_uncoupled <- rbinom(n_data,1,0.5)
x_uncoupled <- rnorm(n_data, mean = ifelse(z == 1, mean(samps_true[,1]), mean(samps_true[,2])), sd = sigma)
hist(x,breaks=50)


# Coupled MCMC

swap_Switch = 1
set.seed(123)
n_data <- 500  
th1 = 6
th2 =-6
sigma =1

z <- rbinom(n_data, size = 1, prob = 0.5)
x <- rnorm(n_data, mean = ifelse(z == 1, th1, th2), sd = sigma)

n = 5
col1 = paste("th1_",seq(1:n),sep="")
col2 = paste("th2_",seq(1:n),sep="")

heats = c(0.1,0.3,0.5,0.7,1)

th1_init = 8    # prior: N(0,10)
th2_init = 8   # prior: N(0,10)
#sigma_init = 5  # prior: Inv-Gamma(0.001,0.001)
tune = 2
iters = 1e5


log_Post <- function(x,t1,t2,sigma,beta){
  
  log_numA = sum(log(dnorm(x,t1,sigma)/2 + dnorm(x,t2,sigma)/2)) + dnorm(t1,0,100,log=TRUE) +
    dnorm(t2,0,100,log=TRUE)   
  
  return (log_numA)
  
}

log_Post_heated <- function(x,t1,t2,sigma,beta){
  
  log_numA = beta*(sum(log(dnorm(x,t1,sigma)/2 + dnorm(x,t2,sigma)/2)) + dnorm(t1,0,100,log=TRUE) +
    dnorm(t2,0,100,log=TRUE))   
  
  return (log_numA)
  
}

chain_coupled=function(x,tune,th1_init,th2_init,sigma,iters,temps,n,heats,swap_Switch)
{
  x_th1 = rep(th1_init,n)
  x_th2 = rep(th2_init,n)
  
  xmat = matrix(0,nrow=iters,ncol=n*2)
  for (i in 1:iters) {
    print(i)
    can_th1= x_th1 + rnorm(n,0,tune)
    can_th2= x_th2 + rnorm(n,0,tune)
    
    
    log_numA = unlist(Map(function(a, b, d) log_Post_heated(x, a, b, sigma, d),
                   can_th1, can_th2,heats))
    
    
    log_denomA = unlist(Map(function(a, b, d) log_Post_heated(x, a, b, sigma, d),
                                  x_th1, x_th2,heats)) 
    
    logA = log_numA - log_denomA
    
    accept = (log(runif(n))<logA) 
    
    x_th1[accept] = can_th1[accept]; x_th2[accept] = can_th2[accept];
    
    if (swap_Switch){
      
      swap = sample(1:n,2)
      
      logA = heats[swap[1]]*log_Post(x,x_th1[swap[2]],x_th2[swap[2]],sigma) + heats[swap[2]]*log_Post(x,x_th1[swap[1]],x_th2[swap[1]],sigma)-
        heats[swap[1]]*log_Post(x,x_th1[swap[1]],x_th2[swap[1]],sigma) - heats[swap[2]]*log_Post(x,x_th1[swap[2]],x_th2[swap[2]],sigma)
      
      if (log(runif(1))<logA){
        x_th1[swap] = rev(x_th1[swap])
        x_th2[swap] = rev(x_th2[swap])
      }
      
    }
    
    xmat[i,] = c(x_th1,x_th2)
  }
  
  x.df <- data.frame(xmat)
  
  names(x.df)<- c(col1,col2)
  
  x.df
}


samps_coupled <- chain_coupled(x,tune,th1_init,th2_init,sigma,iters,temps,n,heats,swap_Switch)



samps_conv <- samps_coupled[(iters*0.8):iters,]



p1 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,5]))+geom_line()+labs(title="th1")+xlab("iters")+ylab("th1")
p2 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,10]))+geom_line()+labs(title="th2")+xlab("iters")+ylab("th2")
p1+p2 + plot_annotation(
  title = "Chain 5"
)


p1 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,3]))+geom_line()+labs(title="th1")+xlab("iters")+ylab("th1")
p2 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,8]))+geom_line()+labs(title="th2")+xlab("iters")+ylab("th2")
p1+p2+ plot_annotation(
    title = "Chain 3"
  )


p1 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,5]))+geom_line()+labs(title="th1")+xlab("iters")+ylab("th1")
p2 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,10]))+geom_line()+labs(title="th2")+xlab("iters")+ylab("th2")
p1+p2


z1 <- rbinom(n_data, size = 1, prob = 0.5)
x1 <- rnorm(n_data, mean = ifelse(z1 == 1, mean(samps_conv[,1]), mean(samps_conv[,6])), sd = sigma)

hist(x1,breaks=50)

z2 <- rbinom(n_data, size = 1, prob = 0.5)
x2 <- rnorm(n_data, mean = ifelse(z2 == 1, mean(samps_conv[,3]), mean(samps_conv[,8])), sd = sigma)

hist(x2,breaks=50)

z3 <- rbinom(n_data, size = 1, prob = 0.5)
x3 <- rnorm(n_data, mean = ifelse(z3 == 1, mean(samps_conv[,5]), mean(samps_conv[,10])), sd = sigma)

hist(x3,breaks=50)

df.dist <- data.frame(ch1 = x1, ch2=x2, ch3=x3)

p1 <- ggplot(data=df.dist,aes(ch1))+geom_histogram()+labs(title="chain 1")
p2 <- ggplot(data=df.dist,aes(ch2))+geom_histogram()+labs(title="chain 2")
p3 <- ggplot(data=df.dist,aes(ch3))+geom_histogram()+labs(title="chain 3")

p1+p2+p3

