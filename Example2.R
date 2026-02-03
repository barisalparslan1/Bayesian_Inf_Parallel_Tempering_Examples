library(ggplot2)
library(patchwork)

# True values
alpha= 3
beta = 2
epsilon = rnorm(1,0,1)

set.seed(123)
n_data <- 500       

x <- rnorm(n_data, mean = (alpha^2)*beta, sd=1)

# Uncoupled MCMC
alpha_init = -2   
beta_init = -5 

tune = 0.1
iters = 1e5


chain=function(x,tune,th1_init,th2_init,iters)
{
  x_alpha = th1_init
  x_beta = th2_init

  xvec=matrix(NA,nrow=iters,ncol=2)
  for (i in 1:iters) {
    print(i)
    can_alpha= x_alpha + rnorm(1,0,tune)
    can_beta= x_beta + rnorm(1,0,tune)
    
    
    log_numA = sum(dnorm(x,mean=(can_alpha^2)*beta,1)) -log(20)-log(10)

    
    log_denomA = sum(dnorm(x,mean=(can_alpha^2)*beta,1)) -log(20)-log(10)
    
    logA = log_numA - log_denomA
    
    if (log(runif(1))<logA){
      x_alpha = can_alpha
      x_beta = can_beta    
    }
    
    xvec[i, 1] = x_alpha; xvec[i,2] = x_beta;
  }
  xvec
}


samps <- chain(x,tune,alpha_init,beta_init,iters)

samps_df <- data.frame(samps)
samps_true <- data.frame(samps[(iters*0.95):iters,])

p1 <- ggplot(data=samps_true,aes(seq((iters*0.95):iters),samps_true[,1]))+geom_line()+labs(title="alpha")+xlab("iters")+ylab("alpha")
p2 <- ggplot(data=samps_true,aes(seq((iters*0.95):iters),samps_true[,2]))+geom_line()+labs(title="beta")+xlab("iters")+ylab("beta")
p1+p2



# Coupled-MCMC

swap_Switch = 1 
alpha= 3
beta = 2
epsilon = rnorm(1,0,1)

set.seed(123)
n_data <- 500       

x <- rnorm(n_data, mean = (alpha^2)*beta, sd=1)

alpha_init = -2   
beta_init = -5 

tune = 0.1
iters = 1e5


n = 5
col1 = paste("alpha_",seq(1:n),sep="")
col2 = paste("beta_",seq(1:n),sep="")

heats = c(0.1,0.3,0.5,0.7,1)
  
tune = 0.1
iters = 1e5


log_Post <- function(x,alpha,beta){
  
  log_numA =  (sum(dnorm(x,mean=(alpha^2)*beta,1)) -log(20)-log(10))

  
  return (log_numA)
  
}

log_Post_heated <- function(x,alpha,beta,heat){
  
  log_numA = heat*(sum(dnorm(x,mean=(alpha^2)*beta,1)) -log(20)-log(10))  
  
  return (log_numA)
  
}

chain_coupled=function(x,tune,alpha_init,beta_init,iters,n,heats,swap_Switch)
{
  x_alpha = rep(alpha_init,n)
  x_beta = rep(beta_init,n)
  
  xmat = matrix(0,nrow=iters,ncol=n*2)
  for (i in 1:iters) {
    print(i)
    can_alpha= x_alpha + rnorm(n,0,tune)
    can_beta= x_beta + rnorm(n,0,tune)
    
    
    log_numA = unlist(Map(function(a, b, c) log_Post_heated(x, a, b,c),
                          can_alpha, can_beta,heats))
    
    
    log_denomA = unlist(Map(function(a, b, c) log_Post_heated(x, a, b, c),
                            x_alpha, x_beta,heats)) 
    
    logA = log_numA - log_denomA
    
    accept = (log(runif(n))<logA) 
    
    x_alpha[accept] = can_alpha[accept]; x_beta[accept] = can_beta[accept];
    
    if (swap_Switch){
      
      swap = sample(1:n,2)
      
      logA = heats[swap[1]]*log_Post(x,x_alpha[swap[2]],x_beta[swap[2]]) + heats[swap[2]]*log_Post(x,x_alpha[swap[1]],x_beta[swap[1]])-
        heats[swap[1]]*log_Post(x,x_alpha[swap[1]],x_beta[swap[1]]) - heats[swap[2]]*log_Post(x,x_alpha[swap[2]],x_beta[swap[2]])
      
      if (log(runif(1))<logA){
        x_alpha[swap] = rev(x_alpha[swap])
        x_beta[swap] = rev(x_beta[swap])
      }
      
    }
    
    xmat[i,] = c(x_alpha,x_beta)
  }
  
  x.df <- data.frame(xmat)
  
  names(x.df)<- c(col1,col2)
  
  x.df
}


samps_coupled <- chain_coupled(x,tune,alpha_init,beta_init,iters,n,heats,swap_Switch)



samps_conv <- samps_coupled[(iters*0.8):iters,]



p1 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,5]))+geom_line()+labs(title="alpha")+xlab("iters")+ylab("alpha")
p2 <- ggplot(data=samps_conv,aes(seq((iters*0.8):iters),samps_conv[,10]))+geom_line()+labs(title="beta")+xlab("iters")+ylab("beta")
print(p1+p2)


