################################################################################################
###################### sim HIV  ##########################################################
################################################################################################
################################################################################################
source('log/ODE_log_growth.R')
require(deSolve)
require(fda)
require(CollocInfer)

log_growth<-function(t,r,K,x0){
  K*x0*exp(r*t)/(K+x0*(exp(r*t)-1))
}


##### CollocInfer - HIV

make.logfn<-function(){
  
  ### set up functions for 
  ### dx = r*x - r/K x^2
  
  logfun<-function(times,x,p,more){
    n=length(times);d=1;npar=2
    r=p[1]
    K=p[2]
    f=matrix(0,n,d)
    f[,1]=r*x[,1]-(r/K)*(x[,1]^2)
    return(f)
  }
  
  
  logdfdx<-function(times,x,p,more){
    n=length(times);d=1;npar=2
    r=p[1]
    K=p[2]
    dfdx=array(0,c(n,d,d))
    dfdx[,1,1]=r-2*(r/K)*x[,1]
    return(dfdx)
  }
  
  ### set up functions for 
  ### dx = r*x - r/K x^2
  
  
  logdfdp<-function(times,x,p,more){
    n=length(times);d=1;npar=2
    r=p[1]
    K=p[2]
    dfdp=array(0,c(n,d,npar))
    dfdp[,1,1]=x[,1]-(x[,1]^2)/K
    dfdp[,1,2]=r*(x[,1]^2)*K^(-2)
    return(dfdp)
  }
  
  
  logd2fdx2<-function(times,x,p,more){
    n=length(times);d=1;npar=2
    r=p[1]
    K=p[2]
    d2fdx2=array(0,c(n,d,d,d))
    
    d2fdx2[,1,1,1]=rep(-2*(r/K),n)
    
    return(d2fdx2)
  }
  
  logd2fdxdp<-function(times,x,p,more){
    n=length(times);d=1;npar=2
    r=p[1]
    K=p[2]
    d2fdxdp=array(0,c(n,d,d,npar))
    
    r-2*(r/K)*x[,1]
    
    
    d2fdxdp[,1,1,1]=rep(1,n)-2*K*x[,1]
    d2fdxdp[,1,1,2]=2*r*x[,1]*K^(-2)
    
    return(d2fdxdp)
  }
  
  return(list(fn=logfun, dfdx=logdfdx,
              dfdp=logdfdp, d2fdx2=logd2fdx2,
              d2fdxdp=logd2fdxdp))
  
}

true.y=matrix(0,100,500)
true.sd1=rep(0,100)
burn.in=5000
nit=10000
noise.coef=0.05

log.fn=make.logfn()
log.knots=seq(0,99,99/1)
log.order=3
log.nbasis=length(log.knots)+log.order-2
log.range=c(0,99)
log.basis=create.bspline.basis(log.range,log.nbasis,log.order,log.knots)

colloc_par=matrix(0,100,2)

medie.teta=list()
dev.teta=list()
medie.beta=list()
dev.beta=list()
media.sigma_d=list()
dev.sigma_d=list()
media.mu=list()
dev.mu=list()
media.post_x=list()
media.post_g1=list()

for(ii in 1:100){
  set.seed(ii)
  
  t_obs=(0:499)/500
  true.y[ii,]=as.matrix(log_growth(t_obs,r=2.5,K=20,x0=0.1))
  
  true.sd1[ii]=mean(true.y[ii,])*noise.coef
  y=true.y[ii,]+rnorm(length(true.y[ii,]),0,sd=true.sd1[ii])
  y=as.matrix(y)
  q=2
  
  mod=ode.mcmc_C(y,t_obs,q=q,nit = nit, burn.in = burn.in, progr=FALSE)
  
  log.xfd=smooth.basis(t_obs,y,log.basis)$fd
  log.coefs=log.xfd$coefs
  log.pars0=c(r=2.5,K=20)
  names(log.pars0)=c("r","K")
  log.lambda=c(0.1)
  resultList=Profile.LS(log.fn,y,t_obs,log.pars0,log.coefs,log.basis,log.lambda)
  
  colloc_par[ii,]=resultList$pars[1:2]
  
  #### SAVE RESULTS TO OBJECTS
  
  medie.teta[[ii]]=mod$mean.teta
  dev.teta[[ii]]=mod$sd.teta
  
  medie.beta[[ii]]=mod$mean.beta
  dev.beta[[ii]]=mod$sd.beta
  
  media.sigma_d[[ii]]=mod$mean.sigma_d
  dev.sigma_d[[ii]]=mod$sd.var_d
  
  media.mu[[ii]]=mod$mean.mu
  dev.mu[[ii]]=mod$sd.mu
  
  media.post_x[[ii]]=mod$post.x[,1]
  
  media.post_g1[[ii]]=mod$post.g1
  
  print(ii)
}

save.image("log_n500_low-noise.RData")


