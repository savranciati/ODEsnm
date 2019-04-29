################################################################################################
###################### sim HIV  ##########################################################
################################################################################################
################################################################################################
source('hiv/ODE_hiv.R')
require(deSolve)
require(fda)
require(CollocInfer)

HIV <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = alpha - 0.108*x - beta*x*z
    dy = beta*x*z - 0.5*y
    dz = 0.5*gamma*y - delta*z
    return(list(c(dx, dy, dz)))
  })
}

##### CollocInfer - HIV

make.hivfn<-function(){
  
  ### set up functions for 
  ## dx = alpha - 0.108*x - beta*x*z
  ## dy = beta*x*z - 0.5*y
  ## dz = 0.5*gamma*y - delta*z
  
  hivfun<-function(times,x,p,more){
    n=length(times);d=3;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    f=matrix(0,n,d)
    f[,1]=alpha-0.108*x[,1]-beta*x[,1]*x[,3]
    f[,2]=beta*x[,1]*x[,3]-0.5*x[,2]
    f[,3]=0.5*gamma*x[,2]-delta*x[,3]
    return(f)
  }
  
  ## dx = alpha - 0.108*x - beta*x*z
  ## dy = beta*x*z - 0.5*y
  ## dz = 0.5*gamma*y - delta*z
  
  
  hivdfdx<-function(times,x,p,more){
    n=length(times);d=3;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    dfdx=array(0,c(n,d,d))
    dfdx[,1,1]=-rep(0.108,n)-beta*x[,3]
    dfdx[,1,3]=-beta*x[,1]
    dfdx[,2,1]=beta*x[,3]
    dfdx[,2,2]=-rep(0.5,n)
    dfdx[,2,3]=beta*x[,1]
    dfdx[,3,2]=gamma*rep(0.5,n)
    dfdx[,3,3]=-rep(delta,n)
    return(dfdx)
  }
  
  hivdfdp<-function(times,x,p,more){
    n=length(times);d=3;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    dfdp=array(0,c(n,d,npar))
    dfdp[,1,1]=rep(1,n)
    dfdp[,1,2]=-x[,1]*x[,3]
    dfdp[,2,2]=x[,1]*x[,3]
    dfdp[,3,3]=0.5*x[,2]
    dfdp[,3,4]=-x[,3]
    return(dfdp)
  }
  
  
  hivd2fdx2<-function(times,x,p,more){
    n=length(times);d=3;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    d2fdx2=array(0,c(n,d,d,d))
    d2fdx2[,1,1,3]=-rep(beta,n)
    d2fdx2[,1,3,1]=-rep(beta,n)
    d2fdx2[,2,1,3]=-rep(beta,n)
    d2fdx2[,2,3,1]=rep(beta,n)
    return(d2fdx2)
  }
  
  hivd2fdxdp<-function(times,x,p,more){
    n=length(times);d=3;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    d2fdxdp=array(0,c(n,d,d,npar))
    d2fdxdp[,1,1,2]=-x[,3]
    d2fdxdp[,1,3,2]=-x[,1]
    d2fdxdp[,2,1,2]=x[,3]
    d2fdxdp[,2,3,2]=x[,1]
    d2fdxdp[,3,2,3]=rep(0.5,n)
    d2fdxdp[,3,3,4]=-rep(1,n)
    return(d2fdxdp)
  }
  
  return(list(fn=hivfun, dfdx=hivdfdx,
              dfdp=hivdfdp, d2fdx2=hivd2fdx2,
              d2fdxdp=hivd2fdxdp))
  
}

true.y=matrix(0,100,25)
true.sd1=true.sd2=true.sd3=rep(0,100)
burn.in=5000
nit=10000
noise.coef=0.04

hiv.fn=make.hivfn()
hiv.knots=seq(0,24,24/4)
hiv.order=3
hiv.nbasis=length(hiv.knots)+hiv.order-2
hiv.range=c(0,24)
hiv.basis=create.bspline.basis(hiv.range,hiv.nbasis,hiv.order,hiv.knots)

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
  Pars <- c(alpha = 20, beta = 9.5*10^(-4), gamma = 10, delta = 1)
  State <- c(x = 60, y = 3, z=10^2)
  Time <- seq(0, 24, by = 1)
  out <- as.data.frame(ode(func = HIV, y = State, parms = Pars, times = Time))
  
  true.y[ii,]=out$x
  
  true.sd1[ii]=mean(out$x)*noise.coef
  true.sd2[ii]=mean(out$y)*noise.coef
  true.sd3[ii]=mean(out$z)*noise.coef
  
  out$x=out$x+rnorm(length(out$x),0,sd=true.sd1[ii])
  out$y=out$y+rnorm(length(out$y),0,sd=true.sd2[ii])
  out$z=out$z+rnorm(length(out$z),0,sd=true.sd3[ii])
  
  t=out$time
  y=as.matrix(cbind(out$x,out$y,out$z))
  q=c(5,5,5)
  
  mod=ode.mcmc_C(y,t,q=q,nit = nit, burn.in = burn.in, progr=FALSE)
  
  hiv.xfd=smooth.basis(t,y,hiv.basis)$fd
  hiv.coefs=hiv.xfd$coefs
  hiv.pars0=c(20,9.5*10^(-4), gamma = 10, delta = 1)
  names(hiv.pars0)=c("alpha","beta","gamma","delta")
  hiv.lambda=1e1*c(1,1,1)
  resultList=Profile.LS(hiv.fn,y,t,hiv.pars0,hiv.coefs,hiv.basis,hiv.lambda)
  
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

save.image("hiv_n25_med-noise.RData")


