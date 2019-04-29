################################################################################################
###################### sim HIV  ##########################################################
################################################################################################
################################################################################################
source('lotka/ODE_lotka.R')
require(deSolve)
require(fda)
require(CollocInfer)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*alpha - beta*x*y
    dy = -y*gamma + delta*y*x
    return(list(c(dx, dy)))
  })
}


##### CollocInfer - HIV

make.lotkafn<-function(){
  
  ### set up functions for 
  ### dx = x*alpha - beta*x*y
  ### dy = -y*gamma + delta*y*x
  
  lotkafun<-function(times,x,p,more){
    n=length(times);d=2;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    f=matrix(0,n,d)
    f[,1]=x[,1]*alpha-beta*x[,1]*x[,2]
    f[,2]=-gamma*x[,2]+delta*x[,2]*x[,1]
    return(f)
  }
  
  
  lotkadfdx<-function(times,x,p,more){
    n=length(times);d=2;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    dfdx=array(0,c(n,d,d))
    dfdx[,1,1]=rep(alpha,n)
    dfdx[,1,2]=-beta*x[,1]
    dfdx[,2,1]=delta*x[,2]
    dfdx[,2,2]=-rep(gamma,n)+delta*x[,1]
    return(dfdx)
  }
  
  ### set up functions for 
  ### dx = x*alpha - beta*x*y
  ### dy = -y*gamma + delta*y*x
  
  
  lotkadfdp<-function(times,x,p,more){
    n=length(times);d=2;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    dfdp=array(0,c(n,d,npar))
    dfdp[,1,1]=x[,1]
    dfdp[,1,2]=-x[,1]*x[,2]
    dfdp[,2,3]=-x[,2]
    dfdp[,2,4]=delta*x[,1]
    return(dfdp)
  }
  
  
  lotkad2fdx2<-function(times,x,p,more){
    n=length(times);d=2;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    d2fdx2=array(0,c(n,d,d,d))
    d2fdx2[,1,2,1]=-rep(beta,n)
    d2fdx2[,2,1,2]=rep(delta,n)
    d2fdx2[,2,2,1]=rep(delta,n)
    return(d2fdx2)
  }
  
  lotkad2fdxdp<-function(times,x,p,more){
    n=length(times);d=2;npar=4
    alpha=p[1]
    beta=p[2]
    gamma=p[3]
    delta=p[4]
    d2fdxdp=array(0,c(n,d,d,npar))
    d2fdxdp[,1,1,1]=rep(1,n)
    d2fdxdp[,1,2,2]=-x[,1]
    d2fdxdp[,2,1,4]=x[,2]
    d2fdxdp[,2,2,3]=-rep(1,n)
    d2fdxdp[,2,2,4]=x[,1]
    
    return(d2fdxdp)
  }
  
  return(list(fn=lotkafun, dfdx=lotkadfdx,
              dfdp=lotkadfdp, d2fdx2=lotkad2fdx2,
              d2fdxdp=lotkad2fdxdp))
  
}

true.y=matrix(0,100,500)
true.sd1=true.sd2=rep(0,100)
burn.in=5000
nit=10000
noise.coef=0.1

lotka.fn=make.lotkafn()
lotka.knots=seq(0,99,99/20)
lotka.order=3
lotka.nbasis=length(lotka.knots)+lotka.order-2
lotka.range=c(0,99)
lotka.basis=create.bspline.basis(lotka.range,lotka.nbasis,lotka.order,lotka.knots)

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
  Pars <- c(alpha = 0.1, beta = 0.2, gamma = 1, delta = 0.3)
  State <- c(x = 2, y = 2)
  Time <- seq(0, 99, by = 99/499)
  out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))
  
  true.y[ii,]=out$x
  
  true.sd1[ii]=mean(out$x)*noise.coef
  true.sd2[ii]=mean(out$y)*noise.coef
  
  out$x=out$x+rnorm(length(out$x),0,sd=true.sd1[ii])
  out$y=out$y+rnorm(length(out$y),0,sd=true.sd2[ii])
  
  t=out$time
  y=as.matrix(cbind(out$x,out$y))
  q=c(20,20)
  
  mod=ode.mcmc_C(y,t,q=q,nit = nit, burn.in = burn.in, progr=FALSE)
  
  lotka.xfd=smooth.basis(t,y,lotka.basis)$fd
  lotka.coefs=lotka.xfd$coefs
  lotka.pars0=c(alpha = 0.1, beta = 0.2, gamma = 1, delta = 0.3)
  names(lotka.pars0)=c("alpha","beta","gamma","delta")
  lotka.lambda=1e1*c(1,1)
  resultList=Profile.LS(lotka.fn,y,t,lotka.pars0,lotka.coefs,lotka.basis,lotka.lambda)
  
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

save.image("lotka_n500_med-noise.RData")


