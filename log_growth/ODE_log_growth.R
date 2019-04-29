ode.mcmc_R<-function(y,t, q,
                     nit=10000, burn.in=round(nit/2),
                     progr=TRUE,
                     prior.lam_teta, ###matrix p x 2
                     prior.lam_beta, ### vector of length 2
                     prior.sigma, ###matrix p x 2
                     init.param=c(0.1,0.1,0.1,0.1)){
  
  time.start=Sys.time()
  require(mvtnorm)
  require(combinat)  
  
  dt=t[2]-t[1]
  t_obs=t/max(t)
  n=length(y)
  p=1
  n.par=1+1 ### x+x^2
  
  prior.lam_teta=c(0.01,0.1)
  prior.lam_beta=c(0.01,0.1)
  
  dep.teta=list()
  dep.lam.beta=rep(0,nit)
  dep.sigma_d=rep(0,nit)
  dep.lam.teta=matrix(0,nit,p)
  dep.sigma_y=matrix(0,nit,p)
  dep.beta=matrix(0,nit,n.par)
  dep.mu=rep(0,nit)
  
  lam.teta=rep(0,p)
  sigma_y=rep(0,p)
  sigma_d=0.1
  t.knots=list()
  Psi=list()
  S_teta=list()
  S_beta=diag(1,n.par)
  
  for(k in 1:p) {
    t.knots[[k]]=1:q[k]/(q[k]+1)
    Psi[[k]]=spl(t_obs,t.knots[[k]])
    S_teta[[k]]=penalty.matrix(t.knots[[k]])
    lam.teta[k]=0.001
    sigma_y[k]=var(y[,k])
    dep.teta[[k]]=matrix(0,nit,q[k]+2)
  }
  lam.beta=1
  beta=rmvnorm(1,mean=rep(0,n.par),sigma = solve(S_beta*lam.beta))
  
  x=y
  int.reg=reg=matrix(0,n,n.par)
  
  
  mu=mean(y[,1])
  
  for(tt in 1:nit){
    
    ### Sample teta1
    var.teta=solve(S_teta[[1]]*lam.teta[1]+(t(Psi[[1]])%*%Psi[[1]])/sigma_d)
    m.teta=var.teta%*%(t(Psi[[1]])%*%(y[,1]/sigma_d))
    teta=rmvnorm(1,mean=m.teta,sigma=var.teta)
    dep.teta[[1]][tt,]=teta
    
    ### Sample lam_teta1
    rt_lam.teta=teta%*%S_teta[[1]]%*%t(teta)
    lam.teta[1]=rgamma(1,shape=((q[1]+2)/2+prior.lam_teta[1]), rate=(rt_lam.teta/2+prior.lam_teta[2]))
    
    ### Compute x1
    x[,1]=Psi[[1]]%*%t(teta)
    ### Compute \tilde{x1}
    
    int.reg[,1]=trap.int(x[,1],dt)
    int.reg[,2]=trap.int(x[,1]^2,dt)
    
    rt_y=(y[,1]-x[,1])%*%as.matrix(y[,1]-x[,1])
    shp_y=n
    sigma_y[1]=1/rgamma(1,shape=shp_y/2,rate=rt_y/2)                                                        
    
    ### Sample beta
    var.beta=solve(S_beta*lam.beta+(t(int.reg)%*%int.reg)/sigma_d)
    m.beta=var.beta%*%(t(int.reg)%*%((y[,1]-rep(mu,n))/sigma_d))
    beta=rmvnorm(1,mean=m.beta,sigma=var.beta)
    
    ### Compute g1 solution of ODE
    g1=int.reg%*%t(beta)
    
    
    ### Sample lam_beta
    rt_b=beta%*%t(beta)
    shp_b=n.par
    lam.beta=rgamma(1,shape=(shp_b/2+prior.lam_beta[1]), rate=(rt_b/2+prior.lam_beta[2]))

    
    rt_d=t(y[,1]-rep(mu,n)-g1)%*%as.matrix(y[,1]-rep(mu,n)-g1)+(y[,1]-x[,1])%*%as.matrix(y[,1]-x[,1])
    
    shp_d=n+n
    sigma_d=1/rgamma(1,shape=shp_d/2,rate=rt_d/2)
    
    mu=rnorm(1,mean=(y[,1]-g1)[1],sd=sqrt(sigma_d))
    
    dep.beta[tt,]=beta
    dep.sigma_y[tt,]=sigma_y
    dep.sigma_d[tt]=sigma_d
    dep.lam.teta[tt,]=lam.teta
    dep.lam.beta[tt]=lam.beta
    dep.mu[tt]=mu
    
    if(progr==TRUE) print(tt)
    
    
  }
  
  mean.teta=lapply(dep.teta,function(x) apply(x[burn.in:nit,],2,mean))
  mean.lam.teta=mean(dep.lam.teta[burn.in:nit,])
  mean.lam.beta=mean(dep.lam.beta[burn.in:nit])
  mean.sigma_y=mean(dep.sigma_y[burn.in:nit,])
  mean.sigma_d=mean(dep.sigma_d[burn.in:nit])
  mean.mu=mean(dep.mu[burn.in:nit])
  
  sd.teta=lapply(dep.teta,function(x) apply(x[burn.in:nit,],2,sd))
  sd.lam.teta=sd(dep.lam.teta[burn.in:nit,])
  sd.lam.beta=sd(dep.lam.beta[burn.in:nit])
  sd.sigma_y=sd(dep.sigma_y[burn.in:nit,])
  sd.sigma_d=sd(dep.sigma_d[burn.in:nit])
  sd.mu=sd(dep.mu[burn.in:nit])
  
  
  post.x=matrix(0,n,p)
  for(k in 1:p){
    
    post.x[,k]=Psi[[k]]%*%mean.teta[[k]]
    
  }
  int.reg[,1]=trap.int(post.x[,1],dt)
  int.reg[,2]=trap.int(post.x[,1]^2,dt)
  
  mean.beta=apply(dep.beta[burn.in:nit,],2,mean)
  sd.beta=apply(dep.beta[burn.in:nit,],2,sd)
  post.g1=int.reg%*%mean.beta
  
  time.end=Sys.time()
  output=list(mean.teta=mean.teta,sd.teta=sd.teta,
              chain.teta=dep.teta,
              mean.lam_teta=mean.lam.teta,sd.lam_teta=sd.lam.teta,
              chain.lam_teta=dep.lam.teta,
              mean.lam_beta=mean.lam.beta,sd.lam_beta=sd.lam.beta,
              chain.lam_beta=dep.lam.beta,
              mean.beta=mean.beta,sd.beta=sd.beta,
              chain.beta=dep.beta,
              mean.sigma_y=mean.sigma_y, sd.var_y=sd.sigma_y,
              chain.sigma_y=dep.sigma_y,
              mean.sigma_d=mean.sigma_d, sd.var_d=sd.sigma_d,
              chain.sigma_d=dep.sigma_d,
              mean.mu=mean.mu,sd.mu=sd.mu,
              chain.mu=dep.mu,
              post.x=post.x,
              post.g1=post.g1,
              int.reg=int.reg,
              comp.time=time.end-time.start)
}
require(compiler)
ode.mcmc_C<-cmpfun(ode.mcmc_R)
### builds spline basis between point "x" and knot "z"
basis<-function(x,z){
  return(((z-0.5)^2-1/12)*((x-0.5)^2-1/12)/4-((abs(x-z)-0.5)^4-(abs(x-z)-0.5)^2/2+7/240)/24)
}

### produces matrix for bayesian spline smoothing [1,t,psi(t)]
spl<-function(x,xk){
  q=length(xk)+2
  n=length(x)
  X=matrix(1,n,q)
  X[,2]=x
  X[,3:q]=outer(x,xk,FUN=basis)
  return(X)
}

### unconstrained penalty matrix S_teta
penalty.matrix<-function(xk){
  q=length(xk)+2
  S=matrix(0,q,q)
  S[3:q,3:q]=outer(xk,xk,FUN=basis)
  diag(S)[1:2]=10^(-10)
  return(S)
}

### Trapezoidal rule for integration
trap.int=function(z,dt){
  intgrl=cumsum(z[1:(length(z)-1)]+z[2:length(z)])*dt*0.5
  intgrl=c(0,intgrl)
  return(as.matrix(intgrl))
}