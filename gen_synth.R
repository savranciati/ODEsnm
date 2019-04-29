################################################################################################
###################### Syntethic-ODE ##########################################################
################################################################################################
################################################################################################

library(deSolve)

SynthMod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    da = a*beta1 - beta2*j*sin(b)
    db = -b*beta3 + beta4*a*sin(c)
    dc = c*beta5 - beta6*b*cos(d)
    dd = -d*beta7 + beta8*c*cos(e)
    de = e*beta9 - gamma1*d*cos(f)
    df = -f*gamma2 + gamma3*e*cos(g)
    dg = g*gamma4 - gamma5*f*cos(h)
    dh = -h*gamma6 + gamma7*g*sin(i)
    di = i*gamma8 - gamma9*h*sin(j)
    dj = -j*delta1 + delta2*i*sin(a)
    return(list(c(da, db, dc, dd, de, df, dg, dh, di, dj)))
  })
}

Pars <- c(beta1 = 0.03, beta2 = 0.08,
          beta3 = 0.05, beta4 = 0.07,
          beta5 = 0.05, beta6 = 0.06,
          beta7 = 0.05, beta8 = 0.08,
          beta9 = 0.05, gamma1 = 0.07,
          gamma2 = 0.05, gamma3 = 0.08,
          gamma4 = 0.02, gamma5 = 0.08,
          gamma6 = 0.04, gamma7 = 0.06,
          gamma8 = 0.05, gamma9 = 0.06,
          delta1 = 0.05, delta2 = 0.08)
State <- c(a = 0.1, b = 0.1, c = 0.15, d = 0.1, e = 0.1,
           f = 0.1, g = 0.1, h = 0.1, i = 0.1, j = 0.1)
Time <- seq(0, 100, by = 0.5)

out <- as.data.frame(ode(func = SynthMod, y = State, parms = Pars, times = Time))
sum(is.na(out))

par(mfrow=c(4,3))
plot(out$time,out$a, type = "l", xlab = "time", ylab = "")
plot(out$time,out$b, type = "l", xlab = "time", ylab = "")
plot(out$time,out$c, type = "l", xlab = "time", ylab = "")
plot(out$time,out$d, type = "l", xlab = "time", ylab = "")
plot(out$time,out$e, type = "l", xlab = "time", ylab = "")
plot(out$time,out$f, type = "l", xlab = "time", ylab = "")
plot(out$time,out$g, type = "l", xlab = "time", ylab = "")
plot(out$time,out$h, type = "l", xlab = "time", ylab = "")
plot(out$time,out$i, type = "l", xlab = "time", ylab = "")
plot(out$time,out$j, type = "l", xlab = "time", ylab = "")
par(mfrow=c(1,1))

t=out$time
y=as.matrix(out[,-1])
set.seed(1530)
true.sd=apply(y,2,mean)
defl <- 0.4
true.sd=true.sd*defl
true.y=y
y <- apply(y,2, function(x,defl) x+rnorm(length(x),0,sd=sd(x)*defl),defl)

# par(mfrow=c(4,3))
# plot(out$time,out$a, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,1],col=2)
# plot(out$time,out$b, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,2],col=2)
# plot(out$time,out$c, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,3],col=2)
# plot(out$time,out$d, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,4],col=2)
# plot(out$time,out$e, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,5],col=2)
# plot(out$time,out$f, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,6],col=2)
# plot(out$time,out$g, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,7],col=2)
# plot(out$time,out$h, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,8],col=2)
# plot(out$time,out$i, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,9],col=2)
# plot(out$time,out$j, type = "l", xlab = "time", ylab = "")
# points(t,pert.y[,10],col=2)
# par(mfrow=c(1,1))


nit=10000
burn.in=5000
q=rep(10,10)
source('ODE_synth1.R')
source('ODE_synth2.R')
source('ODE_synth3.R')
source('ODE_synth4.R')
source('ODE_synth5.R')
source('ODE_synth6.R')
source('ODE_synth7.R')
source('ODE_synth8.R')
source('ODE_synth9.R')
source('ODE_synth10.R')
out1=ode.mcmc_C1(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out2=ode.mcmc_C2(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out3=ode.mcmc_C3(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out4=ode.mcmc_C4(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out5=ode.mcmc_C5(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out6=ode.mcmc_C6(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out7=ode.mcmc_C7(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out8=ode.mcmc_C8(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out9=ode.mcmc_C9(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)
out10=ode.mcmc_C10(y,t,q=q,nit=nit, burn.in=burn.in,progr=FALSE)


########################################################################
######################## Chains and Plots ##############################
########################################################################
# out <- out1
# for(k in 1:dim(y)[2]){
#   plot(t,y[,k],
#        xlab="Time Points",
#        ylab=eval(bquote(expression(y[.(k)]))))
#   lines(t,out$post.x[,k], col=4)
#   # legend("topright", legend="post.x", lty=1,pch=NA,lwd=2,col="blue")
# 
#   plot(out$chain.lam_teta[burn.in:nit,k], type="l",
#        xlab="iteration (after warm-up)",
#        ylab=eval(bquote(expression(theta[.(k)]))))
# 
#   plot(out$chain.sigma_y[burn.in:nit,k], type="l",
#        xlab="iteration (after warm-up)",
#        ylab=expression(sigma["y"]^2))
# 
# }
time.comp <- rep(0,10)
beta.list <- list()
mu.list <- list()
true.params <- matrix(Pars,10,2, byrow=TRUE)
true.params <- cbind(true.params,matrix(0,10,1))
par(mfrow=c(4,3))
for(ii in 1:10){
  out <- eval(parse(text=paste("out",ii,sep="")))
  plot(t,y[,ii],
     main="",
     pch=5,
     cex=1,
     xlab="time", xlim=c(0,99),
     ylab=paste("y",ii))
  points(t,out$post.g1+out$mean.mu,
       type="l",
       pch=19,
       lty="longdash",
       cex=0.8,
       lwd=1.5,
       col=2) ## red
  points(t,true.y[,ii],
       type="l",
       pch=15,
       lty="solid",
       cex=0.8,
       lwd=1.5,
       col=4) ## blue
  points(t,out$post.x[,ii],
       type="l",
       pch=15,
       lty="dotted",
       cex=0.8,
       lwd=1.5,
       col=6) ## purple

  est.beta=matrix(c(round(out$mean.beta,3),round(out$sd.beta,3)),length(out$mean.beta),2,byrow=FALSE)
  colnames(est.beta)=c("mean","sd")
  est.beta
  credib=apply(out$chain.beta[burn.in:nit,],2,quantile,probs=c(0.025,0.5,0.975))
  # colnames(credib)=c("par1","par2")
  beta.list[[ii]] <- round(cbind(t(credib),true=true.params[ii,]),3)

  round(out$mean.mu,3)
  round(out$sd.mu,3)
  mu.list[[ii]] <- round(quantile(out$chain.mu[burn.in:nit],probs=c(0.025,0.5,0.975)),3)
  time.comp[ii] <- out$comp.time
}
par(mfrow=c(1,1))

mean(time.comp)


round(true.sd[1]^2,5)

round(out$mean.sigma_d,5)
round(median(out$chain.sigma_d[burn.in:nit]),5)

round(out$mean.sigma_y,5)[1]
round(median(out$chain.sigma_y[burn.in:nit,1]),5)


round(out$sd.var_d,5)
round(out$sd.var_y[1],5)

# # now the confidence bands 95%
# qtl=abs(qnorm(0.025))
# plot(t,y[,1],ylim=c(0,7.5))
# upp.limit=out$mean.mu+out$post.g1+qtl*sqrt(out$mean.sigma_d)
# low.limit=out$mean.mu+out$post.g1-qtl*sqrt(out$mean.sigma_d)
# polygon(c(t, rev(t)), c(upp.limit, rev(low.limit)), col = "grey80", border = NA)
# points(t,y[,1])
# lines(t,out$mean.mu+out$post.g1,col=4) ###blue
# lines(t, upp.limit , lty = "dotted")
# lines(t, low.limit, lty = "dotted")
# 
# 
# # now the confidence bands 99%
# qtl=abs(qnorm(0.005))
# plot(t,y[,1],ylim=c(0,7.5))
# upp.limit=out$mean.mu+out$post.g1+qtl*sqrt(out$mean.sigma_d)
# low.limit=out$mean.mu+out$post.g1-qtl*sqrt(out$mean.sigma_d)
# polygon(c(t, rev(t)), c(upp.limit, rev(low.limit)), col = "grey80", border = NA)
# points(t,y[,1])
# lines(t,out$mean.mu+out$post.g1,col=4) ###blue
# lines(t, upp.limit , lty = "dotted")
# lines(t, low.limit, lty = "dotted")
# 

tosort <- array(unlist(beta.list),c(3,4,10))
compact <- apply(tosort,2,c)
colnames(compact) <- c("    2.5%","    median","    97.5%","true.par")
rownames(compact) <- c("beta1","beta2",rep("zero",1),
                       "beta3","beta4",rep("zero",1),
                       "beta5","beta6",rep("zero",1),
                       "beta7","beta8",rep("zero",1),
                       "beta9","gamma1",rep("zero",1),
                       "gamma2","gamma3",rep("zero",1),
                       "gamma4","gamma5",rep("zero",1),
                       "gamma6","gamma7",rep("zero",1),
                       "gamma8","gamma9",rep("zero",1),
                       "delta1","delta2",rep("zero",1))
compact

tosortmu <- matrix(unlist(mu.list),10,3,byrow=TRUE)
compactmu <- cbind(tosortmu,State)
colnames(compactmu) <- c("    2.5%","    median","    97.5%","true.par")
rownames(compactmu) <- c("mu1","mu2","mu3","mu4","mu5",
                         "mu6","mu7","mu8","mu9","mu10")
compactmu
