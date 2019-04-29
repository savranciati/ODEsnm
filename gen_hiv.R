################################################################################################
###################### HIV  ##########################################################
################################################################################################
################################################################################################

library(deSolve)

HIV <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = alpha - 0.108*x - beta*x*z
    dy = beta*x*z - 0.5*y
    dz = 0.5*gamma*y - delta*z
    return(list(c(dx, dy, dz)))
  })
}

Pars <- c(alpha = 20, beta = 9.5*10^(-4), gamma = 10, delta = 1)
State <- c(x = 60, y = 3, z=10^2)
Time <- seq(0, 99, by = 0.2)

out <- as.data.frame(ode(func = HIV, y = State, parms = Pars, times = Time))

set.seed(1237)
true.y=out$x

true.sd1=mean(out$x)*0.08
true.sd2=mean(out$y)*0.08
true.sd3=mean(out$z)*0.08

out$x=out$x+rnorm(length(out$x),0,sd=true.sd1)
out$y=out$y+rnorm(length(out$y),0,sd=true.sd2)
out$z=out$z+rnorm(length(out$z),0,sd=true.sd3)

plot(out$time,out$x)
plot(out$time,out$y)
plot(out$time,out$z)


t=out$time
y=as.matrix(cbind(out$x,out$y,out$z))
q=c(5,5,5)
burn.in=5000
nit=10000
mod=ode.mcmc_C(y,t,q=q,nit = nit, burn.in = burn.in, progr=FALSE)


### res


for(k in 1:3){
plot(t,y[,k],
     xlab="Time Points",
     ylab=eval(bquote(expression(y[.(k)]))))
lines(t,mod$post.x[,k], col=4)
# legend("topright", legend="post.x", lty=1,pch=NA,lwd=2,col="blue")
}
for(k in 1:3){
plot(mod$chain.lam_teta[burn.in:nit,k], type="l",
     xlab="iteration (after warm-up)",
     ylab=eval(bquote(expression(theta[.(k)]))))
}
for(k in 1:3){
plot(mod$chain.sigma_y[burn.in:nit,k], type="l",
     xlab="iteration (after warm-up)",
     ylab=expression(sigma["y"]^2))
}
plot(mod$chain.lam_beta[burn.in:nit], type="l")
plot(mod$chain.sigma_d[burn.in:nit], type="l")
plot(mod$chain.mu[burn.in:nit], type="l")

plot(t,y[,1])
lines(t,mod$post.g1+mod$mean.mu,col=4) ###blue
lines(t,mod$post.x[,1],col=2) ###red


est.beta=matrix(c(round(mod$mean.beta,5),round(mod$sd.beta,5)),length(mod$mean.beta),2,byrow=FALSE)
colnames(est.beta)=c("mean","sd")
est.beta

round(mod$mean.mu,3)
round(mod$sd.mu,3)




round(true.sd1^2,5)

round(mod$mean.sigma_d,5)
round(median(mod$chain.sigma_d[burn.in:nit]),5)

round(mod$mean.sigma_y[1],5)
round(median(mod$chain.sigma_y[burn.in:nit,1]),5)

round(mod$sd.var_d,5)
round(mod$sd.var_y[1],5)



plot(t,y[,1],
     main="",
     pch=5,
     cex=1,
     xlab="time", xlim=c(0,99),
     ylab=expression(y["1"]))
points(t,mod$post.g1+mod$mean.mu,
       type="l",
       pch=19,
       lty="longdash",
       cex=0.8,
       lwd=2.5,
       col=2) ## red
points(t,true.y,
       type="l",
       pch=15,
       lty="solid",
       cex=0.8,
       lwd=2.5, 
       col=4) ## blue
points(t,mod$post.x[,1],
       type="l",
       pch=15,
       lty="dotted",
       cex=0.8,
       lwd=2.5, 
       col=6) ## green


