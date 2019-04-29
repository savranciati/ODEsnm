################################################################################################
###################### Lotka-Volterra ##########################################################
################################################################################################
################################################################################################

library(deSolve)

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*alpha - beta*x*y
    dy = -y*gamma + delta*y*x
    return(list(c(dx, dy)))
  })
}

Pars <- c(alpha = 0.1, beta = 0.2, gamma = 1, delta = 0.3)
State <- c(x = 2, y = 2)
Time <- seq(0, 99, by = 1)

out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))

matplot(out[,-1], type = "l", xlab = "time", ylab = "population")
legend("topright", c("Cute bunnies", "Rabid foxes"), lty = c(1,2), col = c(1,2), box.lwd = 0)

t=out$time
y=matrix(c(out$x,out$y),length(out$x),2,byrow=FALSE)
set.seed(1530)
true.sd=apply(y,2,mean)
true.sd=true.sd*0.04
true.y=y[,1]
y[,1]=y[,1]+rnorm(length(y[,1]),0,sd=true.sd[1])
y[,2]=y[,2]+rnorm(length(y[,2]),0,sd=true.sd[2])
plot(y[,1])
plot(y[,2])

q=c(20,20)
nit=10000
burn.in=5000
out=ode.mcmc_C(y,t,q=q,nit=nit, burn.in=burn.in,progr=TRUE)



########################################################################
######################## Chains and Plots ##############################
########################################################################
for(k in 1:dim(y)[2]){
  plot(t,y[,k],
       xlab="Time Points",
       ylab=eval(bquote(expression(y[.(k)]))))
  lines(t,out$post.x[,k], col=4)
  # legend("topright", legend="post.x", lty=1,pch=NA,lwd=2,col="blue")

  plot(out$chain.lam_teta[burn.in:nit,k], type="l",
       xlab="iteration (after warm-up)",
       ylab=eval(bquote(expression(theta[.(k)]))))

    plot(out$chain.sigma_y[burn.in:nit,k], type="l",
     xlab="iteration (after warm-up)",
     ylab=expression(sigma["y"]^2))

}

plot(out$chain.lam_beta[burn.in:nit], type="l")
plot(out$chain.sigma_d[burn.in:nit], type="l")
plot(out$chain.mu[burn.in:nit], type="l")

plot(t,y[,1])
lines(t,out$post.g1+out$mean.mu,col=4) ###blue
lines(t,out$post.x[,1],col=2) ###red


est.beta=matrix(c(round(out$mean.beta,3),round(out$sd.beta,3)),length(out$mean.beta),2,byrow=FALSE)
colnames(est.beta)=c("mean","sd")
est.beta
credib=apply(out$chain.beta[burn.in:nit,],2,quantile,probs=c(0.025,0.5,0.975))
colnames(credib)=c("alpha","beta")
round(t(credib),3)

round(out$mean.mu,3)
round(out$sd.mu,3)
round(quantile(out$chain.mu[burn.in:nit],probs=c(0.025,0.5,0.975)),3)


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



plot(t,y[,1],
     main="",
     pch=5,
     cex=1,
     xlab="time", xlim=c(0,99),
     ylab=expression(y["1"]))
points(t,out$post.g1+out$mean.mu,
       type="l",
       pch=19,
       lty="longdash",
       cex=0.8,
       lwd=1.5,
       col=2) ## red
points(t,true.y,
       type="l",
       pch=15,
       lty="solid",
       cex=0.8,
       lwd=1.5, 
       col=4) ## blue
points(t,out$post.x[,1],
       type="l",
       pch=15,
       lty="dotted",
       cex=0.8,
       lwd=1.5, 
       col=6) ## green

