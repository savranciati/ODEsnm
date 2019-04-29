#### Logistic growth

log_growth<-function(t,r,K,x0){
  K*x0*exp(r*t)/(K+x0*(exp(r*t)-1))
}
t=(0:499)/500
plot(t,log_growth(t,r=2.5,K=20,x0=0.1))

t_obs=(0:499)/500
y=as.matrix(log_growth(t,r=2.5,K=20,x0=0.1))
true.sd=sd(y)*0.1
set.seed(1433)
true.y=y
y=y+rnorm(length(y),mean=0,sd=true.sd)
plot(y)

q=2
nit=10000
burn.in=5000
out=ode.mcmc_C(y,t_obs,q=q,nit=nit, burn.in = burn.in,progr=FALSE)


plot(t_obs,y,
     xlab="Time Points",
     ylab=eval(bquote(expression(y[.(k)]))))
lines(t_obs,out$post.x, col=4) ## blue

plot(out$chain.lam_teta[burn.in:nit], type="l", xlab="iteration (after warm-up)")

plot(out$chain.sigma_d[burn.in:nit], type="l",
     xlab="iteration (after warm-up)",
     ylab=expression(sigma["y"]^2))

plot(out$chain.lam_beta[burn.in:nit], type="l")

plot(t_obs,y[,1])
lines(t_obs,out$post.g1+out$mean.mu,col=4) ##blue line


est.beta=matrix(c(round(out$mean.beta,3),round(out$sd.beta,3)),length(out$mean.beta),2,byrow=FALSE)
colnames(est.beta)=c("mean","sd")
est.beta
credib=apply(out$chain.beta[burn.in:nit,],2,quantile,probs=c(0.025,0.5,0.975))
colnames(credib)=c("r","K")
round(t(credib),3)

round(out$mean.mu,3)
round(out$sd.mu,3)
# plot(out$chain.mu[burn.in:nit],type="l")
round(quantile(out$chain.mu[burn.in:nit],probs=c(0.025,0.5,0.975)),3)


round(out$mean.sigma_d,6)
round(median(out$chain.sigma_d[burn.in:nit]),6)

round(out$mean.sigma_y,6)
round(median(out$chain.sigma_y[burn.in:nit]),6)


round(out$sd.var_d,5)
round(out$sd.var_y,5)


# # now the confidence bands 95%
# qtl=abs(qnorm(0.025))
# plot(t_obs,y[,1])
# upp.limit=out$post.g1+out$mean.mu+qtl*sqrt(out$mean.sigma_d)
# low.limit=out$post.g1+out$mean.mu-qtl*sqrt(out$mean.sigma_d)
# polygon(c(t_obs, rev(t_obs)), c(upp.limit, rev(low.limit)), col = "grey80", border = NA)
# points(t_obs,y[,1])
# lines(t_obs,out$post.g1+out$mean.mu,col=4) ###blue
# lines(t_obs, upp.limit , lty = "dotted")
# lines(t_obs, low.limit, lty = "dotted")
# 
# 
# # now the confidence bands 99%
# qtl=abs(qnorm(0.005))
# plot(t_obs,y[,1])
# upp.limit=out$post.g1+out$mean.mu+qtl*sqrt(out$mean.sigma_d)
# low.limit=out$post.g1+out$mean.mu-qtl*sqrt(out$mean.sigma_d)
# polygon(c(t_obs, rev(t_obs)), c(upp.limit, rev(low.limit)), col = "grey80", border = NA)
# points(t_obs,y[,1])
# lines(t_obs,out$post.g1+out$mean.mu,col=4) ###blue
# lines(t_obs, upp.limit , lty = "dotted")
# lines(t_obs, low.limit, lty = "dotted")
# 


plot(t_obs,y,
  main="",
  pch=5,
  cex=1,
  xlab="time", xlim=c(0,1),
  ylab=expression(y["1"]))
points(t_obs,out$post.g1+out$mean.mu,
       type="l",
       pch=19,
       lty="longdash",
       cex=0.8,
       lwd=1.5,
       col=2) ## red
points(t_obs,true.y,
       type="l",
       pch=15,
       lty="solid",
       cex=0.8,
       lwd=1.5, 
       col=4) ## blue
points(t_obs,out$post.x,
       type="l",
       pch=15,
       lty="dotted",
       cex=0.8,
       lwd=1.5, 
       col=6) ## green
