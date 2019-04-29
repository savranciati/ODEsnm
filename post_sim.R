#### post_sim

################# MU
mat=matrix(unlist(media.mu),ncol=1,byrow=TRUE)
mean.mu=apply(mat,2,mean)
sd.mu=apply(mat,2,sd)


################# BETA
### colloc_infer
mean.colloc=apply(colloc_par,2,mean)
sd.colloc=apply(colloc_par,2,sd)
### ODE our model
mat=matrix(unlist(medie.beta),ncol=2,byrow=TRUE)
mean.ode=apply(mat,2,mean)
sd.ode=apply(mat,2,sd)


################ SIGMA_1  
mat=matrix(unlist(media.sigma_d),ncol=1,byrow=TRUE)
mean.sigma=apply(mat,2,mean)
sd.sigma=apply(mat,2,sd)


#mu
round(mean.mu,3)
round(sd.mu,3)

#beta
round(mean.ode,3)
round(sd.ode,3)
round(mean.colloc,3)
round(sd.colloc,3)

#sigma_1
round((mean(true.sd1)^2)*100,3)
round(mean.sigma*100,3)
round(sd.sigma*100,3)

