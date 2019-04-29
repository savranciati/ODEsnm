require(fda)

ARE<-function(estim,param){
  param=rep(param,length(estim))
  return(mean(abs(estim-param)/abs(estim))*100)
}
MSE<-function(estim,param){
  param=rep(param,length(estim))
  return(mean((estim-param)^2))
}
tradeoff<-function(estim,param){
  bias=mean(estim)-param
  return((bias^2/(bias^2+mean(estim^2)+mean(estim)^2))*100)
}
  
load("~/Documents/Dropbox/Saverio/ODE/script/sim1000/hiv_n500_high-noise.RData")

mat=matrix(unlist(media.mu),ncol=1,byrow=TRUE)

MSE(mat[,1],60)
tradeoff(mat[,1],60)

mat=matrix(unlist(medie.beta),ncol=3,byrow=TRUE)


MSE(mat[,1],20)
tradeoff(mat[,1],20)
MSE(colloc_par[,1],20)
tradeoff(colloc_par[,1],20)


MSE(mat[,2],-0.108)
tradeoff(mat[,2],-0.108)


MSE(mat[,3],-0.00095)*10^5
tradeoff(mat[,3],-0.00095)
MSE(-colloc_par[,2],-0.00095)*10^5
tradeoff(-colloc_par[,2],-0.00095)

