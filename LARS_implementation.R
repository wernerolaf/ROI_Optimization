our.lar = function(X, y, epsilon=1e-8){
  
  n=nrow(X)
  p=ncol(X)
  
  #compute the mean and standard deviation of each column of X
  X.means=apply(X,2,mean)
  X.sds=apply(X,2,sd)
  
  #center and scale the columns of X
  Xstar=X
  for (i in 1:p) Xstar[,i]=(X[,i]-X.means[i])/X.sds[i]
  yc=y-mean(y)
  
  #step up matrix to store the parameters that will be returned from the function
  
  beta.hat=rep(0,p+1)
  alpha.hat=NULL
  
  #initial step
  r=yc
  
  #compute inner product and choose the first variable to enter the model
  
  inner.products=t(X)%*%r
  j.hat=which.max(abs(inner.products))
  
  #algorithm on the ith step
  
  i=1
  
  while ((i==1)||(alpha.hat[i-1]<1)){
    
    beta.hat=rbind(beta.hat,rep(0,p+1))
    alpha.hat=c(alpha.hat,1)
    njhat=length(j.hat)
    j.hat=c(j.hat,0)
    XE=Xstar[,j.hat]
    d=rep(0,p)
    d[j.hat]=solve(t(XE)%*%XE)%*%t(XE)%*%r
    Xd=Xstar%*%d
    
    #find alpha.hat[i]
    
    for (j in 1:p){
      
      alpha=1
      
      if (j%in%j.hat==FALSE){
        if (abs(sum(Xstar[,j.hat[1]]*r)-sum(Xstar[,j]*Xd))>epsilon){
          alpha=(sum(Xstar[,j.hat[1]]*r)-sum(Xstar[,j]*r))/(sum(Xstar[,j.hat[1]]*r)-sum(Xstar[,j]*Xd))
          if ((alpha<epsilon)|(alpha>1-epsilon)){
            alpha=1
            if (abs(sum(Xstar[,j.hat[1]]*r)+sum(Xstar[,j]*Xd))>epsilon){
              alpha2=(sum(Xstar[,j.hat[1]]*r)+sum(Xstar[,j]*r))/(sum(Xstar[,j.hat[1]]*r)+sum(Xstar[,j]*Xd))
              if ((alpha2>epsilon)&(alpha2<1-epsilon))
                alpha=alpha2
            }
          }
          if (alpha+epsilon<alpha.hat[i]){
            alpha.hat[i]=alpha
            j.hat[njhat+1]=j
          }
        }
      }
      
    }
    beta.hat[i+1,2:(p+1)]=beta.hat[i,2:(p+1)]+alpha.hat[i]*d
    r=r-alpha.hat[i]*Xd
    i=i+1
    
  }
  #translate coefficient estimates back to original scale
  beta.hat[,-1]=t(t(beta.hat[,-1])/X.sds)
  beta.hat[,1]=mean(y)-beta.hat[,-1]%*%X.means
  
  #output relevant results
  return(list(beta=beta.hat,alpha=alpha.hat))
}

################################################
#X=rbind(
 # c(1,1,4),
  #c(5,3,5),
  #c(6,4,7),
  #c(6,4,1),
  #c(6,5,4),
  #c(6,7,3))
#y=c(6,8,6,7,5,4)

#print(our.lar(X,y)$beta[, -1]) 
#coefficients(lars::lars(X, y, type = 'lasso')) # tak samo wychodzi
