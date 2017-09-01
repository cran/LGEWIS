
#########################
#   Hypothesis Testing
#########################

GA.prelim<-function(Y,time,X=NULL,corstr="exchangeable"){
  ##Preliminary
  Y<-as.matrix(Y);N<-nrow(Y)
  cluster.id<-unique(time[,1]);m<-length(cluster.id)
  #sort the data by cluster.id
  order.index<-order(match(time[,1],cluster.id))
  if (sum(time[,1]!=time[order.index,1])>0){
    msg<-sprintf("time was not sorted, now data is grouped by subject id")
    warning(msg,call.=F)
    Y<-as.matrix(Y[order.index,])
    if(length(X)!=0){X<-as.matrix(X[order.index])}
    time<-time[order.index,]
  }

  if(length(X)!=0){
    Z.A<-svd(as.matrix(X))$u
    nullgee<-geem(Y~Z.A,family=gaussian,corstr=corstr,id=time[,1])
  }else{
    nullgee<-geem(Y~1,family=gaussian,corstr=corstr,id=time[,1])
  }

  phi<-as.numeric(nullgee$phi)
  rho<-as.numeric(nullgee$alpha);N<-nrow(time) #add the intercept
  cluster.id<-unique(time[,1]);m<-length(cluster.id)

  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  V<-matrix(0,N,N);V.inv<-matrix(0,N,N);V.invsqrt<-matrix(0,N,N)
  for (i in 1:m)
  {
    ni<-n.rep[i]
    index<-n.total:(n.total+ni-1)
    n.total<-n.total + ni
    #working covariance matrix
    if (corstr=="exchangeable"){Vi<-diag(1-rho,ni)+ matrix(rho, ni, ni);Vi<-phi*Vi;Vi.inv<-solve(Vi);Vi.invsqrt<-Get.inv.sqrt(Vi)}
    if (corstr=="ar1"){Vi<-rho^abs(outer(time[index,2], time[index,2],"-"));Vi<-phi*Vi;Vi.inv<-solve(Vi);Vi.invsqrt<-Get.inv.sqrt(Vi)}
    if (corstr=="independence"){Vi<-diag(1,ni);Vi<-phi*Vi;Vi.inv<-solve(Vi);Vi.invsqrt<-Get.inv.sqrt(Vi)}
    V[index,index]<-Vi;V.inv[index,index]<-Vi.inv;V.invsqrt[index,index]<-Vi.invsqrt
    #if (link=="logit"){Ci<-mu[index]*(1-mu[index]);Vi.inv<-outer(sqrt(Ci),sqrt(Ci)^-1,"*")*Vi.inv}
  }
  X0<-cbind(rep(1,N),Z.A)
  #P<-diag(1,N)-X0%*%solve(t(X0)%*%V.inv%*%X0)%*%(t(X0)%*%V.inv)
  P1<-t(X0)%*%V.inv;P2<-X0%*%solve(P1%*%X0)

  #prepare the intermediate results
  result.prelim<-list(Y=Y,time=time,X=X,cluster.id=cluster.id,m=m,corstr=corstr,Z.A=Z.A,nullgee=nullgee,V=V,V.inv=V.inv,V.invsqrt=V.invsqrt,P1=P1,P2=P2)
  return(result.prelim)
}

#G is an m*q matrix, each row coresponds to one subject.
GA.test<-function(result.prelim,G,Gsub.id=NULL,weights='beta',B=5000,B.coef=NULL,impute.method='fixed'){
  ## Load preliminary data
  Y<-result.prelim$Y;corstr<-result.prelim$corstr
  m<-result.prelim$m;X<-result.prelim$X
  time<-result.prelim$time;cluster.id<-result.prelim$cluster.id
  Z.A<-result.prelim$Z.A #Z.A0=svd(cbind(rep(1,N),X, M.E))$u
  nullgee<-result.prelim$nullgee;N<-nrow(time)
  mu<-colSums(nullgee$beta*t(cbind(rep(1,N),Z.A)));Y.res<-Y-mu
  V<-result.prelim$V;V.inv<-result.prelim$V.inv;V.invsqrt<-result.prelim$V.invsqrt
  P1<-result.prelim$P1;P2<-result.prelim$P2

  ## Deal with the genotype
  SNP.list<-colnames(G)
  # match the phenotype and genotype subject ID
  if(length(Gsub.id)==0){Gsub.id<-cluster.id}
  G<-as.matrix(G[match(cluster.id,Gsub.id),])

  # missing genotype imputation
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-Impute(G,impute.method)
  }
  # genotype
  G<-as.matrix(G);center.G<-t(t(G)-colMeans(G))
  MAF<-colMeans(as.matrix(G[,colMeans(center.G^2)!=0]))/2;MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  G<-as.matrix(G[,colMeans(center.G^2)!=0])
  SNP.name<-SNP.list[colMeans(center.G^2)!=0]

  if(length(weights)>1){G<-t(t(G)*weights[colMeans(center.G^2)!=0])}
  if(length(weights)==1){
    if(weights=="rare"){G<-G[,MAF<0.05];SNP.name<-SNP.name[MAF<0.05];MAF<-MAF[MAF<0.05]}
    if(weights=="common"){G<-G[,MAF>0.05];SNP.name<-SNP.name[MAF>0.05];MAF<-MAF[MAF>0.05]}
    if(weights=="beta"){weights<-dbeta(MAF,1,25);G<-t(t(G)*weights)}
  }
  G<-as.matrix(G)

  if(ncol(G) == 0){
    msg<-sprintf("G does not include any heterozygous variant")
    stop(msg)
  }

  # generate long form genotype
  Z<-as.matrix(G[match(time[,1],cluster.id),])
  Z<-Z-P2%*%(P1%*%Z)
  n.G<-ncol(Z)

  GA.fit<-GA.bootstrap(nullgee,Y,time,Z.A,Z,corstr,B,B.coef) #Z.A does not include the intercept
  score<-GA.fit$score;B.score<-t(GA.fit$B.score)
  #V<-GA.fit$V;V.inv<-GA.fit$V.inv

  Q1<-sum(score^2);Q2<-sum(score)^2
  re.Q1<-apply(B.score^2,1,sum);re.Q2<-apply(B.score,1,sum)^2

  #cauculated resampled p-values using resampled test statistics
  #rho.class<-c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1)
  rho.class<-c(0,1)

  all.p<-matrix(0,B+1,length(rho.class))
  for (rho in rho.class){
    Q.rho<-(1-rho)*Q1+rho*Q2#(1-rho)*Q1+rho*Q2
    re.Q.rho<-(1-rho)*re.Q1+rho*re.Q2
    all.p[,which(rho.class==rho)]<-Get.p(as.matrix(c(Q.rho,re.Q.rho)),as.matrix(re.Q.rho))
  }
  re.p<-all.p[-1,] # resampled p-values
  temp.p<-all.p[1,] # p-values

  ###
  #   Combine p-values using Fisher's method
  ###
  p.Fisher<-FCombine.p(temp.p,re.p)
  ###
  #   Combine p-values using MinP method
  ###
  p.MinP<-MCombine.p(temp.p,re.p)

  p.value<-matrix(c(temp.p,p.Fisher,p.MinP),1,length(rho.class)+2)
  colnames(p.value)<-c(rho.class,'Fisher','MinP')

  # implement a single SNP based analysis by GEE score tests
  B.variance.single<-apply(B.score,2,var)
  p.single<-cbind(MAF,as.vector(pchisq(as.vector(score)^2/B.variance.single,df=1,lower.tail=F)))
  colnames(p.single)<-c('MAF','p.value')
  rownames(p.single)<-SNP.name

  return(list(n.marker=n.G,p.value=p.value,p.single=p.single))
}

#################################################
#   Calculate Score Statistics and their Variance
#################################################

GA.bootstrap<-function(nullgee,Y,time,Z.A,Z.I,corstr,B,B.coef) #Z.A does not include the intercept
{
  phi<-as.numeric(nullgee$phi)
  rho<-as.numeric(nullgee$alpha);N<-nrow(time);Z.A<-cbind(rep(1,N),Z.A) #add the intercept
  mu<-colSums(nullgee$beta*t(Z.A));Y.res<-Y-mu
  #rho<-nullgee$geese$alpha;mu<-nullgee$fitted.values;Y.res<-Y-mu
  cluster.id<-unique(time[,1]);m<-length(cluster.id);dimI<-ncol(Z.I)
  ##Generalized score test
  score.array<-matrix(0,dimI,m)
  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  #V<-matrix(0,N,N);V.inv<-matrix(0,N,N)
  for (i in 1:m)
  {
    ni<-n.rep[i]
    index<-n.total:(n.total+ni-1)
    n.total<-n.total + ni
    #working covariance matrix
    if (corstr=="exchangeable"){Vi<-diag(1-rho,ni)+ matrix(rho, ni, ni);Vi<-phi*Vi;Vi.inv<-solve(Vi)}
    if (corstr=="ar1"){Vi<-rho^abs(outer(time[index,2], time[index,2],"-"));Vi<-phi*Vi;Vi.inv<-solve(Vi)}
    if (corstr=="independence"){Vi<-diag(1,ni);Vi<-phi*Vi;Vi.inv<-solve(Vi)}
    #V[index,index]<-Vi;V.inv[index,index]<-Vi.inv
    #if (link=="logit"){Ci<-mu[index]*(1-mu[index]);Vi.inv<-outer(sqrt(Ci),sqrt(Ci)^-1,"*")*Vi.inv}
    #score matrices
    Z.Ii<-Z.I[index,];if (ni==1) {Z.Ii<-t(Z.Ii)}
    score.array[,i]<-t(Z.Ii)%*%(Vi.inv%*%Y.res[index,])/sqrt(m)
  }
  score<-as.matrix(apply(score.array,1,sum))
  # null distribution
  if(length(B.coef)==0){B.coef<-matrix(rbinom(m*B,1,0.5)*2-1,m,B)}
  B.score<-score.array%*%B.coef

  return(list(score=score,B.score=B.score))#Q: test statistic; M: covariance matrix
}

#cauculated p-values using resampled test statistics
Get.p<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
  re.mean<-apply(re.Q,2,mean)
  re.variance<-apply(re.Q,2,var)
  re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  re.p<-t(1-pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df))
  return(re.p)
}

# combine p-values
FCombine.p<-function(p,re.p){ #re.p: a B*b*c array; function to combine multiple p-values
  Fisher.stat<--2*sum(log(p))
  re.Fisher<--2*apply(log(re.p),1,sum);re.Fisher[re.Fisher==Inf]<-NA
  Fisher.mean<-mean(re.Fisher,na.rm=T)
  Fisher.variance<-var(re.Fisher,na.rm=T)
  Fisher.kurtosis<-mean((re.Fisher-Fisher.mean)^4,na.rm=T)/Fisher.variance^2-3
  if (Fisher.kurtosis>0){df<-12/Fisher.kurtosis}else{
    df<-100000
  }
  p.combined<-1-pchisq((Fisher.stat-Fisher.mean)*sqrt(2*df)/sqrt(Fisher.variance)+df,df)
  return(p.combined)
}

MCombine.p<-function(p,re.p){ #re.p: a b*c matrix; function to combine multipla p-values for a minP test
  MinP.stat<-min(p)#;re.p[re.p==1]<-0.99
  re.Normal<-qnorm(re.p);re.Normal[re.Normal==Inf]<-NA
  D<-cor(re.Normal,use='complete.obs')
  #diag(D)<-1
  p.combined<-as.numeric(1-pmvnorm(lower=rep(qnorm(MinP.stat),length(p)),sigma=D))
  return(p.combined)
}

# calculate number of independent tests
effect.n<-function(x,MinP.adjust){
  temp<-0;sum.EV<-sum(x) #summation of eigen values
  for (i in 1:length(x)){
    temp<-temp+x[i];if (temp>sum.EV*MinP.adjust){break}
  }
  return(i)
}

#Data management
Minor.allele<-function(x){if(mean(x)>1){x<-2-x};return(x)}
Standardize<-function(x){return((x-mean(x))/sd(x))}
Center<-function(x){return(x-mean(x))}
Variation<-function(x){x<-apply(x,2,Minor.allele);x<-x[,which(apply(x,2,sd)>0)];return(x)}
Common<-function(x){x<-apply(x,2,Minor.allele);freq<-apply(x,2,mean)/2;x<-x[,which(freq>0.05)];return(x)}
##Matrix calculation
#Get sqrt.root of a matrix
Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
#Get inverse of a matrix; numerically robust
Get.inverse<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix inverse!")}
  a.inverse <- a.eig$vectors[,ID1] %*% diag(a.eig$values[ID1]^-1) %*% t(a.eig$vectors[,ID1])
  return(a.inverse)
}
#Get -1/2 of a matrix; numerically robust
Get.inv.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])^-1,length(ID1)) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
#Time similarity
Check.same<-function(x,y){return(as.numeric(x==y))}
Check.near<-function(x,y){return(as.numeric(abs(x-y)==1))}
D.matrix<-function(X){
  n<-length(X[,1]);temp<-matrix(0,n,n)
  #checking time
  temp<-outer(X[,1],X[,1],Check.same)*outer(X[,2],X[,2],Check.near)
  diag(temp)<-0
  return(temp)
}
# Simple Imputation (from SKAT package)
# Z : an m x p genotype matrix with m samples and p SNPs

Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}




