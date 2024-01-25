Topic_Bayes_Mixture<-function(Y,X,K,T0,T,Opt_Hyper,LML){

if(missing(Opt_Hyper)) { 
             sp1=1;sp2=1;r1=1;r2=0.01; J = 0.01;cutoff=0.3
        } else {
            sp1=Opt_Hyper$sp1; sp2=Opt_Hyper$sp2;
		r1=Opt_Hyper$r1; r2=Opt_Hyper$r2;
		J=Opt_Hyper$J;cutoff=Opt_Hyper$cutoff;
        }   



#### sensitivity hyperparameters
#sp1=1;sp2=1;r1=1;r2=0.01; 

theta=rep(1,K);
liks=NULL;m=1;marls=NULL

#dimensions
n=length(Y);
p=dim(X)[2];

#instore posterior sampling
#postalp=NULL;
postz=array(0,c(p,K,T))
postmu=array(0,c(p,K,T))
posttau2=matrix(0,T,p);
postw=NULL;
postsig2=NULL;
postd=matrix(0,T,K);
postpx=NULL;
postH=matrix(0,T,n);
postnk=matrix(0,T,K);
postp=matrix(0,m,n);

setidT=array(0,c(p,p,K,T))

#initionation
d=theta/sum(theta);
tau2=rep(rinvgamma(1,sp1,sp2),p);
w=0.5;
sig2=rinvgamma(1,r1,r2);

z=matrix(as.numeric(runif(p*K)>w), p,K);
mu=matrix(rnorm(p*K,0,1),p,K);
H=sample(1:K, n, replace=T, prob=d)


#### Reording function from here #####
###############
K=ncol(mu)
p=nrow(mu)
n=length(H)
stop=FALSE;l=1:K
#id=1:K;
for (i in 1:p){
if (stop) break;
lc=l[order(mu[i,l])]
mutemp=mu;Htemp=H;dtemp=d;ztemp=z;
for (t in 1:length(lc)){
mu[,l[t]]=mutemp[,lc[t]];
z[,l[t]]=ztemp[,lc[t]];
H[Htemp==lc[t]]=l[t];
d[l[t]]=dtemp[lc[t]];
}
l=sort(l[mu[i,l]==0])
if (length(l)<=1) stop=TRUE 
}
mu=as.matrix(mu); H=H; z=as.matrix(z); d=d;
#### Reording function by here #####


V=matrix(0,p,K); 

X<-as.matrix(X) #just in case, 
## for computing distance-based neighbors
X1<-X 
X1[X>0]=1

print("running MCMC... the number of iterations are below.")

################## MCMC ##################
for (iter in 1:(T0+T)){
if (iter/100==round(iter/100,dig=0)) print(iter) # print number of itertations

nk=rep(0,K);id=1:n;
for(i in 1:n) {
nk[H[i]]=nk[H[i]]+1
} # n of k segment

for(k in 1:K){
idk=id[H==k]
for(j in 1:p) 
V[j,k]=sum(X[idk,j]^2)
}

setidS=array(0,c(p,p,K))

########## sample sig2 #############

res=0
for(i in 1:n)
res=res+(Y[i]-X[i,]%*%mu[,H[i]])^2;

ap=n/2+r1;
be=res/2+r2
sig2=rinvgamma(1,ap,be);


########## sample z and mu, w and tau_p2 #############
#sample z and mu
for (k in 1:K){
idk=id[H==k]

for (j in 1:p){

if (length(idk)==0) {dis_mm<-rep(0,(p))} # for safety to keep MCMC chain
if (length(idk)==1){ dis_mm<-sqrt((X1[idk,j]-X1[idk,])^2)}
if (length(idk)>1) {dis_mm <- sqrt(colSums((X1[idk,j]-X1[idk,])^2))} #Euclidean distance

thq<-quantile(dis_mm,cutoff) 

setid<-which(dis_mm<thq) ## members of neighbors
set <- z[setid,k] #a set of selected co-occuring words
Zk_1 <- sum(set==1) #how many of selected neighbor words

if(sum(setid)==0) {setidS[j,,k]<-0} else { 
setidS[j,(1:length(setid)),k]<-setid} #save set id

	# sample Z
	eps = Y[idk]-X[idk,-j]%*%mu[-j,k];
	tmpmean = sum(eps*X[idk,j])/sig2;
	tmpvar = V[j,k]/sig2+1/tau2[j]; 
	logr= (-J*Zk_1) + log(V[j,k]*tau2[j]/sig2+1)/2 - tmpmean^2/tmpvar/2
	if (logr>100) logr=100
	z[j,k]=ifelse(runif(1)< 1/(1+exp(logr)),1,0);

	# sample mu
	if (z[j,k]==0) mu[j,k]=0 else {
		bjmu = tmpmean/tmpvar;
		mu[j,k]=rnorm(1,bjmu,1/sqrt(tmpvar))
	}
}#for j
}#for k

#aw=sum(z); bw=p*K-aw;
#w=rbeta(1,aw+a, bw+b);

##taup2
at=rep(0,p)
bt=rep(0,p)

for (j in 1:p){
aw1=sum(z[j,])
bw1=(sum(mu[j,]^2))

at[j]=aw1/2+sp1
bt[j]=bw1*sum(z[j,])/2+sp2
tau2[j]=rinvgamma(1,at[j],bt[j]);
}


########## sample H and d ##################
d=as.vector(rdirichlet(1, nk+theta));

wp=matrix(0,K,n);
for(i in 1:n)
for(k in 1:K)
wp[k,i]=exp(-(Y[i]-X[i,]%*%mu[,k])^2/(2*sig2)) 

for (k in 1:K)
wp[k,]=wp[k,]#*d[k]

for (i in 1:n){
H[i]=sample(1:K,1,replace=TRUE,wp[,i])
}


########## reordering ##################
############### one more time after iteration
K=ncol(mu)
p=nrow(mu)
n=length(H)
stop=FALSE;l=1:K
#id=1:K;
for (i in 1:p){
if (stop) break;
lc=l[order(mu[i,l])]
mutemp=mu;Htemp=H;dtemp=d;ztemp=z;
for (t in 1:length(lc)){
mu[,l[t]]=mutemp[,lc[t]];
z[,l[t]]=ztemp[,lc[t]];
H[Htemp==lc[t]]=l[t];
d[l[t]]=dtemp[lc[t]];
}
l=sort(l[mu[i,l]==0])
if (length(l)<=1) stop=TRUE 
}
mu=as.matrix(mu); H=H; z=as.matrix(z); d=d;
#### Reording function by here #####


########## save posterior samples ########

if(iter>T0){
tmpT=iter-T0;
postz[,,tmpT]=z;
posttau2[tmpT,]=tau2
postmu[,,tmpT]=mu;
postw=c(postw,w);
postsig2=c(postsig2,sig2)
postd[tmpT,]=d
postH[tmpT,]=H;
postnk[tmpT,]=nk;

setidT[,,,tmpT]<-setidS ## selected neighbors per each variable across all MCMC runs

}


} # End of MCMC


## Posterior mean of beta ##

ss=matrix(0,p,K) 
for (i in 1:p)
for (j in 1:K)
{ 
ss[i,j]=mean(postmu[i,j,])
} 



## Posterior mean of Z ##

tt=matrix(0,p,K)
for (i in 1:p)
for (j in 1:K)
{ 
tt[i,j]=mean(postz[i,j,])
}

## membership calucation ##

memb=rep(0,n)
for (i in 1:n)
memb[i]<-as.numeric(names(table(postH[,i]))[which.max(table(postH[,i]))]) 
table(memb)




if(missing(LML))

list("ss"=as.matrix(ss),
		"memb"=memb,"tt"=as.matrix(tt),"setidT"=as.array(setidT))
else{

print("Calulating LML. It takes some time..")


#########################
#Log Marginal Likelihood#
#########################


#############################################Y
# Likelihood of Y|all others (mu, H, sigma2)#
#############################################

sigstar=sum(postsig2)/T

Ylik=0
for(i in 1:n){
meanstar=X[i,]%*%ss[,memb[i]]
temp1=dnorm(Y[i],meanstar,sqrt(sigstar),log=TRUE)
Ylik=Ylik+temp1;
}

Ylik #Y likelihood


###Prior Parts###


## prior of beta
prss=sum(dt(ss, df=sp1, ncp=0, log = TRUE))

## prior of sigma2
prsig=dgamma(sigstar,r1,r2,log=T)

pri=prss+prsig
pri


###Posterior parts###


# H* #

pH=matrix(0,n,T)

for (t in 1:T)
for (i in 1:n){
eps1=exp(-((Y[i]-X[i,]%*%ss[,memb[i]])^2)/(sigstar*2));
eps2=0
for (l in 1:K){
temp3=exp(-((Y[i]-X[i,]%*%ss[,l])^2)/(sigstar*2));
eps2=eps2+temp3
}
pH[i,t]=eps1/eps2
}
pH=log(pH)


pHH=apply(pH,2,sum)
pHH=mean(pHH)
pHH


#sig*#

psig2=rep(0,T)
for (t in 1:T){
res=0
for(i in 1:n){
res=res+(Y[i]-X[i,]%*%postmu[,postH[t,i],t])^2;
}
ap1=n/2+r1;
be1=res/2+r2
psig2[t]=dinvgamma(sigstar,ap1,be1);
}
#psig2=psig2[-which(psig2==0)]

mean(log(psig2))



#mu*#

T1=T

pmu=matrix(0,K,T1)

for (j in 1:p){

for (t in 1:T1){

for (k in 1:K){
idk1=id[postH[t,]==k]
if (sum(postz[j,k,t])==0) pmu[k,t]=0 else
{if (length(idk1)==1)Xidk=X[idk1,postz[j,k,t]==1] else {
Xidk=as.matrix(X[idk1,postz[,k,t]==1])
ss1=ss[(postz[,k,t]==1),k]
#sss=ss[,k]*postz[,k,t]
#ss1=sss[abs(sss)>0]
sp=sigstar/posttau2[t,j]
tempmean=(solve(sp*diag(ncol(Xidk))+t(Xidk)%*%Xidk))%*%(t(Xidk)%*%Y[idk1])
tempvar=solve((diag(ncol(Xidk))/posttau2[t,j])+(t(Xidk)%*%Xidk)/sigstar)
pmu[k,t]=dmvnorm(ss1,tempmean,tempvar)
}}

}#for k

} #for t

}# for j

ppmu=rep(0,T1)
for (t in 1:T1){
ppmu[t]=sum(log(pmu[,t][abs(pmu[,t])>0]))
}
ppmu=mean(ppmu)
ppmu

pos=pHH+ppmu+mean(log(psig2))

#logmarginal likelihood#

lml=Ylik+pri-pos

list("ss"=as.matrix(ss),
		"memb"=memb,"tt"=as.matrix(tt),"setidT"=as.array(setidT),"LML"=lml)}


} # for Topic-based Latent Class Model