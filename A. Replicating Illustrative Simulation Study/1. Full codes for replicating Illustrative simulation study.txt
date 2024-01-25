###################################################################################################################
### This is a full demonstration R code for replicating illustrative simulation study in the Kim et al. (2024). ###
###################################################################################################################
## by runnig the below code from top to bottom, you can see the replication results ##


#Package Loading#

library(mvtnorm)
library(MCMCpack)
library(cluster)


#######################
#Function Reordering###
#######################
reorder<-function(mu,H,z,d){
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
list("mu"=as.matrix(mu),"H"=H,"z"=as.matrix(z),"d"=d)
}#reorder

set.seed(1003)

#######################################################
## Generating Simulation Data p=100, n=1000, 3 topic across 3 segments
#######################################################

## initialization
p=100; ## number of IVs; n number of samples;
n=1000; ## n number of samples
K=3; ## K number of Segments
p1=5  # number of words within each dimensions
G=3 #number of dimensions per segment

d_t=rep(1/K,K); 
X=matrix(sample(c(1,0), p*n, replace = TRUE, prob=c(0.4,0.6)),n,p)
tau2_t=1;
z_t=matrix(0,p,K);
mu_t=matrix(0,p,K);
H_t=sample(1:K,n,replace=TRUE,d_t)

## reordering ##
newp=reorder(mu_t,H_t,z_t,d_t);
mu_t=newp$mu;
H_t=newp$H;
z_t=newp$z;
d_t=newp$d;

Y=rep(0,n);
id=1:n
## generate word matrix ##

idk_s=list()
dm=list()

for (k in 1:K){

group=list()
gr.f <- sample(2:p, (p1*G) , replace = FALSE) #group features

for (g in 1:G){
group[[g]] <- gr.f[(g+(4*(g-1))):(g*p1)] ##a dimension of the selected variables (cooccurring words)
} # for g

dm[[k]]<-group

idk_s[[k]]=id[H_t==k]
idkk <- idk_s[[k]]

idnk<-length(idkk)
tmpX<-X[idkk,]

for (g in 1:G){

dmk <- unlist(dm[[k]][g]) # group variables

coccurringN <- round(idnk*0.7,0)
coccset <- sample(idnk, coccurringN, replace=FALSE)
coccN <- length(coccset)

for (v in 1:idnk){

if (v %in% coccset) {
tmpX[v,dmk]=c(sample(c(1,0), p1, replace = TRUE, prob=c(0.9,0.1))) #co-occur
} else {
tmpX[v,dmk]=sample(c(0,1), p1, replace = TRUE, prob=c(0.7,0.3))
} # else

} ## v
} #for g

X[idkk,]=tmpX

} # for k -- Word Matrix

X1<-X
X1[X>0]=1

X<-(X*sample(1:5, n*p, replace = TRUE, prob=c(0.4,0.3,0.1,0.1,0.1)))

check=rep(0,n)
for (i in 1:n){
check[i]=sum(X[i,])
}

which(check==0) #checking if the word matrix includes any rows with all zeros

ps<-round(p*0.01,0)

mu_t=matrix(0,p,K); 

## set intercept values
mu_t[1,1]=1 
mu_t[1,2]=2
mu_t[1,3]=3

for (k in 1:K){
for (g in 1:G){

mu_t[dm[[k]][[g]],k]=runif(p1,1,3); 

}}

ps<-round(p*0.05,0)

z_t=matrix(0,n,K)
z_t[mu_t>0]=1

mu_t1<-mu_t
z_t1<-z_t

sig2_t=0.01

## generating Y (DV)
for(i in 1:n)
Y[i]=X[i,]%*%mu_t[,H_t[i]]+rnorm(1,0,sqrt(sig2_t))
Y_t=Y

######### Data Generation by here ####


##Pre-specifying############################

K=3 # number of Segments
T0=2000 #N of Burn-in
T=2000 #N of posterior Sampling


####################################################################
###### Latent Class Regression with Group Variable Selection from here 
####################################################################

#### Setting hyperparameters
sp1=1;sp2=1;r1=1;r2=0.01; 

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

#initialization
d=theta/sum(theta);
tau2=rep(rinvgamma(1,sp1,sp2),p);
w=0.5;
sig2=rinvgamma(1,r1,r2);
z=matrix(as.numeric(runif(p*K)>w), p,K);
mu=matrix(rnorm(p*K,0,1),p,K);
H=sample(1:K, n, replace=T, prob=d)

########## reordering ###
newp=reorder(mu,H,z,d);
mu=newp$mu;
H=newp$H;
z=newp$z;
d=newp$d;

V=matrix(0,p,K); 

X1<-X
X1[X>0]=1

error=sample(1:3,1)
cutoff=(p1+error)/p # cutoff point

J = 0.01 # tuning parameter (lambda in paper). As J increases, power for selecting group variables becomes stonger 

##########################################
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
if (length(idk)==1){ dis_mm<-(abs(X1[idk,j]-X1[idk,]))}
if (length(idk)>1) {dis_mm <- (colSums(abs(X1[idk,j]-X1[idk,])))} #distance

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
wp[k,]=wp[k,]*d[k]

for (i in 1:n){
H[i]=sample(1:K,1,replace=TRUE,wp[,i])
}

########## reordering ##################

newp=reorder(mu,H,z,d);
mu=newp$mu;
H=newp$H;
z=newp$z;
d=newp$d;


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
} 

#### End of MCMC ####


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

## finding membership for H 
memb=rep(0,n)
for (i in 1:n)
memb[i]<-as.numeric(names(table(postH[,i]))[which.max(table(postH[,i]))]) 

## print coeffs, variable selection and memberships of segments

round(ss,3)
round(tt, 2)
memb 


#########################################
## post-hoc searching         ###########
## finding topics (neighbor structure)###
#########################################


dim(setidT) ## selected neighbors per each variable

joint_D <- array(0,c(p,p,K)) 

for (k in 1:K){
for (j in 1:p){
for (l in 1:p){
joint_D [j,l,k] <- sum(setidT[j,,k,]==l)
}
}
}

for (k in 1:K)
diag(joint_D[,,k])<-T  ## segment-level neighbors

## Heuristic Search ## 
cut=0.95
P_dim=array(0,c(p,p,K))
#S_sel=array(0,c(p,p,K))

## find segmentment-level neighborsets!

for (k in 1:K){
for(j in 1:p){
pdim = which(joint_D[j,,k] > (cut*T)) # probabilistic segment-level neighbors 
sel_t <-which(tt[,k]>0.95)
pdim0<-intersect(pdim,sel_t) # selected neighbors
ssel<-setdiff(sel_t, pdim0)
if(sum(pdim0)>0){
P_dim[j,1:length(pdim0),k]<-pdim0} #Co-occuring selected variables; as cooccuring variables are active...
}
}

#selected neighbors per each j word

############################################################
#### Generating dendrogram of neighbors for each segment####
############################################################

k=1 # for each segment; please change segment from 1, 2 and 3 (e.g., k=1, k=2 and k=3)
 
dist_R=rep(0,K)

symm_M<-matrix(0,p,p)

for (j in 1:p){

tmp_i<- P_dim[j,,k]
symm_M[j,tmp_i]=1
}
colnames(symm_M)<-as.character(1:p)
rownames(symm_M)<-as.character(1:p)

sel_t <-which(tt[,k]>0.95)

symm_M0<-symm_M[sel_t[-1],sel_t[-1]] #### a matrix of the neighbor structure of selected variables

## each row includes its selected neighbors in columns

group.dist=rep(0,K)

seg.dist1<-as.matrix(dist(symm_M0 ))

dim(as.matrix(seg.dist1)) #distances between 300 members

seg.hc<-hclust(as.dist(seg.dist1), method="complete")
plot(seg.hc) #resulting tree dendrogram

rect.hclust(seg.hc, k=3, border="red")

## Now, let's compare and check the identified dimensions (groups of frequently co-occurring words) with true dimensions

dm[[k]] ## true group words of per dimension k


