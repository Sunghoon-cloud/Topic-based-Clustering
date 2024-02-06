########################
## finding neighbor structure
##############################


## finding dimensions -- post-hoc searching ##

Dendrograms<-function(setidT,tt,K){

p=dim(tt)[1]

joint_D <- array(0,c(p,p,K)) 

for (k in 1:K){
for (j in 1:p){
for (l in 1:p){
joint_D [j,l,k] <- sum(setidT[j,,k,]==l)
}
}
}

for (k in 1:K)
diag(joint_D[,,k])<-T ## segment-level neighbors

## Heuristic Search ## 
cut=0.95
P_dim=array(0,c(p,p,K))
S_sel=array(0,c(p,p,K))

## find segmentment-level neighborsets!


for (k in 1:K){
for(j in 1:p){
pdim = which(joint_D[j,,k] > (cut*T)) # probabilistic segment-level neighbors 
sel_t <-which(tt[,k]>0.95)
pdim0<-intersect(pdim,sel_t) # selected neighbors
ssel<-setdiff(sel_t, pdim0)
if(sum(pdim0)>0){
P_dim[j,1:length(pdim0),k]<-pdim0} #Co-occuring selected variables; as cooccuring variables are active...
if(sum(ssel)>0){
S_sel[j,1:length(ssel),k]<-ssel
} 
}
}

#selected neighbors per each j word


# DV: % of cooccuring variables are successfully recovered.





#### for simulation data ####

dist_R=rep(0,K)

  par(mfrow = c(1, 2))

for (k in 1:K){

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
# complete linkage method evaluates the distance between every member
plot(seg.hc) #resulting tree for all N=300 observations of seg.df.

#rect.hclust(seg.hc, k=2, border="red") ## dendrogram

#p.plot1<-heatmap(as.matrix(seg.dist1), col=grey.colors(20, start=0.2, end=0.8), labCol=FALSE) # heatmap

}

} #end_Dendrogram


