
######################################################
## Extract Word Frequency Matrix from Review Texts ###
######################################################

Data0<-matrix_Brand

n_reviews<-as.numeric(Data0[,2])

sum(as.numeric(Data0[,2]))


n=dim(Data0)[1]

subdata=data.frame(doc_id=as.character(1:n),text=Data0[,4]) #extract only textual reviews


TMD=VCorpus(VectorSource(subdata$text)) #without DataframeSource...

#Representing and computing on corpora. 
#Corpora are collections of documents containing (natural language) text

TMD[[1]]$content
TMD[[1]]$meta ## additional informaiton for this text

##some preprocessing
TMD <- tm_map(TMD, content_transformer(tolower)) ## change to lower case
TMD <- tm_map(TMD, removeWords, stopwords("english")) #remove stopwords
TMD <- tm_map(TMD, stemDocument, language="english") ## stemming words
TMD <- tm_map(TMD, stripWhitespace) #Strip extra whitespace from a text document.
TMD <- tm_map(TMD, removeNumbers)
TMD <- tm_map(TMD, removePunctuation)
TMD <- tm_map(TMD, removeWords, c("does","not","thing")) #remove additional unnecessary words selected by users. 



## creat word matrix ###
TMDM=TermDocumentMatrix(TMD,control = list(removePunctuation = TRUE,stopwords = TRUE))


## Compute frequency distance ##
dtmss <- removeSparseTerms(TMDM, 0.9) ## search frequent words

dtmss0<-t(as.matrix(dtmss))

dim(dtmss0)

range(as.numeric(Data0[,2]))





###########################################################################
## Selecting a set of Frequent words commonly appeared across documents ###
###########################################################################

TMDM_freq=list()

for(i in 1:n_B){

if (i/100==round(i/100,dig=0)) print(i) # print number of itertations

TMDM_freq[[i]]=findFreqTerms(TMDM[,i], round(n_reviews[i]/2,0)) #average 1 appear per 2 review

} #proportionally sampling 10 reviews 5 appearance


length(TMDM_freq) #check how many frequent words

#stemCompletion(TMDM_freq,TMD)

U_Features=NULL

for(i in 1:n){

U_Features = union(U_Features,TMDM_freq[[i]] ) ##Union words

}

length(U_Features)
nf=length(U_Features)

mul.words<-rep(0,nf)

for (j in 1:nf){

tempc<-rep(0,n)

for (i in 1:n){
if (U_Features[j] %in% TMDM_freq[[i]]){tempc[i]=1}
} #appear at ith review

mul.words[j]=sum(tempc)

}

sum(mul.words>10)

wl<-which(mul.words>10) # more than 10 restaurants

m<-length(wl)

words_features<-U_Features[wl]

TMDM0<-matrix(0,n,m)
colnames(TMDM0)<-words_features

for (j in 1:m){

ll<-which(colnames(dtmss0)==(words_features[j]))

if(identical(ll, integer(0))){TMDM0[,j]=0}else{

TMDM0[,j] <- dtmss0[,ll]
}
}

remove_c=rep(0,m)

for (j in 1:m){

remove_c[j]=sum(TMDM0[,j])

}

which(remove_c==0)

remove_r=rep(0,n)

for (i in 1:n){

remove_r[i]=sum(TMDM0[i,])

}

which(remove_c==0)
which(remove_r==0)

TMDM_f <- TMDM0[,-which(remove_c==0)]

######
