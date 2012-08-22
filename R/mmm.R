mmm <-
function(data,nresp,family='gaussian',corstr='independence',coefnames=NULL,tol=0.001,maxiter=25,Mv=1,silent=TRUE){

# Number of coefficients for a univariate response
ncoef<-dim(data)[2]-nresp-1+1
mod<-length(coefnames)%%ncoef
if (mod!=0) stop("Length of coefficient names are not multiple of number of coefficients")

# id of subjects
id2<-unique(data[,1])   

# reconstructing the response matrix
resp<-NULL
for (i in id2){
resp1<-data[data[,1]==i,2:(1+nresp)]   
resp2<-NULL
for (j in 1:nresp){
resp2<-c(resp2,resp1[,j])
}
resp<-c(resp,resp2)
}
resp<-matrix(resp)

# reconstructing the covariate matrix
# additionally 1s in the 1st column
# to have different intercepts for each response
ones<-rep(1,dim(data)[1])
covmat<-cbind(data[,1],ones,data[,(1+nresp+1):dim(data)[2]])
cov3<-NULL
for (k in id2){
cov1<-covmat[covmat[,1]==k,2:dim(covmat)[2]]
cov2<-kronecker(diag(1,nresp),as.matrix(cov1))
cov3<-rbind(cov3,cov2)
}

# reconstructing the id column
id<-NULL
for (t in id2){
id3<-data[data[,1]==t,1]
id4<-rep(id3,nresp)
id<-c(id,id4)
}

# loading gee library
library(gee)

## Fitting Multivariate Marginal Model
gee1<-gee(resp~-1+cov3,id=id,family=family,corstr=corstr,tol=tol,maxiter=maxiter,Mv=Mv,silent=silent)
summary1<-summary(gee1)
# giving coefficient names if the user already defined them
if (length(coefnames)!=0){
row.names(summary1$coefficients)<-coefnames
}
list1<-list("Multivariate Marginal Model",multivout=summary1)


# Fitting Univariate Marginal Models
list2<-list()
for (z in 1:nresp){
gee2<-gee(data[,1+z]~as.matrix(data[,(1+nresp+1):dim(data)[2]]),id=data[,1],family=family,corstr=corstr,tol=tol,maxiter=maxiter,Mv=Mv,silent=silent)
summary2<-summary(gee2)
# giving coefficient names if the user already defined them
if (length(coefnames)!=0){
row.names(summary2$coefficients)<-coefnames[((z-1)*(ncoef)+1):(z*ncoef)]
}
list2[[z]]<-summary2
}
list2<-list("Univariate Marginal Models",univout=list2)

# Ultimate output including outputs of multivariate and univariate models
list3<-list(multiv=list1,univ=list2)
list3
}

