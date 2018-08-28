library(KernSmooth)
library(akima)

RENAME_brim<-function(EDGE,gene){
  BRIM<-paste0(gene[EDGE[1]],"_",gene[EDGE[2]])
  BRIM
}

BKDE_PMI<- function(edge,DataPre,method=c("integers","bandwidth")[2]){
  x <- DataPre[,edge[1]]
  y <- DataPre[,edge[2]]
  N <- max(c(x,y))
  
  Width <- c(bw.nrd0(x),bw.nrd0(y))
  
  gridSize <- switch(method,
                     integers  = c(N, N),
                     bandwidth = ceiling(N / c(min(Width[1]),min(Width[2]))))
  gridSize <- pmax(gridSize,10) # make sure there are at least 100 points in total
  
  BSmooth <- bkde2D(x=cbind(x, y), bandwidth=Width, gridsize=gridSize)
  
  USmoothx<-bkde(x)
  USmoothy<-bkde(y)
  
  BPre <- bicubic(x=BSmooth$x1, y=BSmooth$x2, z=BSmooth$fhat, x0=x,y0=y)
  UPrex <-aspline(USmoothx$x,USmoothx$y,x)
  UPrey <-aspline(USmoothy$x,USmoothy$y,y)
  
  # make sure there are no zeros in the smooth function (since we will take a log of that)
  BPre$z<- pmax(BPre$z, 1e-10)
  UPrex$y<-pmax(UPrex$y, 1e-10)
  UPrey$y<-pmax(UPrey$y, 1e-10)
  
  denpre <- log(BPre$z/(UPrex$y*UPrey$y))
  denpre
}

PMINR <- function(EDGE,ClassLabel,DataSet,ZSet=NULL){
  if(!is.null(ZSet)){
    PMINR.Z(EDGE=EDGE,ClassLabel=ClassLabel,DataSet=DataSet,ZSet=ZSet)
  }else{
    DenPreData <- apply(EDGE,1,BKDE_PMI,DataPre=DataSet)
    name<-colnames(DataSet)
    brim_gene<-apply(EDGE,1,RENAME_brim,name)
    colnames(DenPreData)<-brim_gene
    FitData<-cbind(DataSet,DenPreData)
    
    fit<-glm(ClassLabel~FitData,family = binomial(link = "logit"))
    
    Temp<-summary(fit)
    Result<-Temp$coefficients
    Sign<-matrix((Temp$coefficients[,"Pr(>|z|)"]<0.05)*1,ncol = 1)
    Result<-cbind(Result,Sign)
    rownames(Result)<-c("Intercept",colnames(FitData))
    Result
  }
}


PMINR.Z <- function(EDGE,ClassLabel,DataSet,ZSet){
  DenPreData <- apply(EDGE,1,BKDE_PMI,DataPre=DataSet)
  name<-colnames(DataSet)
  brim_gene<-apply(EDGE,1,RENAME_brim,name)
  colnames(DenPreData)<-brim_gene
  FitData<-cbind(DataSet,DenPreData)
  
  fit<-glm(ClassLabel~ZSet+FitData,family = binomial(link = "logit"))
  
  Temp<-summary(fit)
  Result<-Temp$coefficients
  Sign<-matrix((Temp$coefficients[,"Pr(>|z|)"]<0.05)*1,ncol = 1)
  Result<-cbind(Result,Sign)
  rownames(Result)<-c("Intercept",colnames(ZSet),colnames(FitData))
  Result
}

