# PMINR
Pointwise mutual information-based biological network regression model<br>
## Usage<br>
R code avaiable from https://github.com/LinSmallstrong/PMINR/blob/master/PMINR.R <br>

    PMINR(EDGE,ClassLabel,DataSet,ZSet=NULL)

## Arguments:<br>

    * DataSet: data matrices containing one sample per row, one variable per column.  
    * ClassLabel: must be 0 or 1, e.g. 1 for cases and 0 for controls. 
    * EDGE: array indices. which edge will be estimated.
    * ZSet: covariates in data containing one sample per row, one variable per column.

## Example:
run the R code `PMINR.r` first<br>

    library(mvtnorm)
    sigma<-c( 1 ,0 ,1 ,1 ,0 ,  	
              0 ,1 ,0 ,1 ,0 ,	 
              1 ,0 ,1 ,0 ,1 , 	
              1 ,1 ,0 ,1 ,1 ,	 
              0 ,0 ,1 ,1 ,1 ) 
              
    sigma0<-sigma 
    sigma0<-matrix(sigma0,nrow = 5) 
    sigma<-matrix(sigma,nrow = 5) 
    mu<-c(0 , 0 , 0 , 0 , 0 ) 
    sigma0[lower.tri(sigma0,diag=T)]<-0 
    EDGE <- which(sigma0!=0,arr.ind=T) 
    EDGE 

          row col 
    [1,]   1   3 
    [2,]   1   4 
    [3,]   2   4 
    [4,]   3   5 
    [5,]   4   5 

    set.seed(2018) 
    size0 <- size1 <- 50 
    class0 <- rmvnorm( n = size0, sigma = sigma,method = "svd") 
    class1 <- rmvnorm( n = size1, sigma = sigma,method = "svd") 
    DataSet <- rbind(class0,class1) 
    colnames(DataSet)<-c("a","b","c","d","e") 
    classLabel <- c(rep(0,size0),rep(1,size1)) 

    PMINR(EDGE=EDGE,ClassLabel=classLabel,DataSet=DataSet,ZSet=NULL) 

                 Estimate Std. Error    z value   Pr(>|z|)  
    Intercept -0.29206477  0.2932381 -0.9959987 0.31925077 0
    a         -0.22610544  0.2172268 -1.0408728 0.29793457 0
    b         -0.13009524  0.2663512 -0.4884350 0.62524175 0
    c         -0.33364520  0.2797345 -1.1927208 0.23297874 0
    d         -0.03800546  0.2317497 -0.1639936 0.86973620 0
    e         -0.52132066  0.2620783 -1.9891793 0.04668141 1
    a_c        0.78081787  1.1830310  0.6600147 0.50924438 0
    a_d       -0.93239239  1.0874868 -0.8573828 0.39123339 0
    b_d       -0.26872979  0.8617785 -0.3118316 0.75516849 0
    c_e        0.80331533  0.9525449  0.8433359 0.39904059 0
    d_e        2.73561242  1.3415043  2.0392125 0.04142882 1
    

