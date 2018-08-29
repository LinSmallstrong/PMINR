# PMINR
Pointwise mutual information-based biological network regression model<br>

    PMINR(EDGE,ClassLabel,DataSet,ZSet=NULL)<br>

## Arguments:<br>

    * DataSet: data matrices containing one sample per row, one variable per column. <br>
    * classLabel: must be 0 or 1, e.g. 1 for cases and 0 for controls.<br>
    * EDGE: array indices. which edge will be tested.<br>
    * ZSet:covariate data containing one sample per row, one variable per column.<br>

## Example:

    library(mvtnorm)<br>
    sigma<-c( 1 ,0 ,1 ,1 ,0 ,     <br>		
              0 ,1 ,0 ,1 ,0 ,	<br>	
              1 ,0 ,1 ,0 ,1 ,	<br>	
              1 ,1 ,0 ,1 ,1 ,	<br>	
              0 ,0 ,1 ,1 ,1 )	<br>
    sigma0<-sigma<br>
    sigma0<-matrix(sigma0,nrow = 5)<br>		
    sigma<-matrix(sigma,nrow = 5)<br>
    mu<-c(0 , 0 , 0 , 0 , 0 )<br>
    sigma0[lower.tri(sigma0,diag=T)]<-0<br>
    EDGE <- which(sigma0!=0,arr.ind=T)<br>
    EDGE<br>

          row col<br>
    [1,]   1   3<br>
    [2,]   1   4<br>
    [3,]   2   4<br>
    [4,]   3   5<br>
    [5,]   4   5<br>

    set.seed(2018)<br>
    size0 <- size1 <- 50<br>
    class0 <- rmvnorm( n = size0, sigma = sigma,method = "svd")<br>
    class1 <- rmvnorm( n = size1, sigma = sigma,method = "svd") <br>
    DataSet <- rbind(class0,class1)<br>
    colnames(DataSet)<-c("a","b","c","d","e")<br>
    classLabel <- c(rep(0,size0),rep(1,size1))<br>

    PMINR(EDGE=EDGE,ClassLabel=classLabel,DataSet=DataSet,ZSet=NULL)<br>

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
    

