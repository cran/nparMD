#' @title Nonparametric Test For Multivariate Data With Two-Way Layout Factorial Design - Large Samples
#'
#' @description Analysis of multivariate data with two-way completely randomized
#' factorial design - version for large samples. The analysis is based on fully nonparametric, rank-based methods and
#' uses an F-approximation for 'Dempster's ANOVA' and a chisquare-approximation for the criteria called 'Wilks Lambda', 'Lawley-Hotelling' and
#' 'Bartlett-Nanda-Pillai'.
#' These approximations are given by the asymtotic distribution
#' of these statistics under true null-hypothesis. In contrast to the normal-approximated test (as used in the small sample version) it is designed for
#' data with large samples (see details) while the number of factorial levels is allowed to be small. The multivariate response is allowed to be ordinal, quantitative, binary or a mixture
#' of the different variable types. The test statistics are constructed using nonparametric relative effect estimators.
#' @usage nparml(formula, data)
#' @param formula an object of class "formula" with two explanatory variables (factors), see examples.
#' @param data an object of class "data.frame" containing the variables in the formula
#' @return Returns a list of data frames providing the values of the test statistics, p-values, degrees of freedom, factor levels,
#' and groupsize per factor level combination.
#' @details
#' The data is analysed for main effects
#' and interaction effect of the explanatory factors. In each case the null hypothesis "no effect" is testet. In order to obtain
#' reliable results the considered data should include at least 7 observations per factor level combination.
#' This method is only implemented for complete data sets without missing values.
#' @references 
#' Kiefel M., Freidl J. (2024)
#' Nonparametric Analysis of Multivariate Data in Factorial Designs with Nondetects: A Case Study with Microbiome Data
#' In: Journal of Agricultural, Biological, and Environmental Statistics
#' https://doi.org/10.1007/s13253-024-00671-5.
#' @references
#' Kiefel M., Bathke A.C. (2022)
#' Fully Nonparametric Methods for Multivariate Data in Factorial Designs. Asymptotics, Finite Sample Approximations, and Implementation in R
#' In: Open Statistics, vol. 3, 2022, p. 74
#' https://doi.org/10.1515/stat-2022-0112
#' @references
#' Kiefel M., Bathke A.C. (2020)
#' Rank-Based Analysis of Multivariate Data in Factorial Designs and Its Implementation in R
#' In: Nonparametric Statistics (285-294)
#' Springer Proceedings in Mathematics & Statistics
#' Springer International Publishing, Cham
#' @references
#' Bathke A.C., Harrar S.W. (2016)
#' Rank-Based Inference for Multivariate Data in Factorial Designs.
#' In: Liu R., McKean J. (eds) Robust Rank-Based and Nonparametric Methods.
#' Springer Proceedings in Mathematics & Statistics, vol 168. Springer, Cham
#' @references
#' Harrar S.W., Bathke A.C. (2012)
#' A modified two-factor multivariate analysis of variance:
#' asymptotics and small sample approximations (and erratum).
#' In: Annals of the Institute of Statistical Mathematics, 64(1&5):135-165&1087, 2012.
#' @references
#' Brunner E., Dette H., Munk A. (1997)
#' Box-Type Approximations in Nonparametric Factorial Designs
#' In: Journal of the American Statistical Association, 92(440):1494-1502
#' @examples
#' data(pseudostudy1)
#' nparml(resp1|resp2|resp3~treatment*age, pseudostudy1)
#' @import Formula
#' @import matrixcalc
#' @import matrixStats
#' @import stats
#' @importFrom MASS ginv
#' @importFrom gtools permutations
#' @importFrom methods is
#'
#' @export
nparml<-function(formula, data){

  if(!is.data.frame(data)){stop("data argument is not a data frame")}
  if(!is(formula, "formula")){stop("Please give a formula")}

  formel<-Formula(formula)
  data<-model.frame(formel,data)

  p<-ncol(data)-2
  N<-nrow(data)

  data[,p+1]<-factor(data[,p+1])
  data[,p+2]<-factor(data[,p+2])

  ifelse(length(levels(data[,p+1]))>length(levels(data[,p+2])),
         {NameFactorA<-colnames(data)[p+1]; NameFactorB<-colnames(data)[p+2]},
         {NameFactorA<-colnames(data)[p+2]; NameFactorB<-colnames(data)[p+1]}
  )

  data<-data[order(data[,NameFactorA],data[,NameFactorB]),]

  struc<-table(data[,NameFactorA],data[,NameFactorB])
  a<-nrow(struc)
  b<-ncol(struc)
  groups<-as.vector(t(struc))

  data<-t(data[c(-(p+2),-(p+1))])
  
  #RankTransform matrix pxN
  RT<-(rowRanks(data,ties.method="average")-0.5)/N		

  # Contrast Matrices

  P.a <- diag(1,a)-(rep(1,a) %*% t(rep(1,a)))/a
  P.b <- diag(1,b)-(rep(1,b) %*% t(rep(1,b)))/b

  C.A <- kronecker(P.a, (rep(1,b) %*% t(rep(1,b)))/b)
  C.AB <-kronecker(P.a, P.b)
  C.B <- kronecker((rep(1,a) %*% t(rep(1,a)))/a, P.b)

  T.A<-kronecker( C.A , diag(p))
  T.B<-kronecker( C.B, diag(p))
  T.AB<-kronecker(C.AB, diag(p))

  #### Data grouping and means	####
  
  #group (column) index 
  rcl<-c(0,cumsum(groups))
  
  #Groups as list
  RTlist<-lapply(c(1:(a*b)),function(x){     
    RT[,(rcl[x]+1):rcl[x+1]]                 
  })
  
  #cell (group) means 
  meanRijList<-lapply(RTlist,rowMeans)						
  meanRij<-do.call(cbind,meanRijList)						
  
  #pooled sample mean w.r.t. factor b (a columns)
  Rtilde_i<-t(matrix(colMeans(matrix(t(meanRij),b)),a))	
  #pooled sample mean w.r.t. factor a (b columns)
  Rtilde_j<-matrix(rowMeans(matrix(meanRij,b*p)),p)			
  
  #(total) pooled sample mean
  Rtilde<-rowMeans(meanRij)							





  ####	Dispersion Matrix	####

  #cell wise (group-wise) deviations from mean
  Deviations<-mapply("-",RTlist,meanRijList,SIMPLIFY = FALSE)		
  #squared deviations
  Cproducts<-lapply(Deviations, tcrossprod)					

  nList<-as.list(1/(groups-1))
  
  # factors for sum in G
  nList2<-as.list(1/(groups))       


  Sij<-mapply("*",Cproducts,nList,SIMPLIFY=FALSE)				

  Sij.by.nij<-(mapply("*",Sij,nList2,SIMPLIFY=FALSE))             

  G<-(Reduce("+",Sij.by.nij))/(a*b)	# G

  VhatN<-N*Reduce("+",mapply(kronecker,lapply(c(1:(a*b)),function(i){     
    diag(diag(1,(a*b))[,i])                                               
  }),Sij.by.nij,SIMPLIFY = FALSE))


  ####  H  ####

  HA<-b*(tcrossprod((Rtilde_i)-(Rtilde)))/(a-1)

  HB<-a*(tcrossprod((Rtilde_j)-(Rtilde)))/(b-1)

  HAB<-(1/((a-1)*(b-1)))*(tcrossprod(
    meanRij
    -(matrix(rep(Rtilde_j,a),p))
    -(t(matrix(rep(t(Rtilde_i),each=b),a*b)))
    +(matrix(rep(Rtilde,a*b),p))
  ))

  H<-list(HA,HB,HAB)


  #### estimated d.o.f. ####

  fhat.A <- ((a-1)^2*(matrix.trace(VhatN))^2)/((a*b)^2*matrix.trace(T.A%*%VhatN%*%T.A%*%VhatN))

  fhat.B <- ((b-1)^2*(matrix.trace(VhatN))^2)/((a*b)^2*matrix.trace(T.B%*%VhatN%*%T.B%*%VhatN))

  fhat.AB <- ((a-1)^2*(b-1)^2*(matrix.trace(VhatN))^2)/((a*b)^2*matrix.trace(T.AB%*%VhatN%*%T.AB%*%VhatN))

  fhat <- c(fhat.A,fhat.B,fhat.AB)

  fhat.0 <- (matrix.trace(VhatN))^2/(((sapply(Sij.by.nij,matrix.trace))^2)%*%unlist(nList))


  ####  Test Statistics ####

  TraceG <- matrix.trace(G)
  m <- c((a-1),(b-1),(a-1)*(b-1))

  Results<-lapply(c(MainEffectA=1,MainEffectB=2,Interaction=3),function(i){

    D.Anova <- matrix.trace(H[[i]])/TraceG
    LR <- N*log(det(diag(p)+(m[i]*H[[i]]%*%(ginv(G*N)))))
    LH <- matrix.trace(m[i]*H[[i]]%*%(ginv(G)))
    BNP <- N*matrix.trace( m[i]*H[[i]]%*%(ginv(G*N)) %*% ginv(diag(p)+(m[i]*H[[i]]%*%(ginv(G*N)))))

    test.values<-c(D.Anova, LR, LH, BNP)

      p.values <- 1-c(pf(D.Anova,fhat[i],fhat.0),pchisq(LR,m[i]*p),pchisq(LH,m[i]*p),pchisq(BNP,m[i]*p))
      df <- c(fhat[i],m[i]*p,m[i]*p,m[i]*p)
      Rnames <- c("Dempster's ANOVA (F-approximation with df2='inf')",
                  "Wilks Lambda (Chi^2 approximation)",
                  "Lawley-Hotelling-Criterion (Chi^2 approximation)",
                  "Bartlett-Nanda-Pillai Criterion (Chi^2 approximation)")
      
    data.frame("Test Statistic"=test.values,
               "p value"=p.values,"df1"=df,
               row.names = Rnames
               )

  })

  Results$factorlevels<-data.frame("levels"=c(a,b),row.names = c(NameFactorA,NameFactorB))
  Results$groupsize<-t(struc)

  message(paste0(NameFactorA," is considered as factor 'A'"))
  message(paste0(NameFactorB," is considered as factor 'B'"))
  Results


}
