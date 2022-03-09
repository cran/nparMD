#' @title Nonparametric Test For Multivariate Data With Two-Way Layout Factorial Design - Small Samples
#'
#' @description Analysis of multivariate data with two-way completely randomized
#' factorial design - version for small samples. The analysis is based on fully nonparametric, rank-based methods and
#' uses a N(0,1)-approximation for test statistics based on 'Dempster's ANOVA', 'Wilk's Lambda', 'Lawley-Hotelling' and
#' 'Bartlett-Nanda-Pillai' criteria. This approximation is established by the asymptotic distribution
#' of these four statistics under true null-hypothesis if one of the explanatory factors has a large number
#' of levels. The multivariate response is allowed to be ordinal, quantitative, binary or a mixture
#' of the different variable types. The test statistics are constructed using nonparametric relative effect estimators.
#' @usage nparms(formula,data)
#' @param formula an object of class "formula" with two explanatory variables (factors), see examples.
#' @param data an object of class "data.frame" containing the variables in the formula
#' @return Returns a list of data frames providing the values of the test statistics, p-values, degrees of freedom, factor levels,
#' and groupsize per factor level combination.
#' @details This method is only implemented for complete data sets without missing values. The data is analysed for main effects
#' and interaction effect of the explanatory factors. In each case the null hypothesis "no effect" is testet. The explanatory factor 
#' that has the higher number of levels is automatically designated as factor "A". The covariance matrix
#' estimation requires at least 4 observations (observation vectors) per factor level combination. As the estimation is very time-consuming
#' for large groups it is performed wih a random selection of observations when a group exceeds a size of 6 observation vectors. 
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
#' data(pseudostudy2)
#' nparms(resp1|resp2|resp3~treatment*age, pseudostudy2)
#' @import Formula
#' @import matrixcalc
#' @import matrixStats
#' @import stats
#' @importFrom MASS ginv
#' @importFrom gtools permutations
#' @importFrom methods is
#'
#' @export
nparms<-function(formula, data){

  if(!is.data.frame(data)){stop("data argument is not a data frame")}
  if(!is(formula, "formula")){stop("Please give a formula")}

  formel<-Formula(formula)
  data<-model.frame(Formula(formula),data)

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

  RT<-(rowRanks(data,ties.method="average")-0.5)/N		#RankTransform matrix pxN


  #### Data grouping and means	####

  rcl<-c(0,cumsum(groups))                     #replication-group column-index list

  RTlist<-lapply(c(1:(a*b)),function(x){     #List of a*b matrices while
    #each list-element represents
    RT[,(rcl[x]+1):rcl[x+1]]                 #one replication group

  })

  meanRijList<-lapply(RTlist,rowMeans)						#row means R(bar)_ij
  meanRij<-do.call(cbind,meanRijList)						#matrix of row means

  Rtilde_i<-t(matrix(colMeans(matrix(t(meanRij),b)),a))			#sample mean over all b (a columns)
  Rtilde_j<-matrix(rowMeans(matrix(meanRij,b*p)),p)			#sample mean over all a (b columns)

  Rtilde<-rowMeans(meanRij)							#sample mean over a and b


  # Contrast Matrices

  P.b <- diag(1,b)-(rep(1,b) %*% t(rep(1,b)))/b
  C.B <- kronecker((rep(1,a) %*% t(rep(1,a)))/a, P.b)
  T.B<-kronecker( C.B, diag(p))


  ####	Dispersion Matrix	####

  Deviations<-mapply("-",RTlist,meanRijList,SIMPLIFY = FALSE)					# (R_ijk - Rbar_ij)
  Cproducts<-lapply(Deviations, tcrossprod)					# (R_ijk - Rbar_ij)(R_ijk - Rbar_ij)'

  nList<-as.list(1/(groups-1))			#
  nList2<-as.list(1/(groups))       # factors for sum in G


  Sij<-mapply("*",Cproducts,nList,SIMPLIFY=FALSE)				# S_ij

  Sij.by.nij<-(mapply("*",Sij,nList2,SIMPLIFY=FALSE))             # Sij by nij

  G<-(Reduce("+",Sij.by.nij))/(a*b)	# G

  VhatN<-N*Reduce("+",mapply(kronecker,lapply(c(1:(a*b)),function(i){     #Matrix V^hat_N
    diag(diag(1,(a*b))[,i])                                               #direct sum via kronecker product
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

    H<-list(HA,HAB,HB)


  #### estimated d.o.f. (only for ANOVA type) ####

  fhat.0<-(matrix.trace(VhatN))^2/(((sapply(Sij.by.nij,matrix.trace))^2)%*%unlist(nList))

  fhat.1<-((b-1)^2*(matrix.trace(VhatN))^2)/((a*b)^2*matrix.trace(T.B%*%VhatN%*%T.B%*%VhatN))


  ####  OMEGA ####

  Omega<-diag(p)*1/matrix.trace(G)
  Omega2<-ginv(G)

  ####	Quadruples and PSIijOMEGA	####

  Quadruples<-lapply(c(1:(a*b)),function(x){				#Set of all quadruples

    #Randomization here - in case of large groups
    #check if groupsize is greater than 6 - then take random selection of
    #all quadruples
  if(groups[x]>6){
    GS<-groups[x]
    AllQuads<-permutations(GS,4)
    QuadSample<-sample(1:(GS*(GS-1)*(GS-2)*(GS-3)),360)
    AllQuads[QuadSample,]
  }
    else {permutations(groups[x],4)}

    #List structure according to groups

  })



  PsiOMEGAij<-lapply(c(1:(a*b)),function(rGroup){

    Qparts<-lapply(c(1:(nrow(Quadruples[[rGroup]]))),function(x){	#Summanden eines PSIij(Omega)

      k1<-Quadruples[[rGroup]][x,1]
      k2<-Quadruples[[rGroup]][x,2]
      k3<-Quadruples[[rGroup]][x,3]
      k4<-Quadruples[[rGroup]][x,4]

      Omega%*%(RTlist[[rGroup]][,k1]-RTlist[[rGroup]][,k2])%*%t(RTlist[[rGroup]][,k1]-RTlist[[rGroup]][,k2]) %*% Omega%*%(RTlist[[rGroup]][,k3]-RTlist[[rGroup]][,k4])%*%t(RTlist[[rGroup]][,k3]-RTlist[[rGroup]][,k4])

    })

    if(groups[rGroup]>6){
    (Reduce("+",Qparts))*(1/(4*360))
    }
    else{
      (Reduce("+",Qparts))*(1/(4*(groups[rGroup])*(groups[rGroup]-1)*(groups[rGroup]-2)*(groups[rGroup]-3)))
    }

  })



  PsiOMEGA2ij<-lapply(c(1:(a*b)),function(rGroup){

    Qparts<-lapply(c(1:(nrow(Quadruples[[rGroup]]))),function(x){	#Summanden eines PSIij(Omega)

      k1<-Quadruples[[rGroup]][x,1]
      k2<-Quadruples[[rGroup]][x,2]
      k3<-Quadruples[[rGroup]][x,3]
      k4<-Quadruples[[rGroup]][x,4]

      Omega2%*%(RTlist[[rGroup]][,k1]-RTlist[[rGroup]][,k2])%*%t(RTlist[[rGroup]][,k1]-RTlist[[rGroup]][,k2]) %*% Omega2%*%(RTlist[[rGroup]][,k3]-RTlist[[rGroup]][,k4])%*%t(RTlist[[rGroup]][,k3]-RTlist[[rGroup]][,k4])

    })

    if(groups[rGroup]>6){
      (Reduce("+",Qparts))*(1/(4*360))
    }
    else{
      (Reduce("+",Qparts))*(1/(4*(groups[rGroup])*(groups[rGroup]-1)*(groups[rGroup]-2)*(groups[rGroup]-3)))
    }

  })

  ####  nyOMEGA ####


  ny1OMEGA<-sum(sapply(PsiOMEGAij,matrix.trace)/(groups*(groups-1)))/(a*b)

  ny1OMEGA2<-sum(sapply(PsiOMEGA2ij,matrix.trace)/(groups*(groups-1)))/(a*b)


  Indexj<-permutations(b,2)

  ny2OMEGA<-(1/(a*b))*sum(

    sapply(c(1:a),function(s){

      sum(

        sapply(c(1:nrow(Indexj)),function(x){

          (matrix.trace(Omega %*% Sij[[s-1+Indexj[x,1]]] %*% Omega %*% Sij[[s-1+Indexj[x,2]]]))/
            (groups[s-1+Indexj[x,1]]*groups[s-1+Indexj[x,2]])

        })

      )

    })

  )

  ny2OMEGA2<-(1/(a*b))*sum(

    sapply(c(1:a),function(s){

      sum(

        sapply(c(1:nrow(Indexj)),function(x){

          (matrix.trace(Omega2 %*% Sij[[s-1+Indexj[x,1]]] %*% Omega2 %*% Sij[[s-1+Indexj[x,2]]]))/
            (groups[s-1+Indexj[x,1]]*groups[s-1+Indexj[x,2]])

        })

      )

    })

  )


  ####  tau ####

  tauA<-sqrt((2/b)*(ny1OMEGA+ny2OMEGA))
  tauAB<-sqrt((2/b)*(ny1OMEGA+ny2OMEGA/((b-1)^2)))
  tauA2<-sqrt((2/b)*(ny1OMEGA2+ny2OMEGA2))
  tauAB2<-sqrt((2/b)*(ny1OMEGA2+ny2OMEGA2/((b-1)^2)))

  tau<-list(c(tauA,tauA2),c(tauAB,tauAB2))


  ####  Test Statistics ####

  TraceG<-matrix.trace(G)
  m <- c((a-1),(a-1)*(b-1),(b-1))


  Results<-lapply(c(MainEffectA=1, Interaction=2, MainEffectB=3),function(i){



    if(i == 3){

      D.Anova <- matrix.trace(H[[i]])/TraceG
      LR <- N*log(det(diag(p)+(m[i]*H[[i]]%*%(ginv(G*N)))))
      LH <- matrix.trace(m[i]*H[[i]]%*%(ginv(G)))
      BNP <- N*matrix.trace( m[i]*H[[i]]%*%(ginv(G*N)) %*% ginv(diag(p)+(m[i]*H[[i]]%*%(ginv(G*N)))))

      test.values<-c(D.Anova, LR, LH, BNP)

      p.values <- 1-c(pf(D.Anova,fhat.1,fhat.0),pchisq(LR,m[i]*p),pchisq(LH,m[i]*p),pchisq(BNP,m[i]*p))
      df <- c(fhat.1,m[i]*p,m[i]*p,m[i]*p)
      Rnames <- c("Dempster's ANOVA (F-approximation with df2='inf')",
                  "Wilks Lambda (Chi^2 approximation)",
                  "Lawley-Hotelling-Criterion (Chi^2 approximation)",
                  "Bartlett-Nanda-Pillai Criterion (Chi^2 approximation)")
      
      data.frame("Test Statistic"=test.values,
                 "p value"=p.values,"df1"=df,
                 row.names = Rnames)

    } else {

        D.Anova <- matrix.trace(H[[i]])/TraceG
        LR <- log(det(diag(p)+(H[[i]]%*%(ginv(G)))))
        LH <- matrix.trace(H[[i]]%*%(ginv(G)))
        BNP <- matrix.trace(H[[i]]%*%(ginv(G))%*%ginv(diag(p)+(H[[i]]%*%(ginv(G)))))

        test.values<-c(sqrt(a)*(D.Anova-1)/tau[[i]][1],
                       sqrt(a)*(2*LR-2*p*log(2))/tau[[i]][2],
                       sqrt(a)*(1*LH-p)/tau[[i]][2],
                       sqrt(a)*(4*BNP-2*p)/tau[[i]][2])
        p.values<-pnorm(-test.values)

        Rnames <- c("Dempster's ANOVA",
                    "Wilks Lambda",
                    "Lawley-Hotelling-Criterion",
                    "Bartlett-Nanda-Pillai Criterion")

        data.frame("Test Statistic"=test.values,
                   "p value"=p.values,
                   row.names = Rnames
        )
    }



    

  })






  #names(TestValues)[[1]]<-paste0("main effect of ",NameFactorA)
  Results$factorlevels<-data.frame("levels"=c(a,b),row.names = c(NameFactorA,NameFactorB))
  Results$groupsize<-t(struc)

  # names(TestValues)<-c(paste("Main Effect of", NameFactorA),
  #                      paste("Interaction Effect of",NameFactorA,NameFactorB),
  #                      "Group size per factor level combination, col(row) = factor level of A(B)")


  message(paste0(NameFactorA," is considered as factor 'A'"))
  message(paste0(NameFactorB," is considered as factor 'B'"))
  Results



}
