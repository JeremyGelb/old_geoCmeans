rm(list = ls())
library(ClustGeo)
library(rgeos)
library(sf)
library(ggplot2)
library(sp)
library(tidyverse)
library(spdep)
library(reldist)
library(SDMTools)
data(estuary)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de clalculs algorithmes c-mean#####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#fonction pour calculer la distance euclidienne entre les observations d'un dataframe et un point

#' Calculate the euclidean distance between a numeric matrix n X p and a numeric vector of length p
#'
#' @param m A n X p matrix or dataframe with only numeric columns
#' @param v A numerci vector of length p
#' @return A vector of length n giving the euclidean distance between all matrix row and the vector p
#' @examples
#' mat <- matrix(c(1,4,2,5,3,6),nrow=2,ncol=3)
#' v1 <- c(1,2,3)
#' calcEuclideanDistance(mat,v1)
calcEuclideanDistance <- function(m,v){
  mat1 <- as.matrix(m)
  mat2 <- matrix(as.numeric(v),nrow=nrow(m),ncol=length(v),byrow=TRUE)
  alldistances <- rowSums((mat1-mat2)**2)
  return(alldistances)
}

#' Intermediar step in the calculation of the belonging matric
#'
#' @param centerid An integer representing the column to use as reference in the distance matrix
#' @param alldistances A distance matrix
#' @return m a float representing the fuzzyness degree
#' @examples
#' #no example provided, this is an internal function
#
calcBelongDenom__ <- function(centerid,alldistances,m){
  values <- sapply(1:ncol(alldistances),function(i){
    div <- (alldistances[,centerid] / alldistances[,i])
    return (div**(2/(m-1)))
  })
  return (rowSums(values))
}


#fonction pour calculer la matrice d'appartenance

#' Calculate the belonging matrix according to a set of centroids, the observed data and the fuzzyness degree
#'
#' @param centers A matrix or a dataframe representing the centers of the clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with nrows and p columns
#' @param m An integer representing the fuzzyness degree
#' @return A n X k matrix represening the belongings of each datapoint to each cluster
#' @examples
#' data <- matrix(c(1,4,2,5,3,6,7,8,9,10,11,12),nrow=6,ncol=2)
#' centers <- data[1:2,]
#' calcBelongMatrix(centers,data,1.5)
calcBelongMatrix <- function(centers,data,m){
  centerdistances <- as.data.frame(apply(centers,1,function(x){return(calcEuclideanDistance(data,x))}))
  belongingmatrix <- sapply(1:nrow(centers),function(x){
    return (1/calcBelongDenom__(x,centerdistances,m))
  })
  belongingmatrix[is.na(belongingmatrix)] <- 1
  return(belongingmatrix)
}


#
# calcADWCfNum__ <- function(Data,Center,neighbours){
#   AllDistances <- calcEuclideanDistance(Data,Center)
#   Dists <- apply(neighbours,1,function(x){
#     Dist <- sum(AllDistances[as.logical(x)])
#   })
#   return(Dists)
# }


#fonction pour calculer la matrice d'appartenance avec la ponderation spatiale
#neighbours : une matrice de 0/1 represantant le voisinage

#' Calculate the belonging matrix (spatial version) according to a set of centroids, the observed data, the fuzzyness degree a neighbouring matrix and a spatial ponderation term
#'
#' @param centers A matrix or a dataframe representing the centers of the clusters with p columns and k rows
#' @param data A dataframe or matrix representing the observed data with nrows and p columns
#' @param neighbours A list.w object describing the neighbours typically produced by the spdep package
#' @param m An integer representing the fuzzyness degree
#' @param alpha An integer representing the weight of the space in the analysis (0 is a typical fuzzy-c-mean algorithm)
#' @return A n X k matrix represening the belongings of each datapoint to each cluster
#' @examples
#' data <- matrix(c(1,4,2,5,3,6,7,8,9,10,11,12),nrow=6,ncol=2)
#' centers <- data[1:2,]
#' calcBelongMatrix(centers,data,1.5)

calcSWCMBelongMatrix <- function(centers,data,neighbours,m,alpha){
  #calcul des distances originales
  originaldistances <- apply(centers,1,function(x){return(calcEuclideanDistance(data,x))})
  #calcul des distances lagguees
  nbdistances <- apply(originaldistances,2,function(x){
    return(lag.listw(neighbours,x))
  })
  #calcul de la amtrice d'appartenance

  denom <- rowSums((OrignalDistances + alpha*nbdistances)**(-1/(m-1)))
  belongingmatrix <- sapply(1:ncol(originaldistances),function(x){
    odist <- originaldistances[,x]
    ndist <- nbdistances[,x]
    u <- (odist+alpha*odist)**(-1/(m-1))
    u <- u/denom
  })
  belongingmatrix[is.na(belongingmatrix)] <- 1
  return(belongingmatrix)
}



#fonction pour calculer les nouveaux centroids

#' Calculate the new centroids of the clusters based on the belonging matrix
#'
#' @param data A dataframe or matrix representing the observed data with nrows and p columns
#' @param belongmatrix A n X k matrix giving for each observation n, its probability to belong to the cluster k
#' @param m An integer representing the fuzzyness degree
#' @return A n X k matrix represening the belongings of each datapoint to each cluster
#' @examples
#' data <- matrix(c(1,4,2,5,3,6,7,8,9,10,11,12),nrow=6,ncol=2)
#' centers <- data[1:2,]
#' belongmatrix <- calcBelongMatrix(centers,data,1.5)
#' newcenters <- calcCentroids(data,belongmatrix,1.5)
calcCentroids<- function(data,belongmatrix,m){
  centers <- sapply(1:ncol(belongmatrix),function(i){
    apply(data,2,function(x){
      weighted.mean(x, belongmatrix[,i])
    })
  })
  centers <- as.data.frame(t(Centers))
}

#fonction pour calculer les nouveaux centroids (version spatialle)

#' Calculate the new centroids of the clusters based on the belonging matrix (spatial version)
#'
#' @param data A dataframe or matrix representing the observed data with nrows and p columns
#' @param belongmatrix A n X k matrix giving for each observation n, its probability to belong to the cluster k
#' @param neighbours A list.w object describing the neighbours typically produced by the spdep package
#' @param m An integer representing the fuzzyness degree
#' @param alpha An integer representing the weight of the space in the analysis (0 is a typical fuzzy-c-mean algorithm)
#' @return A n X k matrix represening the belongings of each datapoint to each cluster
#' @examples
#' data <- matrix(c(1,4,2,5,3,6,7,8,9,10,11,12),nrow=6,ncol=2)
#' centers <- data[1:2,]
#' belongmatrix <- calcBelongMatrix(centers,data,1.5)
#' newcenters <- calcCentroids(data,belongmatrix,1.5)
calcSWFCCentroids<- function(data,belongmatrix,neighbours,m,alpha){
 centers <- sapply(1:ncol(belongmatrix),function(i){
   weights <- belongmatrix[,i]**m
   apply(data,2,function(x){
     v <- (x + alpha*lag.listw(neighbours,x))
     return (sum(v*weights)/((1+alpha)*sum(weights)))
   })
 })
 centers <- as.data.frame(t(centers))
}



#fonction pour evaluer la difference entre deux matrices d'appartenance (ici : la moyenne d'appartenance de chaque observation abougee de moins de tol)

#' evaluate if the algorithm converged by comparing two successive belongings matrices. Calculate the absolute difference between the matrices and then calculate the mean of each row. If all the values of the final vector are below the fixed tolerance, then return True, else return False
#'
#' @param mat1 A n X k matrix giving for each observation n, its probability to belong to the cluster k at iteration i
#' @param mat2 A n X k matrix giving for each observation n, its probability to belong to the cluster k at iteration i+1
#' @param tol a float representing the algorithm tolerance
#' @examples
#' mat1 <- rbind(c(0.45,0.55),c(0.35,0.65),c(0.8,0.2),c(0.5,0.5),c(0.9,0.1),c(0.7,0.3))
#' mat2 <- rbind(c(0.45,0.55),c(0.4,0.60),c(0.8,0.2),c(0.45,0.55),c(0.95,0.05),c(0.7,0.3))
#' evaluateMatrices(mat1,mat2,0.01)
evaluateMatrices <- function(mat1,mat2,tol){
  mat1 < -as.matrix(mat1)
  mat2 <- as.matrix(mat2)
  differ <- abs(mat1-mat2)
  diffobs <- rowSums(differ) / ncol(mat1)
  if (length(diffobs[diffobs>=tol])>0){
    return (FALSE)
  }else(
    return(TRUE)
  )
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### algorithme c-mean#####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Fonction appliquant l'algo
#' calssical c-mean algorithm
#'
#' @param data A dataframe with only numerical variable
#' @param k An integer describing the number of cluster to find
#' @param m An integer for the fuzzyness degree
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for convergence assessment
#' @param standardize A boolean to specify if the variable must be centered and reduce (default = True)
#' @examples
#' mat1 <- rbind(c(0.45,0.55),c(0.35,0.65),c(0.8,0.2),c(0.5,0.5),c(0.9,0.1),c(0.7,0.3))
#' mat2 <- rbind(c(0.45,0.55),c(0.4,0.60),c(0.8,0.2),c(0.45,0.55),c(0.95,0.05),c(0.7,0.3))
#' evaluateMatrices(mat1,mat2,0.01)
CFuzzyMeans <- function(data,k,m,maxiter=500,tol=0.01,standardize=TRUE){

  if(standardize){
    for(i in 1:ncol(data)){
      data[,i] <- scale(data[,i])
    }
  }

  #selection aleatoires des centres originaux
  centers <- data[sample(nrow(data), k), ]

  #calcul de la matrice d'appartenance
  belongmatrix <- calcBelongMatrix(centers,data,m)

  CriterioReached <- FALSE
  #lancement de la boucle
  pb <- txtProgressBar(1, maxiter, style=3)
  for(i in 1:maxiter){
    setTxtProgressBar(pb, i)
    newcenters <- calcCentroids(data,belongmatrix,m)
    newbelongmatrix <- calcBelongMatrix(newcenters,data,m)
    if (evaluateMatrices(belongmatrix,newbelongmatrix,tol)==FALSE){
      #si on n'atteint pas encore la convergence
      centers <- newcenters
      belongmatrix <- newbelongmatrix
    }else{
      #si on atteint la convergence
      print("criterion reached")
      CriterioReached <- TRUE
      centers <- newcenters
      belongmatrix <- newbelongmatrix
      break
    }
  }
  #calcul de l'appartenance la plus vraissemblable
  DF <- as.data.frame(NewBelongMatrix)
  groups <- colnames(DF)[max.col(DF,ties.method="first")]

  return(list("Centers"=centers,"Belongings"=belongmatrix,"Groups"=groups,"Data"=data))

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### algorithme c-mean spatial #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#SFCM : c mean with neighbourhood constraint

#' spatial version of the c-mean algorithm
#'
#' @param data A dataframe with only numerical variable
#' @param nblistw A list.w object describing the neighbours typically produced by the spdep package
#' @param k An integer describing the number of cluster to find
#' @param m An integer for the fuzzyness degree
#' @param maxiter An integer for the maximum number of iteration
#' @param tol The tolerance criterion used in the evaluateMatrices function for convergence assessment
#' @param standardize A boolean to specify if the variable must be centered and reduce (default = True)
#' @examples
#' mat1 <- rbind(c(0.45,0.55),c(0.35,0.65),c(0.8,0.2),c(0.5,0.5),c(0.9,0.1),c(0.7,0.3))
#' mat2 <- rbind(c(0.45,0.55),c(0.4,0.60),c(0.8,0.2),c(0.45,0.55),c(0.95,0.05),c(0.7,0.3))
#' evaluateMatrices(mat1,mat2,0.01)
SpatialCFuzzyMeans <- function(data,nblistw,k,m,alpha,maxiter=500,tol=0.01,standardize=TRUE){

  for(i in 1:ncol(data)){
    data[,i] <- scale(data[,i])
  }

  #selection aleatoires des centres originaux
  centers <- data[sample(nrow(data), k), ]

  #calcul de la matrice d'appartenance
  belongmatrix <- calcSWCMBelongMatrix(centers,data,nblistw,m,alpha = alpha)


  CriterioReached <- FALSE
  #lancement de la boucle
  pb <- txtProgressBar(1, maxiter, style=3)
  for(i in 1:maxiter){
    setTxtProgressBar(pb, i)
    newcenters <- calcSWFCCentroids(data,belongmatrix,nblistw,m,alpha)
    newbelongmatrix <- calcSWCMBelongMatrix(newcenters,data,nblistw,m,alpha = alpha)
    if (evaluateMatrices(belongmatrix,newbelongmatrix,tol)==FALSE){
      #si on n'atteint pas encore la convergence
      centers <- newcenters
      belongmatrix <- newbelongmatrix
    }else{
      #si on atteint la convergence
      print("criterion reached")
      CriterioReached <- TRUE
      centers <- newcenters
      belongmatrix <- newbelongmatrix
      break
    }
  }
  #calcul de l'appartenance la plus vraissemblable
  DF <- as.data.frame(newbelongmatrix)
  Groups <- colnames(DF)[max.col(DF,ties.method="first")]

  return(list("Centers"=centers,"Belongings"=belongmatrix,"Groups"=Groups,"Data"=data))

}








#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Tests #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##### classical fuzzy c-means


# #Fonction appliquant l'algo
#
#
# k<-5
# m<-1.5
# alpha <- 0
# Data <- estuary$dat
# SpData <- estuary$map
#
#
# Neighbours <- poly2nb(SpData,queen = TRUE)
# WMat <- nb2listw(Neighbours,style="W",zero.policy = TRUE)
#
#
#
# Values <- sapply(2:10,function(k){
#   Result2 <- SpatialCFuzzyMeans(Data,WMat,k,m,alpha=0,tol=0.0001)
#   idxsf <- SIL.F(Data,Result2$Belongings)#look for maximum
#   idxpe <- PE(Result2$Belongings) #look for minimum
#   idxpc <- PC(Result2$Belongings) #loook for maximum
#   idxmpc <- MPC(Result2$Belongings) #look for maximum
#   return(c(idxsf,idxpe,idxpc,idxmpc,k))
# })
#
# Values <- as.data.frame(t(Values))
# names(Values) <- c("idxsf","idxpe","idxpc","idxmpc","k")
#
# ggplot(data=Values)+
#   geom_point(aes(x=k,y=idxsf))
#
# ggplot(data=Values)+
#   geom_point(aes(x=k,y=idxpe))
#
# ggplot(data=Values)+
#   geom_point(aes(x=k,y=idxpc))
#
# ggplot(data=Values)+
#   geom_point(aes(x=k,y=idxmpc))
#
#
# k<-5
# for(alpha in seq(0,1.2,0.1)){
#   Result2 <- SpatialCFuzzyMeans(Data,WMat,k,m,alpha=alpha,tol=0.0001)
#   Colors <- recode(Result2$Groups,
#                    "V1"="green",
#                    "V2"="red",
#                    "V3"="blue",
#                    "V4"="orange",
#                    "V5"="black")
#   sp::plot(estuary$map,col=Colors,main=paste("5 cluster with alpha=",alpha,sep=""),cex.main=0.8)
#   legend("topleft", legend=paste("cluster",1:5), fill=c("green","red","blue","orange","black"), cex=0.8,bty="n",border="white")
# }
#
# ######test donnees de lyon !
#
# library(rgdal)
# library(rgeos)
#
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#
#
# LyonIris <- readOGR("IRIS_GrandLyonPollutionINSEE_ssHydro_veg2.gpkg")
# ToNumeric <- c("POPTOT","Pop0_14","Pop65plus","PopImg","Chom15_64","Brevet15P","X","Y")
# for (Var in ToNumeric){
#   LyonIris@data[[Var]] <- as.numeric(as.character(LyonIris@data[[Var]]))
# }
#
# Neighbours <- poly2nb(LyonIris,queen = TRUE)
# WMat <- nb2listw(Neighbours,style="W",zero.policy = TRUE)
#
# AnalysisFields <-c("Lden","NO2","PM25","VegHautPrt","Pct0_14","Pct_65","Pct_Img","TxChom1564","Pct_brevet","NivVieMed")
# LyonIris$OID <- as.character(1:nrow(LyonIris))
# FortiData <- broom::tidy(LyonIris,region="OID")
#
# Data <- LyonIris@data[AnalysisFields]
# k <- 4
# m<-1.5
#
# for(alpha in seq(0,2,0.1)){
#   Result2 <- SpatialCFuzzyMeans(Data,WMat,k,m,alpha=alpha,tol=0.0001)
#   Colors <- recode(Result2$Groups,
#                    "V1"="green",
#                    "V2"="red",
#                    "V3"="blue",
#                    "V4"="orange",
#                    "V5"="black")
#   sp::plot(LyonIris,col=Colors,main=paste("4 cluster with alpha=",alpha,sep=""),cex.main=0.8)
#   legend("topleft", legend=paste("cluster",1:4), fill=c("green","red","blue","orange","black"), cex=0.8,bty="n",border="white")
# }
#
# OkAlpha <- 0.5
# Result2 <- SpatialCFuzzyMeans(Data,WMat,k,m,alpha=OkAlpha,tol=0.0001)
#
# Result2$Belongings$OID <- LyonIris$OID
#
# FortiData <- merge(FortiData,Result2$Belongings,by.x="id",by.y="OID")
#
# ggplot(FortiData)+
#   geom_polygon(aes(group=group,fill=V1,x=long,y=lat),color="white")
#
# ggplot(FortiData)+
#   geom_polygon(aes(group=group,fill=V2,x=long,y=lat),color="white")
#
# ggplot(FortiData)+
#   geom_polygon(aes(group=group,fill=V3,x=long,y=lat),color="white")
