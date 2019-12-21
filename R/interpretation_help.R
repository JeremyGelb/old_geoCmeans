library(fclust)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de diagnostic #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#calcul des indices de diagnostic sur les resultats de la classification

#' calculate several clustering quality indexes (most of them come from fclust package)
#'
#' @param data the original dataframe used fot the classification
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param centers The dataframe representing the centers of the clusters
#' @examples
#'
calcqualityIndexes <- function(data,belongmatrix,centers){
  idxsf <- SIL.F(data,belongmatrix)#look for maximum
  idxpe <- PE(belongmatrix) #look for minimum
  idxpc <- PC(belongmatrix) #loook for maximum
  idxmpc <- MPC(belongmatrix) #look for maximum
  #calcul de l'indice d'intertie
  means <- apply(data,2,mean)
  baseinertia <- sum(calcEuclideanDistance(data,means))
  restinertia <- sapply(1:nrow(centers),function(i){
    point <- centers[i,]
    dists <- calcEuclideanDistance(data,point)*belongmatrix[,i]
    return (sum(dists))
  })
  explainedinertia <- 1-(sum(restinertia)/baseinertia)

  return(list("Silhouette.index"=idxsf,
              "Partition.entropy"=idxpe,
              "Partition.coeff"=idxpc,
              "Modified.partition.coeff"=idxmpc,
              "Explained.inertia"=explainedinertia))
}

#calcul des indicateurs spatiaux
#' calculate the moran I of each column of the belonging matrix and the multiple join count test of the most probablt cluster (spdep package)
#'
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param A list.w object describing the neighbours typically produced by the spdep package
#' @param undecided A numeric value giving the minimum value that an observation must get in the belonging matrix to not be considered as uncertain (default = NULL)
#' @examples
#'
spatialdiag <- function(belongmatrix,nblistw,undecided=NULL){
  belongmatrix <- as.data.frame(belongmatrix)
  # calcul des I de Moran pour les appartenances
  Values <- sapply(1:ncol(belongmatrix),function(i){
    x <- belongmatrix[,i]
    morantest <- moran.mc(x,nblistw,nsim=999)
    return(list("MoranI"=morantest$statistic,
                "pvalue"=morantest$p.value,
                "Cluster"=paste("Cluster_",i,sep="")))
  })
  morandf <- as.data.frame(t(Values))
  #attribution des groupes
  groups <- colnames(belongmatrix)[max.col(belongmatrix,ties.method="first")]
  if (is.null(undecided)==FALSE){
    Maximums  <- do.call(pmax, belongmatrix)
    groups[Maximums<undecided]<-"undecided"
  }
  groups <- as.factor(groups)
  #calcul des join count test
  spjctetst <- joincount.multi(groups,nblistw,zero.policy = TRUE)
  return(list("MoranValues"=morandf,
              "JoinCounts"=spjctetst))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##### Fonctions de visualisation #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' build some maps to visualize the results of the clustering
#'
#' @param geodata A spatialpolygondataframe ordered like the original dataset used for the clustering
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param undecided A numeric value giving the minimum value that an observation must get in the belonging matrix to not be considered as uncertain (default = NULL)
#' @examples
#'
mapClusters <- function(geodata,belongmatrix,undecided=NULL){
  geodata@data <- cbind(geodata@data,belongmatrix)
  geodata$OID <- 1:nrow(geodata@data)

  #attribution des groupes
  belongmatrix <- as.data.frame(belongmatrix)
  Groups <- colnames(belongmatrix)[max.col(belongmatrix,ties.method="first")]
  if (is.null(undecided)==FALSE){
    Maximums  <- do.call(pmax, belongmatrix)
    Undecided <- ifelse(Maximums<undecided,"Undecided",paste("Max higher than ",undecided,sep=""))
  }else{
    Undecided <- rep("Ok",length(Groups))
  }
  geodata$Cluster <- Groups
  geodata$Undecided <- Undecided

  FortiData <- broom::tidy(geodata,region="OID")
  FortiData <- merge(FortiData,geodata@data,by.x="id",by.y="OID")
  #realisation des cartes de probabilites
  ProbaPlots <- lapply(names(belongmatrix),function(n){
    Plot <- ggplot(FortiData)+
      geom_polygon(aes_string(x="long",y="lat",group="group",fill=n),colour="black",size=0.1)+
      scale_fill_gradient(low="green", high="blue")+
      coord_fixed(ratio=1)
    return(Plot)
  })
  #realisation de la carte des groupes
  ClusterMap <- ggplot(FortiData)+
    geom_polygon(aes(x=long,y=lat,group=group,fill=Cluster),colour="black",size=0.1)+
    scale_fill_brewer(palette="Set1")+
    geom_polygon(aes(x=long,y=lat,group=group),fill=rgb(1,1,1,0.7),colour="black",size=0.1,data=subset(FortiData,FortiData$Undecided=="Undecided"))+
    coord_fixed(ratio=1)
  return(list("ProbaMaps"=ProbaPlots,"ClusterPlot"=ClusterMap))
}

#fonction pour calculer les statistiques descriptives des differents groupes
#' Calculate some descriptive statistics of each group
#'
#' @param data the original dataframe used fot the classification
#' @param belongmatrix The belonging matrix obtained at the end of the algorithm
#' @param weighted A boolean indicating if the summary statistics must use the belonging matrix columns as weights (TRUE) or simply assign each observation to its most likely cluster and compute the statistics on each subset (default=True)
#' @examples
#'
summarizeClusters<-function(data,belongmatrix,weighted=TRUE,dec=3){
  belongmatrix <- as.data.frame(belongmatrix)
  if(weighted){
    Summaries <- lapply(1:ncol(belongmatrix),function(c){
      W <- belongmatrix[,c]
      Values <- apply(Data,2,function(x){
        Q5 <- wtd.quantile (x, q=0.05, na.rm = TRUE, weight=W)
        Q10 <- wtd.quantile (x, q=0.1, na.rm = TRUE, weight=W)
        Q25 <- wtd.quantile (x, q=0.25, na.rm = TRUE, weight=W)
        Q50 <- wtd.quantile (x, q=0.5, na.rm = TRUE, weight=W)
        Q75 <- wtd.quantile (x, q=0.75, na.rm = TRUE, weight=W)
        Q90 <- wtd.quantile (x, q=0.90, na.rm = TRUE, weight=W)
        Q95 <- wtd.quantile (x, q=0.95, na.rm = TRUE, weight=W)
        Mean <- weighted.mean(x, W)
        Std <- wt.sd(x, W)
        return(list("Q5"=Q5,"Q10"=Q10,"Q25"=Q25,"Q50"=Q50,"Q75"=Q75,"Q90"=Q90,"Q95"=Q95,"Mean"=Mean,"Std"=Std))
      })
      DF <- do.call(cbind,Values)
      print(paste("Statistic summary for cluster ",c,sep=""))
      print(DF)
      return(DF)
    })
    names(Summaries)<-paste("Cluster_",c(1:ncol(belongmatrix)),sep="")
    return(Summaries)

  }else{
    Groups <- colnames(belongmatrix)[max.col(belongmatrix,ties.method="first")]
    Data$Groups <- Groups
    lapply(unique(data$Groups),function(c){
      DF <- subset(data,Data$Groups==c)
      DF$Groups <- NULL
      Values <- apply(DF,2,function(x){
        Q5 <- quantile (x, probs=0.05, na.rm = TRUE)
        Q10 <- quantile (x, probs=0.1, na.rm = TRUE)
        Q25 <- quantile (x, probs=0.25, na.rm = TRUE)
        Q50 <- quantile (x, probs=0.5, na.rm = TRUE)
        Q75 <- quantile (x, probs=0.75, na.rm = TRUE)
        Q90 <- quantile (x, probs=0.90, na.rm = TRUE)
        Q95 <- quantile (x, probs=0.95, na.rm = TRUE)
        Mean <- mean(x)
        Std <- sd(x)
        return(list("Q5"=Q5,"Q10"=Q10,"Q25"=Q25,"Q50"=Q50,"Q75"=Q75,"Q90"=Q90,"Q95"=Q95,"Mean"=Mean,"Std"=Std))
      })
      DF <- do.call(cbind,Values)
      print(paste("Statistic summary for cluster ",c,sep=""))
      print(DF)
      return(DF)
    })
    names(Summaries)<-paste("Cluster_",c(1:ncol(belongmatrix)),sep="")
    return(Summaries)
  }
}
