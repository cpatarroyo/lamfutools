#' Calculates a Phylogenetic Correlation From a Network Using PVR
#'
#' This function calculates the phylogenetic correlation between the genetic distance calculated from a minimal spanning network. This uses the average commute time to calculate the distance between the nodes of the netowrk and it takes phenotypic information from the user. Each individual must be unique (incompatible with various individuals with the same genotype).
#'
#' @param pop Population in a poppr format. Can be imported from a Genalex formatted file.
#' @param traits Vector of the phenotypic characteristics of the members of the populations.
#' @param dismet Metric for distance in the network. Defaults to normalized average commute time.
#' @param pvrmet Metric for the correlation coefficient. Defaults to Moran's I
#' @param gendist Metric for the genetic distance to elaborate the MSN. Defaults to Bruvo distance. Other options include: "nei" for Nei's distance, "bitwise" for Bitwise distance, "absolute" for absolute distance, "rogers" for Roger's distance, and "reynolds" for Reynold's distance.
#' @return Returns the eigenvalue correlation between the traits given and the distance between tne nodes in the network

nPVR<-function(pop,traits,dismet=NULL,pvrmet=NULL,gendist=NULL) {

  #Load required libaries
  require(PVR)
  require(poppr)
  require(linkprediction)

  #Assign default values for distance and PVR methods
  if(is.null(dismet)) {
    dismet<-"act_n"
  }
  if(is.null(pvrmet)) {
    pvrmet<-"moran"
  }

  #Build network
  if(is.null(gendist)) {
    network<-bruvo.msn(pop,vertex.label=NA)
  }
  else {
    if(gendist=="nei") {
      dista<-nei.dist(pop)
    }
    else if(gendist=="bitwise") {
      dista<-bitwise.dist(pop)
    }
    else if(gendist=="absolute") {
      dista<-diss.dist(pop)
    }
    else if(gendist=="rogers") {
      dista<-rogers.dist(pop)
    }
    else if(gendist=="reynolds") {
      dista<-reynolds.dist(pop)
    }
    network<-poppr.msn(pop,dista,vertex.label=NA)
  }

  #Calculate average commute time from the network
  distmat<-proxfun(network$graph,method=dismet)
  netPVR <- new("PVR")
  netPVR@Eigen <- eigen(distmat)
  netPVR@phyDist <- distmat

  #Add names to the PVR element
  numb<-dim(distmat)[1]
  cs<-rep("c",numb)
  cs<-paste(cs,1:numb,sep="")
  names(netPVR@Eigen$values)<-cs
  colnames(netPVR@Eigen$vectors)<-cs

  #Regression
  sizePVRn <- PVR(x=netPVR,phy=NULL,trait=traits,significance=TRUE)

  #Print results
  print(sizePVRn@Selection$Id)
  return(sizePVRn)

}

#' Converts the Results of an Evolutionary Simulation into a Genalex File.
#'
#' Converts the resulting population of an evolutionary simulation from the Simulation function of this same package.
#'
#' @param x Population resulting from the evolutionary simulation.
#' @param fname File name for the resultingg Genalex File

sim2gen <- function(x,fname) {
  #Taking the dimensions of the pop element
  msat <- dim(x@micro)[2]
  num <- length(x@location)
  #Compose the header of the csv file
  header <- c(paste(msat,num,1,paste(rep(",",msat-3),sep="",collapse = ""),sep =",",collapse = ""),paste(rep(",",msat+1),sep = "",collapse = ""))
  temptab <- as.data.frame(x@micro)
  temptab[,3:(msat+2)] <- temptab
  temptab[,1] <- paste("S",1:num,sep = "")
  temptab[,2] <- rep("P1",num)
  colnames(temptab) <- c("Ind","pop",paste("M",1:msat,sep = ""))
  #Write data into file
  write(header,file = fname)
  write.table(temptab, file = fname, sep = ",", quote = FALSE, append = TRUE, row.names = FALSE)
}
