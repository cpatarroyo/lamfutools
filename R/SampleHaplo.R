#' Create a Table With the Possible Haplotypes From a Table of Polyploid Genotypes
#'
#' This function loads a table with a set of polyploid markers and separates its possible haplotypes.
#'
#' @param table Table with the polyploid loci, one organisms per row
#' @param ploidy Ploidy level of the organism
#' @return Table of haploid data from the polyploid table

haplotypes<-function(table,ploidy=3) {
  #Define the position of the first locus
  stlocus<-min(grep("[0-9]",colnames(table)))

  #Define the positions of the loci
  poslist<-vector(mode="list",length = ploidy)
  for(i in 1:length(poslist)) {
    poslist[[i]]<-seq(from=stlocus+(i-1),to=length(table[1,]), by=ploidy)
  }

  #Define the resulting dataframe
  restable<-data.frame()
  for(j in 1:length(table[,1])) {
    if(length(restable)==0) {
      restable<-rbind(restable,c(paste(table[j,1],1,sep="_"),table[j,2],table[j,poslist[[1]]]))
    }
    else {
      restable<-rbind(restable,c(paste(table[j,1],1,sep="_"),as.character(table[j,2]),as.numeric(table[j,poslist[[1]]])))
    }
    restable<-rbind(restable,c(paste(table[j,1],2,sep="_"),as.character(table[j,2]),as.numeric(table[j,poslist[[2]]])))
    if(ploidy==3) {
      restable<-rbind(restable,c(paste(table[j,1],3,sep="_"),as.character(table[j,2]),as.numeric(table[j,poslist[[3]]])))
    }

  }

  #Find the haplotypes
  populations<-levels(as.factor(restable$Pop))
  for(h in populations) {
    indpop<-grep(i,restable$Pop)
    restable[indpop,stlocus:length(restable[1,])]<-apply(restable[indpop,stlocus:length(restable[1,])],MARGIN = 2,FUN = sample)
    restable[indpop,1]<-paste(h,1:length(indpop),sep="_")
  }

  return(restable)
}

#' Randomly Shuffles the Loci From the Haplotypes Within a Population
#'
#' This function takes a table of haplotypes with a population column and shuffles the loci within each population.
#'
#' @param table Table with a column of the population each sample belongs to and with the haploid loci for each marker
#' @return Table of newly generated individuals from the reshuffling of the genetic markers

samhaplo<-function(table) {
  populations<-levels(as.factor(table$Pop))
  stlocus<-min(grep("[0-9]",colnames(table)))
  for(h in populations) {
    indpop<-grep(h,table$Pop)
    table[indpop,stlocus:length(table[1,])] <- apply(table[indpop,stlocus:length(table[1,])],MARGIN = 2,FUN = sample)
    table[indpop,1] <- paste(h,1:length(indpop),sep="_")
  }
  return(table)
}
