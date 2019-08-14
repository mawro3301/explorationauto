library("Rcpp")
#library("corrplot")
require("parallel")
require("digest")
library("blockmodels")
library("stats")
#library("igraph")
library("grDevices")
#require(graphics)
#library(tidygraph)
#library("ggplot2")
#library("ggraph")
#library(RColorBrewer)
library("fastcluster")

library(lsa)

library(ape)

#scales::show_col(few_pal()(8))
#vector=few_pal()(8)
###repertoires
# on commence par donner les chemins des repertoirs
## pour cela, soit il faudrait modifier les repertoir dans les lignes suivantes soit renseigner lors de l'execution







simplifynames<-function(char){
  #arguments :
  # @char : data frame of characters
  #               fichiers des charact?res
  #
  #value :
  # @mat : data frame of characters
  #               fichiers des charact?res simplifi?s
  b=levels(char$Order)
  c=levels(char$Family)
  d=levels(char$Genus)
  e=levels(char$Species)
  n=dim(char)
  mat=matrix(0,n[1],4)
  for (i in 1:n[1]){
    #print("nouvelle boucle")
    x=c(b[char$Order[i]], c[char$Family[i]],d[char$Genus[i]], e[char$Species[i]])
    #print(x)
    liste2=c()
    for (j in 1: 4){
      y=as.character(x[j])
      #print(y)
      toto=unlist(strsplit(y,split=" ", fixed = T))
      #print(toto)
      liste2=append(liste2, toto[[1]])
      for (k in 1: length(liste2)){
        mat[i,k]=liste2[k]
      }
    }
    #liste2
    
  }
  mat=as.data.frame(mat)
  names(mat)=c("order","family","genus","specie")
  mat=as.data.frame(mat)
  return(mat)
}

taxonames<-function(char){
  #arguments :
  # @char : data frame of characters
  #         fichiers des charact?res
  #
  #value :
  # @mat : liste of characters
  #        liste des niveaux taxonomiques avec les noms des taxons simplifi?s
  
  b_order=levels(as.factor(char[,1]))
  b_family=levels(as.factor(char[,2]))
  b_genus=levels(as.factor(char[,3]))
  b_species=levels(as.factor(char[,4]))
  #
  res=list(b_order,b_family,b_genus,b_species)
  return(res)
}


#extraction des noms 
extractfromnames<-function(char,taxoname, taxolevel){
  #arguments :
  # @char : data frame of characters
  #         fichiers des charact?res
  # @taxoname : liste of characters
  #             liste de noms taxonomiques
  # @taxolevel : character
  #              le nom d'un niveau taxonomique
  #value :
  # @mat : vector of integers 
  #               liste des ?l?ments \in  ce  taxon 
  liste=c()
  for (k in 1:length(taxoname)){
    liste2=which(char[,taxolevel]==taxoname[k])
    liste=append(liste2, liste)
  }
  liste3=c()
  for (k in 1 :length(liste)){
    liste3[k]=liste[[k]]
  }
  return(liste3)
}
alltoxosfromtaxo<-function(char, onetaxoname , taxolevel){
  #arguments :
  # @char : data frame of characters
  #         fichiers des charact?res
  # @taxoname : liste of characters
  #             liste de noms taxonomiques
  # @taxolevel : character
  #              le nom d'un niveau taxonomique
  #value :
  # @mat : liste of char 
  #         liste des sous-taxons  \in  ce  taxon 

  l=extractfromnames(char, onetaxoname, taxolevel)
  if (taxolevel=="order"){
    d=c()
    c=char$family[l]
    for (k in 1:length(l)){
      d[k]=as.character(c[[k]]) 
    }
    
    d=as.factor(d)
    
    return(levels(d))
    
  }
  if (taxolevel=="family"){
    d=c()
    c=char$genus[l]
    for (k in 1:length(l)){
      d[k]=as.character(c[[k]]) 
    }
    
    d=as.factor(d)
    return(levels(d))
    
  }
  if (taxolevel=="genus"){
    c=char$specie[l]
    d=c()
    for (k in 1:length(l)){
      d[k]=as.character(c[[k]]) 
    }
    
    d=as.factor(d)
    return(levels(d))
  }
  
}
viewmat<-function(mat){
  #arguments :
  # @mat : matrix of doubles
  #         sous-matrice des distances
  #value :
  # @mat : image of matrix 
  #               heatmap de la sous matrice des distances 
  
  p=dim(mat)[2]
  mat2=t(mat)
  mat3=mat2[p:1,]
  image(mat3)
}


heatmapper<-function(dist, char, taxoname, taxolevel){
  mat=subdismatrix(dist, char, taxoname, taxolevel)
  dist[l,l]
  
  heatmap(mat,  Colv = NA, Rowv = NA, scale="column")
}

histo<-function(mat){
  #arguments :
  # @mat : matrix of doubles
  #         sous-matrice des distances
  #value :
  # @mat : histogramme of matrix 
  #       histogramme de la sous matrice des distances 
  
  hist(mat, col="darkred", main="Titre")
  abline(v=mean(mat), col="darkcyan")
  abline(v=median(mat),col="darkgreen")
  abline(v=quantile(mat)[2], col ="gray")
  abline(v=quantile(mat)[4], col ="gray")
  legend(x= "topleft",legend=c("moyenne", "mediane", "quantile 1","quantile 3"),col=c("darkred", "darkcyan","gray", "gray"), lty=1:2, cex=0.8)
}


extractfromonename<-function(char,onetaxoname, taxolevel){
  #arguments :
  # @char : data frame of characters
  #         fichiers des charact?res
  # @onetaxoname : character
  #              nom taxonomiques
  # @taxolevel : character
  #              le nom d'un niveau taxonomique
  #value :
  # @mat : vector of integers 
  #               liste des ?l?ments \in  ce  taxon 
  taxolevels=c("order","family","genus","specie")
  n=which(taxolevels==taxolevel)
  liste=which(char[,n]==onetaxoname)
  liste2=c()
  for (k in 1 :length(liste)){
    liste2[k]=liste[[k]]
  }
  return(liste2)
}

extractfromlotname<-function(char,taxoname, taxolevel){
  #arguments :
  # @char : data frame of characters
  #         fichiers des charact?res
  # @taxoname : list of character
  #              noms taxonomiques
  # @taxolevel : character
  #              le nom d'un niveau taxonomique
  #value :
  # @liste : vector of characters 
  #           liste des ?l?ments \in  ces  taxons
  liste=c()
  for (k in 1:length(taxoname)){
    liste2=extractfromonename(char,taxoname[k], taxolevel)
    liste=append(liste2, liste)
  }
  return(liste)
}


extractnexttaxosfromonename<-function(char,onetaxoname, taxolevel){
  #arguments :
  # @char : data frame of characters
  #         fichiers des charact?res
  # @onetaxoname : character
  #              nom taxonomiques
  # @taxolevel : character
  #              le nom d'un niveau taxonomique
  #value :
  # @liste3 : vector of integers 
  #       liste des classes des ?l?ments \in  ce  taxon
  liste=extractfromonename(char,onetaxoname, taxolevel)
  liste2=liste
  char2=c("order","family","genus")
  n=which(char2==taxolevel)
  liste2=as.character(char[,n+1][liste]) 
  l=length(unique(liste2))
  l2=unique(liste2)
  liste3=vector("integer",length(liste2))
  for (k in 1:length(liste2)){
    for (i in 1:l){
      if (liste2[k]==l2[i]){
        liste3[k]=i
      }
      
    }
  }
  return(liste3)
}

classifhier<-function(distances,nbrclasses){
  # @distances : matrix of double
  #              sous-matrice des distances
  # @nbrclasses : integer
  #              le nombre de classes
  #value :
  # @vec : vector of integers 
  #       liste des classes des ?l?ments \in  ce  taxon (bas? sur une classification hclust)
  

  vecteurdeclassesg=c()
  hc<-hclust.vector(distances, members = NULL)
  vecteurdeclassesg<- cutree(hc, nbrclasses)
  vec=vector("double", length(vecteurdeclassesg))
  for (i in 1:length(vecteurdeclassesg)){
    vec[i]=vecteurdeclassesg[i]
  }
  plot(hc, hang = -1, cex = 0.6,col="darkred")
  return(vec)
}


classifsbm<-function(distances,nbrclasses){
  # @distances : matrix of double
  #              sous-matrice des distances
  # @nbrclasses : integer
  #              le nombre de classes
  #value :
  # @vec : vector of integers 
  #       liste des classes des ?l?ments \in  ce  taxon (bas? sur une classification sbm)
  my_modelg<- BM_poisson("SBM_sym",distances, explore_min=nbrclasses,explore_max=nbrclasses)
  my_modelg$estimate()
  matricepredg =my_modelg$memberships[[nbrclasses]]$Z
  vecteurdeclassesg=c()
  n=dim(distances)[1]
  for (i in 1:n){
    vecteurdeclassesg[i]=which.max(matricepredg[i,])
  }
  return(vecteurdeclassesg)
}


contingeance<-function(partition1, partition2){
  #arguments :
  # @partition1 : vector of integers
  #               premi?re partition
  # @partition2 : vector of integers
  #               deuxi?me partition
  #
  #value :
  # @conting : matrix of integers 
  #            contigency table betwen two partitions
  conting=matrix(0,max(partition1),max(partition1))
  for (a in 1:length(partition1)){
    i=partition1[a]
    j=partition2[a]
    conting[i,j]=conting[i,j]+1
  }
  
  return(conting)
}


print("end")