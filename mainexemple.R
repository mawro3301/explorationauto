setwd("~/Bureau/rapport/")
source("explorautomatise2.R")

# chargement des fichiers
distfile="atlas_guyane_trnH_dissw.txt"
charfile="atlas_guyane_trnH_char.txt"
#lecture du fichier de carractere
charmat=read.delim(charfile)
#lecture du fichier des distances
dismat=read.delim(distfile, comment.char = "#",row.names=1, header=T)
# on enleve les titres des lignes et colonnes 
rownames(dismat) <- c()
colnames(dismat)<-c()
print(dim(dismat))
dismat=data.matrix(dismat)
print(dim(dismat))

#simplification des noms taxonomiques
charmat=simplifynames(charmat)


#liste des noms taxonomiques  
toto= taxonames(charmat)
ordernames=toto[[1]]
familynames=toto[[2]]
genusnames=toto[[3]]
speciesnames=toto[[4]]

#extract from name
taxolevel="order"
taxoname=c("Laurales","Magnoliales")
whichitem=extractfromlotname(charmat,taxoname, taxolevel)
print(whichitem)

#extraction des laurales et des magnoliales
magno=extractfromonename(charmat,taxoname[2], taxolevel)
laural=extractfromonename(charmat,taxoname[1], taxolevel)

# matrice des distances pour les magnoliales 
mat1=dismat[magno,magno]
histo(mat1)
viewmat(mat1)

# matrice des distances pour les laurales 
mat3=dismat[laural,laural]
histo(mat3)
viewmat(mat3)

#matrices des distances inter laurales et magnoliales
mat2=dismat[magno,laural]
histo(mat2)
viewmat(mat2)

# extraction des distances des deux ordres
mat=dismat[whichitem,whichitem]
# visualisation des distances
viewmat(mat)
histo(mat)

#next taxo level
onetaxoname="Magnoliales"
taxolevel="order"
descent=alltoxosfromtaxo(charmat, onetaxoname , taxolevel)
print(descent)

#partition botanique
onetaxoname="Magnoliales"
taxolevel="order"
botafamilymagno=extractnexttaxosfromonename(charmat,onetaxoname, taxolevel)
botafamilymagno

#partition hclust
nbrclasses=2
dismagno=as.matrix(dismat[magno,magno])
hclustmagno=classifhier(dismagno,nbrclasses)

#partition sbm
nbrclasses=2
sbmmagno=classifsbm(dismagno,nbrclasses)

# matrices de contingeance
#hclust vs bota
botahclust=contingeance(hclustmagno,botafamilymagno)
botahclust
#hclust vs sbm
sbmhclust=contingeance(hclustmagno,sbmmagno)
sbmhclust
#bota vs sbm
sbmhclust=contingeance(botafamilymagno,sbmmagno)
sbmhclust
print("end main")