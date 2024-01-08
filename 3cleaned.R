# new script 06 01 2023

rm(list=ls())

library(ggplot2)
library(picante)
library(ade4)
library(sf)
#library(raster)
library(visreg)

setwd("C:/Users/aleja/OneDrive/tese D - Omar/Dados/analise cap 3/")
load("cap313062023.RData")

dados<-read.table("dados_omar_13.07.2022.txt",sep="\t",header=T)
dadosabu<-dados[,4*(1:59)-1]
rownames(dadosabu)<-dados$Specie
dim(dadosabu) # 93 x 59

dadosmass<-dados[,4*(1:59)+2]
rownames(dadosmass)<-dados$Specie
dim(dadosmass) # 93 x 59

#data<-read.table("dadosfinal11.04.2022.txt",sep="\t",header=T)
data<-read.table("dataset17.10.2022.txt",sep="\t",header=T)

head(data)

env<-data[,c(95:105)]
coord<-data[,106:107]

#library(ggplot2)
#library(GGally) 
#GGally::ggpairs(env) # see Discharge and water velocity

#
#gr<-c(rep(1,6),rep(2,6),rep(3,5),rep(4,6),rep(5,6),rep(6,6),rep(7,6),rep(8,6),
#      rep(9,6),rep(10,6)) -> see the gr object below

coordID<-data.frame(coord,ID=gr)

xys1 = st_as_sf(coordID[grnew], coords=c("Longitude","Latitude"))

polys1 = st_sf(
  aggregate(
    xys1,
    by=list(coordID$ID),
    do_union=FALSE,
    FUN=function(vals){
      vals[1]
    }
  ))

xypoly1 = st_cast(polys1, 'POLYGON')

polyschull1 <- polys1
st_geometry(polyschull1) <- st_convex_hull(polys1$geometry)

plot(polyschull1$geometry)

sf::st_write(polyschull1,dsn="10groups.shp")
#plot(polys$geometry)
#points(coord,pch=16)
#text(coord,labels=1:59,adj=c(1,1))

grnew<-c(1,2,6,3,5,4,7,8,12,9,11,10,13,14,17,15,16,18,19,23,22,21,20,24,25,
         29,26,28,27,31,30,32,35,33,34,37,36,38,41,39,40,43,42,44,47,46,45,
         49,48,50,53,52,51,54,55,59,56,57,58)

coordID<-data.frame(coord[grnew,],ID=gr)

xys2 = st_as_sf(coordID, coords=c("Longitude","Latitude"))

polys2 = st_sf(
  aggregate(
    xys2,
    by=list(coordID$ID),
    do_union=FALSE,
    FUN=function(vals){
      vals[1]
    }
  ))

xypoly2 = st_cast(polys2, 'POLYGON')

par(mfrow=c(1,2))
plot(xypoly1$geometry)
points(coord,pch=16)
text(coord,labels=grnew,adj=c(1,1))
plot(xypoly2$geometry)
points(coord,pch=16)
text(coord,labels=grnew,adj=c(1,1))
dev.off()

#

# body size
indlargelong<- 4*(1:59)+1
long<-dados[,indlargelong]
long[is.na(long)]<-0
bodysize<-apply(long,MARGIN=1,max)

#in the case of Planaltina britskii, we used the standard longitude 
# reported by Deprá et al 2018, equal to 33.6 mm

# data 19.90 taken from FishBase
# https://www.fishbase.se/summary/Bryconamericus-stramineus
# Bryconamericus-stramineus used as synonimous of Piabarchus stramineus
dadosmass[13,37]<-19.90*64 
dadosmass[,37]<-as.numeric(dadosmass[,37])
dadosmass[is.na(dadosmass)]<-0

dadossp<-read.table("comexot24.07.2022.txt",sep="\t",header=T)
head(dadossp)

abunat<-dadosabu[which(dadossp$Exnatgen=="native"),]
dim(abunat) # 47 species x 59 sites
abunat[is.na(abunat)]<-0
prod(rowSums(abunat))

abuexo<-dadosabu[which(dadossp$Exnatgen=="exotic"),]
dim(abuexo) # 28 species x 59 sites
abuexo[is.na(abuexo)]<-0
prod(colSums(abuexo)) # =0 because some sites has not exotic species

dim(abuexo) # 28 species x 59 sites
dim(abunat) # 47 species x 59 sites

#

nums<-c(1,2,3,4,5,6,7,8)
noms<-c("Upper Paranapanema","Lower Paranapanema","Lower Grande","Upper Grande",
        "Peixes","Aguapei","Lower Tiete","Sao Jose dos Dourados")
cods<-c("UP","LP","LG","UG","P","A","LT","SJD")
nsit<-c(8,9,6,12,6,6,6,6) #c(6,12,8,9,6,6,6,6)

zones<-data.frame(rep(nums,nsit),rep(cods,nsit),rep(noms,nsit))
colnames(zones)<-c("nums","cods","noms")

#

# functional traits

traits<-read.table("fishmorph09012023.txt",sep="\t",header=T)
head(traits)
bry<-which(traits$Genus.species=="Bryconamericus stramineus")
# Bryconamericus stramineus -> Piabarchus stramineus
traits$Genus.species[bry]="Piabarchus stramineus"

hyp<-which(traits$Genus.species=="Hyphessobrycon anisitsi")
# Hyphessobrycon anisitsi -> Psalidodon anisitsi  
traits$Genus.species[hyp]="Psalidodon anisitsi"

#which(traits$Genus.species=="Characidium lagosantensis")

#mat1<-match(rownames(dadosabu)[-which(dadossp$Exnatgen=="genus")],
#      traits$Genus.species)
#lea<-which(is.na(mat1))
#rownames(dadosabu)[-which(dadossp$Exnatgen=="genus")][lea]
#View(traits[na.omit(mat1),])

#bodysize[-which(dadossp$Exnatgen=="genus")]
#rownames(dadosabu)[-which(dadossp$Exnatgen=="genus")]

#length(mat1)

tr4ad<-read.table("traits4aditional.txt",sep="\t",header=T)

traitsn<-rbind(traits,tr4ad)
mat2<-match(rownames(dadosabu)[-which(dadossp$Exnatgen=="genus")],
            traitsn$Genus.species)
sum(is.na(mat2))

trait75g<-traitsn[mat2,]
trait75g[which(trait75g$Genus.species=="Apareiodon ibitiensis"),"RMl"]<-0.718
trait75g[which(trait75g$Genus.species=="Aspidoras fuscoguttatus"),"RMl"]<-0.588
trait75g[which(trait75g$Genus.species=="Corydoras aeneus"),"RMl"]<-0.391
trait75g[which(trait75g$Genus.species=="Phalloceros caudimaculatus"),"PFs"]<-0.181
trait75g[which(trait75g$Genus.species=="Parodon nasus"),"RMl"]<-0.435
trait75g[which(trait75g$Genus.species=="Planaltina britskii"),"MBl"]<-33.6
trait75g[which(trait75g$Genus.species=="Rineloricaria latirostris"),"CPt"]<-2.985
trait75g[which(trait75g$Genus.species=="Synbranchus marmoratus"),"PFv"]<-0
sum(is.na(trait75g)) # the 4 NAs are in SPecCode column

colnames(trait75g)
trait75<-trait75g[,6:15]
rownames(trait75)<-trait75g$Genus.species
sum(is.na(trait75)) # it must be equal to zero


## varpart functional MPD

#library(picante)

dadosabun<-dadosabu

dadosabun[is.na(dadosabun)]<-0
sum(dadosabun) # 14941 individuals

abunat<-dadosabun[which(dadossp$Exnatgen=="native"),]
dim(abunat) # 47 x 59
prod(rowSums(abunat))
#distcom<-vegdist(t(abunat),method="bray")
#dadosLCBDnat<-adespatial::LCBD.comp(as.matrix(distcom))

abuexo<-dadosabun[which(dadossp$Exnatgen=="exotic"),]
dim(abuexo) # 28 x 59
#abuexo[is.na(abuexo)]<-0
prod(colSums(abuexo)) # equal to zero
                      # because some sites do not have exotic species

a<-t(abunat)*0
b<-t(abuexo)*0

matcomabu<-rbind(cbind(t(abunat),b),cbind(a,t(abuexo)))

# estimating the functional distances between species

trait75x5<-trait75[,c("MBl","VEp","REs","OGp","RMl")]
trait75MBl<-trait75$MBl

#library(ade4)

wM1<-vegan::decostand(trait75x5,method="standardize")

ktab1<-ktab.list.df(list(as.data.frame(wM1)))

distrait75 <- dist.ktab(ktab1, type= "Q")

# RDA, ANOVA and Kruskal-Wallis

matcomabun<-matcomabu
names(distrait75)<-colnames(matcomabun)


comdistMPDf<-comdist(matcomabun,as.matrix(distrait75),abundance.weighted=T)
natexoMPDf<-diag(as.matrix(comdistMPDf)[1:59,60:118]) # betaMPD

comdistMNTDf<-comdistnt(matcomabun,as.matrix(distrait75),abundance.weighted=T)
natexoMNTDf<-diag(as.matrix(comdistMNTDf)[1:59,60:118]) # betaMNTD
ind<-which(is.na(0==natexoMNTDf))
#natexoMNTDf[ind]<-0
natexoMNTDfres<-natexoMNTDf[-ind]
natexoMPDfres<-natexoMPDf[-ind]

rdaMPD<-rda(natexoMPDf[-ind],env[-ind,],scale=T)
anova(rdaMPD)
rdaMNTD<-rda(natexoMNTDf[-ind],env[-ind,],scale=T)
anova(rdaMNTD)

boxplot(natexoMPDf[-ind] ~ as.factor(zones$cods[-ind]))
boxplot(natexoMNTDf[-ind] ~ as.factor(zones$cods[-ind]))

anova(aov(natexoMPDf[-ind] ~ as.factor(zones$cods[-ind])))
anova(aov(natexoMNTDf[-ind] ~ as.factor(zones$cods[-ind]))) #has effect
kruskal.test(natexoMPDf[-ind] ~ as.factor(zones$cods[-ind]))
kruskal.test(natexoMNTDf[-ind] ~ as.factor(zones$cods[-ind])) # has effect

# group for betadisper()

#envfull<-cbind(env,hetenv,denspoint,denspointpot,gr)[-ind,]#[33:51,]

rdavif<-vegan::rda(natexoMPDfres ~  MeanWidth  + MeanDepth + WaterVelocity + Discharge + 
                     WaterTemprerature + pH + Condutivity + HorizontalTransparency, 
                   data=env[-ind,])
vegan::vif.cca(rdavif)

# without discharge

rdavif<-vegan::rda(natexoMPDfres ~  MeanWidth  + MeanDepth + WaterVelocity  + 
                     WaterTemprerature + pH + Condutivity + HorizontalTransparency, 
                   data=env[-ind,])
vegan::vif.cca(rdavif) # all less than 2

envscl<-scale(env[,c(3,4,5,8,9,10,11)]) # without Discharge

gr<-c(rep(1,6),rep(2,6),rep(3,5),rep(4,6),rep(5,6),rep(6,6),rep(7,6),rep(8,6),
      rep(9,6),rep(10,6))

bdisenv<-betadisper(dist(envscl),gr)
hetenv<-bdisenv$group.distances
hetenv<-rep(hetenv,each=6)
hetenv<-hetenv[-18]

# populational data

#library(sf)

muni<-read_sf(dsn=".",layer="BRA_adm2")

munin<-muni[muni$NAME_1==c("São Paulo")|muni$NAME_1==c("Paraná"),] #
plot(st_geometry(munin),col="gray")

ordpr<-order(munin$NAME_2[1:401])
ordsp<-order(munin$NAME_2[402:1045])

ordprsp<-c(401 + ordsp, ordpr)

muninn<-munin[ordprsp,]

write.csv(as.data.frame(munin[ordprsp,c(5,7)])[,1:2],"namepoly.csv")

pop<-read.csv("popmuni.csv",sep=";")
pop

muniproj<-st_transform(muninn,crs=st_crs("+proj=utm +zone=22 +south=T +ellps=WGS84"))
plot(st_geometry(muniproj),col="gray")


dens<-pop$pop/(st_area(muniproj)/1000000)

#library(raster)

ext<- raster::extent(c(-55,-43,-27,-19))
raster1<- raster::raster(nrows=100,ncols=150,ext=ext)
raster2<-raster::rasterize(muninn,raster1,field=dens)
#plot(raster2)
denspoint<-raster::extract(raster2,coord)
denspointpot<-denspoint^0.25

#number of exotic species

numexo<-28-colSums(abuexo==0)
numnat<-47-colSums(abunat==0)
plot(numnat,numexo)
cor(numnat,numexo)
cor.test(numnat,numexo)#http://127.0.0.1:18433/graphics/plot_zoom_png?width=690&height=508
plot(numnat,numexo,xlab="native richness",ylab="exotic richness",cex=0.5,pch=16)


envall<-cbind(env,hetenv,denspoint,denspointpot,gr,numexo)[-ind,]

cor(envall[,c("Altitude","Order","hetenv","denspointpot","numexo")])

regr1<-lm(natexoMPDfres ~ Altitude+Order+hetenv+denspointpot+numexo,data=envall)
summary(regr1)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
library(visreg)
visreg(regr1)
hist(resid(regr1))

regr2<-lm(natexoMNTDfres ~ Altitude+Order+hetenv+denspointpot+numexo,data=envall)
summary(regr2)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr2)
hist(resid(regr2))

# group 3

envall1<-cbind(env,hetenv,denspoint,denspointpot,gr,numexo)[-ind,][1:14,]

cor(envall1[,c("Altitude","Order","hetenv","denspointpot","numexo")]) # retire Altitude

natexoMPDfres1<-natexoMPDfres[1:14]
natexoMNTDfres1<-natexoMNTDfres[1:14]

regr11<-lm(natexoMPDfres1 ~ Altitude+Order+hetenv+denspointpot+
             numexo,data=envall1)
summary(regr11)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr11)
hist(resid(regr11))

regr21<-lm(natexoMNTDfres1 ~ Altitude+Order+hetenv+denspointpot+
             numexo,data=envall1)
summary(regr21)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr21)
hist(resid(regr21))

# group 1

envall2<-cbind(env,hetenv,denspoint,denspointpot,gr,numexo)[-ind,][15:32,]

cor(envall2[,c("Altitude","Order","hetenv","denspointpot","numexo")]) 

natexoMPDfres2<-natexoMPDfres[15:32]
natexoMNTDfres2<-natexoMNTDfres[15:32]

regr12<-lm(natexoMPDfres2 ~ Altitude + Order+hetenv+denspointpot+
             numexo,data=envall2)
summary(regr12)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr12)
hist(resid(regr12))

regr22<-lm(natexoMNTDfres2 ~ Altitude+Order+hetenv+denspointpot+
             numexo,data=envall2)
summary(regr22)
087par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr22)
hist(resid(regr22))

# group 2

envall3<-cbind(env,hetenv,denspoint,denspointpot,gr,numexo)[-ind,][33:51,]

cor(envall3[,c("Altitude","Order","hetenv","denspointpot","numexo")])

natexoMPDfres3<-natexoMPDfres[33:51]
natexoMNTDfres3<-natexoMNTDfres[33:51]

regr13<-lm(natexoMPDfres3 ~ Altitude+Order+hetenv+denspointpot+
             numexo,data=envall3)
summary(regr13)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr13)
hist(resid(regr13))

regr23<-lm(natexoMNTDfres3 ~ Altitude+Order+hetenv+denspointpot+
             numexo,data=envall3)
summary(regr23)
par(mfrow=c(2,3),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regr23)
hist(resid(regr23))

#

dev.off()

setwd("C:/Users/aleja/OneDrive/tese D - Omar/graficos/cap3/")

matsig<-matrix(ncol=5,nrow=8,
               c(c(0,1,0,0,0),
               c(0,0,0,0,0),
               c(0,0,0,0,1),
               c(0,0,0,0,0),
               c(0,0,0,0,1),
               c(0,0,0,0,0),
               c(0,0,0,0,1),
               c(0,0,0,0,0)))

pdf("Fig3.pdf",width=8,height=12)

xvari<-c("Altitude","Order","hetenv","denspointpot","numexo")

par(mfrow=c(8,5),oma=c(4,4,0.5,0.5),mai=c(0.2,0.2,0.2,0.2))
p1<-visreg(regr1)#,xvar=xvari[1],line="red")
p2<-visreg(regr12)
visreg(regr13)
visreg(regr11)
visreg(regr2)
visreg(regr22)
visreg(regr23)
visreg(regr21)
#text(c(3,1),labels="claro")
dev.off()

setwd("C:/Users/aleja/OneDrive/tese D - Omar/Dados/analise cap 3")


cor.test(numnat[-ind][1:14],numexo[-ind][1:14])
cor.test(numnat[-ind][15:32],numexo[-ind][15:32])
cor.test(numnat[-ind][33:51],numexo[-ind][33:51])

plot(numnat[-ind][33:51],numexo[-ind][33:51])

cor.test(numnat[-ind],numexo[-ind])
lm(numexo[-ind]~numnat[-ind])
#plot(numnat[-ind],numexo[-ind],pch=16,xlab="native richness",
#      ylab="exotic richness")
numnatexowo<-data.frame(numnat=numnat[-ind],numexo=numexo[-ind])
regne<-lm(numexo~numnat,data=numnatexowo)
visreg::visreg(regne,data=numnatexowo,cex=3,xlab="native richness",ylab="exotic richness")
abline abline(a=1.2383,b=0.1022,col="dodgerblue3",lty=1,lwd=3)


#######################################
# with a large data (second data base)

alltraits<-read.table("traits08062023.txt",sep="\t",header=T)
dim(alltraits) #178 x 27 - since BL is functional traits
allcoord<-read.table("coord05062023.txt",sep="\t",header=T)
dim(allcoord) # 253 x 4
allusosolo<-read.table("usosolo05062023.txt",sep="\t",header=T)
dim(allusosolo) # 252 x 11
allenv<-read.table("env05062023.txt",sep="\t",header=T)
dim(allenv) # 252 x 12
allcomabu<-read.table("comabu05062023.txt",sep="\t",header=T)
dim(allcomabu) # 116 x 255
origin<-read.table("origins08062023.txt",sep="\t",header=T)

colnames(allcoord)
plot(allcoord$X,allcoord$Y)

comabug<-allcomabu[,4:255] # the site 242 (row 179) is not in the allcomabu database
rownames(comabug)<-allcomabu$Especies
comabug<-comabug[,-221] # removing the site 268
dim(comabug) # 116 x 251

library(stringi)
rownam<-as.numeric()
for(i in 1:ncol(comabug)) {
  rownam[i]<-as.numeric(stri_sub(colnames(comabug)[i],2))
}
rownam
sort(rownam)==rownam
length(rownam) #251

allcoord1<-allcoord[match(rownam,allcoord$PONTOS[1:252]),][,c(1,3:4)]
rownames(allcoord1)<- allcoord1$PONTOS
dim(allcoord1) # 251 x 3
coordg<-allcoord1[,2:3]
rownames(coordg)<-allcoord1$PONTOS
head(coordg)
tail(coordg)
dim(coordg) # 251 x 2
sum(rownames(coordg)==rownam) # 251

rownames(alltraits)<-alltraits$Species
alltraits<-alltraits[,-1]
traitsg<-alltraits[,16:26]

#myspecies<-gsub(" ","_",rownames(traitsg))
#rownames(traitsg)<-myspecies

match(rownames(comabug),rownames(traitsg))

#install.packages("remotes")
#remotes::install_github("brunobrr/bdc")
library(bdc)

namcl<-bdc_clean_names(rownames(comabug))
namcl$names_clean[17]<-"Characidium lagosantense"
namcl$names_clean[18]<-"Characidium zebra"

namcl1<-bdc_clean_names(rownames(traitsg))
namcl1$names_clean

myspecies<-gsub(" ","_",namcl$names_clean)
rownames(comabug)<-myspecies

myspecies1<-gsub(" ","_",namcl1$names_clean)
rownames(traitsg)<-myspecies1
rownames(traitsg)[82]<-"Poecilia_reticulata"
rownames(traitsg)[83]<-"Poecilia_vivipara"

namori1<-bdc_clean_names(origin$Especies)
origin$specln<-gsub(" ","_",namori1$names_clean)
origin$specln[13:14]<-c("Characidium_lagosantense","Characidium_zebra")

traitg<-traitsg[match(rownames(comabug),rownames(traitsg)),]

sum(is.na(match(rownames(comabug),origin$specln))) # must be zero

originsort<-origin[match(rownames(comabug),origin$specln),2:3]

envg<-allenv[-221,]
rownames(envg)<-envg$Ponto
envg<-envg[,-1]

allusosolo<-allusosolo[-221,]
dim(allusosolo)                       

comabug<-t(comabug)


                       
dim(comabug)
dim(traitg)
dim(coordg)
dim(envg)
dim(originsort)
dim(allusosolo)

# computing the functional distance between species

library(ade4)

wM<-traitg[,c(1,3,4,5)]

distfg<-dist.ktab(ktab.list.df(list(wM)), type="Q")

# computing climate heterogeneity




abunatg<-comabug[,which(originsort$origin=="native")]
abuexog<-comabug[,which(originsort$origin=="exotic")]
ag<-abunatg*0
bg<-abuexog*0

abuartif<-rbind(cbind(abunatg,bg),cbind(ag,abuexog))

library(picante)

comdi<-comdist(abuartif,as.matrix(distfg),abundance.weighted = TRUE)
comdint<-comdistnt(abuartif,as.matrix(distfg),abundance.weighted = FALSE)
gnatexoMPDf<-diag(as.matrix(comdi)[1:251,252:502])
gnatexoMNTDf<-diag(as.matrix(comdint)[1:251,252:502])
indg<-which(is.na(gnatexoMNTDf))

gnatexoMPDfres<-gnatexoMPDf[-indg]
gnatexoMNTDfres<-gnatexoMNTDf[-indg]

plot(envg$Altitude_1km_media_SRTM[-indg],gnatexoMPDfres)
summary(lm(gnatexoMPDfres ~ Altitude_1km_media_SRTM, data=envg[-indg,]))

plot(envg$Altitude_1km_media_SRTM[-indg],gnatexoMNTDfres)
summary(lm(gnatexoMNTDfres ~ Altitude_1km_media_SRTM, data=envg[-indg,]))

plot(allusosolo$VREM.[-indg],gnatexoMPDfres)
summary(lm(gnatexoMPDfres ~ VREM., data=allusosolo[-indg,]))

plot(allusosolo$VREM.[-indg],gnatexoMNTDfres)
summary(lm(gnatexoMNTDfres ~ VREM., data=allusosolo[-indg,]))

par(mfrow=c(2,2))
plot(envg$Altitude_1km_media_SRTM[-indg],gnatexoMPDfres,pch=16,cex=0.5,col=8)
plot(envg$Altitude_1km_media_SRTM[-indg],gnatexoMNTDfres,pch=16,cex=0.5,col=8)
plot(allusosolo$VREM.[-indg],gnatexoMPDfres,pch=16,cex=0.5,col=8)
plot(allusosolo$VREM.[-indg],gnatexoMNTDfres,pch=16,cex=0.5,col=8)

dev.off()

# populational data

library(sf)

muni<-read_sf(dsn=".",layer="BRA_adm2")

muning<-muni[muni$NAME_1==c("São Paulo")|muni$NAME_1==c("Paraná")|muni$NAME_1==c("Goiás")|muni$NAME_1==c("Minas Gerais"),] #
plot(st_geometry(muning),col="gray")

gordgo<-order(muning$NAME_2[1:242])
gordmg<-order(muning$NAME_2[243:1096])
gordpr<-order(muning$NAME_2[1097:1497])
gordsp<-order(muning$NAME_2[1498:2141])


gordgomgprsp<-c(gordgo,242+gordmg,1096+gordpr,1497+gordsp)

muninng<-muning[gordgomgprsp,]

write.csv(as.data.frame(muning[gordgomgprsp,c(5,7)])[,1:2],"namepolyg.csv")

popg<-read.csv("popmunig.csv",sep=";")

muniprojg<-st_transform(muninng,crs=st_crs("+proj=utm +zone=22 +south=T +ellps=WGS84"))
plot(st_geometry(muniprojg),col="gray")


densg<-popg$pop/(st_area(muniprojg)/1000000)

#library(raster)

extg<- raster::extent(c(-55,-43,-27,-13))
rasterg1<- raster::raster(nrows=150,ncols=200,ext=extg)
rasterg2<-raster::rasterize(muninng,rasterg1,field=densg)
#plot(rasterg2)
denspointg<-raster::extract(rasterg2,coordg)
denspointpotg<-denspointg^0.25

# number of exotic species (46 = total exotic species)

numexog<-(46-rowSums(abuexog==0))
numnatg<-70-rowSums(abunatg==0)
plot(numnatg,numexog)
cor(numnatg,numexog)
cor.test(numnatg,numexog)

envallg<-cbind(envg,allusosolo,denspointg,denspointpotg,numexog)[-indg,] #,hetenv ,gr

cor(envallg[,c("Altitude_1km_media_SRTM","VREM.","denspointpotg","numexog")])

regrg1<-lm(gnatexoMPDfres ~ Altitude_1km_media_SRTM+VREM.+denspointpotg+
             numexog,data=envallg)
regrg1b<-lm(gnatexoMPDfres ~ Altitude_1km_media_SRTM+VREM.+denspointpotg,data=envallg)
summary(regrg1)
par(mfrow=c(2,2),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regrg1)
hist(resid(regrg1))

plot(regrg1$fitted.values,resid(regrg1))

#regrg2<-glm(gnatexoMNTDfres ~ Altitude_1km_media_SRTM+VREM.+
#              denspointpotg, data=envallg, family="Gamma")
regrg2<-lm(gnatexoMNTDfres ~ Altitude_1km_media_SRTM+VREM.+
              denspointpotg+numexog, data=envallg)
summary(regrg2)
par(mfrow=c(2,2),oma=c(1,1.5,0.5,0.5),mai=c(0.7,0.2,0.2,0.2))
visreg(regrg2)
hist(resid(regrg2))

plot(regrg2$fitted.values,resid(regrg2))

dev.off()

################################################################
# Matheus

load("cap313062023.RData")

abunati<-t(abunat)
abuexot<-t(abuexo)
abutota<-cbind(abunati,abuexot)

colnames(abutota)==colnames(as.matrix(distrait75))
colnames(matcomabu)==colnames(abutota)

mpdtot<-picante::mpd(abutota,as.matrix(distrait75))
mntdtot<-picante::mntd(abutota,as.matrix(distrait75))

nritot<-picante::ses.mpd(abutota,as.matrix(distrait75))
ntitot<-picante::ses.mntd(abutota,as.matrix(distrait75))

mpdnat<-picante::mpd(abunati,as.matrix(distrait75)[1:47,1:47])
mntdnat<-picante::mntd(abunati,as.matrix(distrait75)[1:47,1:47])

nrinat<-picante::ses.mpd(abunati,as.matrix(distrait75)[1:47,1:47])
ntinat<-picante::ses.mntd(abunati,as.matrix(distrait75)[1:47,1:47])

matrizglobal<-data.frame(natexoMNTDfres,natexoMPDfres,
                         mpdtotres=mpdtot[-ind],mntdtotres=mntdtot[-ind],
                         nritotres=nritot$mpd.obs.z[-ind],ntitotres=ntitot$mntd.obs.z[-ind],
                         mpdnatres=mpdnat[-ind],mntdnatres=mntdnat[-ind],
                         nrinatres=nrinat$mpd.obs.z[-ind],ntinatres=ntinat$mntd.obs.z[-ind],
                         envall)

plot(matrizglobal$mpdtotres,matrizglobal$natexoMPDfres)

plot(matrizglobal$mpdnatres,matrizglobal$natexoMPDfres)

plot(matrizglobal$nrinatres,matrizglobal$natexoMPDfres)

cor.test(matrizglobal$nrinatres,matrizglobal$natexoMPDfres)

plot(matrizglobal$numexo,matrizglobal$nrinatres)

################################################################

modelm<-lm(matrizglobal$mpdtotres ~ nritotres + numexo, data=matrizglobal)
summary(modelm)
par(mfrow=c(1,2))
visreg(modelm)

modelm<-lm(matrizglobal$mntdtotres ~ nritotres + numexo + hetenv + Order , data=matrizglobal)
summary(modelm)
par(mfrow=c(1,3))
visreg(modelm)

