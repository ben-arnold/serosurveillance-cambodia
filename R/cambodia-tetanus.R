

#----------------------------------
# cambodia-tetanus.R
#
# estimate age-dependent antibody
# curves for tetanus in cambodia
#----------------------------------

#----------------------------------
# input files:
#   cambodia_serology_public.rds
#
# output files:
#  cambodia-tetanus-parity-quant-seroprev.pdf
#----------------------------------

#----------------------------------
# preamble
#----------------------------------

rm(list=ls())
library(SuperLearner)
library(tmle)
library(tmleAb)  # this is a development package on GitHub: https://github.com/ben-arnold/tmleAb
library(scales)

#----------------------------------
# load the dataset
#----------------------------------
d <- readRDS("~/dropbox/cambodia/data/final/cambodia_serology_public.rds")



# recode 4 observations with negative
# values to 1
table(d$ttmb<=0)
d$ttmb[d$ttmb<=0] <- 1

# create an indicator of positive based on
# the multiplex bead assay
d$ttpos <- ifelse(d$ttmb>100,1,0)

# create an indicator for parous vs nulliparus
d$parous <- ifelse(d$parity=='1:Parous',1,0)

# create an age category variable
d$agecat <- cut(d$age,breaks=c(14,19,24,29,34,39))

#----------------------------------
# estimate the quantiative
# age-dependent antibody
# curve for tetanus toxiod, 
# stratified by parity
#----------------------------------

SL.library <- c("SL.mean","SL.glm","SL.gam","SL.loess","SL.randomForest","SL.polymars")
set.seed(987325)
pacurve <- agecurveAb(Y=log10(d$ttmb[d$parity=='1:Parous']),Age=d$age[d$parity=='1:Parous'],
                      W=subset(d,parity=='1:Parous',select=c("region")),
                      id=d$psuid[d$parity=='1:Parous'],
                      SL.library=SL.library)

set.seed(987325)
npcurve <- agecurveAb(Y=log10(d$ttmb[d$parity=='2:Nulliparous']),Age=d$age[d$parity=='2:Nulliparous'],
                      W=subset(d,parity=='2:Nulliparous',select=c("region")),
                      id=d$psuid[d$parity=='2:Nulliparous'],
                      SL.library=SL.library)


#----------------------------------
# estimate overall mean MFI 
# by parity and age group
#
# not used in the mansucript, due to
# slight redundancy with the
# seroprevalence estimates (below)
# (albeit with higher precision on
#  on the means here)
#----------------------------------
# count of obs in strata and crude means
Nobs <- tapply(log10(d$ttmb),list(d$agecat,d$parity),function(x) length(x))
Nobs
mus <- tapply(log10(d$ttmb),list(d$agecat,d$parity),function(x) mean(x,na.rm=T))
mus


# reduced library to just estimate simple means flexibly adjusted for age and region 
# w/ correct CIs that account for clustering
SL.library <- c("SL.mean","SL.glm","SL.gam")
agecats <- levels(d$agecat)
set.seed(987325)
EYx_p <- sapply(agecats,function (x) tmleAb(Y=log10(d$ttmb[d$parity=='1:Parous' & d$agecat==x]),
                                            W=subset(d,parity=='1:Parous' & agecat==x,select=c("age","region")),
                                            id=d$psuid[d$parity=='1:Parous' & d$agecat==x],
                                            SL.library=SL.library)[c("psi","se","lb","ub")])

set.seed(987325)
EYx_n <- sapply(agecats,function (x) tmleAb(Y=log10(d$ttmb[d$parity=='2:Nulliparous' & d$agecat==x]),
                                                 W=subset(d,parity=='2:Nulliparous' & agecat==x,select=c("age","region")),
                                                 id=d$psuid[d$parity=='2:Nulliparous' & d$agecat==x],
                                                 SL.library=SL.library)[c("psi","se","lb","ub")])

set.seed(987325)
EYx_diff <- sapply(agecats,function (x) tmleAb(Y=log10(d$ttmb[d$agecat==x]),
                                               X=d$parous[d$agecat==x],
                                               W=subset(d,agecat==x,select=c("age","region")),
                                               id=d$psuid[d$agecat==x],
                                               SL.library=SL.library)[c("psi","se","lb","ub","p")])


#----------------------------------
# estimate seroprevalence 
# by parity and age group
#----------------------------------

# count of obs in strata and crude means
Nobs <- tapply(log10(d$ttmb),list(d$agecat,d$parity),function(x) length(x))
Nobs
prevs <- tapply(log10(d$ttpos),list(d$agecat,d$parity),function(x) mean(x,na.rm=T))
prevs


# reduced library to just estimate simple means flexibly adjusted for age and region 
# w/ correct CIs that account for clustering
SL.library <- c("SL.mean","SL.glm","SL.gam")
agecats <- levels(d$agecat)
set.seed(987325)
EYxp_p <- sapply(agecats,function (x) tmleAb(Y=d$ttpos[d$parity=='1:Parous' & d$agecat==x],
                                            W=subset(d,parity=='1:Parous' & agecat==x,select=c("age","region")),
                                            id=d$psuid[d$parity=='1:Parous' & d$agecat==x],
                                            SL.library=SL.library)[c("psi","se","lb","ub")])

set.seed(987325)
EYxp_n <- sapply(agecats,function (x) tmleAb(Y=d$ttpos[d$parity=='2:Nulliparous' & d$agecat==x],
                                            W=subset(d,parity=='2:Nulliparous' & agecat==x,select=c("age","region")),
                                            id=d$psuid[d$parity=='2:Nulliparous' & d$agecat==x],
                                            SL.library=SL.library)[c("psi","se","lb","ub")])

set.seed(987325)
EYxp_diff <- sapply(agecats,function (x) tmleAb(Y=d$ttpos[d$agecat==x],
                                               X=d$parous[d$agecat==x],
                                               W=subset(d,agecat==x,select=c("age","region")),
                                               id=d$psuid[d$agecat==x],
                                               SL.library=SL.library)[c("psi","se","lb","ub","p")])


#----------------------------------
# Create a figure of results
#----------------------------------

# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"


pdf("~/dropbox/cambodia/results/figs/cambodia-tetanus-parity-quant-seroprev.pdf",width=9,height=4)
lo <- layout(mat=matrix(1:2,nrow=1,ncol=2),widths=c(1,1))

cols=c(cblue,cgreen)


## Age dependent mean response
# jitter the ages just a little to improve visualization
ytics <- 0:5
xtics <- seq(15,40,by=5)
op <-par(mar=c(4,4,2,4)+0.1,xpd=T)
plot(1,1,type="n",yaxt="n",xaxt="n",bty="n",ylim=c(0,5),xlim=range(xtics),xlab="",ylab="")
axis(2,at=ytics,labels=c(
  expression(10^0),
  expression(10^1),
  expression(10^2),
  expression(10^3),
  expression(10^4),
  expression(10^5)), las=1,cex.axis=1.25)
axis(1,at=xtics,las=1,cex.axis=1.25)
mtext("Age (years)",side=1,line=2.5,las=1)
mtext("Luminex response (MFI-bg)",side=2,line=3)
mtext("Tetanus age-dependent mean response",side=3,line=0,adj=0)

# parous
points(jitter(pacurve$Age),pacurve$Y,pch=16,cex=0.5,col=alpha(cols[1],alpha=0.3))
lines(pacurve$Age,pacurve$pY,col=cols[1])

# nulliparous
points(jitter(npcurve$Age),npcurve$Y,pch=16,cex=0.5,col=alpha(cols[2],alpha=0.3))
lines(npcurve$Age,npcurve$pY,col=cols[2])

# group labels
text(40,pacurve$pY[length(pacurve$pY)],"Parous", col=cols[1],cex=1,adj=0)
text(40,npcurve$pY[length(npcurve$pY)],"Nulliparous", col=cols[2],cex=1,adj=0)


## Seroprotection by age category
ytics <- seq(0,1,by=0.2)
op <-par(mar=c(4,4,2,1)+0.1)
xc <- 1:5
x1s <- xc-0.1
x2s <- xc+0.1
plot(xc,xc,type="n",yaxt="n",xaxt="n",bty="n",ylim=c(0,1),xlim=c(0.5,5.5),xlab="",ylab="")
axis(2,at=ytics,labels=sprintf("%1.0f",ytics*100),las=1)
mtext(c("15-19","20-24","25-29","30-34","35-39"),at=xc,side=1,line=1,las=1,cex=1.1)
mtext("Age category (years)",side=1,line=2.5,las=1)
mtext("Seroprotected (%)",side=2,line=2.5)
mtext("Tetanus seroprotected (%)",side=3,line=0,adj=0)

arrows(x0=x1s,y0=unlist(EYxp_p[3,]),y1=unlist(EYxp_p[4,]),angle=90,length=0.05,code=3,col=cols[1])
points(x1s,unlist(EYxp_p[1,]),pch=19,bg="white",col=cols[1])

arrows(x0=x2s,y0=unlist(EYxp_n[3,]),y1=unlist(EYxp_n[4,]),angle=90,length=0.05,code=3,col=cols[2])

points(x2s,unlist(EYxp_n[1,]),pch=21,bg="white",col=cols[2])


# group labels
text(x=x1s[2],y=unlist(EYxp_p[1,2]),"Parous", col=cols[1],cex=1,pos=2,adj=1)
text(x=x2s[1],y=unlist(EYxp_n[1,1]),"Nulliparous",col=cols[2],cex=1,pos=4,adj=0)


par(op)
lo <- layout(mat=matrix(1))
dev.off()






