setwd("CHANGE TO YOUR DIRECTORY")
library(astrochron)
library(quantmod)
rm(list=ls())

Basin_correction=read.csv("Basin_corrections_Cramer_carbon.csv")
Pacific_d13C_trend=read.csv("Pacific_d13C_trend_carbon.csv",header = T)


Site926=read.csv("Site_926.csv")
colnames(Site926) <- c("Age (Ma)", "Benthic d13C")
#Site982_Andersson=read.csv("Site_982_Andersson2003_ETP.csv")
#colnames(Site982_Andersson) <- c("Age (Ma)", "Benthic d13C")
Site982=read.csv("Site_982_Drury_wuellstorfi.csv",header = T)
colnames(Site982) <- c("Age (Ma)", "Benthic d13C")
Site1090=read.csv("Site_1090.csv")
colnames(Site1090) <- c("Age (Ma)", "Benthic d13C")
Site1146=read.csv("Site_1146_Holbourn2018_AJD.csv")
colnames(Site1146) <- c("Age (Ma)", "Benthic d13C")
# Site1146=read.csv("Site_1146_Holbourn2018.csv")
# colnames(Site1146) <- c("Age (Ma)", "Benthic d13C")
Site1218=read.csv("Site_1218.csv")
colnames(Site1218) <- c("Age (Ma)", "Benthic d13C")
idx=which(Site1218$`Age (Ma)`<35)
Site1218=Site1218[idx,]
Site1337=read.csv("Site_1337.csv")
colnames(Site1337) <- c("Age (Ma)", "Benthic d13C")
Site1338=read.csv("Site_1338.csv")
colnames(Site1338) <- c("Age (Ma)", "Benthic d13C")
Site1264=read.csv("Site_1264.csv")
colnames(Site1264) <- c("Age (Ma)", "Benthic d13C")
Site1267=read.csv("Site_1267.csv")
colnames(Site1267)<- c("Age (Ma)", "Benthic d13C")
Site1143=read.csv("Site_1143.csv")
colnames(Site1143)<- c("Age (Ma)", "Benthic d13C")
Site1337_Drury=read.csv("Site_1337_Drury.csv")
colnames(Site1337_Drury)<- c("Age (Ma)", "Benthic d13C")

###########################
# Bring to Pacific trend
###########################

for (name in c("Site926", "Site982", "Site1090","Site1146","Site1218","Site1337","Site1338","Site1264","Site1267","Site1143","Site1337_Drury")){

b=c()
eval(parse(text=paste('output=mwStats(sortNave(',name,'),win=0.5,conv = 1,genplot=F)')))
for (i in 1:nrow(output)){
  a=which(abs(Pacific_d13C_trend$Age..GTS2004.-output[i,1])==min(abs(Pacific_d13C_trend$Age..GTS2004.-output[i,1])))  
  b[i]=output[i,2]-Pacific_d13C_trend$Pacific.d13Ctrend[a[1]]
}
b=cbind(output$Center_win,b)
colnames(b)<- c("Center_win","Offset_record_Pacific_trend")
eval(parse(text=paste('N=nrow(',name,')')))
for (j in 1:N) {
  eval(parse(text=paste('a=which(abs(b[,1]-',name,'[j,1])==min(abs(b[,1]-',name,'[j,1])))')))  
  eval(parse(text=paste(name,'[j,3]=',name,'[j,2]-b[a[1],2]')))
  }
eval(parse(text=paste('colnames(',name,') <- c("Age (Ma)", "Benthic d13C","Benthic d13C Pacific")')))
eval(parse(text=paste0('write.csv(',name,',"Pacific_corrected_data/',name,'_Pacific_corr.csv",row.names = F)')))
}

####################################
# Concatenation of records
####################################

BP=c(3.2844,5.174,6.677,7.891,12.834,15.8,20.025,20.215,22.56,23.9)

pos <- which(Site1267[,1]>=0 & Site1267[,1]<=BP[1])
pos1 <- which(Site1264[,1]>=BP[1] & Site1264[,1]<=BP[2])
pos2 <- which(Site982[,1]>=BP[2] & Site982[,1]<=BP[3])
pos3 <- which(Site1337_Drury[,1]>=BP[3] & Site1337_Drury[,1]<=BP[4])
#pos4 <- which(Site1143[,1]>=BP[4] & Site1143[,1]<=BP[5])
pos5 <- which(Site1146[,1]>=BP[4] & Site1146[,1]<=BP[5])
pos6 <- which(Site1338[,1]>=BP[5] & Site1338[,1]<=BP[6])
pos7 <- which(Site1337[,1]>=BP[6] & Site1337[,1]<=BP[7])
pos8 <- which(Site1090[,1]>=BP[7] & Site1090[,1]<=BP[8])
pos9 <- which(Site926[,1]>=BP[8] & Site926[,1]<=BP[9])
pos10 <- which(Site1090[,1]>=BP[9] & Site1090[,1]<=BP[10])
pos11 <- which(Site1218[,1]>=BP[10] & Site1218[,1]<=35)

Splice=rbind(Site1267[pos,c(1,3)],Site1264[pos1,c(1,3)],Site982[pos2,c(1,3)], Site1337_Drury[pos3,c(1,3)],Site1146[pos5,c(1,3)], Site1338[pos6,c(1,3)], Site1337[pos7,c(1,3)],Site1090[pos8,c(1,3)],Site926[pos9,c(1,3)],Site1090[pos10,c(1,3)],Site1218[pos11,c(1,3)])
Splice=Splice[is.finite(Splice[,2]), ]
plot(Splice,type="l")

write.csv(Splice,"d13C_megasplice.csv",row.names=F)
