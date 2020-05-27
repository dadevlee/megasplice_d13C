setwd("CHANGE TO YOUR DIRECTORY")
library(astrochron)
library(quantmod)
library(IRISSeismic)
rm(list=ls())

#### Here, we define the frequency ranges for the phase-analysis (in kyr-1). Not that in the manuscript, we focus on the 100-kyr eccentricity frequency range
f_400_low=1/0.435 
f_400_high=1/0.370
f_obl_low=1/0.045
f_obl_high=1/0.038
f_ecc_low=1/0.135
f_ecc_high=1/0.090


# Site 926 Pälike 2006a ----
Site926_C=read.csv("Site_926.csv") #Reading data
Site926_O=read.csv("Site926.csv") #Reading data
Site926_CO=cbind(Site926_C,Site926_O[,2]) #Combining d13C and d18O data in one matrix

end=length(Site926_CO[,1])

T1 = round(Site926_CO[1,1], digits = 2) # T1 is the younger limit of the FIRST analysis window (sliding window approach)
if (T1 - Site926_CO[1,1] < 0) {
    if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
  }
if (T1 - Site926_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
  }

T2 = round(Site926_CO[end,1], digits = 2) # T2 is the older limit of the LAST analysis window (sliding window approach)
if (T2 - Site926_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site926_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2) # T1_all is an array that contains the younger limits of ALL analysis windows that will be applied to this site
T2_all = seq(T1+1.2, T2, by = 0.2) # T2_all is an array that contains the older limits of ALL analysis windows that will be applied to this site

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2) # This array contains the midpoints (in Ma) of all analyses windows. 
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1)) # These matrixes will be filled with phase and coherence output in the "Loop" that follows
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc
#coh_400=phase_400

# Loop 
for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site926_CO[,1]<T2 & Site926_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site926_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site926_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site926_CO[,1], Site926_CO[,2], xout = time)$y
  d18O_win1=approx(Site926_CO[,1], Site926_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1) # Making time series out of the dataset in question
  
  DF <- crossSpectrum(win1, spans=c(3,5)) # Calculating the cross spectrum
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low) # Finding the frequency window that correspond to 100-kyr eccentricity
  coh_ecc[i,2]=max(DF$coh[idx_ecc]) # Finding the frequency with maximum coherence within that frequency window
  idx_coh=which(DF$coh == coh_ecc[i,2]) 
  phase_ecc[i,2]=DF$phase[idx_coh] # Getting the phase result for the frequency with maximum coherence
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
   # # 405 kyr eccentricity ----
   # idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
   # coh_400=max(DF$coh[idx_400])
   # idx_coh=which(DF$coh == coh_400)
   # phase_400[i,2]=DF$phase[idx_coh]
}

plot(coh_400, ylim = c(0, 1), type="l")

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

# phase_400=cbind(phase_400, phase_400[,1])
# for (i in 1:length(coh_400[,1])){
#   if (coh_400[i,2]<0.3) {phase_400[i,3]=3}
#   else if (coh_400[i,2]<0.6) {phase_400[i,3]=2} 
#   else {phase_400[i,3]=1}
# }  

Site <- c("926_Pälike")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F) # Saving results in the form of csv files
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 926 Wilkens 2017 4 - 8 Ma----
Site926_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/926_Wilkens.csv")
Site926_CO[,1]=Site926_CO[,1]/1000
plot(Site926_CO[,1],Site926_CO[,2], type = "l")
idx2=which(Site926_CO[,1]>4 & Site926_CO[,1]<8)
Site926_CO=Site926_CO[idx2,]

end=length(Site926_CO[,1])

T1 = round(Site926_CO[1,1], digits = 2)
if (T1 - Site926_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site926_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site926_CO[end,1], digits = 2)
if (T2 - Site926_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site926_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site926_CO[,1]<T2 & Site926_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site926_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site926_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site926_CO[,1], Site926_CO[,2], xout = time)$y
  d18O_win1=approx(Site926_CO[,1], Site926_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

plot(phase_ecc, ylim = c(-pi, pi))

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("926_Wilkens")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 926 Wilkens 2017 11 - 14 Ma----
Site926_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/926_Wilkens.csv")
Site926_CO[,1]=Site926_CO[,1]/1000
plot(Site926_CO[,1],Site926_CO[,2], type = "l")
idx2=which(Site926_CO[,1]>11 & Site926_CO[,1]<14)
Site926_CO=Site926_CO[idx2,]

end=length(Site926_CO[,1])

T1 = round(Site926_CO[1,1], digits = 2)
if (T1 - Site926_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site926_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site926_CO[end,1], digits = 2)
if (T2 - Site926_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site926_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site926_CO[,1]<T2 & Site926_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site926_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site926_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site926_CO[,1], Site926_CO[,2], xout = time)$y
  d18O_win1=approx(Site926_CO[,1], Site926_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

plot(phase_ecc, ylim = c(-pi, pi))

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("926_Wilkens_11-14Ma")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site 982 Drury 2018 ----
Site982_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/982_Drury_wuellerstorfi.csv")
Site982_CO=Site982_CO[,c(1,4,5)]
Site982_CO[,1]=Site982_CO[,1]/1000

end=length(Site982_CO[,1])

T1 = round(Site982_CO[1,1], digits = 2)
if (T1 - Site982_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site982_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site982_CO[end,1], digits = 2)
if (T2 - Site982_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site982_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site982_CO[,1]<T2 & Site982_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site982_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site982_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site982_CO[,1], Site982_CO[,2], xout = time)$y
  d18O_win1=approx(Site982_CO[,1], Site982_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

plot(phase_ecc, ylim = c(-pi, pi))

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("982_Drury")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 982 Andersson 2003 ??????----
Site982_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/982_Andersson2003.csv")
Site982_CO=Site982_CO[,c(2:4)]

end=length(Site982_CO[,1])

T1 = round(Site982_CO[1,1], digits = 2)
if (T1 - Site982_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site982_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site982_CO[end,1], digits = 2)
if (T2 - Site982_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site982_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site982_CO[,1]<T2 & Site982_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site982_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site982_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site982_CO[,1], Site982_CO[,2], xout = time)$y
  d18O_win1=approx(Site982_CO[,1], Site982_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("982_Andersson")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 1267 # Bell 2014 ----
Site1267_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/Site_1267/Site1267_all.csv")
Site1267_CO=Site1267_CO[,c(2:4)]
Site1267_CO[,1]=Site1267_CO[,1]/1000

end=length(Site1267_CO[,1])

T1 = round(Site1267_CO[1,1], digits = 2)
if (T1 - Site1267_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1267_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1267_CO[end,1], digits = 2)
if (T2 - Site1267_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1267_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1267_CO[,1]<T2 & Site1267_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1267_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1267_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1267_CO[,1], Site1267_CO[,2], xout = time)$y
  d18O_win1=approx(Site1267_CO[,1], Site1267_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1267_Bell")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 1264 # Bell 2014 ----
Site1264_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1264_Bell.csv")
Site1264_CO=Site1264_CO[,c(2:4)]
Site1264_CO[,1]=Site1264_CO[,1]/1000

end=length(Site1264_CO[,1])

T1 = round(Site1264_CO[1,1], digits = 2)
if (T1 - Site1264_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1264_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1264_CO[end,1], digits = 2)
if (T2 - Site1264_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1264_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1264_CO[,1]<T2 & Site1264_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1264_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1264_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1264_CO[,1], Site1264_CO[,2], xout = time)$y
  d18O_win1=approx(Site1264_CO[,1], Site1264_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1264_Bell")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site 1264-1265 Liebrand 2016 ----
Site1264_1265_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1264_1265_Liebrand.csv")
Site1264_1265_CO=Site1264_1265_CO[,c(2:4)]
Site1264_1265_CO[,1]=Site1264_1265_CO[,1]/1000

end=length(Site1264_1265_CO[,1])

T1 = round(Site1264_1265_CO[1,1], digits = 2)
if (T1 - Site1264_1265_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1264_1265_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1264_1265_CO[end,1], digits = 2)
if (T2 - Site1264_1265_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1264_1265_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1264_1265_CO[,1]<T2 & Site1264_1265_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1264_1265_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1264_1265_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1264_1265_CO[,1], Site1264_1265_CO[,2], xout = time)$y
  d18O_win1=approx(Site1264_1265_CO[,1], Site1264_1265_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1264_1265_Liebrand")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site U1337 Drury 2017 ----
Site1337_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1337_Drury.csv")
Site1337_CO=Site1337_CO[,c(2:4)]
Site1337_CO[,1]=Site1337_CO[,1]/1000

end=length(Site1337_CO[,1])

T1 = round(Site1337_CO[1,1], digits = 2)
if (T1 - Site1337_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1337_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1337_CO[end,1], digits = 2)
if (T2 - Site1337_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1337_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1337_CO[,1]<T2 & Site1337_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1337_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1337_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1337_CO[,1], Site1337_CO[,2], xout = time)$y
  d18O_win1=approx(Site1337_CO[,1], Site1337_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1337_Drury")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site U1337 Tian ----
Site1337_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1337_Tian.csv")
Site1337_CO=Site1337_CO[,c(2:4)]

end=length(Site1337_CO[,1])

T1 = round(Site1337_CO[1,1], digits = 2)
if (T1 - Site1337_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1337_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1337_CO[end,1], digits = 2)
if (T2 - Site1337_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1337_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T2 = 5.92 # Not to create an overlap with Drury on the same Site!

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1337_CO[,1]<T2 & Site1337_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1337_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1337_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1337_CO[,1], Site1337_CO[,2], xout = time)$y
  d18O_win1=approx(Site1337_CO[,1], Site1337_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1337_Tian")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site U1337 Holbourn ----
Site1337_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1337_Holbourn.csv")
Site1337_CO=Site1337_CO[,c(2:4)]
Site1337_CO[,1]=Site1337_CO[,1]/1000

end=length(Site1337_CO[,1])

T1 = round(Site1337_CO[1,1], digits = 2)
if (T1 - Site1337_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1337_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1337_CO[end,1], digits = 2)
if (T2 - Site1337_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1337_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1337_CO[,1]<T2 & Site1337_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1337_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1337_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1337_CO[,1], Site1337_CO[,2], xout = time)$y
  d18O_win1=approx(Site1337_CO[,1], Site1337_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1337_Holbourn")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site 1146 Holbourn 2018 ----
Site1146_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/Site 1146/Re-scaled U1146 for phase analysis.csv")

end=length(Site1146_CO[,1])

T1 = round(Site1146_CO[1,1], digits = 2)
if (T1 - Site1146_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1146_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1146_CO[end,1], digits = 2)
if (T2 - Site1146_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1146_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1146_CO[,1]<T2 & Site1146_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1146_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1146_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1146_CO[,1], Site1146_CO[,2], xout = time)$y
  d18O_win1=approx(Site1146_CO[,1], Site1146_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1146_Holbourn")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site U1338 Holbourn 2014  ----
Site1338_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1338_Holbourn.csv")
Site1338_CO=Site1338_CO[,c(2:4)]

end=length(Site1338_CO[,1])

T1 = round(Site1338_CO[1,1], digits = 2)
if (T1 - Site1338_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1338_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1338_CO[end,1], digits = 2)
if (T2 - Site1338_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1338_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1338_CO[,1]<T2 & Site1338_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1338_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1338_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1338_CO[,1], Site1338_CO[,2], xout = time)$y
  d18O_win1=approx(Site1338_CO[,1], Site1338_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1338_Holbourn")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)
# Site U1338 Drury 2016 ----
Site1338_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1338_Drury_2018agemodel.csv")
Site1338_CO=Site1338_CO[,c(2:4)]
Site1338_CO[,1]=Site1338_CO[,1]/1000

end=length(Site1338_CO[,1])

T1 = round(Site1338_CO[1,1], digits = 2)
if (T1 - Site1338_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1338_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1338_CO[end,1], digits = 2)
if (T2 - Site1338_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1338_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1338_CO[,1]<T2 & Site1338_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1338_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1338_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1338_CO[,1], Site1338_CO[,2], xout = time)$y
  d18O_win1=approx(Site1338_CO[,1], Site1338_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1338_Drury_2018agemodel")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 1090 Billups 2004 ----
Site1090_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1090_Billups.csv")
Site1090_CO=Site1090_CO[,c(2:4)]
Site1090_CO[,1]=Site1090_CO[,1]/1000

end=length(Site1090_CO[,1])

T1 = round(Site1090_CO[1,1], digits = 2)
if (T1 - Site1090_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1090_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1090_CO[end,1], digits = 2)
if (T2 - Site1090_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1090_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1090_CO[,1]<T2 & Site1090_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1090_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1090_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1090_CO[,1], Site1090_CO[,2], xout = time)$y
  d18O_win1=approx(Site1090_CO[,1], Site1090_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1090_Billups")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 1218 Pälike 2006b ----
Site1218_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/1218_Pälike.csv")
Site1218_CO=Site1218_CO[,c(2:4)]

end=length(Site1218_CO[,1])

T1 = round(Site1218_CO[1,1], digits = 2)
if (T1 - Site1218_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1218_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1218_CO[end,1], digits = 2)
if (T2 - Site1218_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1218_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1 =22
T2 = 35

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1218_CO[,1]<T2 & Site1218_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1218_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1218_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1218_CO[,1], Site1218_CO[,2], xout = time)$y
  d18O_win1=approx(Site1218_CO[,1], Site1218_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1218_Pälike")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)

# Site 1334 Beddow 2018 ----
Site1334_CO=read.csv("/Users/pdoc3/Documents/Cenozoic_stack/MegaSplice carbon/April 2019/Site1334_Beddow.csv")
Site1334_CO=Site1334_CO[,c(2:4)]
Site1334_CO[,1]=Site1334_CO[,1]/1000

end=length(Site1334_CO[,1])

T1 = round(Site1334_CO[1,1], digits = 2)
if (T1 - Site1334_CO[1,1] < 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1 + 0.02}
  else {T1 = T1 + 0.01}
}
if (T1 - Site1334_CO[1,1] >= 0) {
  if (round(T1*100) %% 2 == 0) {T1 = T1}
  else {T1 = T1 + 0.01}
}

T2 = round(Site1334_CO[end,1], digits = 2)
if (T2 - Site1334_CO[end,1] < 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2}
  else {T2 = T2 - 0.01}
}
if (T2 - Site1334_CO[end,1] >= 0) {
  if (round(T2*100) %% 2 == 0) {T2 = T2 - 0.02}
  else {T2 = T2 - 0.01}
}

T1_all = seq(T1, T2-1.2, by = 0.2)
T2_all = seq(T1+1.2, T2, by = 0.2)

phase_ecc=seq(T1+0.6,T2-0.6, by = 0.2)
phase_ecc=cbind(phase_ecc, matrix(nrow = length(phase_ecc), ncol = 1))
phase_obl=phase_ecc
phase_400=phase_ecc
coh_ecc=phase_ecc

# Loop 

for (i in 1:length(T1_all)){
  
  T1=T1_all[i]
  T2=T2_all[i]
  
  idx1=which(Site1334_CO[,1]<T2 & Site1334_CO[,1]>T1)
  d13C_win1=astrochron::detrend(Site1334_CO[idx1,c(1,2)], genplot = F)
  d18O_win1=astrochron::detrend(Site1334_CO[idx1,c(1,3)], genplot = F)
  
  res_d13C=length(d13C_win1[,1])/(T2-T1)
  dt = signif(1/res_d13C,1)
  time = seq(from = T1, to = T2, by = dt)
  
  d13C_win1=approx(Site1334_CO[,1], Site1334_CO[,2], xout = time)$y
  d18O_win1=approx(Site1334_CO[,1], Site1334_CO[,3], xout = time)$y
  
  d13C_win1=ts(d13C_win1, frequency = 1/dt, start = T1)
  d18O_win1=ts(d18O_win1, frequency = 1/dt, start = T1)
  win1=ts.union(d13C_win1,d18O_win1)
  
  DF <- crossSpectrum(win1, spans=c(3,5))
  
  # Eccentricity ----
  idx_ecc=which(DF$freq < f_ecc_high & DF$freq > f_ecc_low)
  coh_ecc[i,2]=max(DF$coh[idx_ecc])
  idx_coh=which(DF$coh == coh_ecc[i,2])
  phase_ecc[i,2]=DF$phase[idx_coh]
  
  # Obliquity ----
  idx_obl=which(DF$freq < f_obl_high & DF$freq > f_obl_low)
  coh_obl=max(DF$coh[idx_obl])
  idx_coh=which(DF$coh == coh_obl)
  phase_obl[i,2]=DF$phase[idx_coh]
  
  # 405 kyr eccentricity ----
  idx_400=which(DF$freq < f_400_high & DF$freq > f_400_low)
  coh_400=max(DF$coh[idx_400])
  idx_coh=which(DF$coh == coh_400)
  phase_400[i,2]=DF$phase[idx_coh]
}

# Make an extra colomn in phase_ecc that labels results with low coherency (<0.3) with "3" and high coherency (>0.6) with "1", and medium coherency with "2"
phase_ecc=cbind(phase_ecc, phase_ecc[,1])
for (i in 1:length(coh_ecc[,1])){
  if (coh_ecc[i,2]<0.3) {phase_ecc[i,3]=3}
  else if (coh_ecc[i,2]<0.6) {phase_ecc[i,3]=2} 
  else {phase_ecc[i,3]=1}
}  

Site <- c("1334_Beddow")
write.csv(coh_ecc, paste0("Phase time-series/", Site, "_ecc_coh.csv"), row.names = F)
write.csv(phase_ecc, paste0("Phase time-series/", Site, "_ecc.csv"), row.names = F)
write.csv(phase_obl, paste0("Phase time-series/", Site, "_obl.csv"), row.names = F)
write.csv(phase_400, paste0("Phase time-series/", Site, "_400.csv"), row.names = F)


