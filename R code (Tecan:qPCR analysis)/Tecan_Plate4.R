
#Libraries
require(pracma)

#Set working directory
setwd("~/OneDrive - Leland Stanford Junior University/Github/Tecan-Analysis/")

#Import well IDs
well.id.file = "SAP30_ananalysis_Plate1_Well_IDs.txt"
x = as.matrix(read.csv(well.id.file, header = F, strip.white = T))
ids = x[2,]

#Import Tecan Data
tecan.file = "2018-05-21_18-46_180521_Adam_Mia_4_96_T7.txt"
x = read.csv(tecan.file, skip = 19)
x = x[1:(nrow(x) - 2), ]

#Matix of ODs
od = as.matrix(x[,3:ncol(x)])

#Vector of temperatures
temp = x[,2]

#Vector of Times (s)
y = as.vector(x[,1])
z = strsplit(y, "\t")
time = 1:length(z)
for(i in 1:length(time)){
  time[i] = as.numeric(z[[i]][2])
}

#Find the area under the curve
auc = 1:ncol(od)
for(i in 1:ncol(od)){
  auc[i] = trapz(time, od[,i])
}

#Compare 
cog8_sap30_P = c(1:2, 49:50)
cog8_neg_P = c(4:6, 52:54)
reb1_sap30_P = c(7:9, 55:57)
reb1_neg_P = c(10:12, 58:60)
ret2_sap30_P = c(13:15, 61:63)
ret2_neg_P = c(16:18, 64:66)
neg_sap30_P = cog8_sap30_P + 18
neg_neg_P = cog8_sap30_P + 21
cog8_sap30_N = c(25:26, 73:74)
cog8_neg_N = cog8_sap30_N + 3
reb1_sap30_N = cog8_sap30_N + 6
reb1_neg_N = cog8_sap30_N + 9
ret2_sap30_N = cog8_sap30_N + 12
ret2_neg_N = cog8_sap30_N + 15
neg_sap30_N = cog8_sap30_N + 18
neg_neg_N = cog8_sap30_N + 21


condition = c("cog8_sap30_P", "cog8_neg_P", "reb1_sap30_P", "reb1_neg_P", "ret2_sap30_P", "ret2_neg_P", "neg_sap30_P", "neg_neg_P", "cog8_sap30_N", "cog8_neg_N", "reb1_sap30_N", "reb1_neg_N", "ret2_sap30_N", "ret2_neg_N", "neg_sap30_N", "neg_neg_N" )



#plot

#General
pdf(file = "P4_Rep.pdf", width = 14, height = 14)
par(mfrow = c(4,4))
for(j in 1:length(condition)){
  x = eval(parse(text = condition[j])) 
  plot(time/60, od[,x[1]], type = "l" , xlab = "Time (minutes)", ylab = "OD595", main = condition[j], ylim = c(0,0.6))
  for(i in 2:length(x)){
    points(time/60, od[,x[i]], type = "l")
  }
}
dev.off()

#+ATC vs -ATC
pdf(file = "P4_ATC.pdf", width = 16, height = 8)
par(mfrow = c(2,4))
for(j in 1:(length(condition)/2)){
  x = eval(parse(text = condition[j])) 
  plot(time/60, od[,x[1]], type = "l" , xlab = "Time (minutes)", ylab = "OD595",  ylim = c(0,0.6))
  for(i in 2:length(x)){
    points(time/60, od[,x[i]], type = "l")
  }
  x = eval(parse(text = condition[j+8])) 
  for(i in 1:length(x)){
    points(time/60, od[,x[i]], type = "l", col = "red")
  }
  text(400, 0.6, labels = condition[j])
  text(400, 0.55, labels = condition[j+8], col = "red")
}
dev.off()


#+/-SAP30
pdf(file = "P4_SAP30.pdf", width = 16, height = 8)
par(mfrow = c(2,4))
for(j in (1:8)*2 -1){
  x = eval(parse(text = condition[j])) 
  plot(time/60, od[,x[1]], type = "l" , xlab = "Time (minutes)", ylab = "OD595", ylim = c(0,0.6))
  for(i in 2:length(x)){
    points(time/60, od[,x[i]], type = "l")
  }
  x = eval(parse(text = condition[j+1])) 
  for(i in 1:length(x)){
    points(time/60, od[,x[i]], type = "l", col = "red")
  }
  text(400, 0.6, labels = condition[j])
  text(400, 0.55, labels = condition[j+1], col = "red")
}
dev.off()

