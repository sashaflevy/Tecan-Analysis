
#Set working directory
setwd("~/OneDrive - Leland Stanford Junior University/Github/Tecan-Analysis/")

#Import well IDs
well.id.file = "SAP30_ananalysis_Plate1_Well_IDs.txt"
x = read.csv(well.id.file, header = F, row.names = NULL, col.names = NULL)

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
z = strsplit(z, "\t")[[1]][2]
time = 1:length(z)
for(i in 1:length(time)){
  time[i] = as.numeric(z[[i]][2])
}