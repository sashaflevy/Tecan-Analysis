require("readxl")
#goal: analyze qPCR results - plot 
getwd()
#setwd("/Users/adam/Documents/Levy Lab/Projects/Mia's double CRISPRi screening review/qPCR results")
setwd("~/OneDrive - Leland Stanford Junior University/Github/Tecan-Analysis/qPCR analysis 2/")
#Import qPCR Amplification Data (txt file)
Amp_Data_Ct = read.table("Tecan_Mia_qPCR_RQmanager_1.sdm-Amplification Data_2.txt",
                         skip = 4,
                         header = TRUE,
                         sep = "\t",
                         fill = TRUE)

#Delete unnecessary columns
Amp_Data_Ct = Amp_Data_Ct[,-c(1,2,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)]                         

#split matrix (gene test & endogenous control)
gene_test = Amp_Data_Ct[seq(1,nrow(Amp_Data_Ct), by = 2),]
endogenous_control = Amp_Data_Ct[seq(2,nrow(Amp_Data_Ct), by = 2),]

#calculate delta Ct
delta_ct = as.vector(gene_test[,2] - endogenous_control[,2])
x = c("COG8_SAP30", "COG8_control", "REB1_SAP30", "REB1_control", "RET2_SAP30", "RET2_control", "control_SAP30", "control_control")
y = 1:48
for(i in 1:length(x)){
  y[(i*3-2):(i*3)] = paste(x[i], "+ATC", sep = "")
  y[(i*3+22):(i*3+24)] = paste(x[i], "-ATC", sep = "")
}
names(delta_ct) = y



mean_delta_ct = 1:16
names(mean_delta_ct) = c(paste(x, "+ATC", sep = ""), paste(x, "-ATC", sep = ""))
for(i in 1:length(mean_delta_ct)){
  mean_delta_ct[i] = mean(delta_ct[(i*3-2):(i*3)])
}

delta_delta_ct = mean_delta_ct[1:8] - mean_delta_ct[9:16]
names(delta_delta_ct) = x



par(mar = c(10, 4, 4, 2) + 0.1)
i = 1
plot(c(i,i,i), delta_ct[(i*3-2):(i*3)], xlim = c(0,16), ylim = c(-1, 5), xlab = "", xaxt = 'n', ylab = "deltaCt")
for(i in 2:16){
  points(c(i,i,i), delta_ct[(i*3-2):(i*3)], xlim = c(0,16), ylim = c(-1, 5))
}
axis(1, at = 1:16, labels = names(mean_delta_ct), las = 2)

