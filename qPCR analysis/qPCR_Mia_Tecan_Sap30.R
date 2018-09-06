require("readxl")
#goal: analyze qPCR results - plot 
getwd()
setwd("/Users/adam/Documents/Levy Lab/Projects/Mia's double CRISPRi screening review/qPCR results")
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
delta_ct = matrix(NA,48,1)
for (i in 1:nrow(gene_test)){
  delta_ct[i,] = gene_test[i,2] - endogenous_control[i,2]
}


