require("readxl")
#goal: analyze qPCR results - plot 
getwd()
setwd("/Users/adam/Documents/Levy Lab/Projects/Mia's double CRISPRi screening review/qPCR results")
#setwd("/Users/adam/Documents/Levy Lab/Projects/Mia's double CRISPRi screening review/qPCR results")
#Import qPCR Amplification Data (txt file)
Amp_Data_Ct = read.table("Tecan_Mia_qPCR_RQmanager_2.sdm-Amplification Data_2.txt",
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
x = c("COG1_SAP30", "COG1_control", "REB1_SAP30", "REB1_control", "RET2_SAP30", "RET2_control", "control_SAP30", "control_control")
y = 1:96
for(i in 1:length(x)){
  y[(i*3-2):(i*3)] = paste(x[i], "+ATC", sep = "")
  y[(i*3+22):(i*3+24)] = paste(x[i], "-ATC", sep = "")
  y[(i*3+46):(i*3+48)] = paste(x[i], "+ATC", sep = "")
  y[(i*3+70):(i*3+72)] = paste(x[i], "-ATC", sep = "")
}
names(delta_ct) = y
z = delta_ct[1:96]


#CAUTION: control_control and control_SAP30 use different primers

mean_delta_ct = 1:16
names(mean_delta_ct) = c(paste(x, "+ATC", sep = ""), paste(x, "-ATC", sep = ""))
for(i in 1:length(mean_delta_ct)){
  mean_delta_ct[i] = mean(delta_ct[c((i*3-2):(i*3), (i*3+46):(i*3+48))])
}

delta_delta_ct = mean_delta_ct[1:8] - mean_delta_ct[9:16]
names(delta_delta_ct) = x

#test
for(i in 1){
  yyy = delta_ct[(i*3-2):(i*3) + (i*3+46):(i*3+48)]
}



par(mar = c(5, 10, 4, 2) + 0.1)
i = 1
plot( delta_ct[c((i*3-2):(i*3), (i*3+46):(i*3+48))], c(i,i,i,i,i,i), ylim = c(0,16), xlim = c(-2, 3), ylab = "", yaxt = 'n', xlab = "deltaCt")
for(i in 2:16){
  points(delta_ct[c((i*3-2):(i*3), (i*3+46):(i*3+48))], c(i,i,i,i,i,i))
}
axis(2, at = 1:16, labels = names(mean_delta_ct), las = 2)

#controls +ATC
cc_cog = delta_ct[22]
cc_reb = delta_ct[23:24]
cc_ret = delta_ct[70]



#Make a barplot comparing +/- SAP30 in +ATC

pdf(file = "SAP30 effect in ATC.pdf")
mean_matrix = matrix(NA, 2, 3)
colnames(mean_matrix) = c("RET2", "REB1", "COG1")
rownames(mean_matrix) = c("-SAP30","+SAP30")
sem_matrix = mean_matrix

n = delta_ct[which(names(delta_ct) == 'RET2_control+ATC')]
p = delta_ct[which(names(delta_ct) == 'RET2_SAP30+ATC')]
mean_matrix[2,1] = mean(mean(n) - p + 1)
mean_matrix[1,1] = mean(mean(n) - n + 1)
sem_matrix[2,1] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,1] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'REB1_control+ATC')]
p = delta_ct[which(names(delta_ct) == 'REB1_SAP30+ATC')]
mean_matrix[2,2] = mean(mean(n) - p + 1)
mean_matrix[1,2] = mean(mean(n) - n + 1)
sem_matrix[2,2] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,2] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'COG1_control+ATC')]
p = delta_ct[which(names(delta_ct) == 'COG1_SAP30+ATC')]
mean_matrix[2,3] = mean(mean(n) - p + 1)
mean_matrix[1,3] = mean(mean(n) - n + 1)
sem_matrix[2,3] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,3] = sd(mean(n) - n + 1)/sqrt(6)

x = barplot(mean_matrix, main="with ATC",
            ylab="Relative expression", col=c("white","grey"),
            legend = rownames(mean_matrix), beside=TRUE, ylim = c(0,2.5))
segments(x, mean_matrix - sem_matrix, x,
         mean_matrix + sem_matrix, lwd = 1.5)
arrows(x, mean_matrix - sem_matrix, x,
       mean_matrix + sem_matrix,, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()


#Make a barplot comparing +/- SAP30 in -ATC

pdf(file = "SAP30 effect without ATC.pdf")
mean_matrix = matrix(NA, 2, 3)
colnames(mean_matrix) = c("RET2", "REB1", "COG1")
rownames(mean_matrix) = c("-SAP30","+SAP30")
sem_matrix = mean_matrix

n = delta_ct[which(names(delta_ct) == 'RET2_control-ATC')]
p = delta_ct[which(names(delta_ct) == 'RET2_SAP30-ATC')]
mean_matrix[2,1] = mean(mean(n) - p + 1)
mean_matrix[1,1] = mean(mean(n) - n + 1)
sem_matrix[2,1] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,1] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'REB1_control-ATC')]
p = delta_ct[which(names(delta_ct) == 'REB1_SAP30-ATC')]
mean_matrix[2,2] = mean(mean(n) - p + 1)
mean_matrix[1,2] = mean(mean(n) - n + 1)
sem_matrix[2,2] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,2] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'COG1_control-ATC')]
p = delta_ct[which(names(delta_ct) == 'COG1_SAP30-ATC')]
mean_matrix[2,3] = mean(mean(n) - p + 1)
mean_matrix[1,3] = mean(mean(n) - n + 1)
sem_matrix[2,3] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,3] = sd(mean(n) - n + 1)/sqrt(6)

x = barplot(mean_matrix, main="without ATC",
            ylab="Relative expression", col=c("white","grey"),
            legend = rownames(mean_matrix), beside=TRUE, ylim = c(0,2.5))
segments(x, mean_matrix - sem_matrix, x,
         mean_matrix + sem_matrix, lwd = 1.5)
arrows(x, mean_matrix - sem_matrix, x,
       mean_matrix + sem_matrix,, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()



#Make a barplot comparing +/- ATC (no SAP30)
pdf(file = "ATC effect without SAP30 guide.pdf")
mean_matrix = matrix(NA, 2, 3)
colnames(mean_matrix) = c("RET2", "REB1", "COG1")
rownames(mean_matrix) = c("+ATC","-ATC")
sem_matrix = mean_matrix

n = delta_ct[which(names(delta_ct) == 'RET2_control+ATC')]
p = delta_ct[which(names(delta_ct) == 'RET2_control-ATC')]
mean_matrix[2,1] = mean(mean(n) - p + 1)
mean_matrix[1,1] = mean(mean(n) - n + 1)
sem_matrix[2,1] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,1] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'REB1_control+ATC')]
p = delta_ct[which(names(delta_ct) == 'REB1_control-ATC')]
mean_matrix[2,2] = mean(mean(n) - p + 1)
mean_matrix[1,2] = mean(mean(n) - n + 1)
sem_matrix[2,2] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,2] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'COG1_control+ATC')]
p = delta_ct[which(names(delta_ct) == 'COG1_control-ATC')]
mean_matrix[2,3] = mean(mean(n) - p + 1)
mean_matrix[1,3] = mean(mean(n) - n + 1)
sem_matrix[2,3] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,3] = sd(mean(n) - n + 1)/sqrt(6)

x = barplot(mean_matrix, main="No SAP30 guide",
            ylab="Relative expression", col=c("white","grey"),
            legend = rownames(mean_matrix), beside=TRUE, ylim = c(0,4))
segments(x, mean_matrix - sem_matrix, x,
         mean_matrix + sem_matrix, lwd = 1.5)
arrows(x, mean_matrix - sem_matrix, x,
       mean_matrix + sem_matrix, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()


#Make a barplot comparing +/- ATC (w/ SAP30 guide)
pdf(file = "ATC effect with SAP30 guide.pdf")
mean_matrix = matrix(NA, 2, 3)
colnames(mean_matrix) = c("RET2", "REB1", "COG1")
rownames(mean_matrix) = c("+ATC","-ATC")
sem_matrix = mean_matrix

n = delta_ct[which(names(delta_ct) == 'RET2_SAP30+ATC')]
p = delta_ct[which(names(delta_ct) == 'RET2_SAP30-ATC')]
mean_matrix[2,1] = mean(mean(n) - p + 1)
mean_matrix[1,1] = mean(mean(n) - n + 1)
sem_matrix[2,1] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,1] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'REB1_SAP30+ATC')]
p = delta_ct[which(names(delta_ct) == 'REB1_SAP30-ATC')]
mean_matrix[2,2] = mean(mean(n) - p + 1)
mean_matrix[1,2] = mean(mean(n) - n + 1)
sem_matrix[2,2] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,2] = sd(mean(n) - n + 1)/sqrt(6)
n = delta_ct[which(names(delta_ct) == 'COG1_SAP30+ATC')]
p = delta_ct[which(names(delta_ct) == 'COG1_SAP30-ATC')]
mean_matrix[2,3] = mean(mean(n) - p + 1)
mean_matrix[1,3] = mean(mean(n) - n + 1)
sem_matrix[2,3] = sd(mean(n) - p + 1)/sqrt(6)
sem_matrix[1,3] = sd(mean(n) - n + 1)/sqrt(6)

x = barplot(mean_matrix, main="With SAP30 guide",
            ylab="Relative expression", col=c("white","grey"),
            legend = rownames(mean_matrix), beside=TRUE, ylim = c(0,4))
segments(x, mean_matrix - sem_matrix, x,
         mean_matrix + sem_matrix, lwd = 1.5)
arrows(x, mean_matrix - sem_matrix, x,
       mean_matrix + sem_matrix, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
dev.off()