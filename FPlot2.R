F.PLOT <- function(baseline="N", reference="r-1", color="Y") {
print(paste("reference = ", reference))
#color="Y"
#reference="r-1"
#baseline="N"

setwd("/Users/justinsmith/Dropbox/CRISPR/TECAN_Plots")
require("RColorBrewer")
reds <- brewer.pal(5, "Reds")
greens <- brewer.pal(5, "Greens")
gradient <- c(rev(reds),"white",greens)
data <- readLines(file.choose())
# where does the data start
start <- grep("^Cycle_Header", data)

# get file name
exp_name <- grep("^Exp_Name", data)
exp_name <- strsplit(data[exp_name], "\t")[[1]][2]
 
# get experiment date
exp_date <- grep("^Exp_StartDate", data)
exp_date <- strsplit(data[exp_date], "\t")[[1]][2]
 
# SOME MORE INFORMATION
plate <- "not assigned"
name <- exp_name
user <- grep("^User_Name", data)
user <- strsplit(data[user], "\t")[[1]][2]
exp_time <- grep("^Exp_StartTime", data)
exp_time <- strsplit(data[exp_time], "\t")[[1]][2]
 
reader <- grep("^Reader_Name", data)
reader <- strsplit(data[reader], "\t")[[1]][2]
gsub
 
# number of reads of the complete run
run.length <- length(data)-start
D <- data[(start+1):length(data)]
od <- c()
for (i in 1:run.length) {
if (length(grep("^Cycle_", D[i])) == 1) {
od <- cbind(od, as.numeric(strsplit(strsplit(D[i], "\t")[[1]][2], ",")[[1]]))
}
}
OD <- od[3:dim(od)[1],]
 
setwd("/Users/justinsmith/Dropbox/CRISPR/TECAN_Plots")
 
# if user wants to have a baseline correction
# use first 10 reads in every well to normalize the baseline to 0
if (baseline == "Y") {
od.norm <- sweep(OD, 1, rowMeans(OD[,1:10]), "-")
OD <- round(od.norm, digit=4)
}
 
# variable definition
coo.384 <- paste(rep(LETTERS[1:16], each=24), c("01","02","03","04","05","06","07","08", "09", 10:24), sep="")
coo.96 <- paste(rep(LETTERS[1:8], each=12), c("01","02","03","04","05","06","07","08", "09", 10:12), sep="")
coo.48 <- paste(rep(LETTERS[1:6], each=8), c("01","02","03","04","05","06","07","08"), sep="")
 
# create a table of reference wells to use for plotting and RelG-calculation
# let's parse the input
ref <- strsplit(reference, "-")[[1]]
# create reference table for 48-well plate
if (dim(OD)[1] == 48) {
  if (ref[1] == "w") { colref <- cbind(as.numeric(ref[2]), 1:48) }
  if (ref[1] == "c") { colref <- cbind(rep(seq(as.numeric(ref[2]), 48, 8), each=8), 1:48) }
  rr <- matrix(1:48, ncol=8, byrow=T)
  if (ref[1] == "r") { colref <- cbind(rep(rr[as.numeric(ref[2]),], 6), 1:48) }
}
#if (dim(OD)[1] == 48) {
#if (ref[1] == "w") { colref <- cbind(as.numeric(ref[2]), 1:48) }
#if (ref[1] == "c") { colref <- cbind(rep(seq(as.numeric(ref[2]), 48, 6), each=6), 1:48) }
#rr <- matrix(1:48, ncol=8, byrow=T)
#if (ref[1] == "r") { colref <- cbind(rep(rr[as.numeric(ref[2]),], 8), 1:48) }
#  
#}
# create reference table for 96-well plate
if (dim(OD)[1] == 96) {
if (ref[1] == "w") { colref <- cbind(as.numeric(ref[2]), 1:96) }
if (ref[1] == "c") { colref <- cbind(rep(seq(as.numeric(ref[2]), 96, 12), each=12), 1:96) }
rr <- matrix(1:96, ncol=12, byrow=T)
if (ref[1] == "r") { colref <- cbind(rep(rr[as.numeric(ref[2]),], 8), 1:96) }
}
# create reference table for 384-well plate
if (dim(OD)[1] == 384) {
if (ref[1] == "w") { colref <- cbind(as.numeric(ref[2]), 1:384) }
}
 
#############################
# Create growth curves here
#############################
 
# calculate the Area Under Curve (AUC) for the OD reads in the matrix
od.sum <- round(apply(OD, 1, sum), digit=4)
 
# prepare an object for RelG
relg <- c()
 
# plot of the curves (48-well plate)
if (dim(OD)[1] == 48) {
pdf(paste("Tecan.Plate", exp_name, exp_date, "pdf", sep="."), width=6.5, height=8)
par(mar=c(0,0,0,0))
par(oma=c(3,3,0,1))
# t(matrix(1:48, ncol=8, byrow=T))[,c(6:1)]
layout(matrix(c(rep(49,6), as.vector(t(t(matrix(1:48, ncol=8, byrow=T))[,c(6:1)]))), ncol=6, byrow=T))
for (i in 1:48) {
y.upper.limit <- max(OD)*1.1
len <- dim(OD)[2]
plot(OD[i,seq(1,len,1)], pch = 1, axes = F, ylim = c(-0.1, y.upper.limit), lab=c(3,3,3), type="n")
if (color == "Y") {
col.ratio <- cut(od.sum[colref[i,2]]/od.sum[colref[i,1]], breaks=c(-100, 0.2, 0.4, 0.6, 0.8, 0.9, 1.1, 2, 3, 5, 10, 100), labels=F)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = gradient[col.ratio])
} # color boxes
lines(OD[c(colref[i,2]),seq(1,len,1)], pch = 1, col = 1, lwd=2) # black line
lines(OD[c(colref[i,1]),seq(1,len,1)], pch = 1, col = 2, lwd=2) # Red reference lines
 
#mtext(drugs[i], side = 3, line = -1.5, adj = 0.1, cex = 0.5)
#mtext(strain[i], side = 3, line = -2.5, adj = 0.1, cex = 0.5)
atx <- axTicks(1)
aty <- axTicks(2)
if (i == 48) { axis(2, las=1, cex.axis=1); axis(1, at=atx, labels=atx/4, las=1, cex.axis=1); }

if (i %in% unique(colref[,1])) { box(lwd = 3) } else { box() }
relg <- c(relg, round((od.sum[colref[i,2]]-od.sum[colref[i,1]]) / od.sum[colref[i,1]], digit=4))
}
plot(x=1:10, y=seq(0.1,1,0.1),type="n", axes=F)
text(1,0.7,paste("ID: ", plate, "     Name: ", exp_name, "     User: ", user, sep=""), cex=1.7, adj=c(0,0))
text(1,0.3,paste("Started on ", paste(exp_date, " at ", exp_time, paste=""), " on ", reader, ".          Plotted all ", run.length," reads.", sep=""), cex=1.3, adj=c(0,0))
dev.off()
}
 
# plot of the curves (96-well plate)
if (dim(OD)[1] == 96) {
pdf(paste("Tecan.Plate", exp_name, exp_date, "pdf", sep="."), width=11, height=8.5)
par(mar=c(0,0,0,0))
par(oma=c(3,3,0,1))
layout(matrix(c(rep(97,12), 1:96), ncol=12, byrow=T))
for (i in 1:96) {
y.upper.limit <- max(OD)*1.1
len <- dim(OD)[2]
plot(OD[i,seq(1,len,1)], pch = 1, axes = F, ylim = c(-0.1, y.upper.limit), lab=c(3,3,3), type="n")
if (color == "Y") {
col.ratio <- cut(od.sum[colref[i,2]]/od.sum[colref[i,1]], breaks=c(-100, 0.2, 0.4, 0.6, 0.8, 0.9, 1.1, 2, 3, 5, 10, 100), labels=F)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = gradient[col.ratio])
}

lines(OD[c(colref[i,2]),seq(1,len,1)], pch = 1, col = 1, lwd=2)
lines(OD[c(colref[i,1]),seq(1,len,1)], pch = 1, col = 2, lwd=2)
 
mtext(coo.96[i], side = 3, line = -1.5, adj = 0.1, cex = 0.7)
#mtext(drugs[i], side = 3, line = -1.5, adj = 0.1, cex = 0.5)
#mtext(strain[i], side = 3, line = -2.5, adj = 0.1, cex = 0.5)
atx <- axTicks(1)
aty <- axTicks(2)
if (i == 85) { axis(2, las=1, cex.axis=1.2); axis(1, at=atx, labels=atx/4, cex.axis=1.2) }
if (i %in% unique(colref[,1])) { box(lwd = 3) } else { box() }
relg <- c(relg, round((od.sum[colref[i,2]]-od.sum[colref[i,1]]) / od.sum[colref[i,1]], digit=4))
}
plot(x=1:10, y=seq(0.1,1,0.1),type="n", axes=F)
text(1,0.7,paste("ID: ", plate, "     Name: ", exp_name, "     User: ", user, sep=""), cex=1.7, adj=c(0,0))
text(1,0.3,paste("Started on ", paste(exp_date, " at ", exp_time, paste=""), " on ", reader, ".          Plotted all ", run.length," reads.", sep=""), cex=1.3, adj=c(0,0))
dev.off()
}
 
# plot of the curves (384-well plate)
if (dim(OD)[1] == 384) {
pdf(paste("Tecan.Plate", exp_name, exp_date, "pdf", sep="."), width=11, height=8.5)
par(mar=c(0,0,0,0))
par(oma=c(3,3,0,1))
layout(matrix(c(rep(385,24), 1:384), ncol=24, byrow=T))
for (i in 1:384) {
y.upper.limit <- max(OD)*1.1
len <- dim(OD)[2]
plot(OD[i,seq(1,len,1)], pch = 1, axes = F, ylim = c(-0.1, y.upper.limit), lab=c(3,3,3), type="n")
if (color == "Y") {
col.ratio <- cut(od.sum[colref[i,2]]/od.sum[colref[i,1]], breaks=c(-100, 0.2, 0.4, 0.6, 0.8, 0.9, 1.1, 2, 3, 5, 10, 100), labels=F)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = gradient[col.ratio])
}
lines(OD[c(colref[i,2]),seq(1,len,1)], pch = 1, col = 1, lwd=2)
lines(OD[c(colref[i,1]),seq(1,len,1)], pch = 1, col = 2, lwd=2)
 
mtext(coo.384[i], side = 3, line = -1, adj = 0.1, cex = 0.4)
atx <- axTicks(1)
aty <- axTicks(2)
if (i == 361) { axis(2, las=1, cex.axis=1); axis(1, at=atx, labels=atx/4, las=1, cex.axis=1); }
if (i == as.numeric(ref[2])) { box(lwd = 3) } else { box() }
}
plot(x=1:10, y=seq(0.1,1,0.1),type="n", axes=F)
text(1,0.7,paste("ID: ", plate, "     Name: ", exp_name, "     User: ", user, sep=""), cex=1.7, adj=c(0,0))
text(1,0.3,paste("Started on ", paste(exp_date, " at ", exp_time, paste=""), " on ", reader, ".          Plotted all ", run.length," reads.", sep=""), cex=1.3, adj=c(0,0))
dev.off()
}
}

