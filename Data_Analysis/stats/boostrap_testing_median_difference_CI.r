library(boot)
library(simpleboot)
d1 <- read.csv(file=paste("Merged_Rate_dependence_1.csv", sep=''), header=TRUE, sep=",")
d3 <- read.csv(file=paste("Merged_Rate_dependence_3.csv", sep=''), header=TRUE, sep=",")
cols = c('"IK2p"', 'ISK', 'IKur', 'IK2p_ISK', 'IK2p_IKur', 'ISK_IKur','IK2P_ISK_IKur')

cols <- c("IK2p", "ISK","IKur","IK2p_ISK","IK2p_IKur","ISK_IKur","IK2P_ISK_IKur");



model = "CZ %"
filename = "RD_figure_CI.dat"
col=cols[[1]]
for (col in cols) {
res=two.boot(na.omit(d3[,col])*100, na.omit(d1[,col])*100, median, 10000)
print(col)
# bootcorr=boot.ci(boot.out = res, type = c("norm", "basic"), conf=1-0.05) #, "perc", "bca"
ci=boot.ci(boot.out = res, type = c("norm", "basic"), conf=1-0.05/7) #, "perc", "bca"
cat (model, file=filename, append=TRUE, sep=" ")
cat(' ',file=filename,append=TRUE)
cat (ci$t0, file=filename, append=TRUE, sep=" ")
cat(' ',file=filename,append=TRUE)
cat (ci$basic, file=filename, append=TRUE, sep=" ")
cat(' ',file=filename,append=TRUE)
cat (ci$normal, file=filename, append=TRUE, sep=" ")
cat('\n',file=filename,append=TRUE)
}



model = "CZ %"
filename = "RD_figure_CI_ABS.dat"
col=cols[[1]]
for (col in cols) {
res=two.boot(na.omit((d3[,col]-1)*d3[,'APD_baseline']), na.omit((d1[,col]-1)*d1[,'APD_baseline']), median, 10000)
print(col)
# bootcorr=boot.ci(boot.out = res, type = c("norm", "basic"), conf=1-0.05) #, "perc", "bca"
ci=boot.ci(boot.out = res, type = c("norm", "basic"), conf=1-0.05/7) #, "perc", "bca"
cat (model, file=filename, append=TRUE, sep=" ")
cat(' ',file=filename,append=TRUE)
cat (ci$t0, file=filename, append=TRUE, sep=" ")
cat(' ',file=filename,append=TRUE)
cat (ci$basic, file=filename, append=TRUE, sep=" ")
cat(' ',file=filename,append=TRUE)
cat (ci$normal, file=filename, append=TRUE, sep=" ")
cat('\n',file=filename,append=TRUE)
}



