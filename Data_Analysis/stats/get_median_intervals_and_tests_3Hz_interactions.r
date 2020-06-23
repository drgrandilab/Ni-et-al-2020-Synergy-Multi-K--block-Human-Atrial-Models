for (i in 0:11) {
model = toString(i)

# Data <-read.table(paste("Ratio_data_", model, ".csv", sep=''))
d1 <- read.csv(file=paste("Rate_dependence_3_", model, ".csv", sep=''), header=TRUE, sep=",")

# print (summary(d1))
# d2 <- read.csv(file=paste("../../100_percent_1000ms_1000ms/3Hz/Ratio_data_", model, ".csv", sep=''), header=TRUE, sep=",")
# d2 <- read.csv(file=paste("Rate_dependence_3_", model, ".csv", sep=''), header=TRUE, sep=",")

# print (summary(d2))


IK2p_ISK = d1['IK2p_ISK'] + d1['Baseline'] - d1['IK2p'] -d1['ISK']
IK2p_IKur = d1['IK2p_IKur'] + d1['Baseline'] - d1['IK2p'] -d1['IKur']
ISK_IKur = d1['ISK_IKur'] + d1['Baseline'] - d1['ISK'] -d1['IKur']
IK2P_ISK_IKur = (d1['IK2P_ISK_IKur'] + 2*d1['Baseline'] - d1['ISK'] -d1['IKur'] -d1['IK2p']) #/ ( d1['ISK'] +d1['IKur'] +d1['IK2p']-3)


a<-na.omit(IK2p_ISK)
b<-na.omit(IK2p_IKur)
c<-na.omit(ISK_IKur)
d<-na.omit(IK2P_ISK_IKur)


output_res <- function(model, test1,filename){
	
	cat (model, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (test1$statistic, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (test1$p.value, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (test1$conf.int, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (test1$estimate, file=filename, append=TRUE, sep=" ")
	cat('\n',file=filename,append=TRUE)

}


test1=wilcox.test(a$IK2p_ISK,conf.level=0.9875,conf.int = TRUE)
print(test1)
file<-"out_IK2p_ISK_wilcox.dat.3Hz"
output_res(model, test1, file)

test2=wilcox.test(b$IK2p_IKur,conf.level=0.9875,conf.int = TRUE)
print(test2)
file<-"out_IK2p_IKur_wilcox.dat.3Hz"
output_res(model, test2, file)


test2=wilcox.test(c$ISK_IKur,conf.level=0.9875,conf.int = TRUE)
print(test2)
file<-"out_ISK_IKur_wilcox.dat.3Hz"
output_res(model, test2, file)


test2=wilcox.test(d$IK2P_ISK_IKur,conf.level=0.9875,conf.int = TRUE)
print(test2)
file<-"out_IK2P_ISK_IKur_wilcox.dat.3Hz"
output_res(model, test2, file)

# print (test1$p.value)



# print (wilcox.test(a$IK2p_ISK,conf.int = TRUE))
# print (wilcox.test(b$IK2p_IKur,conf.int = TRUE))
# print (median(b$ISK_IKur))
# print(boxplot.stats(b$IK2P_ISK_IKur))




### bootstrap analysis here

fc <- function(d, i){
	d2 <- d[i]
	median(d2);  
}
library(boot)
set.seed(626)  # set seed in order to be reproducible




bootstrap_model <- function(model, data, filename) {
	bootcorr <- boot(data, fc, R=10000)
	print(bootcorr)
	ci=boot.ci(boot.out = bootcorr, type = c("norm", "basic", "perc", "bca"), conf=1-0.05/4)
	cat (model, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (ci$t0, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (ci$basic, file=filename, append=TRUE, sep=" ")
	cat(' ',file=filename,append=TRUE)
	cat (ci$normal, file=filename, append=TRUE, sep=" ")
	cat('\n',file=filename,append=TRUE)
	# print(ci)
	# print(model)
}

bootstrap_model(model, a$IK2p_ISK, 'out_boot_IK2p_ISK.dat.3Hz')
bootstrap_model(model, b$IK2p_IKur, 'out_boot_IK2p_IKur.dat.3Hz')
bootstrap_model(model, c$ISK_IKur, 'out_boot_ISK_IKur.dat.3Hz')
bootstrap_model(model, d$IK2P_ISK_IKur, 'out_boot_IK2P_ISK_IKur.dat.3Hz')

# # IK2p_IKur = (d1['IK2p_IKur'] -1) / ( d1['IK2p'] +d1['IKur'] - 2*d1['Baseline'])
# # IK2p_IKur_500 = (d2['IK2p_IKur'] -1) / ( d2['IK2p'] +d2['IKur'] - 2*d2['Baseline'])

# # IK2P_ISK_IKur = (d1['IK2P_ISK_IKur']) / (  d1['ISK'] +d1['IKur'] +d1['IK2p']- 2*d1['Baseline']) 
# # IK2P_ISK_IKur_500 = (d2['IK2P_ISK_IKur']) / (  d2['ISK'] +d2['IKur'] +d2['IK2p']- 2*d2['Baseline']) 

# # a<-na.omit(IK2P_ISK_IKur)
# # b<-na.omit(IK2P_ISK_IKur_500)

# a<-na.omit(IK2p_IKur)
# b<-na.omit(IK2p_IKur_500)
# # boxplot(c(IK2p_ISK, IK2p_IKur, ISK_IKur, IK2P_ISK_IKur))
# boxplot(c(a,b ))
# # print(boxplot.stats(IK2P_ISK_IKur))
# # print(boxplot.stats(IK2P_ISK_IKur_500))
# # library(fitdistrplus)



# # # descdist(a$IK2P_ISK_IKur, discrete=FALSE)
# # # descdist(a$IK2P_ISK_IKur, discrete=FALSE)
# # # summary statistics
# # # ------
# # # min:  -0.04439865   max:  0.1336406 
# # # median:  0.007239819 
# # # mean:  0.01363605 
# # # estimated sd:  0.02165565 
# # # estimated skewness:  1.276983 
# # # estimated kurtosis:  5.662521 
# # # hist()
# # # Error in hist.default() : argument "x" is missing, with no default
# # # hist(a$IK2P_ISK_IKur)
# # # hist(a$IK2P_ISK_IKur,10)
# # # hist(a$IK2P_ISK_IKur,100)
# # # normal_dist <- fitdist(a$IK2P_ISK_IKur, "norm")
# # # plot(normal_dist)


# test1=wilcox.test(a$IK2p_IKur,conf.level=0.9875,conf.int = TRUE)
# test2=wilcox.test(b$IK2p_IKur,conf.level=0.9875,conf.int = TRUE)






# print (median(b$IK2p_IKur))
# print(boxplot.stats(b$IK2p_IKur))


# fc <- function(d, i){
# 	d2 <- d[i]
# 	median(d2);  #get median instead
# }
# library(boot)
# set.seed(626)
# bootcorr <- boot(a$IK2p_IKur, fc, R=10000)
# print(bootcorr)
# ci=boot.ci(boot.out = bootcorr, type = c("norm", "basic", "perc", "bca"), conf=1-0.05/4)
# print(ci)
# print(model)


# cat (summary(G, test=adjusted("bonferroni"))$test$pvalues, file="out_pvalues.dat", append=TRUE, sep="\t")

print (">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

}