# #oldwd=getwd()
# setwd('~/Downloads/TPS_code')
# a=read.table("input.data.txt", sep='\t')
# shuffled=sample(length(a[,1])-1)+1
# barriers=round(seq(1,length(shuffled),length.out=11))
# 
# sapply(c(2:length(barriers)), function(i){
#     mat=c(as.character(a[1,]),
#           as.character(a[shuffled[barriers[i-1]:barriers[i]],]))
#     print(mat)
#     write.table(mat, file=paste('tenthTPS', (i-1), sep=''), quote=F, row.names = F, col.names = F)
# })
# i=1
# mat=read.table(paste('tenthTPS', i, sep=''), sep=',', header=T)
# tp=as.numeric(sapply(as.character(colnames(mat)), function(j){substring(j, 2, nchar(j))}))
# 
# b_TPSData=lapply(c(1:10), function(i){
#     mat=read.table(paste('tenthTPS', i, sep=''), sep=',', header=T)
#     tp=as.numeric(sapply(as.character(colnames(mat)), function(j){substring(j, 2, nchar(j))}))
#     print(tp)
#     print(as.matrix(mat))
#     findPath(tp, rep(0,length(colnames(mat))),t(as.matrix(mat)/100000), 8, 1, multiple=F, type=1, numPerts=40, resampleTraining = T)
# 
# })
# save(b_TPSData, file='b_TPSData.RData')
# 
# ##read in TPS results:
# tps=apply(read.table('timepointsTPS_8', sep=','), c(1,2), function(i){as.numeric(i)})
# tps_asB=apply(tps, 1, function(i){
#     c(1, which(tp %in% i)+1, length(tp))
# })
# 
# #collapse NITPicker table
# nit_asB=sapply(b_TPSData, function(i){i})
# 
# #read in data to use as test set
# splitTable=lapply(c(1:10), function(i){
#     read.table(paste('tenthTPS', i, sep=''), header=T, sep=',')
# })
# #combine into test sets (i.e. all but i)
# combineTestSets=lapply(c(1:10), function(i){
#        index=c(1:10)[-i]
#        print(index)
#        l=sapply(index, function(j){
#            print(dim(splitTable[[j]]))
#            as.matrix(splitTable[[j]])
#        })
#        do.call("rbind", l)
# })
# 
# #calculate L2 score for test set:
# nitL2=sapply(c(1:10), function(i){
# 
#     inds=nit_asB[,i]
# 
#     testMat=combineTestSets[[i]]/10000
# 
# sum(apply(testMat, 1, function(gene){
#     L2(c(tp[1], tp, tp[length(tp)]),
#        as.numeric(c(gene[inds[1]],gene, gene[inds[length(inds)]])),
#        rep(0, length(gene)+2), min(tp), max(tp), c(1, 1+inds, length(tp)+2))
# 
# }))*10000
# 
# })
# 
# 
# tpsL2=sapply(c(1:10), function(i){
# 
#     inds=tps_asB[,i]
# 
#     testMat=combineTestSets[[i]]/10000
# 
#     sum(apply(testMat, 1, function(gene){
#         L2(c(tp[1], tp, tp[length(tp)]),
#            as.numeric(c(gene[inds[1]],gene, gene[inds[length(inds)]])),
#            rep(0, length(gene)+2), min(tp), max(tp), c(1, 1+inds, length(tp)+2))
# 
#     }))*10000
# 
# })
# 
# 
# #draw gene expression curves and corresponding perturbations
# 
# pdf(paste("Figure_ComparisonToTPS.pdf", sep=""), width=9, height=4,pointsize = 14)
# par(mfrow=c(1,2));
# par(mar=c(4, 4, 1, 1)+0.1);
# par(cex=0.9)
# 
# 
# 
# i=1
# mat=read.table(paste('tenthTPS', i, sep=''), sep=',', header=T)
# tp=as.numeric(sapply(as.character(colnames(mat)), function(j){substring(j, 2, nchar(j))}))
# 
# 
# set.seed(123)
# perts=generatePerturbations(t(mat), tp, spline=1, numPert=1000)
# plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(mat)-4, max(mat)+2), ylab='gene expression', xlab='time')
# apply(perts$ft, 2, function(i){lines(perts$time, i, col=rgb(0.5, 0.5, 0.5, 0.2))})
# apply(mat, 1, function(i){lines(tp, i)})
# legend(16, -7, c('data', 'sampled'), lty=1, col=c('black', 'grey'), bty='n')
# 
# pval=paste('p-value <', round(t.test(nitL2, tpsL2)$p.value, 5))
# boxplot(list(nitL2, tpsL2), names=c("NITPicker", "TPS"), ylab="L2", xlab='method for time point selection', main=pval)
# points(1+jitter(rep(0, length(nitL2)), amount=0.1), nitL2, pch=19, col='orange')
# points(2+jitter(rep(0, length(nitL2)), amount=0.1), tpsL2, pch=19, col='orange')
# dev.off()
