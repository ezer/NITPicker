#oldwd=getwd()

setwd('NITPicker')

library(vioplot)
scores=c()
namesScores=c()
set.seed(123)
dataGenes=read.table("input.data.txt", sep=',')
tp=as.numeric(sapply(as.character(colnames(dataGenes)), function(j){substring(j, 2, nchar(j))}))
#     
# split into training and testing sets 50 sets
# sapply(c(1:50), function(i){
#     trainIDs=sample.int(dim(dataGenes)[1], ceiling(0.5*dim(dataGenes)[1]))
#     testIDs=which(! (c(1:dim(dataGenes)[1]) %in% trainIDs ))
#     tp=as.numeric(sapply(as.character(colnames(dataGenes)), function(j){substring(j, 2, nchar(j))}))
#     
#    # rownames(dataGenes)=paste("A", rownames(dataGenes), sep="")
#     write.table(rbind(tp, dataGenes[trainIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("tenth_tpsTrain_", i, ".txt", sep=""))
#     write.table(rbind(tp, dataGenes[testIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("tenth_tpsTest_", i, ".txt", sep=""))
# })
#output_tenthTPS.csv

#read in test set
readOrder=read.table("output_fileNames_tenth_final.txt")
testSet_tenth=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][3]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
    print(number)
    name=paste('tenth_tpsTest_', number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})

trainSet_tenth=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][3]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
    print(number)
    name=paste('tenth_tpsTrain_', number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})


#run NITPicker on this
set.seed(123)
b_NITData=sapply(c(1:50), function(i){
    
    
    
    b=findPathF1(tp, t(trainSet_tenth[[i]]/10), 8, spline=1, numPerts=10000, resampleTraining = T, fast=T, iter=40)
    #dataGenes=generatePerturbations(trainSet_elf4[[i]], tp)$f0#
    dataGenes=trainSet_tenth[[i]]
    plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
    apply(dataGenes, 1, function(j){
        lines(tp, j)
    })
    print(b)
    print(tp[b])
    abline(v=tp[b])
    tp[b]
})

b_TPSData=read.csv("output_tenthTPS.csv", header=F)
b_TPSData=apply(b_TPSData, 1, function(i){sort(i)})

#get random
set.seed(123)
b_randomData=sapply(c(1:50), function(i){
    sort(sample(tp, 8))
})
#get even
b_evenData=sapply(c(1:50), function(i){
    c(0.5, 4.0, 8.0, 12.0, 16, 19, 23, 27)
})

#test sets


l2_TPSData= sapply(c(1:50), function(j){
    i=b_TPSData[,j]
    sum(apply(testSet_tenth[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_NITData= sapply(c(1:50), function(j){
    i=b_NITData[,j]
    sum(apply(testSet_tenth[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    

l2_random= sapply(c(1:50), function(j){
    i=b_randomData[,j]
    sum(apply(testSet_tenth[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})  

l2_even= sapply(c(1:50), function(j){
    i=b_evenData[,j]
    sum(apply(testSet_tenth[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})  
pdf(paste("Figure_testTPS_tenth.pdf", sep=""), width=8, height=5,pointsize = 12)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

vioplot(log(l2_NITData), log(l2_TPSData), log(l2_even), log(l2_random), col="grey", names=c('NITPicker', 'TPS', 'Even', 'Random'))
points(1+jitter(rep(0, length(l2_NITData)), amount=0.2), log(l2_NITData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(2+jitter(rep(0, length(l2_TPSData)), amount=0.2), log(l2_TPSData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(3+jitter(rep(0, length(l2_even)), amount=0.2), log(l2_even), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(4+jitter(rep(0, length(l2_random)), amount=0.2), log(l2_random), pch=19, col=rgb(0, 0, 0.5, 0.5))
dev.off()
namesScores=c(namesScores, 'tenth')
scores=c(scores, t.test(l2_NITData, l2_TPSData)$p.value)
save(l2_NITData, l2_TPSData, l2_even, l2_random, file='tenth.RData')

###########Also compare it to something that is a little less random...
geneExpr=read.csv("col0GeneExpression.csv", header=T)
genesOfInterest=read.table("elf4Genes.txt")
dataGenes=geneExpr[which(as.character(geneExpr[,1]) %in% as.character(genesOfInterest[,1])),4:19]
tp=sapply(colnames(geneExpr)[4:19], function(i){
    a=substring(i, 2, nchar(i))
    a=gsub('\\.', '-', a)
    as.numeric(a)
})
#plot these
plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
apply(dataGenes, 1, function(i){
    lines(tp, i)
})

#split into training and testing sets 50 sets
# sapply(c(1:50), function(i){
#     trainIDs=sample.int(dim(dataGenes)[1], ceiling(0.5*dim(dataGenes)[1]))
#     testIDs=which(! (c(1:dim(dataGenes)[1]) %in% trainIDs ))
#     rownames(dataGenes)=paste("A", rownames(dataGenes), sep="")
#     write.table(rbind(tp, dataGenes[trainIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elf4_tpsTrain_", i, ".txt", sep=""))
#     write.table(rbind(tp, dataGenes[testIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elf4_tpsTest_", i, ".txt", sep=""))
# })

#read in test set
readOrder=read.table("tps_elf4_fileOrder.txt")
testSet_elf4=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][3]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
    print(number)
    name=paste('elf4_tpsTest_', number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})

trainSet_elf4=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][3]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
    print(number)
    name=paste('elf4_tpsTrain_', number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})


#run NITPicker on this
set.seed(123)
b_NITData=sapply(c(1:50), function(i){
    
    
    
    b=findPathF1(4+tp, t(trainSet_elf4[[i]]/10), 8, spline=1, numPerts=10000, resampleTraining = T, fast=T, iter=40)
    #dataGenes=generatePerturbations(trainSet_elf4[[i]], tp)$f0#
    dataGenes=trainSet_elf4[[i]]
    plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
    apply(dataGenes, 1, function(j){
        lines(tp, j)
    })
    print(b)
    print(tp[b])
    abline(v=tp[b])
    tp[b]
})

b_TPSData=read.csv("tps_elf4_output.csv", header=F)
b_TPSData=apply(b_TPSData, 1, function(i){sort(i)})

#test sets
  

l2_TPSData= sapply(c(1:50), function(j){
    i=b_TPSData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_NITData= sapply(c(1:50), function(j){
    i=b_NITData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    
#get random
set.seed(123)
b_randomData=sapply(c(1:50), function(i){
    sort(sample(tp, 8))
})
#get even
b_evenData=sapply(c(1:50), function(i){
    c(-4, 1, 8, 16, 22, 28, 32, 40)
})

l2_even= sapply(c(1:50), function(j){
    i=b_evenData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_random= sapply(c(1:50), function(j){
    i=b_randomData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    

pdf(paste("Figure_testTPS_elf4.pdf", sep=""), width=8, height=6,pointsize = 12)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

vioplot(log(l2_NITData), log(l2_TPSData), log(l2_even), log(l2_random), col="grey", names=c('NITPicker', 'TPS', 'Even', 'Random'))
points(1+jitter(rep(0, length(l2_NITData)), amount=0.2), log(l2_NITData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(2+jitter(rep(0, length(l2_TPSData)), amount=0.2), log(l2_TPSData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(3+jitter(rep(0, length(l2_even)), amount=0.2), log(l2_even), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(4+jitter(rep(0, length(l2_random)), amount=0.2), log(l2_random), pch=19, col=rgb(0, 0, 0.5, 0.5))
dev.off()

t.test(l2_NITData, l2_TPSData)
namesScores=c(namesScores, 'elf4')
scores=c(scores, t.test(l2_NITData, l2_TPSData)$p.value)
save(l2_NITData, l2_TPSData, l2_even, l2_random, file='elf4.RData')
#################


#Do one more example: C-elegans

#read in C. elegans table
genes_elegans=read.csv('Celegans_Table_S2.csv', header=T)
tp=sapply(colnames(genes_elegans), function(i){
    a=strsplit(as.character(i), "_")[[1]][2]
    as.numeric(substring(a, 5, nchar(a)))
})

#read in interesting C. elegans genes
interesting_genes=as.character(read.table("CelegansInterestingGenes")[,1])
dataGenes=genes_elegans[interesting_genes,]


plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
apply(dataGenes, 1, function(i){
    lines(tp, i)
})
#set.seed(123)
# #split into training and testing sets 50 sets
# sapply(c(1:50), function(i){
#     trainIDs=sample.int(dim(dataGenes)[1], ceiling(0.5*dim(dataGenes)[1]))
#     testIDs=which(! (c(1:dim(dataGenes)[1]) %in% trainIDs ))
#    # rownames(dataGenes)=paste("A", rownames(dataGenes), sep="")
#     write.table(rbind(tp, dataGenes[trainIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elegans_tpsTrain_", i, ".txt", sep=""))
#     write.table(rbind(tp, dataGenes[testIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elegans_tpsTest_", i, ".txt", sep=""))
# })
#output_elegansTPS.csv
#output_fileNames_elegans_final.txt

#read in test set
readOrder=read.table("output_fileNames_elegans_final.txt")
testSet_elegans=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][3]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
    print(number)
    name=paste('elegans_tpsTest_', number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})

trainSet_elegans=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][3]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
    print(number)
    name=paste('elegans_tpsTrain_', number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})


#run NITPicker on this
set.seed(123)
b_NITData=sapply(c(1:50), function(i){
    
    
    
    b=findPathF1(4+tp, t(trainSet_elegans[[i]]), 8, spline=1, numPerts=10000, resampleTraining = T, fast=T, iter=40)
    #dataGenes=generatePerturbations(trainSet_elf4[[i]], tp)$f0#
    dataGenes=trainSet_elegans[[i]]
    plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
    apply(dataGenes, 1, function(j){
        lines(tp, j)
    })
    print(b)
    print(tp[b])
    abline(v=tp[b])
    tp[b]
})

b_TPSData=read.csv("output_elegansTPS.csv", header=F)
b_TPSData=apply(b_TPSData, 1, function(i){sort(i)})

#test sets


l2_TPSData= sapply(c(1:50), function(j){
    i=b_TPSData[,j]
    sum(apply(testSet_elegans[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_NITData= sapply(c(1:50), function(j){
    i=b_NITData[,j]
    sum(apply(testSet_elegans[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    
#get random
set.seed(123)
b_randomData=sapply(c(1:50), function(i){
    sort(sample(tp, 6))
})
#get even
b_evenData=sapply(c(1:50), function(i){
    c(0, 100, 200, 300, 400, 500, 600)
})

l2_even= sapply(c(1:50), function(j){
    i=b_evenData[,j]
    sum(apply(testSet_elegans[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_random= sapply(c(1:50), function(j){
    i=b_randomData[,j]
    sum(apply(testSet_elegans[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    

pdf(paste("Figure_testTPS_elegans.pdf", sep=""), width=8, height=6,pointsize = 12)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

vioplot(log(l2_NITData), log(l2_TPSData), log(l2_even), log(l2_random), col="grey", names=c('NITPicker', 'TPS', 'Even', 'Random'))
points(1+jitter(rep(0, length(l2_NITData)), amount=0.2), log(l2_NITData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(2+jitter(rep(0, length(l2_TPSData)), amount=0.2), log(l2_TPSData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(3+jitter(rep(0, length(l2_even)), amount=0.2), log(l2_even), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(4+jitter(rep(0, length(l2_random)), amount=0.2), log(l2_random), pch=19, col=rgb(0, 0, 0.5, 0.5))
dev.off()
namesScores=c(namesScores, 'elf4_elegans')
scores=c(scores, t.test(l2_NITData, l2_TPSData)$p.value)
save(l2_NITData, l2_TPSData, l2_even, l2_random, file='elegans.RData')
t.test(l2_NITData, l2_TPSData)

#####I want to do the two patterned ones again, but with smaller training sets
genes_elegans=read.csv('Celegans_Table_S2.csv', header=T)
tp=sapply(colnames(genes_elegans), function(i){
    a=strsplit(as.character(i), "_")[[1]][2]
    as.numeric(substring(a, 5, nchar(a)))
})
interesting_genes=as.character(read.table("CelegansInterestingGenes")[,1])
dataGenes=genes_elegans[interesting_genes,]
set.seed(123)
# #split into training and testing sets 50 sets
# sapply(c(1:50), function(i){
#     trainIDs=sample.int(dim(dataGenes)[1], ceiling(0.33*dim(dataGenes)[1]))
#     testIDs=which(! (c(1:dim(dataGenes)[1]) %in% trainIDs ))
#    # rownames(dataGenes)=paste("A", rownames(dataGenes), sep="")
#     write.table(rbind(tp, dataGenes[trainIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elegans_st_tpsTrain_", i, ".txt", sep=""))
#     write.table(rbind(tp, dataGenes[testIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elegans_st_tpsTest_", i, ".txt", sep=""))
# })
fileCSV="output_elegans_stTPS.csv"
fileNames="output_fileNames_elegans_st_final.txt"
rootTest="elegans_st_tpsTest_"
rootTrain="elegans_st_tpsTrain_"
#read in test set
readOrder=read.table(fileNames)
testSet_elf4=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][4]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][4], "\\.")[[1]][1]
    print(number)
    name=paste(rootTest, number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})

trainSet_elf4=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][4]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][4], "\\.")[[1]][1]
    print(number)
    name=paste(rootTrain, number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})


#run NITPicker on this
set.seed(123)
b_NITData=sapply(c(1:50), function(i){
    
    
    
    b=findPathF1(tp, t(trainSet_elf4[[i]]/10), 8, spline=1, numPerts=10000, resampleTraining = T, fast=T, iter=40)
    #dataGenes=generatePerturbations(trainSet_elf4[[i]], tp)$f0#
    dataGenes=trainSet_elf4[[i]]
    plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
    apply(dataGenes, 1, function(j){
        lines(tp, j)
    })
    print(b)
    print(tp[b])
    abline(v=tp[b])
    tp[b]
})

b_TPSData=read.csv(fileCSV, header=F)
b_TPSData=apply(b_TPSData, 1, function(i){sort(i)})

#test sets


l2_TPSData= sapply(c(1:50), function(j){
    i=b_TPSData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_NITData= sapply(c(1:50), function(j){
    i=b_NITData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    
#get random
set.seed(123)
b_randomData=sapply(c(1:50), function(i){
    sort(sample(tp, 6))
})
#get even
b_evenData=sapply(c(1:50), function(i){
    c(0, 90, 180, 270, 330, 420, 510, 600) 
})

l2_even= sapply(c(1:50), function(j){
    i=b_evenData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    

l2_random= sapply(c(1:50), function(j){
    i=b_randomData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    

pdf(paste("Figure_testTPS_elegans_smallerTrainingSet.pdf", sep=""), width=8, height=6,pointsize = 12)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

vioplot(log(l2_NITData), log(l2_TPSData), log(l2_even), log(l2_random), col="grey", names=c('NITPicker', 'TPS', 'Even', 'Random'))
points(1+jitter(rep(0, length(l2_NITData)), amount=0.2), log(l2_NITData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(2+jitter(rep(0, length(l2_TPSData)), amount=0.2), log(l2_TPSData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(3+jitter(rep(0, length(l2_even)), amount=0.2), log(l2_even), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(4+jitter(rep(0, length(l2_random)), amount=0.2), log(l2_random), pch=19, col=rgb(0, 0, 0.5, 0.5))
dev.off()
namesScores=c(namesScores, 'elegans_smallerTrain')
scores=c(scores, t.test(l2_NITData, l2_TPSData)$p.value)
save(l2_NITData, l2_TPSData, l2_even, l2_random, file='elegans_smaller.RData')





geneExpr=read.csv("col0GeneExpression.csv", header=T)
genesOfInterest=read.table("elf4Genes.txt")
dataGenes=geneExpr[which(as.character(geneExpr[,1]) %in% as.character(genesOfInterest[,1])),4:19]
tp=sapply(colnames(geneExpr)[4:19], function(i){
    a=substring(i, 2, nchar(i))
    a=gsub('\\.', '-', a)
    as.numeric(a)
})
set.seed(123)
#split into training and testing sets 50 sets
# sapply(c(1:50), function(i){
#     trainIDs=sample.int(dim(dataGenes)[1], ceiling(0.33*dim(dataGenes)[1]))
#     testIDs=which(! (c(1:dim(dataGenes)[1]) %in% trainIDs ))
#    # rownames(dataGenes)=paste("A", rownames(dataGenes), sep="")
#     write.table(rbind(tp, dataGenes[trainIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elf4_st_tpsTrain_", i, ".txt", sep=""))
#     write.table(rbind(tp, dataGenes[testIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elf4s_st_tpsTest_", i, ".txt", sep=""))
# })
fileCSV="output_elf4_stTPS.csv"
fileNames="output_fileNames_elf4_st_final.txt"
rootTest="elf4s_st_tpsTest_"
rootTrain="elf4_st_tpsTrain_"
#read in test set
readOrder=read.table(fileNames)
testSet_elf4=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][4]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][4], "\\.")[[1]][1]
    print(number)
    name=paste(rootTest, number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})

trainSet_elf4=lapply(readOrder[,1], function(i){
    i=as.character(i)
    print(i)
    number=strsplit(i, "_")[[1]][4]
    print(number)
    
    number=strsplit(strsplit(i, "_")[[1]][4], "\\.")[[1]][1]
    print(number)
    name=paste(rootTrain, number, '.txt', sep="")
    print(name)
    a=read.csv(name, header=T)
    rownames(a)=a[,1]
    a=a[,2:dim(a)[2]]
    apply(a, c(1,2), function(j){as.numeric(as.character(j))})
})


#run NITPicker on this
set.seed(123)
b_NITData=sapply(c(1:50), function(i){
    
    
    
    b=findPathF1(tp, t(trainSet_elf4[[i]]/10), 8, spline=1, numPerts=10000, resampleTraining = T, fast=T, iter=40)
    #dataGenes=generatePerturbations(trainSet_elf4[[i]], tp)$f0#
    dataGenes=trainSet_elf4[[i]]
    plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
    apply(dataGenes, 1, function(j){
        lines(tp, j)
    })
    print(b)
    print(tp[b])
    abline(v=tp[b])
    tp[b]
})

b_TPSData=read.csv(fileCSV, header=F)
b_TPSData=apply(b_TPSData, 1, function(i){sort(i)})

#test sets


l2_TPSData= sapply(c(1:50), function(j){
    i=b_TPSData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_NITData= sapply(c(1:50), function(j){
    i=b_NITData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    
#get random
set.seed(123)
b_randomData=sapply(c(1:50), function(i){
    sort(sample(tp, 8))
})
#get even
b_evenData=sapply(c(1:50), function(i){
    c(-4, 1, 8, 16, 22, 28, 32, 40)
    #-4
})

l2_even= sapply(c(1:50), function(j){
    i=b_evenData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_random= sapply(c(1:50), function(j){
    i=b_randomData[,j]
    sum(apply(testSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    

pdf(paste("Figure_testTPS_elf4_smallerTrainingSet.pdf", sep=""), width=8, height=6,pointsize = 12)
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=1)

vioplot(log(l2_NITData), log(l2_TPSData), log(l2_even), log(l2_random), col="grey", names=c('NITPicker', 'TPS', 'Even', 'Random'))
points(1+jitter(rep(0, length(l2_NITData)), amount=0.2), log(l2_NITData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(2+jitter(rep(0, length(l2_TPSData)), amount=0.2), log(l2_TPSData), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(3+jitter(rep(0, length(l2_even)), amount=0.2), log(l2_even), pch=19, col=rgb(0, 0, 0.5, 0.5))
points(4+jitter(rep(0, length(l2_random)), amount=0.2), log(l2_random), pch=19, col=rgb(0, 0, 0.5, 0.5))
dev.off()
namesScores=c(namesScores, 'elf4_smallerTrain')
scores=c(scores, t.test(l2_NITData, l2_TPSData)$p.value)
save(l2_NITData, l2_TPSData, l2_even, l2_random, file='elf4_smaller.RData')



















# 
# set.seed(123)
# b_NITData=lapply(c(1:10), function(i){
#     mat=read.table(paste('tenthTPS', i, sep=''), sep=',', header=T)
#    
#     tp=as.numeric(sapply(as.character(colnames(mat)), function(j){substring(j, 2, nchar(j))}))
#     plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(as.matrix(mat)), max(as.matrix(mat))))
#     
#     apply(as.matrix(mat), 1, function(temp){lines(tp, temp)})
#     print(tp)
#     print(as.matrix(mat))
#     b=findPathF1(tp, t(as.matrix(mat)), 8, spline=3, numPerts=10000, resampleTraining = T, fast=T, iter=100)
#     write.table(b, file=paste('TPS_compare_', i, '.txt', sep=''))
#     abline(v=tp[b])
#     b
# })
# 
# 
# 
# 
# #save(b_TPSData, file='b_TPSData.RData')
# 
# ##read in TPS results:
# tps=read.csv('tps_tenth_output_inorder.csv', header=F)
# tps_asB=apply(tps, 1, function(i){
#     which(tp %in% i)
# })
# 
# #collapse NITPicker table
# nit_asB=sapply(b_NITData, function(i){i})
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
#     inds=as.numeric(nit_asB[,i])
# 
#     testMat=combineTestSets[[i]]
# 
# sum(apply(testMat, 1, function(gene){
#    print(inds)
#    L2(tp,
#       gene,
#       rep(0, length(gene)), min(tp), max(tp), inds)
# 
# }))
# 
# })
# 
# 
# tpsL2=sapply(c(1:10), function(i){
# print(i)
#     inds=tps_asB[,i]
# print(inds)
# print(tp[inds])
# print(length(tp))
#     testMat=combineTestSets[[i]]
# 
#     sum(apply(testMat, 1, function(gene){
#         
#         L2(tp,
#            gene,
#            rep(0, length(gene)), min(tp), max(tp), inds)
#         
#     }))
# 
# })
# 
# 
# #draw gene expression curves and corresponding perturbations
# 
# pdf(paste("Figure_ComparisonToTPS6.pdf", sep=""), width=9, height=4,pointsize = 14)
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
# points(2+jitter(rep(0, length(tpsL2)), amount=0.1), tpsL2, pch=19, col='orange')
# dev.off()
# 
# 
# 
# 
# 
# 
# ###########Also compare it to something that is a little less random...
# geneExpr=read.csv("col0GeneExpression.csv", header=T)
# genesOfInterest=read.table("elf4Genes.txt")
# dataGenes=geneExpr[which(as.character(geneExpr[,1]) %in% as.character(genesOfInterest[,1])),4:19]
# tp=sapply(colnames(geneExpr)[4:19], function(i){
#     a=substring(i, 2, nchar(i))
#     a=gsub('\\.', '-', a)
#     as.numeric(a)
# })
# #plot these
# plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
# apply(dataGenes, 1, function(i){
#     lines(tp, i)
# })
# 
# #split into training and testing sets 50 sets
# # sapply(c(1:50), function(i){
# #     trainIDs=sample.int(dim(dataGenes)[1], ceiling(0.5*dim(dataGenes)[1]))
# #     testIDs=which(! (c(1:dim(dataGenes)[1]) %in% trainIDs ))
# #     rownames(dataGenes)=paste("A", rownames(dataGenes), sep="")
# #     write.table(rbind(tp, dataGenes[trainIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elf4_tpsTrain_", i, ".txt", sep=""))
# #     write.table(rbind(tp, dataGenes[testIDs, ]), quote=F, col.names = FALSE, sep=",", file=paste("elf4_tpsTest_", i, ".txt", sep=""))
# # })
# 
# #read in test set
# readOrder=read.table("tps_elf4_fileOrder.txt")
# testSet_elf4=lapply(readOrder[,1], function(i){
#    i=as.character(i)
#    print(i)
#    number=strsplit(i, "_")[[1]][3]
#    print(number)
#    
#     number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
#     print(number)
#     name=paste('elf4_tpsTest_', number, '.txt', sep="")
#     print(name)
#     a=read.csv(name, header=T)
#     rownames(a)=a[,1]
#     a=a[,2:dim(a)[2]]
#     apply(a, c(1,2), function(j){as.numeric(as.character(j))})
# })
# 
# trainSet_elf4=lapply(readOrder[,1], function(i){
#     i=as.character(i)
#     print(i)
#     number=strsplit(i, "_")[[1]][3]
#     print(number)
#     
#     number=strsplit(strsplit(i, "_")[[1]][3], "\\.")[[1]][1]
#     print(number)
#     name=paste('elf4_tpsTrain_', number, '.txt', sep="")
#     print(name)
#     a=read.csv(name, header=T)
#     rownames(a)=a[,1]
#     a=a[,2:dim(a)[2]]
#     apply(a, c(1,2), function(j){as.numeric(as.character(j))})
# })
# 
# 
# #run NITPicker on this
# set.seed(123)
# b_NITData=sapply(c(1:50), function(i){
#     
#     
#     
#     b=findPathF1(4+tp, t(trainSet_elf4[[i]]/10), 8, spline=3, numPerts=10000, resampleTraining = T, fast=T)
#     #dataGenes=generatePerturbations(trainSet_elf4[[i]], tp)$f0#
#     dataGenes=trainSet_elf4[[i]]
#     plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
#     apply(dataGenes, 1, function(j){
#         lines(tp, j)
#     })
#     print(b)
#     print(tp[b])
#     abline(v=tp[b])
#     tp[b]
# })
# 
# b_TPSData=read.csv("tps_elf4_output.csv", header=F)
# b_TPSData=apply(b_TPSData, 1, function(i){sort(i)})
# 
# #test sets
# 
# 
# l2_TPSData= sapply(c(1:50), function(j){
#     i=b_TPSData[,j]
#     sum(apply(testSet_elf4[[j]], 1, function(k){
#         L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
#     }))
#     
# })    
#     
#     
#     
#     
# l2_NITData= sapply(c(1:50), function(j){
#     i=b_NITData[,j]
#     sum(apply(testSet_elf4[[j]], 1, function(k){
#         L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
#     }))
#     
# })    
# 
# vioplot(l2_TPSData, l2_NITData, col="grey")
# 
# sapply(c(1:length(l2_NITData)), function(i){
#     
#         if(l2_TPSData[i]<l2_NITData[i]){
#             print(paste("TPS better", i))
#             dataGenes=trainSet_elf4[[i]]
#             plot(c(), ylim=c(min(dataGenes), max(dataGenes)), xlim=c(min(tp), max(tp)))
#             apply(dataGenes, 1, function(j){
#                 lines(tp, j)
#             })
# 
#             abline(v=b_NITData[,i], lwd=5)
#             abline(v=b_TPSData[,i], col='red')
#         }
#    
# })
# 
# 
#    



l2_TPSData= sapply(c(1:50), function(j){
    i=b_TPSData[,j]
    sum(apply(trainSet_elf4[[j]], 1, function(k){
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_NITData= sapply(c(1:50), function(j){
    i=b_NITData[,j]
    sum(apply(trainSet_elf4[[j]], 1, function(k){
        print(length(which(tp %in% i)))
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    
#get random
set.seed(123)
b_randomData=sapply(c(1:50), function(i){
    sort(sample(tp, 8))
})
#get even
b_evenData=sapply(c(1:50), function(i){
    c(-4, 1, 8, 16, 22, 28, 32, 40)
    #-4
})

l2_even= sapply(c(1:50), function(j){
    i=b_evenData[,j]
    sum(apply(trainSet_elf4[[j]], 1, function(k){
        print(length(which(tp %in% i)))
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})    




l2_random= sapply(c(1:50), function(j){
    i=b_randomData[,j]
    sum(apply(trainSet_elf4[[j]], 1, function(k){
        print(length(which(tp %in% i)))
        L2(tp, rep(0, length(tp)), k, min(tp), max(tp), which(tp %in% i))
    }))
    
})

####Generate specific images for paper, utilising this data:

###Draw graphs of all the raw data:
pdf(paste("Figure_testTPS_data.pdf", sep=""), width=8, height=3,pointsize = 8)
par(mfrow=c(1,3));
par(mar=c(4, 4, 1, 1));
par(cex=1)

dataGenes=read.table("input.data.txt", sep=',')
tp=as.numeric(sapply(as.character(colnames(dataGenes)), function(j){substring(j, 2, nchar(j))}))
plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(dataGenes), max(dataGenes)), xlab="time", ylab="gene expression")
apply(dataGenes, 1, function(i){lines(tp, i, lwd=0.2)})

geneExpr=read.csv("col0GeneExpression.csv", header=T)
genesOfInterest=read.table("elf4Genes.txt")
dataGenes=geneExpr[which(as.character(geneExpr[,1]) %in% as.character(genesOfInterest[,1])),4:19]
tp=sapply(colnames(geneExpr)[4:19], function(i){
    a=substring(i, 2, nchar(i))
    a=gsub('\\.', '-', a)
    as.numeric(a)
})
plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(dataGenes), max(dataGenes)), xlab="time", ylab="gene expression")
apply(dataGenes, 1, function(i){lines(tp, i, lwd=0.5)})

genes_elegans=read.csv('Celegans_Table_S2.csv', header=T)
tp=sapply(colnames(genes_elegans), function(i){
    a=strsplit(as.character(i), "_")[[1]][2]
    as.numeric(substring(a, 5, nchar(a)))
})

#read in interesting C. elegans genes
interesting_genes=as.character(read.table("CelegansInterestingGenes")[,1])
dataGenes=genes_elegans[interesting_genes,]
plot(c(), xlim=c(min(tp), max(tp)), ylim=c(min(dataGenes), max(dataGenes)), xlab="time", ylab="gene expression")
apply(dataGenes, 1, function(i){lines(tp, i, lwd=0.5)})

dev.off()
###Draw comparative results for each pair of values:
pdf(paste("Figure_TPSvsNITPicker_logL2error.pdf", sep=""), width=12, height=3,pointsize = 10)
par(mfrow=c(1,5));
par(mar=c(3, 4, 1, 1));
par(cex=1)


files=c('tenth.RData', 'elegans.RData', 'elegans_smaller.RData', 'elf4.RData', 'elf4_smaller.RData')

TPSResults=list()
NITResults=list()
pval=c()
for(f in files){
    print(f)
    load(f)
    TPSResults[[f]]=l2_TPSData
    NITResults[[f]]=l2_NITData
    pval=ceiling(t.test(log(TPSResults[[f]]), log(NITResults[[f]]))$p.value*1000)/1000
    vioplot(log(TPSResults[[f]]), log(NITResults[[f]]), col="grey", names=c("TPS", "NITPicker"))
    title(paste("p-value<", pval, sep=""), ylab="log(L2-error)")
    points(1+jitter(rep(0, length(l2_TPSData)), amount=0.2), log(l2_TPSData), pch=19, col=rgb(0, 0, 0.5, 0.5))
    points(2+jitter(rep(0, length(l2_NITData)), amount=0.2), log(l2_NITData), pch=19, col=rgb(0, 0, 0.5, 0.5))
    # pval=c(pval, t.test(TPSResults[f], NITResults[f])$p.value)
}

dev.off()