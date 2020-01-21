library(magrittr)
library(Rcpi)
library(eiR)
library(ChemmineR)
library(ChemmineOB)
library(fmcsR)
library(e1071)
library(caret)
library(RSNNS)
library(randomForest)


source('lib.R')

type <- c(1,2,3,4,5)
name <- c('Muscle','Ganglion','Heteromeric CNS','Further CNS','Homomeric CNS')
symbol <- c('(α1)2β1δε[24] or (α1)2β1δγ', '(α3)2(β4)3', '(α4)2(β2)3', '(α3)2(β4)3', '(α7)5')
## 1/Kd ~ anatoxin > epibatidine > acetylcholine > DMPP >> cytisine > pyrantel > nicotine > coniine > tubocurare > lobeline

n1A.agonists <- n1AChAg.smi <- read.SMIset('data/smi/n1AChAg.smi')
n1A.antagonists <- n1AChAg.smi <- read.SMIset('data/smi/n1AChAn.smi')

m1A.agonists <- m1AChAg.smi <- read.SMIset('data/smi/m1AChAg.smi')
m1A.antagonists <- m1AChAn.smi <- read.SMIset('data/smi/m1AChAn.smi')

# m1A.antagonists
ap <- smi2ap(m1A.antagonists)
cmp.similarity(ap[1], ap[2])
cmp.similarity(ap[1], ap[3])
cmp.similarity(ap[1], ap[10])

cluster <- cmp.cluster(ap, 2)
model.knn <- nearestNeighbors(ap, cutoff=0.5)
cluster.visualize(ap, cluster, size.cutoff=1)
cluster.jp <- jarvisPatrick(model.knn, k=2, mode="a1b")

#####
####
###
##
# fmcsBatch(m1A.agonists[[1]])
# data(sdfsample)
# sdfset <- sdfsample
# 
# m1A.agonists[[1]]
# fmcsBatch(sdfset[[1]], sdfset[1:3], au=2)
# fmcsBatch(sdfset[[1]], sdfset[1:3], bu=1)
# fmcsBatch(sdfset[[1]], sdfset[1:3], matching.mode="aromatic", au=1, bu=1)

musc.agonists <- c()
musc.antagonists <- c()
nico.agonists <- c()
nico.antagonists <- c()
inactive <- c()

for (i in 1:5){
  musc.agonists[[i]] <- readMolFromSDF(paste('data/smi', sprintf('m%iAChAg.smi',i), sep='/')); 
  musc.antagonists[[i]] <- readMolFromSDF(paste('data/smi', sprintf('m%iAChAn.smi',i), sep='/'));
  nico.agonists[[i]] <- readMolFromSDF(paste('data/smi', sprintf('n%iAChAg.smi',i), sep='/'));
  nico.antagonists[[i]] <- readMolFromSDF(paste('data/smi', sprintf('n%iAChAn.smi',i), sep='/'));
  if(length(inactive) > 0) inactive[[i]] <- readMolFromSDF(paste('data/smi', 'inactives.smi', sep='/'));
}

dats <- list(musc.agonists, musc.antagonists, nico.agonists, nico.antagonists)
descriptors <- getDescriptorsByBatch(dats)
descriptors.clean <- removeNADescriptors(descriptors)
names(descriptors.clean) <- c('mAg', 'mAn', 'nAg', 'nAn')

x.musc.agonists <- matrix(unlist(descriptors.clean$mAg), ncol=5)
x.musc.antagonists <- matrix(unlist(descriptors.clean$mAn), ncol=5)
x.nico.agonists <- matrix(unlist(descriptors.clean$nAg), ncol=5)
x.nico.antagonists <- matrix(unlist(descriptors.clean$nAn), ncol=5)
y.musc.agonists <- rep(1, nrow(x.musc.agonists))
y.musc.antagonists <- rep(2, nrow(x.musc.antagonists))
y.nico.agonists <- rep(3, nrow(x.nico.agonists))
y.nico.antagonists <- rep(4, nrow(x.nico.antagonists))

y <- c(y.musc.agonists, y.musc.antagonists, y.nico.agonists, y.nico.antagonists)
x <- rbind(x.musc.agonists, x.musc.antagonists, x.nico.agonists, x.nico.antagonists)
dset <- splitForTrainingAndTest(x=x, y=y)

model.glm <- glm(y~., data=data.frame(x=x,y=y))
model.rft <- randomForest(x, y, ntree=500) 
model.mlp <- mlp(dset$inputsTrain, dset$targetsTrain, size=2, learnFunc="Rprop", learnFuncParams=c(.1, .2, .5, 1.5), maxit=2000, inputsTest=dset$inputsTest, targetsTest=dset$targetsTest)
confusionMatrix(dset$targetsTrain, encodeClassLabels(fitted.values(model.mlp), method="402040", l=0.4, h=0.6))

summary(model.glm)
summary(model.rft)
summary(model.mlp)

par(mfrow=c(2,1))
plot(model.glm)
predictions.mlp <- predict(model.mlp, dset$inputsTest)
plotIterativeError(model.mlp)
plotRegressionError(dset$targetsTest, predictions.mlp)

## ##
## eiR Similarity Searching
## ##
eiInit(m1AChAg.smi)
validSDF(read.SDFset(paste('data/smi', sprintf('m%iAChAg.smi',1), sep='/')))
sapply(1:5, function(i){
  musc.agonists[[i]] <- eiInit(paste('data/smi', sprintf('m%iAChAg.smi',i), sep='/')) 
  musc.antagonists[[i]] <- eiInit(paste('data/smi', sprintf('m%iAChAn.smi',i), sep='/'))
  nico.agonists[[i]] <- readMolFromSDF(paste('data/smi', sprintf('n%iAChAg.smi',i), sep='/'))
  nico.antagonists[[i]] <- readMolFromSDF(paste('data/smi', sprintf('n%iAChAn.smi',i), sep='/'))
})
