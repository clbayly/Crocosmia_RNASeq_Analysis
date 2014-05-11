#5x LINEAR MODELLING WITH CONTRAST SUMS PARAMETERIZATION------------

#This script assumes you have already prepared the following objects using PrepareData.R:
#1. the log2cpm counts data matrix: dat
#2. the factor: tissue
#3. the design matrix (*Which is actually a data.frame*) containing sample information: design
#4. the numeric vector: norm.factors

#it also assumes you have already loaded the packeges listed in PrepareData.R.

#finally, this will be less heavily annotated with information that also applies to doing a single
#fit with a contrast sums matrix.
#instead, annotations will relate more to 'new stuff' being done when doing 5 of everything.

#Step 1: Create 5 contr.sums model matrices, switching out each of the groups once.-----------------
#NOTE: whereas Ref+Tr. Eff. makes the first level of the factor the X intercept,
#contrast sums dummies out the last factor in the level. 
#original order in 'tissue' is corm, flower, leaf, stem, stolon.

#NOTES for releveling factors. for ex, see the 1st line of #Dummying stem:
#this first line identifies a new order for the levels. The numbers correspond to the indices of the levels:
#in 'tissue' 1 = corm, 2 = flower, etc. the numbers in [c(5,1:4)] denote the leves in 'tissue', 
#but the order in which they are given in [c(5,1:4)] denotes thier order in the 'tissue_stem' 
#(that is, the 5th level of 'tissue' is the 1st listed in 'tissue.stem', the 1st of 'tissue' = the 2nd in 'tissue.stem', etc...

#Dummying stolon
tissue.stol <- tissue
modelMat.CS.stol <- model.matrix(~tissue.stol, contrasts = list(tissue.stol = "contr.sum"))
show(modelMat.CS.stem)  #you can turn this off if you want.

#Dummying stem
tissue.stem <- factor(tissue,levels(tissue)[c(5,1:4)])
modelMat.CS.stem <- model.matrix(~tissue.stem, contrasts = list(tissue.stem = "contr.sum"))
show(modelMat.CS.stem) #you can turn this off if you want.

#Dummying leaf
tissue.leaf <- factor(tissue,levels(tissue)[c(4:5,1:3)])
modelMat.CS.leaf <- model.matrix(~tissue.leaf, contrasts = list(tissue.leaf = "contr.sum"))
show(modelMat.CS.leaf) #you can turn this off if you want.

#Dummying flower
tissue.flow <- factor(tissue,levels(tissue)[c(3:5,1:2)])
modelMat.CS.flow <- model.matrix(~tissue.flow, contrasts = list(tissue.flow = "contr.sum"))
show(modelMat.CS.flow)  #you can turn this off if you want.

#Dummying corm
tissue.corm  <-  factor(tissue,levels(tissue)[c(2:5,1)])
modelMat.CS.corm <- model.matrix(~tissue.corm, contrasts = list(tissue.corm = "contr.sum"))
show(modelMat.CS.corm)  #you can turn this off if you want.

#Step 2: Voom, Limma, and eBayes -ing all of these------------
#Note: Voom calculates different weights (apparently) for different model matrices.
#checked the voom weights from 2 different contrast sum model matrices, and the numbers differed by <2e-15.
#just in case, I will continue making separate voom values for each version, but it might not be important.

#VoomLimBays stolon
dat.voomed.stol <- voom(dat, modelMat.CS.stol, plot = TRUE, lib.size = colSums(dat) * norm.factor) #this calculates more accurate variances.
fit.CS.stol <- lmFit(dat.voomed.stol, modelMat.CS.stol) #for each genes, this fits coefficients to observed log2(cpm) as described by design.
fit.CS.stol <- eBayes(fit.CS.stol) #I will make notes later on what this does.

#VoomLimBays stem
dat.voomed.stem <- voom(dat, modelMat.CS.stem, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.CS.stem <- lmFit(dat.voomed.stem, modelMat.CS.stem)
fit.CS.stem <- eBayes(fit.CS.stem)

#VoomLimBays leaf
dat.voomed.leaf <- voom(dat, modelMat.CS.leaf, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.CS.leaf <- lmFit(dat.voomed.leaf, modelMat.CS.leaf)
fit.CS.leaf <- eBayes(fit.CS.leaf)

#VoomLimBays flower
dat.voomed.flow <- voom(dat, modelMat.CS.flow, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.CS.flow <- lmFit(dat.voomed.flow, modelMat.CS.flow)
fit.CS.flow <- eBayes(fit.CS.flow)

#VoomLimBays corm
dat.voomed.corm <- voom(dat, modelMat.CS.corm, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.CS.corm <- lmFit(dat.voomed.corm, modelMat.CS.corm)
fit.CS.corm <- eBayes(fit.CS.corm)

#Step 3: topTabling the results---------------
#This is just to see the coefficients. 
#The F value and p values are useless, because coefficient 1 means something different from all the others.
ttfit.CS.stol.all <- topTable(fit.CS.stol, n = Inf)
ttfit.CS.stem.all <- topTable(fit.CS.stem, n = Inf)
ttfit.CS.leaf.all <- topTable(fit.CS.leaf, n = Inf)
ttfit.CS.flow.all <- topTable(fit.CS.flow, n = Inf)
ttfit.CS.corm.all <- topTable(fit.CS.corm, n = Inf)

#Step 4: topTabling the results minus the grand average---------------
#The F values from these (and associated adj. P Values) will indicate the "statistically significant trendedness"
#of the data - the extent to which gene expression behaves differently in different groups. 
#However, this value will be biased by whatever the reference was chosen to be (in ref+tx).

#recall that while F is supposed to represent 'explained differences'(ie - variation of group means from grand mean) 
#divided by 'unexplained differences' (namely, variation of group samples from group mean)
#in topTable, F is only the likelihood that all the groups are equal to zero. 
#so if there is another group with a similar expression level to the reference group, with ref+tx parameterization,
#it would contribute mostly to 'unexplained variance' (the denominator in F) wheras with contrast sums parameterization, 
#only a sample with a mean equal to the grand mean would contribute just to 'unexplained variance'. 
ttfit.CS.stol <- topTable(fit.CS.stol, coef = c(2,3,4,5), n = Inf)
ttfit.CS.stem <- topTable(fit.CS.stem, coef = c(2,3,4,5), n = Inf)
ttfit.CS.leaf <- topTable(fit.CS.leaf, coef = c(2,3,4,5), n = Inf)
ttfit.CS.flow <- topTable(fit.CS.flow, coef = c(2,3,4,5), n = Inf)
ttfit.CS.corm <- topTable(fit.CS.corm, coef = c(2,3,4,5), n = Inf)

#MORE ANALYSIS TO BE DONE IN LATER SCRIPTS-----------
#the main "Important objects" to be passed on are the 'fit' objects:
#fit.CS.stol
#fit.CS.stem
#fit.CS.leaf
#fit.CS.flower
#fit.CS.corm

...