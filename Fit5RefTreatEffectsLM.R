#5x LINEAR MODELLING WITH REFERENCE+TREATMENT EFFECTS PARAMETERIZATION------------

#This script assumes you have already prepared the following objects using PrepareData.R:
#1. the log2cpm counts data matrix: dat
#2. the factor: tissue
#3. the design matrix (*Which is actually a data.frame*) containing sample information: design
#4. the numeric vector: norm.factors

#it also assumes you have already loaded the packeges listed in PrepareData.R.

#finally, this will be less heavily annotated with information that also applies to doing a single
#fit with a contrast sums matrix.
#instead, annotations will relate more to 'new stuff' being done when doing 5 of everything.

#Step 1: Create 5 reference + treatment effects model matrices, setting each of the groups as reference once.-----------------
#NOTE: whereas Ref+Tr. Eff. makes the first level of the factor the X intercept,
#contrast sums dummies out the last factor in the level. 
#original order in 'tissue' is corm, flower, leaf, stem, stolon.

#NOTES for releveling factors. for ex, see the 1st line of #Dummying stem:
#this first line identifies a new order for the levels. The numbers correspond to the indices of the levels:
#in 'tissue' 1 = corm, 2 = flower, etc. the numbers in [c(5,1:4)] denote the leves in 'tissue', 
#but the order in which they are given in [c(5,1:4)] denotes thier order in the 'tissue_stem' 
#(that is, the 5th level of 'tissue' is the 1st listed in 'tissue.stem', the 1st of 'tissue' = the 2nd in 'tissue.stem', etc...f

#Reference = stolon
tissue.stol <- factor(tissue,levels(tissue)[c(5,1:4)])
modelMat.RT.stol <- model.matrix(~tissue.stol)
show(modelMat.RT.stol)  #you can turn this off if you want.

#Reference = stem
tissue.stem <- factor(tissue,levels(tissue)[c(4:5,1:3)])
modelMat.RT.stem <- model.matrix(~tissue.stem)
show(modelMat.RT.stem) #you can turn this off if you want.

#Reference = leaf
tissue.leaf <- factor(tissue,levels(tissue)[c(3:5,1:2)])
modelMat.RT.leaf <- model.matrix(~tissue.leaf)
show(modelMat.RT.leaf) #you can turn this off if you want.

#Reference = flower
tissue.flow <- factor(tissue,levels(tissue)[c(2:5,1)])
modelMat.RT.flow <- model.matrix(~tissue.flow)
show(modelMat.RT.flow)  #you can turn this off if you want.

#Reference = corm
tissue.corm  <-  tissue
modelMat.RT.corm <- model.matrix(~tissue.corm)
show(modelMat.RT.corm)  #you can turn this off if you want.

#Step 2: Voom, Limma, and eBayes -ing all of these------------
#Note: Voom calculates different weights (apparently) for different model matrices.
#checked the voom weights from 2 different contrast sum model matrices, and the numbers differed by <2e-15.
#just in case, I will continue making separate voom values for each version, but it might not be important.

#VoomLimBays stolon
dat.voomed.stol <- voom(dat, modelMat.RT.stol, plot = TRUE, lib.size = colSums(dat) * norm.factor) #this calculates more accurate variances.
fit.RT.stol <- lmFit(dat.voomed.stol, modelMat.RT.stol) #for each genes, this fits coefficients to observed log2(cpm) as described by design.
fit.RT.stol <- eBayes(fit.RT.stol) #I will make notes later on what this does.

#VoomLimBays stem
dat.voomed.stem <- voom(dat, modelMat.RT.stem, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.RT.stem <- lmFit(dat.voomed.stem, modelMat.RT.stem)
fit.RT.stem <- eBayes(fit.RT.stem)

#VoomLimBays leaf
dat.voomed.leaf <- voom(dat, modelMat.RT.leaf, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.RT.leaf <- lmFit(dat.voomed.leaf, modelMat.RT.leaf)
fit.RT.leaf <- eBayes(fit.RT.leaf)

#VoomLimBays flower
dat.voomed.flow <- voom(dat, modelMat.RT.flow, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.RT.flow <- lmFit(dat.voomed.flow, modelMat.RT.flow)
fit.RT.flow <- eBayes(fit.RT.flow)

#VoomLimBays corm
dat.voomed.corm <- voom(dat, modelMat.RT.corm, plot = TRUE, lib.size = colSums(dat) * norm.factor)
fit.RT.corm <- lmFit(dat.voomed.corm, modelMat.RT.corm)
fit.RT.corm <- eBayes(fit.RT.corm)

#Step 3: topTabling the results---------------
#This is just to see the coefficients. 
#The F value and p values are useless, because coefficient 1 means something different from all the others.
ttfit.RT.stol.all <- topTable(fit.RT.stol, n = Inf)
ttfit.RT.stem.all <- topTable(fit.RT.stem, n = Inf)
ttfit.RT.leaf.all <- topTable(fit.RT.leaf, n = Inf)
ttfit.RT.flow.all <- topTable(fit.RT.flow, n = Inf)
ttfit.RT.corm.all <- topTable(fit.RT.corm, n = Inf)

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

ttfit.RT.stol <- topTable(fit.RT.stol, coef = c(2,3,4,5), n = Inf)
ttfit.RT.stem <- topTable(fit.RT.stem, coef = c(2,3,4,5), n = Inf)
ttfit.RT.leaf <- topTable(fit.RT.leaf, coef = c(2,3,4,5), n = Inf)
ttfit.RT.flow <- topTable(fit.RT.flow, coef = c(2,3,4,5), n = Inf)
ttfit.RT.corm <- topTable(fit.RT.corm, coef = c(2,3,4,5), n = Inf)

#MORE ANALYSIS TO BE DONE IN LATER SCRIPTS-----------
#the main "Important objects" to be passed on are the 'fit' objects:
#fit.RT.stol
#fit.RT.stem
#fit.RT.leaf
#fit.RT.flower
#fit.RT.corm
