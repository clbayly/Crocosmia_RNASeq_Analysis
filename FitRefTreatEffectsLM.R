#LINEAR MODELLING WITH REFERENCE+TREATMENT EFFECTS PARAMETERIZATION------------

#This script assumes you have already prepared the following objects using PrepareData.R:
#1. the log2cpm counts data matrix: dat
#2. the factor: tissue
#3. the design matrix (*Which is actually a data.frame*) containing sample information: des
#4. the numeric vector: norm.factors

#it also assumes you have already loaded the packeges listed in PrepareData.R.

#Step 1. create model matrix.-------

#WARNING: In reference + treatment effects, you must select one group to be the reference (and its coefficient is called 'Intercept'.) 
#This means: 
#To the extent that the coefficients from cell means represent the average gene expression of each group,
#the coefficients from Ref+TrEff represent the average gene expression of the intercept group (column/coefficient 1) and
#the differences between (the average gene expression in all other groups) and (that of the intercept group) (columns/coefs 2-5)

modelMat.RT.corm <- model.matrix(~tissue) #Apologies! this is actually a model matrix. calling it 'design' is misleading.
show(modelMat.RT.corm) #so you can see that the first coef must fit all columns, but is uniquely fitted to the intercept group
                #and the other coefs have to be fitted to both the intercept group and one other group 
                #(accounting for systematic differences between the intercept and that other group)

#To change the reference group:
#you have to change the order of the levels ("corm","flower","leaf","stem","stolon" of the factor "tissue".
#I currently use relevel() but there are probably better ways (see code in the constrast sums script)
#NOTE: the row order does not change, but the level order does.
show(tissue)
tissue.ex <- relevel(tissue, "stolon")
show(modelMat.RT.ex <- model.matrix(~tissue.ex)) #note that the column titles have moved, and so have the 1's representing each group.

#Step 2. Voom and Limma------------------
dat.voomed.RT.corm <- voom(dat, modelMat.RT.corm, plot = TRUE, lib.size = colSums(dat) * norm.factor) #this calculates more accurate variances
fit.RT.corm <- lmFit(dat.voomed.RT.corm, modelMat.RT.corm) #for each genes, this fits coefficients to observed cpm as described by design.
fit.RT.corm <- eBayes(fit.RT.corm) #I will make notes later on what this does.

#Step 3. topTabling the results-----------
head(ttfit.RT.corm <- topTable(fit.RT.corm, n = Inf))

#Step 4. Play around with different comparisons--------------

#by adjusting the coefficient numbers, you can look at just one coefficent (get t values)
#look at various combinations (getting F, **WHICH represents the likelihood that any coefficient is zero**)
#p.values and adj.P.values (corrected for multiple hypothesis testing with with BH method) are always given.
#I have no idea what B means right now (This column shows up when fitting just 1 coefficient)
#the list is organized in order of ascending adj.P.

#WARNING: Recall that the Intercept (coef 1) is just the average expression of one (reference) tissue group, whereas all the other 
#coefficients represent differences between other group averages and the reference. Since (in topTable) F is just the probability that all 
#of these coefficients are 0, if you include the coef 1 (the reference) when you calculate F for your combination of coefficients,
#your F will no longer represent how much group averages vary from a constant value (the reference average, which is set to zero for
#nonreference groups) but a combination of that along with how much the reference average itself varies from zero (Not very informative!!!)

#so, put different numbers in the 'coef' option
#and adjust the output table accordingly, ex:
#NAMEofTT <- topTable(fit, coef = c(NUMBER,NUMBER,etc...), n = Inf)
ttfit.RT.corm.c1 <- topTable(fit.RT.corm, coef = c(1), n = Inf)

#use head() to look at just the top 6 lines
head(ttfit.RT.corm.c2 <- topTable(fit.RT.corm, coef = c(2), n = Inf))
head(ttfit.RT.corm.c23 <- topTable(fit.RT.corm, coef = c(2, 3), n = Inf))
head(ttfit.RT.corm.c523 <- topTable(fit.RT.corm, coef = c(5, 2, 3), n = Inf))

#Step 5. get the list length you want -----------------
#Ex: ttfit.RTall_trimmed <- topTable(fit.RT, coef = c(2), n = PUT NUMBER OF LINES HERE)
trimmed_ttfit.RT.c2 <- topTable(fit.RT.corm, coef = c(2), n = 1000)

#Step 6. get a list of genes with the adj.P.Val cutoffs you want--------------
#NOTE: You can use genelist with PlotListOfGenes.R to see expression plots. 
#to adapt this code to look at any table and column of interest:
#cutoff <- whateverNumber
#listOfRownames <- rownames(tableWithColumnOfInterest[which(tableWithColumnOfInterest$columnOfInterest), ])

cutoff <- 1e-5
genelist<- rownames(ttfit.RT.corm[which(ttfit.RT.corm$adj.P.Val < cutoff), ])
genelist_segment <- genelist[100:150]

#Step 7. See the TopTable results for these genes----------
#to adapt this code to look at any table and column of interest:

trimmed_ttfit.RT.corm <- (ttfit.RT.corm[genelist, ])

#FOR ATTEMPTS AT "UNBIASED" ANALYSIS------------
#Still working on this. You'd have to do it 5 times and find some way of combining the values. 
#will put out a "larger analysis" script soon. 