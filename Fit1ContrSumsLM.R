#LINEAR MODELLING WITH CONTRAST SUMS PARAMETERIZATION------------

#This script assumes you have already prepared the following objects using PrepareData.R:
#1. the log2cpm counts data matrix: dat
#2. the factor: tissue
#3. the design matrix (*Which is actually a data.frame*) containing sample information: des
#4. the numeric vector: norm.factors

#it also assumes you have already loaded the packeges listed in PrepareData.R.

#Step 1. create model matrix.-------

#WARNING: In contrast sums parameterization, you must select one group to be "dummied out". 
#its coefficient is still called 'Intercept' (as in Ref.+treat. eff.) but it is fitted to all expression values for that gene.
#in effect, it represents the "grand average" of all expression values (before you break them up into any group).
#Therefore, you will always have to sacrifice one of the tissue groups in order to compare other groups to the grand mean.
#for this fit, I have selected corm.

#Dummying corm
tissue_corm  <-  factor(tissue,levels(tissue)[c(2:5,1)]) #this was 'the way better than relevel' I was talking about.
modelMat.CS.corm <- model.matrix(~tissue_corm, contrasts = list(tissue_corm = "contr.sum"))
show(modelMat.CS.corm)

#To change the reference group:
#you have to change the order of the levels ("corm","flower","leaf","stem","stolon" of the factor "tissue".
#I currently use relevel() but there are probably better ways (now see above!)
#NOTE: the row order does not change, but the level order does.
show(tissue)
tissue_ex <- relevel(tissue, "stolon")
show(modelMat.CS_ex <- model.matrix(~tissue_corm, contrasts = list(tissue_corm = "contr.sum")))
#note that the column titles have moved, and so have the 1's representing each group.

#Step 2. Voom and Limma------------------#
dat.voomed.CS.corm <- voom(dat, modelMat.CS.corm, plot = TRUE, lib.size = colSums(dat) * norm.factor) #this calculates more accurate variances
fit.CS.corm  <- lmFit(dat.voomed.CS.corm, modelMat.CS.corm) #for each genes, this fits coefficients to observed cpm as described by design.
fit.CS.corm <- eBayes(fit.CS.corm) #I will make notes later on what this does.

#Step 3. topTabling the results-----------
ttfit.CS.corm <- topTable(fit.CS.corm, n = Inf)

#Step 4. Play around with different comparisons--------------

#by adjusting the coefficient numbers, you can look at just one coefficent (get t values)
#look at various combinations (getting F, **WHICH = the likelihood that any coefficient is zero**)
#p.values and adj.P.values (corrected for multiple hypothesis testing with with BH method) are always given.
#I have no idea what B means right now (This column shows up when fitting just 1 coefficient)
#the list is organized in order of ascending adj.P.

#WARNING: Recall that the Intercept (coef 1) is just the overall average gene expression, whereas all the other coefficients
#represent differences in averages. Since in topTable, F is just the probability that any of these coefficients are 0,
#if you include the coef 1 (the intercept), your F will no longer represent how much group averages vary from 
#the grand average, but a combination of that along with how much the grand average itself varies from zero (Not very informative!!!)

#so, put different numbers in the 'coef' option
#and adjust the output table accordingly, ex:
#NAMEofTT <- topTable(fit, coef = c(NUMBER,NUMBER,etc...), n = Inf)
ttfit.CS.corm.c1 <- topTable(fit.CS.corm, coef = c(1), n = Inf)

#use head() to look at just the top 6 lines
head(ttfit.CS.corm.c2 <- topTable(fit.CS.corm, coef = c(2), n = Inf))
head(ttfit.CS.corm.c23 <- topTable(fit.CS.corm, coef = c(2, 3), n = Inf))
head(ttfit.CS.corm.c523 <- topTable(fit.CS.corm, coef = c(5, 2, 3), n = Inf))

#Step 5. get the list length you want -----------------
#Ex:
#ttfit.CS.cormall_trimmed <- topTable(fit.CS.corm, coef = c(2), n = PUT NUMBER OF LINES HERE)
trimmed_ttfit.CS.corm.c2 <- topTable(fit.CS.corm, coef = c(2), n = 1000)

#Step 6. Get a list of genes with the adj.P.Val cutoffs you want--------------
#NOTE: You can use genelist with PlotListOfGenes.R to see expression plots. 
#to adapt this code to look at any table and column of interest:
#cutoff <- whateverNumber
#listOfRownames <- rownames(tableWithColumnOfInterest[which(tableWithColumnOfInterest$columnOfInterest), ])

cutoff <- 1e-5
genelist<- rownames(ttfit.CS.corm[which(ttfit.CS.corm$adj.P.Val < cutoff), ])
genelist_segment <- genelist[100:150]

#Step 7. See the TopTable results for these genes----------
#to adapt this code to look at any table and column of interest:
#tableWithOnlyRowsInList <- (tableWithColumnOfInterest[listOfRownames, ])
trimmed_ttfit.CS.corm <- (ttfit.CS.corm[genelist, ])

#FOR ATTEMPTS AT "UNBIASED" ANALYSIS------------
#Still working on this. You'd have to do it 5 times and find some way of combining the values. 
#will put out a "larger analysis" script soon. 