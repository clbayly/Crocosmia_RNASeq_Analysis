#LINEAR MODELLING WITH CELL MEANS PARAMETERIZATION-----------

#This script assumes you have already prepared the following objects using PrepareData.R:
#1. the log2cpm counts data matrix: dat
#2. the factor: tissue
#3. the design matrix (*Which is actually a data.frame*) containing sample information: des
#4. the numeric vector: norm.factors

#it also assumes you have already loaded the packeges listed in PrepareData.R.

#Step 1. create model matrix.-------
modelMat <- model.matrix(~0 + tissue)
show(modelMat) #this is optional.

#Step 2. Voom and Limma------------------
dat.voomed <- voom(dat, modelMat, plot = TRUE, lib.size = colSums(dat) * norm.factor) #this calculates more accurate variances
fit <- lmFit(dat.voomed, modelMat) #for each genes, this fits coefficients to observed cpm as described by design.
fit <- eBayes(fit)

#Step 3. topTabling the results-----------
head(ttfitall <- topTable(fit, n = Inf))

#Step 4. Play around with different comparisons--------------
#by adjusting the coefficient numbers, you can look at just one coefficent (get t values) or
#look at various combinations (getting F, **WHICH represents the likelihood that any coefficient is zero**)
#for cell means, this is not a particularly helpful number! genes with high overall expression will give the highest F values.
#p.values and adj.P.values (corrected for multiple hypothesis testing with with BH method) are always given.
#I have no idea what B means right now (This column shows up when fitting just 1 coefficient)
#the list is organized in order of ascending adj.P.

#so, put different numbers in the 'coef' option
#and adjust the output table accordingly, ex:
#NAMEofTT <- topTable(fit, coef = c(NUMBER,NUMBER,etc...), n = Inf)
ttfit.c1 <- topTable(fit, coef = c(1), n = Inf)

#use head() to look at just the top 6 lines
head(ttfit.c2 <- topTable(fit, coef = c(2), n = Inf))
head(ttfit.c23 <- topTable(fit, coef = c(2, 3), n = Inf))
head(ttfit.c523 <- topTable(fit, coef = c(5, 2, 3), n = Inf))

#Step 5. get the list length you want -----------------
#Ex:
#trimmed_ttfit.c2 <- topTable(fit, coef = c(2), n = PUT NUMBER OF LINES HERE)
trimmed_ttfit.c2 <- topTable(fit, coef = c(2), n = 1000)

#Step 6. get a list of genes with the adj.P.Val cutoffs you want--------------
cutoff <- 1e-5
genelist<- which(ttfitall$adj.P.Val < cutoff)
trimmed_ttfitall <- (ttfitall[genelist, ])

#to adapt this code to look at any table and column of interest:
#cutoff <- whateverNumber
#listOfRownames <- which(tableWithColumnOfInterest$columnOfInterest)
#tableWithOnlyRowsInList <- (tableWithColumnOfInterest[listOfRownames, ])

  