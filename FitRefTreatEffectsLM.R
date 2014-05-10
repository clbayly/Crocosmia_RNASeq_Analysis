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

designRT <- model.matrix(~tissue) #Apologies! this is actually a model matrix. calling it 'design' is misleading.
show(designRT) #so you can see that the first coef must fit all columns, but is uniquely fitted to the intercept group
                #and the other coefs have to be fitted to both the intercept group and one other group 
                #(accounting for systematic differences between the intercept and that other group)

#To change reference group:
#I currently use Relevel() but there are probably better ways.
#trick: the row order does not change, but the level order does.
show(tissue)
tissue_ex <- relevel(tissue, "stolon")
show(designRT_ex <- model.matrix(~tissue_ex)) #note that the column titles have moved, and so have the 1's representing each group.

#Step 2. Voom and Limma------------------
dat.voomed.RT <- voom(dat, designRT, plot = TRUE, lib.size = colSums(dat) * norm.factor) #this calculates more accurate variances
fit.RT <- lmFit(dat.voomed.RT, designRT) #for each genes, this fits coefficients to observed cpm as described by design.
fit.RT <- eBayes(fit.RT) #I will make notes later on what this does.

#Step 3. topTabling the results-----------
head(ttfitall.RT <- topTable(fit.RT, n = Inf))

#Step 4. Play around with different comparisons--------------

#WARNING: Recall that in reference

#by adjusting the coefficient numbers, you can look at just one coefficent (get t values)
#look at various combinations (getting F, **WHICH represents the likelihood that any coefficient is zero**)
#p.values and adj.P.values (corrected for multiple hypothesis testing with with BH method) are always given.
#I have no idea what B means right now (This column shows up when fitting just 1 coefficient)
#the list is organized in order of ascending adj.P.

#so, put different numbers in the 'coef' option
#and adjust the output table accordingly, ex:
#NAMEofTT <- topTable(fit, coef = c(NUMBER,NUMBER,etc...), n = Inf)
ttfit.RT_c1 <- topTable(fit.RT, coef = c(1), n = Inf)

#use head() to look at just the top 6 lines
head(ttfit.RT_c2 <- topTable(fit.RT, coef = c(2), n = Inf))
head(ttfit.RT_c23 <- topTable(fit.RT, coef = c(2, 3), n = Inf))
head(ttfit.RT_c523 <- topTable(fit.RT, coef = c(5, 2, 3), n = Inf))

#Step 5. get the list length you want -----------------
#Ex:
#ttfit.RTall_trimmed <- topTable(fit.RT, coef = c(2), n = PUT NUMBER OF LINES HERE)
ttfit.RTall_trimmed <- topTable(fit.RT, coef = c(2), n = 1000)

#Step 6. get a list of genes with the adj.P.Val cutoffs you want--------------
cutoff <- 1e-5
genelist<- which(ttfit.RTall$adj.P.Val < cutoff)
trimmed_ttfit.RTall <- (ttfit.RTall[genelist, ])

#to adapt this code to look at any table and column of interest:
#cutoff <- whateverNumber
#listOfRownames <- which(tableWithColumnOfInterest$columnOfInterest)
#tableWithOnlyRowsInList <- (tableWithColumnOfInterest[listOfRownames, ])

#FOR ATTEMPTS AT "UNBIASED" ANALYSIS------------
#Still working on this. You'd have to do it 5 times and find some way of combining the values. 
#will put out a "larger analysis" script soon. 