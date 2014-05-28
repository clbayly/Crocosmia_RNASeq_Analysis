#gathering Fs for RefTx.

#i. Prepare a table of F-values from the toptables of the five fits---------------------------------------
#First: make mini-dataframes with F-values and rownames from each contrast sum fit.
column <- 6 #this is for the F.value.
ttfit.RT.stol.F <- data.frame(rownames = rownames(ttfit.RT.stol), RT.stol.F = ttfit.RT.stol[,column])
ttfit.RT.stem.F <- data.frame(rownames = rownames(ttfit.RT.stem), RT.stem.F = ttfit.RT.stem[,column])
ttfit.RT.leaf.F <- data.frame(rownames = rownames(ttfit.RT.leaf), RT.leaf.F = ttfit.RT.leaf[,column])
ttfit.RT.flow.F <- data.frame(rownames = rownames(ttfit.RT.flow), RT.flow.F = ttfit.RT.flow[,column])
ttfit.RT.corm.F <- data.frame(rownames = rownames(ttfit.RT.corm), RT.corm.F = ttfit.RT.corm[,column])

#Second: merge mini-dataframes by rownames.
ttfit.RT.all.F <- merge(ttfit.RT.stol.F, ttfit.RT.stem.F, by = "rownames")
ttfit.RT.all.F <- merge(ttfit.RT.all.F, ttfit.RT.leaf.F, by = "rownames")
ttfit.RT.all.F <- merge(ttfit.RT.all.F, ttfit.RT.flow.F, by = "rownames")
ttfit.RT.all.F <- merge(ttfit.RT.all.F, ttfit.RT.corm.F, by = "rownames")
head(ttfit.RT.all.F)

#ii. Prepare a table of F adj.P.Val's (one column for each contrast sum fit)----------------
#First: make mini-dataframes with F adj.P.Val's and rownames from each contrast sum fit.
column <- 8 #this is for the F.value.
ttfit.RT.stol.FadjP <- data.frame(rownames = rownames(ttfit.RT.stol), RT.stol.FadjP =  ttfit.RT.stol[,column])
ttfit.RT.stem.FadjP <- data.frame(rownames = rownames(ttfit.RT.stem), RT.stem.FadjP = ttfit.RT.stem[,column])
ttfit.RT.leaf.FadjP <- data.frame(rownames = rownames(ttfit.RT.leaf), RT.leaf.FadjP = ttfit.RT.leaf[,column])
ttfit.RT.flow.FadjP <- data.frame(rownames = rownames(ttfit.RT.flow), RT.flow.FadjP = ttfit.RT.flow[,column])
ttfit.RT.corm.FadjP <- data.frame(rownames = rownames(ttfit.RT.corm), RT.corm.FadjP = ttfit.RT.corm[,column])

#Second: merge mini-dataframes by rownames.
ttfit.RT.all.FadjP <- merge(ttfit.RT.stol.FadjP, ttfit.RT.stem.FadjP, by = "rownames")
ttfit.RT.all.FadjP <- merge(ttfit.RT.all.FadjP, ttfit.RT.leaf.FadjP, by = "rownames")
ttfit.RT.all.FadjP <- merge(ttfit.RT.all.FadjP, ttfit.RT.flow.FadjP, by = "rownames")
ttfit.RT.all.FadjP <- merge(ttfit.RT.all.FadjP, ttfit.RT.corm.FadjP, by = "rownames")
head(ttfit.RT.all.FadjP)

#iii. Using the FadjP matrix, create a "decision matrix" yielding 1's or 0's based on whether 
cutoff  <- 0.05
ttfit.RT.all.FadjP.cut <- apply(ttfit.RT.all.FadjP[2:ncol(ttfit.RT.all.FadjP)], 
                                2, function(x) as.numeric(x<cutoff)) #these two lines encode the 'decision-making' part
rownames(ttfit.RT.all.FadjP.cut) <- ttfit.RT.all.FadjP$rownames #We might nees these. I'm not actually sure.

ttfit.RT.all.FadjP.cut.Fsum <- as.data.frame(apply(ttfit.RT.all.FadjP.cut, 1, sum)) #this adds all the 1's in each row.
ttfit.RT.all.FadjP.cut.Fsum <- cbind(rownames(ttfit.RT.all.FadjP.cut), ttfit.RT.all.FadjP.cut.Fsum)
colnames(ttfit.RT.all.FadjP.cut.Fsum) <- c("comps", "Fsum")
head(ttfit.RT.all.FadjP.cut.Fsum) 

length(genesWith5F <- which(ttfit.RT.all.FadjP.cut.Fsum$Fsum == 5))
length(genesWith4F <- which(ttfit.RT.all.FadjP.cut.Fsum$Fsum == 4))
length(genesWith3F <- which(ttfit.RT.all.FadjP.cut.Fsum$Fsum == 3))
length(genesWith2F <- which(ttfit.RT.all.FadjP.cut.Fsum$Fsum == 2))
length(genesWith1F <- which(ttfit.RT.all.FadjP.cut.Fsum$Fsum == 1)) #dummyCheck:
length(genesWith0F <- which(ttfit.RT.all.FadjP.cut.Fsum$Fsum == 0)) #if you try this with Fsum == 6, it should 0.

#Now we have lists of genes.
#we must now prepare a table of tadjP values for all possible combinations,
#to determine which relationships are responsible for F values below cutoff.

