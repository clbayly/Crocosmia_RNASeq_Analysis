
#previously...
ttfit.RT.stol.all <- topTable(fit.RT.stol, n = Inf)
ttfit.RT.stem.all <- topTable(fit.RT.stem, n = Inf)
ttfit.RT.leaf.all <- topTable(fit.RT.leaf, n = Inf)
ttfit.RT.flow.all <- topTable(fit.RT.flow, n = Inf)
ttfit.RT.corm.all <- topTable(fit.RT.corm, n = Inf)

#now, we must get all t relationships. 
head(ttfit.RT.stol.all)
<- topTable(fit.RT.stol, coef = 2 n = Inf) #hmmm I did this before...

#(i-a). make ttfit tables for corm (only) from each Ref+treatEff fit (except where corm was dummied out)
#fit.RT.stol coef order: "gndAvg", "corm", "flow", "leaf", "stem"
ttfit.RT.stol.corm <- topTable(fit.RT.stol, coef = c(2), n = Inf)
#fit.RT.stem coef order: "gndAvg", "stol", "corm", "flow", "leaf"
ttfit.RT.stem.corm <- topTable(fit.RT.stem, coef = c(3), n = Inf)
#fit.RT.leaf coef order: "gndAvg", "stem", "stol", "corm", "flow"
ttfit.RT.leaf.corm <- topTable(fit.RT.leaf, coef = c(4), n = Inf)
#fit.RT.flow coef order:"gndAvg", "leaf", "stem", "stol", "corm"
ttfit.RT.flow.corm <- topTable(fit.RT.flow, coef = c(5), n = Inf)

#now, there might be a better way to do this...
#this should account for all relationships between points.--------------------
head(ttfit.RT.corm)
head(ttfit.RT.corm2flow <- topTable(fit.RT.corm, coef = c(2), n = Inf))
head(ttfit.RT.corm2leaf <- topTable(fit.RT.corm, coef = c(3), n = Inf))
head(ttfit.RT.corm2stem <- topTable(fit.RT.corm, coef = c(4), n = Inf))
head(ttfit.RT.corm2stol <- topTable(fit.RT.corm, coef = c(5), n = Inf))

head(ttfit.RT.flow)
head(ttfit.RT.flow2leaf <- topTable(fit.RT.flow, coef = c(2), n = Inf))
head(ttfit.RT.flow2stem <- topTable(fit.RT.flow, coef = c(3), n = Inf))
head(ttfit.RT.flow2stol <- topTable(fit.RT.flow, coef = c(4), n = Inf))

head(ttfit.RT.leaf)
head(ttfit.RT.leaf2stem <- topTable(fit.RT.leaf, coef = c(2), n = Inf))
head(ttfit.RT.leaf2stol <- topTable(fit.RT.leaf, coef = c(3), n = Inf))

head(ttfit.RT.stem)
head(ttfit.RT.stem2stol <- topTable(fit.RT.stem, coef = c(2), n = Inf))

#plot t hitlists for each.
ttplotManyGenes <- function(toptable, title) {
  plotManyGenes(head(rownames(toptable), 16), title)
}

ttplotManyGenes(ttfit.RT.corm2flow, "ttfit.RT.corm2flow")
ttplotManyGenes(ttfit.RT.corm2leaf, "ttfit.RT.corm2leaf")
ttplotManyGenes(ttfit.RT.corm2stem, "ttfit.RT.corm2stem")
ttplotManyGenes(ttfit.RT.corm2stol, "ttfit.RT.corm2stol")
ttplotManyGenes(ttfit.RT.flow2leaf, "ttfit.RT.flow2leaf")
ttplotManyGenes(ttfit.RT.flow2stem, "ttfit.RT.flow2stem")
ttplotManyGenes(ttfit.RT.flow2stol, "ttfit.RT.flow2stol")
ttplotManyGenes(ttfit.RT.leaf2stem, "ttfit.RT.leaf2stem")
ttplotManyGenes(ttfit.RT.leaf2stol, "ttfit.RT.leaf2stol")
ttplotManyGenes(ttfit.RT.stem2stol, "ttfit.RT.stem2stol")

#yes - the plots reflect the t-values heres.

#Column 5 from each of these should be tadjP--------------------
#This is to check to see if the function below works---------------------------
head(ttfit.RT.corm)
head(ttfit.RT.corm2flow <- topTable(fit.RT.corm, coef = c(2), n = Inf))

get_tadjP_Table <- function(inputTable, adjPcolName) {
  outputTable  <- data.frame(rownames = rownames(inputTable),
                             blarg1  = inputTable[,5])
  colnames(outputTable) <- c("rownames", adjPcolName)
  return(outputTable)
}

tt.RT.corm2flow.tadjP <- get_tadjP_Table(ttfit.RT.corm2flow, "RT.corm2flow.tadjP")
dim(tt.RT.corm2flow.tadjP)

#OK the function works. now, using it in the bigger system...---------------

head(ttfit.RT.corm)
head(coefficients(fit.RT.corm))
head(ttfit.RT.corm2flow <- topTable(fit.RT.corm, coef = c(2), n = Inf))

get_tadjP_Table2 <- function(inputFitTable, coefNum, adjPcolName) {
  tvalue_ttfit <- topTable(inputFitTable, coef = c(coefNum), n = Inf)
  outputTable  <- data.frame(rownames = rownames(tvalue_ttfit),
                             blarg1  = tvalue_ttfit[,5])
  colnames(outputTable) <- c("rownames", adjPcolName)
  return(outputTable)
}

tt.RT.corm2flow.tadjP2 <- get_tadjP_Table2(fit.RT.corm, 2, "RT.corm2flow.tadjP")
identical(tt.RT.corm2flow.tadjP2, tt.RT.corm2flow.tadjP) #this also seems to work---------


#Continuing on below------------------
#have to check later is the coef comparison is actually the same as measuring diff tween 2 grps (vs diff to known mean)

#this function uses Fit to create a topTable comparing only 1 tissue to zero at a time
get_tadjP_Table2 <- function(inputFitTable, coefNum, adjPcolName) {
  tvalue_ttfit <- topTable(inputFitTable, coef = c(coefNum), n = Inf)
  outputTable  <- data.frame(rownames = rownames(tvalue_ttfit),
                             blarg1  = tvalue_ttfit[,5])
  colnames(outputTable) <- c("rownames", adjPcolName)
  return(outputTable)
}

#now, using our function
head(tt.RT.corm2flow.tadjP <- get_tadjP_Table2(fit.RT.corm, 2, "corm2flow"))
tt.RT.corm2leaf.tadjP <- get_tadjP_Table2(fit.RT.corm, 3, "corm2leaf")
tt.RT.corm2stem.tadjP <- get_tadjP_Table2(fit.RT.corm, 4, "corm2stem")
tt.RT.corm2stol.tadjP <- get_tadjP_Table2(fit.RT.corm, 5, "corm2stol")

tt.RT.flow2leaf.tadjP <- get_tadjP_Table2(fit.RT.flow, 2, "flow2leaf")
tt.RT.flow2stem.tadjP <- get_tadjP_Table2(fit.RT.flow, 3, "flow2stem")
tt.RT.flow2stol.tadjP <- get_tadjP_Table2(fit.RT.flow, 4, "flow2stol")

tt.RT.leaf2stem.tadjP <- get_tadjP_Table2(fit.RT.leaf, 2, "leaf2stem")
tt.RT.leaf2stol.tadjP <- get_tadjP_Table2(fit.RT.leaf, 3, "leaf2stol")

tt.RT.stem2stol.tadjP <- get_tadjP_Table2(fit.RT.stem, 2, "stem2stol")

#Then merge mini-dataframes by rowname--------------------
ttfit.RT.tadjP <- merge(tt.RT.corm2flow.tadjP, tt.RT.corm2leaf.tadjP, by = "rownames")
ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.corm2stem.tadjP, by = "rownames")
ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.corm2stol.tadjP, by = "rownames")

ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.flow2leaf.tadjP, by = "rownames")
ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.flow2stem.tadjP, by = "rownames")
ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.flow2stol.tadjP, by = "rownames")

ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.leaf2stem.tadjP, by = "rownames")
ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.leaf2stol.tadjP, by = "rownames")

ttfit.RT.tadjP <- merge(ttfit.RT.tadjP, tt.RT.stem2stol.tadjP, by = "rownames")

head(ttfit.RT.tadjP) # table of t values. #As you can see, they turn out all the same!

#then apply decidetests again.--------------
cutoff  <- 5e-7
ttfit.RT.tadjP.cut <- apply(ttfit.RT.tadjP[2:ncol(ttfit.RT.tadjP)], 
                                2, function(x) as.numeric(x<cutoff)) #these two lines encode the 'decision-making' part
rownames(ttfit.RT.tadjP.cut) <- ttfit.RT.tadjP$rownames #We might needs these. I'm not actually sure.

ttfit.RT.tadjP.cut.tsum <- as.data.frame(apply(ttfit.RT.tadjP.cut, 1, sum)) #this adds all the 1's in each row.
ttfit.RT.tadjP.cut.tsum <- cbind(rownames(ttfit.RT.tadjP.cut), ttfit.RT.tadjP.cut.tsum)
colnames(ttfit.RT.tadjP.cut.tsum) <- c("comps", "tsum")
head(ttfit.RT.tadjP.cut.tsum)
head(ttfit.RT.tadjP.cut, 50)

#NOTE: had to fix work because "Which" ONLY RETURNS INDICES!!! so the 'Comps" column may be unnecesary/
#Figure out amounts/types of each group---------------

length(genesWith10t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 10), ]))
length(genesWith9t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 9), ]))
length(genesWith8t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 8), ]))
length(genesWith7t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 7), ]))
length(genesWith6t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 6), ]))
length(genesWith5t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 5), ]))
length(genesWith4t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 4), ]))
length(genesWith3t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 3), ]))
length(genesWith2t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 2), ]))
length(genesWith1t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 1), ]))
length(genesWith0t <- rownames(ttfit.RT.tadjP.cut.tsum[which(ttfit.RT.tadjP.cut.tsum$tsum == 0), ]))

plotManyGenes(head(genesWith4t,30), "genesWith4t") #NOW IT WORKS!!!!!!!!!!!!!!
plotManyGenes(head(genesWith3t,30), "genesWith3t")
plotManyGenes(head(genesWith1t,30), "genesWith1t") #these expression patterns should have more than 1 t... Ohhh, it's the cutoff. 
#that is, expression patters with additional t's JUST below the cutoff would still end up in the "1t" section while really looking
#like a 2-t pattern.

#in short, to get well-defined patterns, I would have to pick define second threshold above which no other t's could go.
#but then, I SHOULD be able to search for well-defined expression patterns.
#next up: distinguishing between t's caused by proximity of mean to 0 (that is, to reference) and t's caused by large error bars.
#Also next up: define expression patterns using additional threshold requirements (as described above)
#Final next up: Input tissue and expression pattern, get back hit list.

#--------Below is more confused ramblings to help me pick up where I;ve left off-------------
#Now I'm getting the right hits, but sometimes things still don't seem to make sense. -----------
#sometimes expression patterns show up in the 4t section where I'm like "how did THAT get in there??"
#and other times I saw patterns in the 3t section and thought 'how is it that's NOT a 4th t??'
#I think the issue is that, in measuring t's, you get a small value for something precise but close to zero,
#and an equally small value for something farther but imprecise. So where I'm losing out... is where there's 
#a really disperse value - I'd like to be able to bin all 'tight' graphs in one section, and ones with loose genes in another.
#that would actually be rather helpful!
#what do I go by? t/FC? this value should be large if the t's small but precise, and large if the t's farther but disperse.
#or, what am I doing? just take the variance, silly. Or SD.
#you're basically going to end up imposing a 'characterization grid' on these genes (oh, references to philosophy of science abound...)
#one axis = thresholds of variance, other = # genes per tissue group past each threshold.
#one axis = thrsholds of t values, another = # genes per relation group past each threshold.

#so we could first sort the genes w/r to tightness (dump the genes that vary all over the place within tissue group)
#then we could sort those w/r to # of t relationships
#and then we could hope to "better bin" genes behaving similarly and differently.
#now, why do I care about this? does this achieve what CR wanted? 
