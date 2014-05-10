#INSTRUCTIONS:----------
#0. SET YOUR WORKING DIRECTORY. Go to the 'Session' tab, go to 'Choose directory', 
              #and select the directory that contains your raw counts data matrix.
#1.Go to the bottom of the page. 
#2.Replace the genes in the list with your own genes of interest.
#3.Replace the plot title with something informative.
#4.Run the code!
#your plots should appear in the plot window. You can export them if you like.

#PREPARATION--------
#Step 1. Load packages---------
library(edgeR)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#Step 2. Load data matrix---------
predat <- read.table("C_x_c_C100_gene.count.matrix_20130314",
                     header = TRUE, row.names = 1)

#Step 3. making tissue groups levels of a factor.------------
#this is not essential for this script, but I'm trying not to chnage everything between scripts.
tissue <- factor(c(rep("corm", 4), rep("flower", 3,), rep("leaf", 3), rep("stem", 3), rep("stolon", 3)))
str(tissue)

#Step 4. using tissue groups to make a design matrix, used here by voom.----------
#Not actually important for our goal in using voom, but just avoiding changing between scripts.
design <- model.matrix(~0 + tissue)
show(design)

#Step 5. Making a design matrix. It will be used by the plotting function below.-----------
des <- as.data.frame(list(tissue = c(rep("corm", 4), rep("flower", 3,), rep("leaf", 3), rep("stem", 3), rep("stolon", 3)),
                          sample = c("corm2", "corm3", "corm4", "corm7", "flower2", "flower4", "flower7", "leaf1", "leaf2", "leaf3", 
                                     "stem2","stem4", "stem7", "stolon2", "stolon4", "stolon7")))

#Step 6. Tidy(make smaller) column names and row names in dat--------
str(predat)
dat <- predat
colnames(dat) <- str_split_fixed(colnames(predat), "_", 6)[, 6]
rownames(dat) <- str_split_fixed(rownames(predat), "_", 6)[, 6]
str(dat)

#Step 7. Remove genes in dat which have too little counts in all genes---------
#(Making a DGElist pretty much so I can use Mac's code to trim)
#NOTE: This removes a LOT of genes! but with such low counts across tissues, they will only add noise.
y <- DGEList(counts=dat, group = tissue)  #This makes your data group into a DGE list. 
z <- y[(rowSums(cpm(y) > 1) >= 3), ]      #This means: for a given gene (row), only keep rows which have at least 3 entries with >1 count. 
z$samples$lib.size <- colSums(z$counts)   #Resets the library size (adjusting for removed rows). Possibly unnecessary for fitting with limma.
dat <- z$counts                           #This takes your data out of the DGE list (I think)
str(dat)

#Step 8. Calculate the normalization factors for dat. This calculates 1 value per library, 
norm.factor <- calcNormFactors(dat)

#Step 9. Vooming dat, but for this script we only do it change from counts to logCpm values.------------
dat.voomed <- voom(dat, lib.size = colSums(dat) * norm.factor)
head(dat.voomed$E) #this contains the log2(counts) expression data that gets plotted

#Step 10. plotManyGenes is a function plots multiple genes ----------
#NOTE: For now, if you want to change the plot scale, you must alter the line:
#     scale_x_continuous(limits=c(-5, 13)) +
#The -5 is what I have put for the minimum, and 13 is the max. you can change those as you need.

plotManyGenes <- function(GENELIST, TITLE_OF_PLOT) {
  toplist <- cbind(GENELIST)
  str(reDat <- as.data.frame(dat.voomed$E[toplist,]))
  cDat <- cbind(des, t(reDat))
  
  prMDat <- melt(cDat,
                 id.vars=c('tissue', 'sample'),
                 variable.name='gene', value.name='gExp')
  
  (p <- ggplot(prMDat, aes(gExp, tissue)) + 
     geom_point() +
     geom_point(position = position_jitter(height = 0.1)) +
     facet_wrap(~ gene) +
     scale_x_continuous(limits=c(-5, 13)) +
     ggtitle(TITLE_OF_PLOT))
}

#Step 11. plotOneGene is a function that plots a single gene------------------
#NOTE: It differs from the above function in that it doesnt' use cbind() to make a toplist
# and it doesnt' transpose reDat when making cDat.
#   toplist <- cbind(GENELIST)
#   cDat <- cbind(des, t(reDat))
#NOTE: For now, if you want to change the plot scale, you must alter the line:
#     scale_x_continuous(limits=c(-5, 13)) +
#The -5 is what I have put for the minimum, and 13 is the max. you can change those as you need.

plotOneGene <- function(GENE, TITLE_OF_PLOT) {
  toplist <- cbind(GENE)
  reDat <- as.data.frame(dat.voomed$E[toplist,])
  colnames(reDat) <- GENE
  cDat <- cbind(des, reDat)
  
  prMDat <- melt(cDat,
                 id.vars=c('tissue', 'sample'),
                 variable.name='gene', value.name='gExp')
  
  (p <- ggplot(prMDat, aes(gExp, tissue)) + 
     geom_point() +
     geom_point(position = position_jitter(height = 0.1)) +
     facet_wrap(~ gene) +
     scale_x_continuous(limits=c(-5, 13)) +
     ggtitle(TITLE_OF_PLOT))
}

#INPUT YOUR GENELIST---------------
GENELIST <- c("comp103267_c0_seq1", "comp102969_c0_seq7","comp103431_c0_seq1", 
                 "comp85415_c0_seq1", "comp99171_c0_seq1", "comp105809_c1_seq1", 
                 "comp90891_c0_seq2", "comp103135_c0_seq8", "comp87326_c0_seq1", 
                 "comp87326_c0_seq3", "comp107131_c0_seq9")
#or
GENE <- "comp103267_c0_seq1"

#ARE ALL YOUR GENES PRESENT IN THE MATRIX?-----
#If they don't, the plotter won't work! 
#If you do have some extra genes, the setdiff line should give you thier names.
length(setdiff(GENELIST, rownames(dat.voomed))) #gives number of gene names not contained in the filtered counts matrix
setdiff(GENELIST, rownames(dat.voomed))
#or
length(setdiff(GENE, rownames(dat.voomed)))
setdiff(GENE, rownames(dat.voomed))

#PLOT YOUR GENES---------------
plotManyGenes(GENELIST, "TITLE_OF_PLOT")
#or
plotOneGene(GENE, "TITLE_OF_PLOT")

