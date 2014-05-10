#DESCRIPTION-----------
#REMEMBER TO SET YOUR WORKING DIRECTORY--------
#this prepares the data matrix for whatever method of parameterization you wish to then use (separate scripts)
#this uses the same code as Plot_list_of_genes, but some of the notes are different,
#and preparation of the design matrix (step 4 in Plot_list_of_genes.R) is not done here.
#(as that is parameterization-specific)

#Step 1. Load packages---------
#possibly not all of this is necessary.
library(edgeR)
library(limma)
library(lattice)
library(ggplot2)
library(data.table) #for renaming my function
library(DESeq)
library(reshape2)
library(stringr)
library(RColorBrewer)

#Step 2. Load data matrix (which will be referrec to as dat after cleaning)---------
predat <- read.table("C_x_c_C100_gene.count.matrix_20130314",
                     header = TRUE, row.names = 1)

#Step 3. Make a factor with tissue groups as levels. Used in Voom and optionally in making a DGElist------------
tissue <- factor(c(rep("corm", 4), rep("flower", 3,), rep("leaf", 3), rep("stem", 3), rep("stolon", 3)))
str(tissue)

#Step 4. Make a design matrix. It will be used by the plotting function below.-----------
des <- as.data.frame(list(tissue = c(rep("corm", 4), rep("flower", 3,), rep("leaf", 3), rep("stem", 3), rep("stolon", 3)),
                          sample = c("corm2", "corm3", "corm4", "corm7", "flower2", "flower4", "flower7", "leaf1", "leaf2", "leaf3", 
                                     "stem2","stem4", "stem7", "stolon2", "stolon4", "stolon7")))

#Step 5. Tidy (make smaller) column names and row names in dat--------
str(predat)
dat <- predat
colnames(dat) <- str_split_fixed(colnames(predat), "_", 6)[, 6]
rownames(dat) <- str_split_fixed(rownames(predat), "_", 6)[, 6]
str(dat)

#Step 6. Remove genes in dat which have too little counts in all genes---------
#(Making a DGElist pretty much so I can use Mac's code to trim)
#NOTE: This removes a LOT of genes! but with such low counts across tissues, they will only add noise.
y <- DGEList(counts=dat, group = tissue)  #This makes your data group into a DGE list. 
z <- y[(rowSums(cpm(y) > 1) >= 3), ]      #This means: for a given gene (row), only keep rows which have at least 3 entries with >1 count. 
z$samples$lib.size <- colSums(z$counts)   #Resets the library size (adjusting for removed rows). Possibly unnecessary for fitting with limma.
dat <- z$counts                           #This takes your data out of the DGE list (I think)
str(dat)

#Step 7. Calculate the normalization factors for dat-------------
#This calculates 1 value per library,
#and more details on what they actually do are in A_scaling_normalization_method_for_differential_expression_(TMM)
#in the protocol package.
norm.factor <- calcNormFactors(dat)

#OPEN SCRIPT WITH DESIRED PARAMETERIZATION AND CONTINUE-----------
#DO NOT clear your environment, or R will 'forget' all the things it's just calculated.
#Choices are:
    #cell means
    #reference + treatment effects
    #contr.sums

