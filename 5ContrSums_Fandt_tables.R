#F_and_t tables of 5FitContrSumsLM.R--------------
#DESCRIPTION--------
#in this script, I will generate numerous summary tables using the fit.CS objects.

#QUESTIONS TO BE ANSWERED:
#First set: Grasping the data
#1. Are the t-values and associated adj.P.Values the same between all contrast sums parameterization fits? (answer: yes)
#2. Are the F and associated adj.P.Values the same between all contrast sums parameterization fits? (answer: No)

#Step_1: for each tissue, get the t values for each coefficient-----------------------------------------------
#i. Getting all t-values and adj.P.Val 's for corm ----------------------------------------------------------
#(i-a). make ttfit tables for corm (only) from each contrast sum fit (expect where corm was dummied out)
#fit.CS.stol coef order: "gndAvg", "corm", "flow", "leaf", "stem"
ttfit.CS.stol.corm <- topTable(fit.CS.stol, coef = c(2), n = Inf)
#fit.CS.stem coef order: "gndAvg", "stol", "corm", "flow", "leaf"
ttfit.CS.stem.corm <- topTable(fit.CS.stem, coef = c(3), n = Inf)
#fit.CS.leaf coef order: "gndAvg", "stem", "stol", "corm", "flow"
ttfit.CS.leaf.corm <- topTable(fit.CS.leaf, coef = c(4), n = Inf)
#fit.CS.flow coef order:"gndAvg", "leaf", "stem", "stol", "corm"
ttfit.CS.flow.corm <- topTable(fit.CS.flow, coef = c(5), n = Inf)

#(i-b-1): create dataframe with corm t-values
#First create mini-dataframes with t-values for corm and rownames from each contrast sum fit
column <- 3 #this is for the corm t.value (likelihood that corm average is not zero 
            #(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.stol.corm.t <- data.frame(rownames = rownames(ttfit.CS.stol.corm), CS.stol.corm.t=  ttfit.CS.stol.corm[,column])
ttfit.CS.stem.corm.t <- data.frame(rownames = rownames(ttfit.CS.stem.corm), CS.stem.corm.t = ttfit.CS.stem.corm[,column])
ttfit.CS.leaf.corm.t <- data.frame(rownames = rownames(ttfit.CS.leaf.corm), CS.leaf.corm.t = ttfit.CS.leaf.corm[,column])
ttfit.CS.flow.corm.t <- data.frame(rownames = rownames(ttfit.CS.flow.corm), CS.flow.corm.t = ttfit.CS.flow.corm[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.corm.t <- merge(ttfit.CS.stol.corm.t, ttfit.CS.stem.corm.t, by = "rownames")
ttfit.CS.corm.t <- merge(ttfit.CS.corm.t, ttfit.CS.leaf.corm.t, by = "rownames")
ttfit.CS.corm.t <- merge(ttfit.CS.corm.t, ttfit.CS.flow.corm.t, by = "rownames")
head(ttfit.CS.corm.t) # table of t values for corm. #As you can see, they turn out all the same!

#(i-b-2): create dataframe with adj.P.Val 's for corm t adj.P.Val's and rownames from each contrast sum fit
#First create mini-dataframes with adj.P.Vals for corm t-values and rownames from each contrast sum fit
column <- 5 #this is for the t adj. P val
ttfit.CS.stol.corm.tadjP <- data.frame(rownames = rownames(ttfit.CS.stol.corm), CS.stol.corm.tadjP = ttfit.CS.stol.corm[,column])
ttfit.CS.stem.corm.tadjP <- data.frame(rownames = rownames(ttfit.CS.stem.corm), CS.stem.corm.tadjP = ttfit.CS.stem.corm[,column])
ttfit.CS.leaf.corm.tadjP <- data.frame(rownames = rownames(ttfit.CS.leaf.corm), CS.leaf.corm.tadjP = ttfit.CS.leaf.corm[,column])
ttfit.CS.flow.corm.tadjP <- data.frame(rownames = rownames(ttfit.CS.flow.corm), CS.flow.corm.tadjP = ttfit.CS.flow.corm[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.corm.tadjP <- merge(ttfit.CS.stol.corm.tadjP, ttfit.CS.stem.corm.tadjP, by = "rownames")
ttfit.CS.corm.tadjP <- merge(ttfit.CS.corm.tadjP, ttfit.CS.leaf.corm.tadjP, by = "rownames")
ttfit.CS.corm.tadjP <- merge(ttfit.CS.corm.tadjP, ttfit.CS.flow.corm.tadjP, by = "rownames")
head(ttfit.CS.corm.tadjP) # table of t adj.P.Val's for corm. #as you can see, these also turned out all the same!

#ii. Getting all t-values and adj.P.Val 's for stolon ---------------------------------------------------------------

#(ii-a): Make ttfit tables for stolon only from each contrast sum fit (expect where stolom was dummied out)

#fit.CS.stem coef order: "gndAvg", "stol", "corm", "flow", "leaf"
ttfit.CS.stem.stol <- topTable(fit.CS.stem, coef = c(2), n = Inf)
#fit.CS.leaf coef order: "gndAvg", "stem", "stol", "corm", "flow"
ttfit.CS.leaf.stol <- topTable(fit.CS.leaf, coef = c(3), n = Inf)
#fit.CS.flow coef order: "gndAvg", "leaf", "stem", "stol", "corm"
ttfit.CS.flow.stol <- topTable(fit.CS.flow, coef = c(4), n = Inf)
#fit.CS.corm coef order: "gndAvg", "flow", "leaf", "stem", "stol"
ttfit.CS.corm.stol <- topTable(fit.CS.corm, coef = c(5), n = Inf)

#ii-b. Pull desired info form ttfit and associated adj.P.val tables for each stolon
#(ii-b-1): create dataframe with stolon t-values
#First create mini-dataframes with t-values for corm and rownames from each contrast sum fit
column <- 3 #this is for the stol t.value (likelihood that stol average is not zero 
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.stem.stol.t <- data.frame(rownames = rownames(ttfit.CS.stem.stol), CS.stem.stol.t = ttfit.CS.stem.stol[,column])
ttfit.CS.leaf.stol.t <- data.frame(rownames = rownames(ttfit.CS.leaf.stol), CS.leaf.stol.t = ttfit.CS.leaf.stol[,column])
ttfit.CS.flow.stol.t <- data.frame(rownames = rownames(ttfit.CS.flow.stol), CS.flow.stol.t = ttfit.CS.flow.stol[,column])
ttfit.CS.corm.stol.t <- data.frame(rownames = rownames(ttfit.CS.corm.stol), CS.corm.stol.t = ttfit.CS.corm.stol[,column])

#iThen merge mini-dataframes by rowname
ttfit.CS.stol.t <- merge(ttfit.CS.stem.stol.t, ttfit.CS.leaf.stol.t, by = "rownames")
ttfit.CS.stol.t <- merge(ttfit.CS.stol.t, ttfit.CS.flow.stol.t, by = "rownames")
ttfit.CS.stol.t <- merge(ttfit.CS.stol.t, ttfit.CS.corm.stol.t, by = "rownames")
head(ttfit.CS.stol.t)

#(ii-b-2): create dataframe with adj.P.Val 's for stol t-values and rownames from each contrast sum fit
#First create mini-dataframes with adj.P.Vals for stol t-values and rownames from each contrast sum fit
column <- 5 #this is for the t.value
ttfit.CS.stem.stol.tadjP <- data.frame(rownames = rownames(ttfit.CS.stem.stol), CS.stem.stol.tadjP = ttfit.CS.stem.stol[,column])
ttfit.CS.leaf.stol.tadjP <- data.frame(rownames = rownames(ttfit.CS.leaf.stol), CS.leaf.stol.tadjP = ttfit.CS.leaf.stol[,column])
ttfit.CS.flow.stol.tadjP <- data.frame(rownames = rownames(ttfit.CS.flow.stol), CS.flow.stol.tadjP = ttfit.CS.flow.stol[,column])
ttfit.CS.corm.stol.tadjP <- data.frame(rownames = rownames(ttfit.CS.corm.stol), CS.corm.stol.tadjP = ttfit.CS.corm.stol[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.stol.tadjP <- merge(ttfit.CS.stem.stol.tadjP, ttfit.CS.leaf.stol.tadjP, by = "rownames")
ttfit.CS.stol.tadjP <- merge(ttfit.CS.stol.tadjP, ttfit.CS.flow.stol.tadjP, by = "rownames")
ttfit.CS.stol.tadjP <- merge(ttfit.CS.stol.tadjP, ttfit.CS.corm.stol.tadjP, by = "rownames")
head(ttfit.CS.stol.tadjP)

#iii. Getting all t-values and adj.P.Val 's for stem -----------------------------------------------------------

#(iii-a). make ttfit tables for stem only from each contrast sum fit (expect where stolom was dummied out)
#fit.CS.leaf coef order: "gndAvg", "stem", "stol", "corm", "flow"
ttfit.CS.leaf.stem <- topTable(fit.CS.leaf, coef = c(2), n = Inf)
#fit.CS.flow coef order: "gndAvg", "leaf", "stem", "stol", "corm"
ttfit.CS.flow.stem <- topTable(fit.CS.flow, coef = c(3), n = Inf)
#fit.CS.corm coef order: "gndAvg", "flow", "leaf", "stem", "stol"
ttfit.CS.corm.stem <- topTable(fit.CS.corm, coef = c(4), n = Inf)
#fit.CS.stol coef order: "gndAvg", "corm", "flow", "leaf", "stem"
ttfit.CS.stol.stem <- topTable(fit.CS.stol, coef = c(5), n = Inf)

#iii-b. Pull desired info form ttfit and associated adj.P.val tables for each stem

#(iii-b-1): create dataframe with stem t-values
#First create mini-dataframes with t-values for stem and rownames from each contrast sum fit
column <- 3 #this is for the stem t.value (likelihood that stem average is not zero 
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.leaf.stem.t <- data.frame(rownames = rownames(ttfit.CS.leaf.stem), CS.leaf.stem.t = ttfit.CS.leaf.stem[,column])
ttfit.CS.flow.stem.t <- data.frame(rownames = rownames(ttfit.CS.flow.stem), CS.flow.stem.t = ttfit.CS.flow.stem[,column])
ttfit.CS.corm.stem.t <- data.frame(rownames = rownames(ttfit.CS.corm.stem), CS.corm.stem.t = ttfit.CS.corm.stem[,column])
ttfit.CS.stol.stem.t <- data.frame(rownames = rownames(ttfit.CS.stol.stem), CS.stol.stem.t = ttfit.CS.stol.stem[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.stem.t <- merge(ttfit.CS.leaf.stem.t, ttfit.CS.flow.stem.t, by = "rownames")
ttfit.CS.stem.t <- merge(ttfit.CS.stem.t, ttfit.CS.corm.stem.t, by = "rownames")
ttfit.CS.stem.t <- merge(ttfit.CS.stem.t, ttfit.CS.stol.stem.t, by = "rownames")
head(ttfit.CS.stem.t)

#(ii-b-2): create dataframe with adj.P.Val 's for stem t-values from each contrast sum fit

#First create mini-dataframes with adj.P.Vals for corm t-values and rownames from each contrast sum fit
column <- 5 #this is for the stem t.value adj. P value
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.leaf.stem.tadjP <- data.frame(rownames = rownames(ttfit.CS.leaf.stem), CS.leaf.stem.tadjP = ttfit.CS.leaf.stem[,column])
ttfit.CS.flow.stem.tadjP <- data.frame(rownames = rownames(ttfit.CS.flow.stem), CS.flow.stem.tadjP = ttfit.CS.flow.stem[,column])
ttfit.CS.corm.stem.tadjP <- data.frame(rownames = rownames(ttfit.CS.corm.stem), CS.corm.stem.tadjP = ttfit.CS.corm.stem[,column])
ttfit.CS.stol.stem.tadjP <- data.frame(rownames = rownames(ttfit.CS.stol.stem), CS.stol.stem.tadjP = ttfit.CS.stol.stem[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.stem.tadjP <- merge(ttfit.CS.leaf.stem.tadjP, ttfit.CS.flow.stem.tadjP, by = "rownames")
ttfit.CS.stem.tadjP <- merge(ttfit.CS.stem.tadjP, ttfit.CS.corm.stem.tadjP, by = "rownames")
ttfit.CS.stem.tadjP <- merge(ttfit.CS.stem.tadjP, ttfit.CS.stol.stem.tadjP, by = "rownames")
head(ttfit.CS.stem.tadjP)

#iv. Getting all t-values and adj.P.Val 's for leaf -----------------------------------------------------------

#(iv-a). make ttfit tables for leaf only from each contrast sum fit (expect where stolom was dummied out)
#fit.CS.flow coef order: "gndAvg", "leaf", "stem", "stol", "corm"
ttfit.CS.flow.leaf <- topTable(fit.CS.flow, coef = c(2), n = Inf)
#fit.CS.corm coef order: "gndAvg", "flow", "leaf", "stem", "stol"
ttfit.CS.corm.leaf <- topTable(fit.CS.corm, coef = c(3), n = Inf)
#fit.CS.stol coef order: "gndAvg", "corm", "flow", "leaf", "stem"
ttfit.CS.stol.leaf <- topTable(fit.CS.stol, coef = c(4), n = Inf)
#fit.CS.leaf coef order: "gndAvg", "stol", "corm", "flow", "leaf"
ttfit.CS.stem.leaf <- topTable(fit.CS.stem, coef = c(5), n = Inf)

#iv-b. Pull desired info form ttfit and associated adj.P.val tables for each leaf

#(iv-b-1): create dataframe with leaf t-values
#First create mini-dataframes with t-values for leaf and rownames from each contrast sum fit
column <- 3 #this is for the leaf t.value (likelihood that stem average is not zero 
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.flow.leaf.t <- data.frame(rownames = rownames(ttfit.CS.flow.leaf), CS.flow.leaf.t = ttfit.CS.flow.leaf[,column])
ttfit.CS.corm.leaf.t <- data.frame(rownames = rownames(ttfit.CS.corm.leaf), CS.corm.leaf.t = ttfit.CS.corm.leaf[,column])
ttfit.CS.stol.leaf.t <- data.frame(rownames = rownames(ttfit.CS.stol.leaf), CS.stol.leaf.t = ttfit.CS.stol.leaf[,column])
ttfit.CS.stem.leaf.t <- data.frame(rownames = rownames(ttfit.CS.stem.leaf), CS.stem.leaf.t = ttfit.CS.stem.leaf[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.leaf.t <- merge(ttfit.CS.flow.leaf.t, ttfit.CS.corm.leaf.t, by = "rownames")
ttfit.CS.leaf.t <- merge(ttfit.CS.leaf.t, ttfit.CS.stol.leaf.t, by = "rownames")
ttfit.CS.leaf.t <- merge(ttfit.CS.leaf.t, ttfit.CS.stem.leaf.t, by = "rownames")
head(ttfit.CS.leaf.t)

#(iv-b-2): create dataframe with adj.P.Val 's for leaf t-values from each contrast sum fit

#First create mini-dataframes with adj.P.Vals for leaf t-values and rownames from each contrast sum fit
column <- 5 #this is for the leaf t.value adj. P value
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.flow.leaf.tadjP <- data.frame(rownames = rownames(ttfit.CS.flow.leaf), CS.flow.leaf.tadjP = ttfit.CS.flow.leaf[,column])
ttfit.CS.corm.leaf.tadjP <- data.frame(rownames = rownames(ttfit.CS.corm.leaf), CS.corm.leaf.tadjP = ttfit.CS.corm.leaf[,column])
ttfit.CS.stol.leaf.tadjP <- data.frame(rownames = rownames(ttfit.CS.stol.leaf), CS.stol.leaf.tadjP = ttfit.CS.stol.leaf[,column])
ttfit.CS.stem.leaf.tadjP <- data.frame(rownames = rownames(ttfit.CS.stem.leaf), CS.stem.leaf.tadjP = ttfit.CS.stem.leaf[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.leaf.tadjP <- merge(ttfit.CS.flow.leaf.tadjP, ttfit.CS.corm.leaf.tadjP, by = "rownames")
ttfit.CS.leaf.tadjP <- merge(ttfit.CS.leaf.tadjP, ttfit.CS.stol.leaf.tadjP, by = "rownames")
ttfit.CS.leaf.tadjP <- merge(ttfit.CS.leaf.tadjP, ttfit.CS.stem.leaf.tadjP, by = "rownames")
head(ttfit.CS.leaf.tadjP)

#v. Getting all t-values and adj.P.Val 's for flow -----------------------------------------------------------

#(v-a). make ttfit tables for flower only from each contrast sum fit (expect where stolom was dummied out)

#fit.CS.corm coef order: "gndAvg", "flow", "leaf", "stem", "stol"
ttfit.CS.corm.flow <- topTable(fit.CS.corm, coef = c(2), n = Inf)
#fit.CS.stol coef order: "gndAvg", "corm", "flow", "leaf", "stem"
ttfit.CS.stol.flow <- topTable(fit.CS.stol, coef = c(3), n = Inf)
#fit.CS.leaf coef order: "gndAvg", "stol", "corm", "flow", "leaf"
ttfit.CS.stem.flow <- topTable(fit.CS.stem, coef = c(4), n = Inf)
#fit.CS.flow coef order: "gndAvg", "stem", "stol", "corm", "flow"
ttfit.CS.leaf.flow <- topTable(fit.CS.leaf, coef = c(5), n = Inf)

#v-b. Pull desired info form ttfit and associated adj.P.val tables for each flower

#(v-b-1): create dataframe with flower t-values
#First create mini-dataframes with t-values for flower and rownames from each contrast sum fit
column <- 3 #this is for the flower t.value (likelihood that stem average is not zero 
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.corm.flow.t <- data.frame(rownames = rownames(ttfit.CS.corm.flow), CS.corm.flow.t = ttfit.CS.corm.flow[,column])
ttfit.CS.stol.flow.t <- data.frame(rownames = rownames(ttfit.CS.stol.flow), CS.stol.flow.t = ttfit.CS.stol.flow[,column])
ttfit.CS.stem.flow.t <- data.frame(rownames = rownames(ttfit.CS.stem.flow), CS.stem.flow.t = ttfit.CS.stem.flow[,column])
ttfit.CS.leaf.flow.t <- data.frame(rownames = rownames(ttfit.CS.leaf.flow), CS.leaf.flow.t = ttfit.CS.leaf.flow[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.flow.t <- merge(ttfit.CS.corm.flow.t, ttfit.CS.stol.flow.t, by = "rownames")
ttfit.CS.flow.t <- merge(ttfit.CS.flow.t, ttfit.CS.stem.flow.t, by = "rownames")
ttfit.CS.flow.t <- merge(ttfit.CS.flow.t, ttfit.CS.leaf.flow.t, by = "rownames")
head(ttfit.CS.flow.t)

#(v-b-2): create dataframe with adj.P.Val 's for leaf t-values from each contrast sum fit

#First create mini-dataframes with adj.P.Vals for leaf t-values and rownames from each contrast sum fit
column <- 5 #this is for the flower t.value adj. P value
#(that is, in contr.sums not the grand average, and in ref+tx, not the same average as the reference))
ttfit.CS.corm.flow.tadjP <- data.frame(rownames = rownames(ttfit.CS.corm.flow), CS.corm.flow.tadjP = ttfit.CS.corm.flow[,column])
ttfit.CS.stol.flow.tadjP <- data.frame(rownames = rownames(ttfit.CS.stol.flow), CS.stol.flow.tadjP = ttfit.CS.stol.flow[,column])
ttfit.CS.stem.flow.tadjP <- data.frame(rownames = rownames(ttfit.CS.stem.flow), CS.stem.flow.tadjP = ttfit.CS.stem.flow[,column])
ttfit.CS.leaf.flow.tadjP <- data.frame(rownames = rownames(ttfit.CS.leaf.flow), CS.leaf.flow.tadjP = ttfit.CS.leaf.flow[,column])

#Then merge mini-dataframes by rowname
ttfit.CS.flow.tadjP <- merge(ttfit.CS.corm.flow.tadjP, ttfit.CS.stol.flow.tadjP, by = "rownames")
ttfit.CS.flow.tadjP <- merge(ttfit.CS.flow.tadjP, ttfit.CS.stem.flow.tadjP, by = "rownames")
ttfit.CS.flow.tadjP <- merge(ttfit.CS.flow.tadjP, ttfit.CS.leaf.flow.tadjP, by = "rownames")
head(ttfit.CS.flow.tadjP)

#Conclusion: 
#this was fairly useless, as these are all the same. T values and associated adj.P.values DO NOT CHANGE with dummying out tissues in contrasts.

#Step_2: Prepare a table of F-values of the five fits, along with a table of associated adj.P.Val's.------------------

#i. Prepare toptables for each tissue which exclude the 'grand average' coefficient----------------
#First: Giving useful labels to the coefficients, so I don't mess myself up later
ttfit.CS.stol <- topTable(fit.CS.stol, coef = c(2,3,4,5), n = Inf)
colnames(ttfit.CS.stol) <- c("corm", "flow", "leaf", "stem", "AveExpr", "F", "P.Value", "adj.P.Val")
ttfit.CS.stem <- topTable(fit.CS.stem, coef = c(2,3,4,5), n = Inf)
colnames(ttfit.CS.stem) <- c("stol", "corm", "flow", "leaf", "AveExpr", "F", "P.Value", "adj.P.Val")
ttfit.CS.leaf <- topTable(fit.CS.leaf, coef = c(2,3,4,5), n = Inf)
colnames(ttfit.CS.leaf) <- c("stem", "stol", "corm", "flow", "AveExpr", "F", "P.Value", "adj.P.Val")
ttfit.CS.flow <- topTable(fit.CS.flow, coef = c(2,3,4,5), n = Inf)
colnames(ttfit.CS.flow) <- c("leaf", "stem", "stol", "corm", "AveExpr", "F", "P.Value", "adj.P.Val")
ttfit.CS.corm <- topTable(fit.CS.corm, coef = c(2,3,4,5), n = Inf)
colnames(ttfit.CS.corm) <- c("flow", "leaf", "stem", "stol", "AveExpr", "F", "P.Value", "adj.P.Val")

#ii. Prepare a table of F-values from the toptables of the five fits----------------------------------------
#First: make mini-dataframes with F-values and rownames from each contrast sum fit.
column <- 6 #this is for the F.value.
ttfit.CS.stol.F <- data.frame(rownames = rownames(ttfit.CS.stol), CS.stol.F =  ttfit.CS.stol[,column])
ttfit.CS.stem.F <- data.frame(rownames = rownames(ttfit.CS.stem), CS.stem.F = ttfit.CS.stem[,column])
ttfit.CS.leaf.F <- data.frame(rownames = rownames(ttfit.CS.leaf), CS.leaf.F = ttfit.CS.leaf[,column])
ttfit.CS.flow.F <- data.frame(rownames = rownames(ttfit.CS.flow), CS.flow.F = ttfit.CS.flow[,column])
ttfit.CS.corm.F <- data.frame(rownames = rownames(ttfit.CS.corm), CS.corm.F = ttfit.CS.corm[,column])

#Third: merge mini-dataframes by rownames.
ttfit.CS.all.F <- merge(ttfit.CS.stol.F, ttfit.CS.stem.F, by = "rownames")
ttfit.CS.all.F <- merge(ttfit.CS.all.F, ttfit.CS.leaf.F, by = "rownames")
ttfit.CS.all.F <- merge(ttfit.CS.all.F, ttfit.CS.flow.F, by = "rownames")
ttfit.CS.all.F <- merge(ttfit.CS.all.F, ttfit.CS.corm.F, by = "rownames")
head(ttfit.CS.all.F)

#iii. Prepare a table of F adj.P.Val's (one column for each contrast sum fit)----------------
#First: make mini-dataframes with F adj.P.Val's and rownames from each contrast sum fit.
column <- 8 #this is for the F.value.
ttfit.CS.stol.FadjP <- data.frame(rownames = rownames(ttfit.CS.stol), CS.stol.FadjP =  ttfit.CS.stol[,column])
ttfit.CS.stem.FadjP <- data.frame(rownames = rownames(ttfit.CS.stem), CS.stem.FadjP = ttfit.CS.stem[,column])
ttfit.CS.leaf.FadjP <- data.frame(rownames = rownames(ttfit.CS.leaf), CS.leaf.FadjP = ttfit.CS.leaf[,column])
ttfit.CS.flow.FadjP <- data.frame(rownames = rownames(ttfit.CS.flow), CS.flow.FadjP = ttfit.CS.flow[,column])
ttfit.CS.corm.FadjP <- data.frame(rownames = rownames(ttfit.CS.corm), CS.corm.FadjP = ttfit.CS.corm[,column])

#Second: merge mini-dataframes by rownames.
ttfit.CS.all.FadjP <- merge(ttfit.CS.stol.FadjP, ttfit.CS.stem.FadjP, by = "rownames")
ttfit.CS.all.FadjP <- merge(ttfit.CS.all.FadjP, ttfit.CS.leaf.FadjP, by = "rownames")
ttfit.CS.all.FadjP <- merge(ttfit.CS.all.FadjP, ttfit.CS.flow.FadjP, by = "rownames")
ttfit.CS.all.FadjP <- merge(ttfit.CS.all.FadjP, ttfit.CS.corm.FadjP, by = "rownames")
head(ttfit.CS.all.FadjP)

#Step 3: Make a single summary table each of t-values and t adj.P.Val's---------------
##i. Prepare a table of t-values from the toptables of the five fits-----------------
#First: make mini-dataframes with t-values from a column in each table
column <- 3
ttfit.CS.stol.corm.t <- data.frame(rownames = rownames(ttfit.CS.stol.corm), CS.corm.t = ttfit.CS.stol.corm[,column])
ttfit.CS.stem.stol.t <- data.frame(rownames = rownames(ttfit.CS.stem.stol), CS.stol.t = ttfit.CS.stem.stol[,column])
ttfit.CS.leaf.stem.t <- data.frame(rownames = rownames(ttfit.CS.leaf.stem), CS.stem.t = ttfit.CS.leaf.stem[,column])
ttfit.CS.flow.leaf.t <- data.frame(rownames = rownames(ttfit.CS.flow.leaf), CS.leaf.t = ttfit.CS.flow.leaf[,column])
ttfit.CS.corm.flow.t <- data.frame(rownames = rownames(ttfit.CS.corm.flow), CS.flow.t = ttfit.CS.corm.flow[,column])

#second: merge mini-dataframes by rownames.
ttfit.CS.all.t <- merge(ttfit.CS.stol.corm.t, ttfit.CS.stem.stol.t, by = "rownames")
ttfit.CS.all.t <- merge(ttfit.CS.all.t, ttfit.CS.leaf.stem.t, by = "rownames")
ttfit.CS.all.t <- merge(ttfit.CS.all.t, ttfit.CS.flow.leaf.t, by = "rownames")
ttfit.CS.all.t <- merge(ttfit.CS.all.t, ttfit.CS.corm.flow.t, by = "rownames")
head(ttfit.CS.all.t)

#ii. Prepare a table of t adj.P.Val's from the toptables of the five fits--------------------------------------
#First: make mini-dataframes with t adj.P.Val's from a column in each table
column <- 5
ttfit.CS.stol.corm.tadjP <- data.frame(rownames = rownames(ttfit.CS.stol.corm), CS.corm.tadjP = ttfit.CS.stol.corm[,column])
ttfit.CS.stem.stol.tadjP <- data.frame(rownames = rownames(ttfit.CS.stem.stol), CS.stol.tadjP = ttfit.CS.stem.stol[,column])
ttfit.CS.leaf.stem.tadjP <- data.frame(rownames = rownames(ttfit.CS.leaf.stem), CS.stem.tadjP = ttfit.CS.leaf.stem[,column])
ttfit.CS.flow.leaf.tadjP <- data.frame(rownames = rownames(ttfit.CS.flow.leaf), CS.leaf.tadjP = ttfit.CS.flow.leaf[,column])
ttfit.CS.corm.flow.tadjP <- data.frame(rownames = rownames(ttfit.CS.corm.flow), CS.flow.tadjP = ttfit.CS.corm.flow[,column])

ttfit.CS.all.tadjP <- merge(ttfit.CS.stol.corm.tadjP, ttfit.CS.stem.stol.tadjP, by = "rownames")
ttfit.CS.all.tadjP <- merge(ttfit.CS.all.tadjP, ttfit.CS.leaf.stem.tadjP, by = "rownames")
ttfit.CS.all.tadjP <- merge(ttfit.CS.all.tadjP, ttfit.CS.flow.leaf.tadjP, by = "rownames")
ttfit.CS.all.tadjP <- merge(ttfit.CS.all.tadjP, ttfit.CS.corm.flow.tadjP, by = "rownames")
head(ttfit.CS.all.tadjP)

#OBJECTS TO BE PASSED ON TO THE NEXT ANALYSIS:-------------------
#ttfit.CS.all.F
#ttfit.CS.all.FadjP
#ttfit.CS.all.t
#ttfit.CS.all.tadjP

#and of course, dat.voomed, design, tissue 
