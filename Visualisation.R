## A principal component approach to analysing the 1000 Genome Project data through R
## data analysis and visualisations
devtools::install_github
library(devtools)             ## needs running bfore vqv/ggbiplot - 
install.packages("installr")
library(installr)
install.packages("stringr")
library(stringr)
install.packages("htmltab")  ## not needed after reinstallation of Rstudio/R
library(htmltab)            ## not needed after reinstallation of Rstudio/R
install.Rtools(TRUE)
install.packages("usethis") ## not needed after reinstallation of Rstudio/R// again required next day
library(usethis)
install.packages("processx")      ## needed to remove error for devtools /## not needed after reinstallation of Rstudio/R
library(processx)                ## needs running bfore vqv/ggbiplot
install.packages("ggplot2")
install.packages("reshape2")
install.packages("R.utils")

install_github("vqv/ggbiplot", force = TRUE)
install.packages("backports")   ## needed for ggbplot - all, but no compilation needed
library(backports) ##pass
library(usethis)
library(devtools)   ##run library usethis, backports once again before calling devtools

install_github("vqv/ggbiplot")  ##pass, option 1, yes - all incl compile
library(ggbiplot)   ##plyr , scales, grid required
install.packages("plyr")
install.packages("scales")
install.packages("grid")
library(scales)
library(plyr)
library(grid)
library(ggbiplot) ##needed library load for ggplot2,plyr,scales
library(ggplot2)
library(ggbiplot) ##successful, after doing all installations and library loading
library(reshape2)  ##pass
library(vcfR)   ##pass

install.packages(rCharts)
install.packages("devtools")  ## To install vcfR package, needs compilation - yes
##load usethis first and then run this command
library(usethis)
library(devtools)            ##pass after loading usethis
#install_github("vqv/ggbiplot")  ##error, keep erroring for backports/htmltab,usethis etc
#library(ggbiplot)    # error
library('R.utils')    ##pass
library(vcfR)         ##pass
### ALL REQD PACKAGES LOADED

##----------------------------------------------------------------------------------------------------------
#### INSTALL RTOOLS AND CHECKING THE PATH AND MAKE
# download and install Rtools from https://cran.rstudio.com/bin/windows/Rtools
#----
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron") #creates env file
# restart Rstudio at this point
Sys.which("make") ## to verify if Rstudio is able to identify make application
# if this works, install jsonlite
install.packages("jsonlite", type = "source") # to verify if rtools/make is working 

## assign working directory - default - [1] "C:/Users/Owner/Documents"
getwd()
setwd("C:/Users/Owner/Documents")
getwd()
#setwd("D:/GENOMICS/1000genome/ftp download")

find_rtools(T) ## true
Sys.getenv()['PATH']  ## output values
Sys.getenv("PATH")    ## output values
Sys.which("ls.exe")   ## output values
Sys.which("gcc.exe")  ## output values

memory.size(max=FALSE)  ## check memory size - success!!
memory.limit()   ##check memory limit - success!!
memory.limit(size=56000)  ## increase memory limit - success!!

install.Rtools(choose_version = TRUE)
#-------------
# https://github.com/r-lib/devtools/issues/1772
library(devtools)
## the below statement gave error- but till find_rtools - came as true
assignInNamespace("version_info", c(devtools:::version_info, list("4.0" = list(version_min = "3.3.0", version_max = "99.99.99", path = "bin"))), "devtools")
find_rtools() # is TRUE on 21 sep 2020 13:03


##-------------------------------------------------------------------------------------------------------
# To read panel file --  PASS!!
Popdata <- read.table("D:/GENOMICS/1000genome/DataSlicerFiles/New folder/integrated_call_samples_v3.20130502.ALL.panel", header = TRUE)
colnames(Popdata)
dimension(Popdata)               ## 2504 X 4
Popdata[1:10, ]
typeof(Popdata)                  ## list

##--------------------------------------------------
#TO CLEAR ALL DATA OBJECTS CREATED FROM GLOBAL ENV
rm(list = ls(all.names = TRUE)) ## to clear all data objects created
##----------------------------------------------------------------------------
##------------------------------------------------------------------------------------------
##dataset         Chromosome Snips samples  Base positions     population details
##data set 1 DS1- CHR22 - 78104   975   17500001-20000000   10 populations (ACB,CHS,ESN,GBR,GIH,IBS,ITU,JPT,MXL,PUR)
##data set 2 DS2- CHR22 - 48784  2504   25000001 - 26500000  all 26 populations 
##data set 3 DS3- CHR20 - 29323  1165   2500001 - 3500000   AFR, EUR populations 
##data set 4 DS4- CHR20 - 29323  1497   2500001 - 3500000   AFR, AMR, SAS populations 
##data set 5 DS5- CHR20 - 29323  1496   2500001 - 3500000   EUR, SAS, EAS populations 
##data set 6 DS6- CHR20 - 29323  836    2500001 - 3500000   AMR, SAS populations 
##data set 7 DS7- CHR20 - 29323  494    2500001 - 3500000   Yoruba, Kenya, Nigeria, BEB, ITU populations (pure dataset) 
##data set 8 DS9- CHR22 - 48784  393    25000001 - 26500000  ACB, ESN, GBR and IBS (pure dataset) 

## NAMING convention followed
## vcf file = starts with DSx_population_details , first file is DS1_allpop as it contains all populations
## list file created after genotype (GT) extract = vcf filename appended with m at the start = mDS1_allpop
##  name of pca file = pcaDSx , first file name = pcaDS1
##  combined file with panel file name = DSx , first file name = DS1

## *********************************************************************************************
##   DS1 :: dataset 1 = 10 population file DS_10pop chrom 22 ----------------------------------------
library(vcfR)
DS_10pop <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/10 pop/10pop_22.17500001-20000000.ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS_10pop)), decreasing = TRUE)[1:10]    ## NO duplicates
dimension(DS_10pop)   ##78104 X 976 

##===Data Analysis ----- VCF file ------------------
getFIX(DS_10pop) %>% head      # output the fixed region of the data
getPOS(DS_10pop) %>% head      # output the position of reference genome
getID(DS_10pop) %>% head       # output the reference ID
getREF(DS_10pop) %>% head      # output the Reference allele information
getALT(DS_10pop) %>% head      # output the alternate allele information
getQUAL(DS_10pop) %>% head     # output the quality score of the data sample
getFILTER(DS_10pop) %>% head   # output the filter value
getINFO(DS_10pop) %>% head     # output the info, and this file consists only the AA (ancestral allele info)
dimension(DS_10pop) %>% head   # output the dimention of the list
typeof(DS_10pop)               # type of variable created = S4
head(DS_10pop)                # Provides info from meta, fix and GT region of the VCF file


## CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS22_10pop <- t(extract.gt(DS_10pop, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))
dimension(mDS22_10pop)     ## 975 X 78104

## PCA -----------------
pcaDS22_10pop <- prcomp(mDS22_10pop) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS22_10pop)                    ## 5
summary(pcaDS22_10pop$x[,1:5])             ## Lists min,max, 1,3 Q, mean and median
str(pcaDS22_10pop)                          ## 78104 X 975
sapply(pcaDS22_10pop, mode)

## DS1 :: Let us find the BARPLOT of percent variation represented by each PCs
pca1.var     <- pcaDS22_10pop$sdev^2 ##Finds how much variation in the orig. data each PC accounts for!!
pca1.var.per <- round(pca1.var/sum(pca1.var)*100,1) ##Finds % of variation each PC accounts for!!
pca1.var[1:10]
dimension(pca3.var.per) ##1497
pca1.var.per[1:200]
pca1.var.per[1:20]  ## decide the limit of PC selection based on <1 criteria

##------------------------------------------------------------------------------------------
# another method to find the proportional covariance
DS22_10pop_pov <- summary(pcaDS22_10pop)$importance[2,] ## The same pov can be calculated through this formula also
DS22_10pop_pov[1:5]                                 ## displays the first 300 pov values
typeof(DS22_10pop_pov)                              ## double
dimension(DS22_10pop_pov)                           ## 975
## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca1.var.per[1:20], main = "DS1 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca1.var.per[1:20], main = "DS1 :: Bar plot to show the coverage of \n each principal component, first 20 ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(15,5,"PC1 = 5.2%, PC2 = 3.4%", col = "blue", font = 2))

##------------------------------------------------------------------------------------------
# DS1:: Plotting for '10 pop dataset' - 10 pop population file - variance explained for each principal component
##-----------------------------------------------------------------------
## DS1 :: PLOT - PCA coverage ---------------------------
plot(pcaDS22_10pop, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, ylim=c(1,80),
     main="DS1:: PCA variant coverage graph   \n in DS1 dataset, all 10 population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))
## ~80% of variants is covered by first 2 Principal components

plot(DS22_10pop_pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS1:: PoV - Proportion of Variance Explained graph \n in DS1 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2)

# DS1 :: Plot - cumulative proportion of variance explained
DS22_10pop_y = sum(DS22_10pop_pov[1:20])
sum(DS22_10pop_y)
plotds22_10pop <- plot((cumsum(DS22_10pop_pov[1:20])), xlab = "Principal Components ",
                       ylab = "Cumulative PoV Explained  ",
                       main="DS1:: CPV - Cumulative Proportion of Variance Explained \n in DS1 dataset for first 20 samples ", 
                       col.main="blue", font.main=3, cex.main=1.2,
                       ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)


## FACTOMINER for scree plot ---------------------------------------------------

install.packages("FactoMineR") ##pass
library(FactoMineR)            ## pass
pcaDS22_10pop_FM <- PCA(mDS22_10pop, scale.unit = TRUE, graph = FALSE)  ##scale=true added
#res.pca <- PCA(mDS2, graph = FALSE)                      ## no scale parameter
print(pcaDS22_10pop_FM)               ##pass

install.packages("factoextra")
library("factoextra")
DS22_10pop_eigval <- get_eigenvalue(pcaDS22_10pop_FM)
DS22_10pop_eigval  ## 3.2% variations are explained by first Dim1 eigenvalue and 2.2% by Dim2 eigenvalue
typeof(DS22_10pop_eigval)        ## double
dimension(DS22_10pop_eigval)     ## 974 X 3

##scree plot
fviz_eig(pcaDS22_10pop_FM, addlabels = TRUE, ylim = c(0, 3))  #3% covered by first 2 PCs
##Total 0f 8% variations are represented by first 10 principal components
##------------------------------------------------------------------------------------------
## create a combined file DS1 with GT data and population panel file--------------------
library(R.utils)
pcaDS22_10pop$x <- cbind(as.data.frame(pcaDS22_10pop$x), "sample" = rownames(pcaDS22_10pop$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS22_10pop$x[1:5,]
dimension(pcaDS22_10pop$x)                 ## 975 X 975+1

DS22_10pop <- merge(pcaDS22_10pop$x[,c(1:6,976)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS22_10pop)               ## Output matrix 975 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS22_10pop)                 ## list
DS22_10pop[1:5,1:10]

str(pcaDS1_all)
pcaDS1_all$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS1_all$center                   ##provides mean of the variables
pcaDS1_all$scale                    ## provides Std Dev of the variables
dimension(pcaDS22_10pop$rotation)      ## 78104 X 975

#------------ PLOTTING - 10 population file DS22_10pop -------------------------------------------
library(ggplot2)
##------------------------------------------------------
## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = African populations-------------
DS22_10pop <- DS22_10pop[DS22_10pop$super_pop %in% c("AFR", "AMR", "EUR", "SAS", "EAS"),]
library(ggplot2)
ggplot(DS22_10pop, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.1 : Population distribution for all 10 population
             from 10 population dataset ")

##-------------------------------------------------------------------
##------------------------------------------------------
## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = African populations-------------
DS22_10pop_AFR <- DS22_10pop[DS22_10pop$super_pop %in% c("AFR"),]
library(ggplot2)
ggplot(DS22_10pop_AFR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.2 : Population distribution for African population
             from 10 population dataset ")

##------------------------------------------------------------------

## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = Europe populations-------------
DS22_10pop_EUR <- DS22_10pop[DS22_10pop$super_pop %in% c("EUR"),]
library(ggplot2)
ggplot(DS22_10pop_EUR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.3 : Population distribution for Europan population
             from 10 population dataset ")

#--------------------------------------------------------------------
## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = East Asian populations-------------
##------------------------------------------------------
DS22_10pop_EAS <- DS22_10pop[DS22_10pop$super_pop %in% c("EAS"),]
library(ggplot2)
ggplot(DS22_10pop_EAS, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.4 : Population distribution for East Asian population
             from 10 population dataset ")

##------------------------------------------------------
## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = American populations-------------
DS22_10pop_AMR <- DS22_10pop[DS22_10pop$super_pop %in% c("AMR"),]
library(ggplot2)
ggplot(DS22_10pop_AMR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.5 : Population distribution for American population
             from 10 population dataset ")

##------------------------------------------------------
## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = South Asian populations-------------
DS22_10pop_SAS <- DS22_10pop[DS22_10pop$super_pop %in% c("SAS"),]
library(ggplot2)
ggplot(DS22_10pop_SAS, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.6 : Population distribution for South Asian population
             from 10 population dataset ")

DS22_10pop_AFREUR <- DS22_10pop[DS22_10pop$super_pop %in% c("AFR", "EUR"),]
##------------------------------------------------------
## PLOT FOR POPULATION ; data set 1 = 10 population, sub set = East Asian populations-------------
library(ggplot2)
ggplot(DS22_10pop_AFREUR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle(" Chr22:: DS1.7 : Population distribution for African/European population
             from 10 population dataset ")

##------------------------------------------------------
## PLOT FOR POPULATION-----10 population file --- DS22_10pop ---------
library(ggplot2)
ggplot(DS22_10pop_AFREUR, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS1.8 : Super Population distribution for African/European population
             from 10 population dataset ")

##---------------------------------------------------------------------------------------
##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS1 :: Let us find which genes have the largest effect , finding top 10 genes 
##DS1 :: loading score calc for PC1
pca1loading <- pcaDS22_10pop$rotation[,1] #as PC1 accounts for max variation of data
##---------------------------------------------------------------------------------------
# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca1loading)   #29323
DS1genes <- abs(pca1loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS1genes[1:25] ## displays loading scores of all genes
DS1gene_rank <- sort(DS1genes, decreasing = TRUE) ## sorts the genes based on score
DS1gene_rank[1:25]  ##displays the ranked genes
DS1top_10 <- names(DS1gene_rank[1:10]) ## gets the names of top 10 genes and display
DS1top_10 ## displays top 10 genes affecting the PC1 for variants in DS1
pcaDS2$rotation[DS2top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##------------------------------------------------------------------------------
##DS1 :: loading score calc for PC2
pca1loading1 <- pcaDS22_10pop$rotation[,2] #as PC2 accounts for next max variation of data
##---------------------------------------------------------------------------------------
# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca1loading1)   #29323
DS1genes1 <- abs(pca1loading1) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS1genes1[1:25] ## displays loading scores of all genes
DS1gene1_rank <- sort(DS1genes1, decreasing = TRUE) ## sorts the genes based on score
DS1gene1_rank[1:25]  ##displays the ranked genes
DS1top1_10 <- names(DS1gene1_rank[1:10]) ## gets the names of top 10 genes and display
DS1top1_10 ## displays top 10 genes affecting the PC1 for variants in DS1
pcaDS2$rotation[DS1top1_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left
##------------------------------------------------------------------------------
##******* end of 10 population  ************  file analysis --------------------------------
##------------------------------------------------------------------------------------------
## DS2 :: DATASET 2 - all population data --------line 277 - --------------------------------------
##------------------------------------------------------------------------------------------
library(vcfR)
DS2_allpop <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/all pop/22.25000001-26500000.ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS2_allpop)), decreasing = TRUE)[1:10]    ##  NO duplicates found

dimension(DS2_allpop)   ##48784 X 2505

## CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS2_allpop <- t(extract.gt(DS2_allpop, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS2_allpop)         ## 2504 X 48784

##   PCA = ALL population file ----------------------------------------

pcaDS2 <- prcomp(mDS2_allpop) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS2)                    ## 5
summary(pcaDS2$x[,1:5])             ## Lists min,max, 1,3 Q, mean and median
head(pcaDS2)
str(pcaDS2)                          ## 48784 X 2504
sapply(pcaDS2, mode)

## DS2 :: Let us find the BARPLOT of percent variation represented by each PCs
pca2.var     <- pcaDS2$sdev^2 ##Finds how much variation in the orig. data each PC accounts for!!
pca2.var.per <- round(pca2.var/sum(pca2.var)*100,1) ##Finds % of variation each PC accounts for!!
pca2.var[1:10]
dimension(pca2.var.per) ##2504
pca2.var.per[1:20]  ## decide the limit of PC selection based on <1 criteria
##------------------------------------------------------------------------------------------
# another method to find the proportional covariance
pca2_pov <- summary(pcaDS2)$importance[2,] ## The same pov can be calculated through this formula also
pca2_pov[1:5]                                 ## displays the first 300 pov values

## Barplot for % of variations represented by each pc,
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca2.var.per[1:20], main = "DS2 :: Bar plot to show the coverage of \n each principal component, first 20
        for 'all population' dataset from chr22 ", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca2.var.per[1:20], main = "DS2 :: Bar plot to show the coverage of \n each principal component, first 20 
        for 'all population' dataset from chr22 ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(15,5,"PC1 = 6.6%, PC2 = 3.8%", col = "blue", font = 2))

##------------------------------------------------------------------------------------------
# DS2:: Plotting for 'all population dataset' - variance explained for each principal component
##-----------------------------------------------------------------------
## DS2 :: PLOT - PCA coverage ---------------------------
plot(pcaDS2, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, ylim=c(1,80),
     main="DS2:: PCA variant coverage graph   \n in DS2 dataset, all 26 population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))
## ~66% of variants is covered by first Principal component

plot(pca2_pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS2:: PoV - Proportion of Variance Explained graph \n in DS2 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2)

# DS2 :: Plot - cumulative proportion of variance explained
DS2_y = sum(pca2_pov[1:20])
sum(DS2_y)
plotds2_all <- plot((cumsum(pca2_pov[1:20])), xlab = "Principal Components ",
                       ylab = "Cumulative PoV Explained  ",
                       main="DS2:: CPV - Cumulative Proportion of Variance Explained \n in DS2 dataset for first 20 samples ", 
                       col.main="blue", font.main=3, cex.main=1.2,
                       ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)


## FACTOMINER for scree plot ---------------------------------------------------

install.packages("FactoMineR") ##pass
library(FactoMineR)            ## pass
pcaDS2_FM <- PCA(mDS2_allpop, scale.unit = TRUE, graph = FALSE)  ##scale=true added
#res.pca <- PCA(mDS2, graph = FALSE)                      ## no scale parameter
print(pcaDS2_FM)               ##pass

install.packages("factoextra")
library("factoextra")
DS2_eigval <- get_eigenvalue(pcaDS2_FM)
DS2_eigval  ## 3.2% variations are explained by first Dim1 eigenvalue and 2.2% by Dim2 eigenvalue
typeof(DS2_eigval)        ## double
dimension(DS2_eigval)     ## 2504 X 3

##scree plot
fviz_eig(pcaDS2_FM, addlabels = TRUE, ylim = c(0, 3))  #3% covered by first 2 PCs
##Total 0f 8% variations are represented by first 10 principal components
##------------------------------------------------------------------------------------------

## create a combined file DS2 with GT data and population panel file--------------------

pcaDS2$x <- cbind(as.data.frame(pcaDS2$x), "sample" = rownames(pcaDS2$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS2$x[1:5,1:10]
dimension(pcaDS2$x)                 ## 2504 X 2504+1 = last column is added by statement above (code at 261)

DS2 <- merge(pcaDS2$x[,c(1:6,2505)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS2)               ## Output matrix 2504 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
DS2[1:5,]
str(pcaDS2)
pcaDS2$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS2$center                   ##provides mean of the variables
pcaDS2$scale                    ## provides Std Dev of the variables
dimension(pcaDS2$rotation)      ## 48784 X 2504

#------------ PLOTTING dataset 2 all pop file --------------------------------------------
##-----------------------
## data = DS2
##-----------------------
## PLOT FOR POPULATION ; data set 2.1 = all population, Plotting all populations-----------
library(ggplot2)
ggplot(DS2, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all
           26 population from all population dataset DS2 ")
##-----------------------
DS2_allviz <- DS2[DS2$super_pop %in% c("EUR", "AFR", "EAS", "SAS", "AMR"),]
## PLOT FOR POPULATION ; data set 2.2 = all population, Plotting all populations-----------
ggplot(DS2_allviz, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all 5
           super population from all population dataset DS2 ")

##-----------------------
DS2_allviz <- DS2[DS2$super_pop %in% c("EUR", "AFR", "EAS", "SAS", "AMR"),]
## PLOT FOR POPULATION ; data set 2.3 = all population, Plotting all populations-----------
ggplot(DS2_allviz, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20:: DS2 - all 5 super pop : Population distribution for all 5
           super population from all population DS2 dataset ")

##-----------------------
DS2_AfrEur <- DS2[DS2$pop %in% c("ACB", "ESN", "GBR","IBS"),]
## PLOT FOR POPULATION ; data set 2.2 = all population, Plotting all populations-----------
ggplot(DS2_AfrEur, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for ACB, ESN, GBR and IBS 
            population from all population dataset DS2 ")

##-----------------------
DS2_AfEu <- DS2[DS2$super_pop %in% c("AFR", "EUR"),]
## PLOT FOR POPULATION ; data set 2.2 = all population, Plotting all populations-----------
ggplot(DS2_AfEu, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for African and European 
            population from all population dataset DS2 ")

##-----------------------

DS2_AFR <- DS2[DS2$super_pop %in% c("AFR"),]
## PLOT FOR POPULATION ; data set 2.4 = all population, Plotting all populations-----------
ggplot(DS2_AFR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all African
           super population from all population DS2 dataset ")

##-----------------------
DS2_AMR <- DS2[DS2$super_pop %in% c("AMR"),]
## PLOT FOR POPULATION ; data set 2.5 = all population, Plotting all populations-----------
ggplot(DS2_AMR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all American
           super population from all population DS2 dataset ")

##-----------------------
DS2_EUR <- DS2[DS2$super_pop %in% c("EUR"),]
## PLOT FOR POPULATION ; data set 2.6 = all population, Plotting all populations-----------
ggplot(DS2_EUR, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all European
           super population from all population DS2 dataset ")

##-----------------------
DS2_EAS <- DS2[DS2$super_pop %in% c("EAS"),]
## PLOT FOR POPULATION ; data set 2.7 = all population, Plotting all populations-----------
ggplot(DS2_EAS, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all East Asian
           super population from all population DS2 dataset ")

##-----------------------
DS2_SAS <- DS2[DS2$super_pop %in% c("SAS"),]
## PLOT FOR POPULATION ; data set 2.8 = all population, Plotting all populations-----------
ggplot(DS2_SAS, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS2 - all 26 pop : Population distribution for all South Asian
           super population from all population DS2 dataset ")
##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS2 :: Let us find which genes have the largest effect , finding top 10 genes 
##DS2 :: loading score calc for PC1
pca2loading <- pcaDS2$rotation[,1] #as PC1 accounts for max variation of data
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca2loading)   #29323
DS2genes <- abs(pca2loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS2genes[1:25] ## displays loading scores of all genes
DS2gene_rank <- sort(DS2genes, decreasing = TRUE) ## sorts the genes based on score
DS2gene_rank[1:25]  ##displays the ranked genes
DS2_top_10 <- names(DS2gene_rank[1:10]) ## gets the names of top 10 genes and display
DS2_top_10 ## displays top 10 genes affecting the PC1 for variants in DS2
pcaDS2$rotation[DS2_top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left
##------------------------------------------------------------------------------
##DS2 :: loading score calc for PC2
pca2loading1 <- pcaDS2$rotation[,2] #as PC2 accounts for next max variation of data
##---------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------
# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca2loading1)   #29323
DS2genes1 <- abs(pca2loading1) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS2genes1[1:25] ## displays loading scores of all genes
DS2gene1_rank <- sort(DS2genes1, decreasing = TRUE) ## sorts the genes based on score
DS2gene1_rank[1:25]  ##displays the ranked genes
DS2_top1_10 <- names(DS2gene1_rank[1:10]) ## gets the names of top 10 genes and display
DS2_top1_10 ## displays top 10 genes affecting the PC1 for variants in DS2
pcaDS2$rotation[DS2_top1_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left
##------------------------------------------------------------------------------
###************************************************************************************************
###*##--------------------------------------------------------------------------------------------
#------------ POPULATION dataset 3 : African/European population -------------------------------------------
##  DS3 :: dataset - 3 Africa, Europe (12 pop)
library(vcfR)
library(R.utils)
#DS3_AfEu <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/1-afr-eur/20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS3_AfEu)), decreasing = TRUE)[1:5]    ## 4,2,1,1 duplicates
#rs369887745 rs111364731  esv3644998  esv3645015   
#4           2           1           1 
DS3_AfEu <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/1-afr-eur/cleanAfEu_20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
dimension(DS3_AfEu)                     ## 29323 X 1165

# CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS3 <- t(extract.gt(DS3_AfEu, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS3)         ## 1164 X 29323

##-----  PCA -------------------------------------------------------------------------------------
pcaDS3 <- prcomp(mDS3) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS3) 
str(pcaDS3)                          ## 29323 X 1164

# another method to find the proportional covariance
DS3pov <- summary(pcaDS3)$importance[2,] ## The same pov can be calculated through this formula also
DS3pov[1:5]                                 ## displays the first 5 pov values
typeof(DS3pov)                              ## double
dimension(DS3pov)                           ## 1164
##------------------------------------------------------------------------------------------
# Plot dataset 3 - variance explained for each principal component
## DS3 :: PLOT - PCA coverage ---------------------------
plot(pcaDS3, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, ylim=c(1,80),
     main="DS3:: PCA variant coverage graph   \n in DS3 dataset, AFR and EUR population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))

## 4.4% of variants is covered by first 2 Principal components

plot(DS3pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS3:: PoV - Proportion of Variance Explained graph \n in DS3 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2)

# Plot - cumulative proportion of variance explained

DS3y = sum(DS3pov[1:20])
sum(DS3y)
plotds2 <- plot((cumsum(DS3pov[1:20])), xlab = "Principal Components ",
                ylab = "Cumulative PoV Explained  ",
                main="DS3:: CPV - Cumulative Proportion of Variance Explained \n in DS3 dataset for first 20 samples ", 
                col.main="blue", font.main=3, cex.main=1.2,
                ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)

## DS3 :: Let us find the BARPLOT of percent variation represented by each PCs
pca3.var     <- pcaDS3$sdev^2 ##Finds how much variation in the orig. data each PC accounts for!!
pca3.var.per <- round(pca3.var/sum(pca3.var)*100,1) ##Finds % of variation each PC accounts for!!
pca3.var[1:10]
dimension(pca3.var.per) ##1164
pca3.var.per[1:200]
pca3.var.per[1:20]  ## decide the limit based on <1 criteria

## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca3.var.per[1:20], main = "DS3 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca3.var.per[1:20], main = "DS3 :: Bar plot to show the coverage of \n each principal component, first 20,
        AFR , EUR dataset ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(20,8,"PC1 = ~10%, PC2 = ~7%", col = "blue", font = 2))

##-----------------------------------------------------------------------------------
## 17% of variations are covered by first 2 Principal components
##------------------------------------------------------------------------------------------
## create a combined file DS1 with GT data and population panel file--------------------
pcaDS3$x <- cbind(as.data.frame(pcaDS3$x), "sample" = rownames(pcaDS3$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS3$x[,1:10]
dimension(pcaDS3$x)                 ## 1164X 1165

DS3 <- merge(pcaDS3$x[,c(1:6,1165)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS3)               ## Output matrix 1164 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS3)                 ## list
DS3[1:5,]
str(pcaDS3)
pcaDS3$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS3$center                   ##provides mean of the variables
pcaDS3$scale                    ## provides Std Dev of the variables
dimension(pcaDS3$rotation)      ## 29323 X 1164

#------------ PLOTTING dataset 3  --------------------------------------------
DS3_EurAfr <- DS3[DS3$super_pop %in% c("EUR", "AFR"),]  ## 1164 obs
##------------------------------------------------------------------------------------------
## DS3 :: PLOT FOR POPULATION ; data set 3.1 = EUR/AFR population, populations-------------
library(ggplot2)
ggplot(DS3_EurAfr, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS3 - : Population distribution for African and European 
            population from Afr/Eur super population dataset DS3 ")

## Scree plot ------------------------------------
library("FactoMineR")
DS3.pca <- PCA(mDS3, graph = FALSE)
DS3_ev <- get_eigenvalue(DS3.pca)
DS3_ev
fviz_eig(DS3.pca, addlabels = TRUE, ylim = c(0, 50))


##------------------------------------------------------------------------------------------
DS3_EurAfr1 <- DS3[DS3$pop %in% c("ACB", "ESN", "GBR", "IBS"),]  ## 1164 obs
## --------------------
## PLOT FOR POPULATION- pure data set 3.1(a)- ACB, ESN, GBR and IBS population ONLY -------------

library(ggplot2)
ggplot(DS3_EurAfr1, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS3 - : Population distribution for ACB, ESN, GBR and IBS 
            population from Afr/Eur super population dataset DS3 ")

##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS3 :: Let us find which genes have the largest effect , finding top 10 genes 
##DS3 :: loading score calc for PC1
pca3loading <- pcaDS3$rotation[,1] #as PC1 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca3loading)   #29323
DS3genes <- abs(pca3loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS3genes[1:25] ## displays loading scores of all genes
DS3gene_rank <- sort(DS3genes, decreasing = TRUE) ## sorts the genes based on score
DS3gene_rank[1:25]  ##displays the ranked genes
DS3_top_10 <- names(DS3gene_rank[1:10]) ## gets the names of top 10 genes and display
DS3_top_10 ## displays top 10 genes affecting the PC1 for variants in DS3
pcaDS3$rotation[DS3_top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##DS3 :: loading score calc for PC2
pca3loading1 <- pcaDS3$rotation[,2] #as PC2 accounts for next max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca3loading1)   #29323
DS3genes1 <- abs(pca3loading1) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS3genes1[1:25] ## displays loading scores of all genes
DS3gene1_rank <- sort(DS3genes1, decreasing = TRUE) ## sorts the genes based on score
DS3gene1_rank[1:25]  ##displays the ranked genes
DS3_top1_10 <- names(DS3gene1_rank[1:10]) ## gets the names of top 10 genes and display
DS3_top1_10 ## displays top 10 genes affecting the PC1 for variants in DS3
pcaDS3$rotation[DS3_top1_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##********************end of DS3 analysis ----------------------------------------------------
##------------------------------------------------------------------------------------------
## DS4 :: DATASET 4 - AFRICA, AMERICA, SOUTH ASIAN population data -------------------------------------------------------------
##----------------------------------------------------------------------------------------------
#dataset - 4 Africa, America, South Asian (16 pop)
#DS4_AfAmEa <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/2-afr-amr-sas/20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
DS4_AfAmSa <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/2-afr-amr-sas/clean_20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS4_AfAmSa)), decreasing = TRUE)[1:5]    ## 4,2,1,1 duplicates
#rs369887745 rs111364731  esv3644998  esv3645015   
#4           2           1           1
dimension(DS4_AfAmSa)                    ## 29323 X 1498 

## CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS4_AfAmSa <- t(extract.gt(DS4_AfAmSa, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS4_AfAmSa)         ## 1497 X 29323

##-----  PCA  -------------------------------------------------------------------------------------
pcaDS4_AfAmSa <- prcomp(mDS4_AfAmSa) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS4_AfAmSa) 
str(pcaDS4_AfAmSa)                          ## 29323 X 1497

pca4.var     <- pcaDS4_AfAmSa$sdev^2
pca4.var.per <- round(pca4.var/sum(pca4.var)*100,1)
pca4.var[1:10]
dimension(pca4.var.per) ##1497
pca4.var.per[1:200]
# another method to find the proportional covariance
DS4pov <- summary(pcaDS4_AfAmSa)$importance[2,] ## The same pov can be calculated through this formula also
DS4pov[1:25]                                 ## displays the first 25 pov values

## ---------------------------------------------------------------------------------
## DS4 :: PLOT - PCA coverage ---------------------------
plot(pcaDS4_AfAmSa, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, 
     main="DS4:: PCA variant coverage graph   \n in DS4 dataset, AFR,AMR,EUR population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))
## ~58 % of variants is covered by PC1 and ~40% covered by PC2

## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca4.var.per[1:20], main = "DS4 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca4.var.per[1:20], main = "DS4 :: Bar plot to show the coverage of \n each principal component, first 20,
        AFR , AMR , EUR dataset ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(20,8,"PC1 = ~10%, PC2 = ~7%", col = "blue", font = 2))

##-----------------------------------------------------------------------------------
#biplot(pcaDS4_AfAmSa, scale = 0)          ## plots biplot graph, takes time
##-----------------------------------------------------------------------------------
# DS4 :: Plot dataset 4 - variance explained for each principal component
plot(DS4pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS4:: PoV - Proportion of Variance Explained graph \n in DS4 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2 )
#----------------------------------
# DS4 :: Plot - cumulative proportion of variance explained
DS4y = sum(DS4pov[1:20])
sum(DS4y)
plotds2 <- plot((cumsum(DS4pov[1:20])), xlab = "Principal Components ",
                ylab = "Cumulative PoV Explained  ",
                main="DS4:: CPV - Cumulative Proportion of Variance Explained \n in DS4 dataset for first 20 samples ", 
                col.main="blue", font.main=3, cex.main=1.2,
                ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)
##------------------------------------------------------------------------------------------
## create a combined file DS4 with GT data and population panel file--------------------
pcaDS4_AfAmSa$x <- cbind(as.data.frame(pcaDS4_AfAmSa$x), "sample" = rownames(pcaDS4_AfAmSa$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS4_AfAmSa$x[,1:10]
dimension(pcaDS4_AfAmSa$x)                 ##1497 X 1498

DS4 <- merge(pcaDS4_AfAmSa$x[,c(1:6,1498)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS4)               ## Output matrix 1497 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS4)                 ## list
DS4[1:5,]
str(pcaDS4_AfAmSa)
pcaDS4_AfAmSa$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS4_AfAmSa$center                   ##provides mean of the variables
pcaDS4_AfAmSa$scale                    ## provides Std Dev of the variables
dimension(pcaDS4_AfAmSa$rotation)      ## 29323 X 1497

#- DS4 :: PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## DS4 :: PLOT FOR POPULATION ; data set 4 = AFR/AMR/SAS population, populations-------------
library(ggplot2)
ggplot(DS4, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20: DS4 : Population distribution for AFR, AMR & SAS
           from DS4 population dataset ")

#-- DS4 ---------- PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## PLOT FOR POPULATION ; data set 4 = AFR/AMR/SAS population, populations-------------
library(ggplot2)
ggplot(DS4, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20:: DS4 : Population distribution for AFR, AMR & SAS
            from DS4 population dataset ")

#DS4_CYNT <- DS4[DS4$pop %in% c("CEU", "YRI", "TSI"),]
## PLOT FOR POPULATION ; data set 4 = AFR/AMR/SAS population, populations-------------
## Utah, Yorbuba in Nigeria, Utah and TSI analysis
library(ggplot2)
ggplot(DS4_CYNT, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20::DS4 : Population distribution for Yoruba in Nigeria,
           Utah and West Europe TSI from DS4 population dataset ")

##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS4 :: Let us find top 10 measurements (genes) that contribute to PC1 
##DS4 :: loading score calc for PC1
pca4loading <- pcaDS4_AfAmSa$rotation[,1] #as PC1 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca4loading)   #29323
DS4genes <- abs(pca4loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS4genes[1:25] ## displays loading scores of first 25 genes
DS4gene_rank <- sort(DS4genes, decreasing = TRUE) ## sorts the genes based on score
DS4gene_rank[1:25]  ##displays the ranked genes, first 25
DS4_top_10 <- names(DS4gene_rank[1:10]) ## gets the names of top 10 genes and display
DS4_top_10 ## displays top 10 genes affecting the PC1 for variants in DS4
pcaDS4_AfAmSa$rotation[DS4_top_10,1] ## see which genes has + and which has _ loading score pushing
                                     ## the samples to right or left

## DS4 :: Let us find top 10 measurements (genes) that contribute to PC2
##DS4 :: loading score calc for PC2
pca4loading_2 <- pcaDS4_AfAmSa$rotation[,2] #as PC2 accounts for next max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca4loading_2)   #29323
DS4genes_2 <- abs(pca4loading_2) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS4genes_2[1:25] ## displays loading scores of first 25 genes
DS4gene_rank_2 <- sort(DS4genes_2, decreasing = TRUE) ## sorts the genes based on score
DS4gene_rank_2[1:25]  ##displays the ranked genes, first 25
DS4_top_10_2 <- names(DS4gene_rank_2[1:10]) ## gets the names of top 10 genes and display
DS4_top_10_2 ## displays top 10 genes affecting the PC1 for variants in DS4
pcaDS4_AfAmSa$rotation[DS4_top_10_2,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##scree plot-----------------------------------------------------------
library("FactoMineR")
#DS2.pca <- PCA(mDS2, graph = FALSE) ## pass, takes time
#DS2_ev <- get_eigenvalue(DS2.pca)
#DS2_ev
#fviz_eig(DS2.pca, addlabels = TRUE, ylim = c(0, 50))
### --------------------------------------------------------------------------------------
##-------END OF DS4 DATA ANALYSIS --------------------------------------------------------
##******************************************************************************************
library(vcfR)
library(R.utils)
#dataset - 5 East Asian , 5 South Asian & , 5 East Asian ,5 Europe (15 pop)
DS5_ <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/3-eas-eur-sas/20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
DS5_EES <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/3-eas-eur-sas/clean_20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS5_EES)), decreasing = TRUE)[1:5]    ## 4,2,1,1 duplicates
#rs369887745 rs111364731  esv3644998  esv3645015   
#4           2           1           1
dimension(DS5_EES)                      ## 29323 X 1497

## CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS5_EES <- t(extract.gt(DS5_EES, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS5_EES)         ## 1496 X 29323

##-----  PCA  -------------------------------------------------------------------------------------
pcaDS5_EES <- prcomp(mDS5_EES) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS5_EES) 
str(pcaDS5_EES)                          ## 29323 X 1497

# another method to find the proportional covariance
DS5pov <- summary(pcaDS5_EES)$importance[2,] ## The same pov can be calculated through this formula also
DS5pov[1:25]                                 ## displays the first 25 pov values

## ---------------------------------------------------------------------------------

## DS5 :: PLOT - PCA coverage ---------------------------
plot(pcaDS5_EES, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, 
     main="DS5:: PCA variant coverage graph   \n in DS5 dataset, EUR, EAS AND SAS population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))
## 47 % of variants is covered by PC1 and ~37% covered by PC2

pca5.var     <- pcaDS5_EES$sdev^2
pca5.var.per <- round(pca5.var/sum(pca5.var)*100,1)
pca5.var[1:10]
dimension(pca5.var.per) ##1496
pca5.var.per[1:200]

## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca5.var.per[1:22], main = "DS5 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca5.var.per[1:22], main = "DS5 :: Bar plot to show the coverage of \n each principal component, first 20,
        AFR , AMR , EUR dataset ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(20,8,"PC1 = ~10%, PC2 = ~8%", col = "blue", font = 2))

##-----------------------------------------------------------------------------------
#biplot(pcaDS5_AfAmSa, scale = 0)          ## plots biplot graph, takes time
##-----------------------------------------------------------------------------------
# DS5 :: Plot dataset 5 - variance explained for each principal component
plot(DS5pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS5:: PoV - Proportion of Variance Explained graph \n in DS5 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2 )
#----------------------------------
# DS5 :: Plot - cumulative proportion of variance explained
DS5y = sum(DS5pov[1:20])
sum(DS5y)
plotds2 <- plot((cumsum(DS5pov[1:20])), xlab = "Principal Components ",
                ylab = "Cumulative PoV Explained  ",
                main="DS5:: CPV - Cumulative Proportion of Variance Explained \n in DS5 dataset for first 20 samples ", 
                col.main="blue", font.main=3, cex.main=1.2,
                ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)
##------------------------------------------------------------------------------------------
## create a combined file DS5 with GT data and population panel file--------------------
pcaDS5_EES$x <- cbind(as.data.frame(pcaDS5_EES$x), "sample" = rownames(pcaDS5_EES$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS5_EES$x[,1:10]
dimension(pcaDS5_EES$x)                 ## 1496 X 1497

DS5 <- merge(pcaDS5_EES$x[,c(1:6,1497)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS5)               ## Output matrix 1496 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS5)                 ## list
DS5[1:5,]
str(pcaDS5_EES)
pcaDS5_EES$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS5_EES$center                   ##provides mean of the variables
pcaDS5_EES$scale                    ## provides Std Dev of the variables
dimension(pcaDS5_EES$rotation)      ## 29323 X 1496

#- DS5 :: PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## DS5 :: PLOT FOR POPULATION ; data set 5 = EAS/EUR/SAS population, populations-------------
library(ggplot2)
ggplot(DS5, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20 DS5 : Population distribution for EUR, EAS and SAS
           super population from pure DS5 population dataset ")

#-- DS5 ---------- PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## PLOT FOR POPULATION ; data set set 5 = EAS/EUR/SAS population, populations-------------
library(ggplot2)
ggplot(DS5, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20 DS5 : Population distribution for EAS, EUR & SAS
           super population from pure DS5 population dataset ")

## find population dist between Utah, Yoruba in Nigeria and WEst Eur Toscani in Italia
DS5_ChsChbFin <- DS5[DS5$pop %in% c("CHS", "CHB", "FIN"),]

## PLOT FOR POPULATION ; data set set 5 = EAS/EUR/SAS population, populations-------------
## Utah, Yorbuba in Nigeria, Utah and TSI analysis
library(ggplot2)
ggplot(DS5_ChsChbFin, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20 DS5 : Population distribution for CHS, CHB, FIN,
            from DS5 population dataset ")

##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS5 :: Let us find top 10 measurements (genes) that contribute to PC1 
##DS5 :: loading score calc for PC1
pca5loading <- pcaDS5_AfAmSa$rotation[,1] #as PC1 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca5loading)   #29323
DS5genes <- abs(pca5loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS5genes[1:25] ## displays loading scores of first 25 genes
DS5gene_rank <- sort(DS5genes, decreasing = TRUE) ## sorts the genes based on score
DS5gene_rank[1:25]  ##displays the ranked genes, first 25
DS5_top_10 <- names(DS5gene_rank[1:10]) ## gets the names of top 10 genes and display
DS5_top_10 ## displays top 10 genes affecting the PC1 for variants in DS5
pcaDS5_AfAmSa$rotation[DS5_top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

## DS5 :: Let us find top 10 measurements (genes) that contribute to PC2
##DS5 :: loading score calc for PC2
pca5loading_2 <- pcaDS5_AfAmSa$rotation[,2] #as PC2 accounts second max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca5loading_2)   #29323
DS5genes_2 <- abs(pca5loading_2) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS5genes_2[1:25] ## displays loading scores of first 25 genes
DS5gene_rank_2 <- sort(DS5genes_2, decreasing = TRUE) ## sorts the genes based on score
DS5gene_rank_2[1:25]  ##displays the ranked genes, first 25
DS5_top_10_2 <- names(DS5gene_rank_2[1:10]) ## gets the names of top 10 genes and display
DS5_top_10_2 ## displays top 10 genes affecting the PC1 for variants in DS5
pcaDS5_AfAmSa$rotation[DS5_top_10_2,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##scree plot-----------------------------------------------------------
library("FactoMineR")
#DS5.pca <- PCA(mDS5, graph = FALSE)
#DS5_ev <- get_eigenvalue(DS5.pca)
#DS5_ev
#fviz_eig(DS5.pca, addlabels = TRUE, ylim = c(0, 50))
### --------------------------------------------------------------------------------------
##---------END of DS5 data analysis ------------------------------------------------------
##******************************************************************************************

#dataset - 6 America South Asian (9 pop)
#DS6_AmSa <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/4-amr-sas/20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
DS6_AmSa <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/4-amr-sas/AmrSas_20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS6_AmSa)), decreasing = TRUE)[1:5]    ## 4,2,1,1  duplicates
dimension(DS6_AmSa)  ## 29323 X 837

## CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS6_AmSa <- t(extract.gt(DS6_AmSa, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS6_AmSa)         ## 836 X 29323

##-----  PCA  -------------------------------------------------------------------------------------
pcaDS6_AmSa <- prcomp(mDS6_AmSa) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS6_AmSa) 
str(pcaDS6_AmSa)                          ## 29323 X 837

# another method to find the proportional covariance
DS6pov <- summary(pcaDS6_AmSa)$importance[2,] ## The same pov can be calculated through this formula also
DS6pov[1:25]                                 ## displays the first 25 pov values

## ---------------------------------------------------------------------------------

## DS6 :: PLOT - PCA coverage ---------------------------
plot(pcaDS6_AmSa, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, 
     main="DS6:: PCA variant coverage graph   \n in DS6 dataset, AMR and SAS population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))
## ~47 % of variants is covered by PC1 and ~39% covered by PC2; from pca6.var[1:10]

pca6.var     <- pcaDS6_AmSa$sdev^2
pca6.var.per <- round(pca6.var/sum(pca6.var)*100,1)
pca6.var[1:10]
dimension(pca6.var.per) ##836
pca6.var.per[1:200]

## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca6.var.per[1:20], main = "DS6 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca6.var.per[1:20], main = "DS6 :: Bar plot to show the coverage of \n each principal component, first 20,
        AFR , AMR , EUR dataset ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(20,8,"PC1 = ~9.7%, PC2 = 8%", col = "blue", font = 2))

##-----------------------------------------------------------------------------------
#biplot(pcaDS6_AfAmSa, scale = 0)          ## plots biplot graph, takes time
##-----------------------------------------------------------------------------------
# DS6 :: Plot dataset 6 - variance explained for each principal component
plot(DS6pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS6:: PoV - Proportion of Variance Explained graph \n in DS6 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2 )
#----------------------------------
# DS6 :: Plot - cumulative proportion of variance explained
DS6y = sum(DS6pov[1:20])
sum(DS6y)
plotds2 <- plot((cumsum(DS6pov[1:20])), xlab = "Principal Components ",
                ylab = "Cumulative PoV Explained  ",
                main="DS6:: CPV - Cumulative Proportion of Variance Explained \n in DS6 dataset for first 20 samples ", 
                col.main="blue", font.main=3, cex.main=1.2,
                ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)
##------------------------------------------------------------------------------------------
## create a combined file DS6 with GT data and population panel file--------------------
pcaDS6_AmSa$x <- cbind(as.data.frame(pcaDS6_AmSa$x), "sample" = rownames(pcaDS6_AmSa$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS6_AmSa$x[,1:10]
dimension(pcaDS6_AmSa$x)                 ## 836 X 837

DS6 <- merge(pcaDS6_AmSa$x[,c(1:6,837)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS6)               ## Output matrix 836 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS6)                 ## list
DS6[1:5,]
str(pcaDS6_AmSa)
pcaDS6_AmSa$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS6_AmSa$center                   ##provides mean of the variables
pcaDS6_AmSa$scale                    ## provides Std Dev of the variables
dimension(pcaDS6_AmSa$rotation)      ## 29323 X 836

#- DS6 :: PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## DS6 :: PLOT FOR POPULATION ; data set 6 = AMR/SAS population, populations-------------
library(ggplot2)
ggplot(DS6, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20: DS6 : Population distribution for  AMR & SAS
           super population from pure DS6 population dataset ")

#-- DS6 ---------- PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## PLOT FOR POPULATION ; data set 6 = AMR/SAS population, populations-------------
library(ggplot2)
ggplot(DS6, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20: DS6 : Population distribution for  AMR & SAS
           super population from pure DS6 population dataset ")

##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS6 :: Let us find top 10 measurements (genes) that contribute to PC1 
##DS6 :: loading score calc for PC1
pca6loading <- pcaDS6_AmSa$rotation[,1] #as PC1 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca6loading)   #29323
DS6genes <- abs(pca6loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS6genes[1:25] ## displays loading scores of first 25 genes
DS6gene_rank <- sort(DS6genes, decreasing = TRUE) ## sorts the genes based on score
DS6gene_rank[1:25]  ##displays the ranked genes, first 25
DS6_top_10 <- names(DS6gene_rank[1:10]) ## gets the names of top 10 genes and display
DS6_top_10 ## displays top 10 genes affecting the PC1 for variants in DS6
pcaDS6_AmSa$rotation[DS6_top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

## DS6 :: Let us find top 10 measurements (genes) that contribute to PC2
##DS6 :: loading score calc for PC2
pca6loading_2 <- pcaDS6_AmSa$rotation[,2] #as PC2 accounts for next max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca6loading_2)   #29323
DS6genes_2 <- abs(pca6loading_2) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS6genes_2[1:25] ## displays loading scores of first 25 genes
DS6gene_rank_2 <- sort(DS6genes_2, decreasing = TRUE) ## sorts the genes based on score
DS6gene_rank_2[1:25]  ##displays the ranked genes, first 25
DS6_top_10_2 <- names(DS6gene_rank_2[1:10]) ## gets the names of top 10 genes and display
DS6_top_10_2 ## displays top 10 genes affecting the PC1 for variants in DS6
pcaDS6_mSa$rotation[DS6_top_10_2,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##scree plot-----------------------------------------------------------
library("FactoMineR")
#DS6.pca <- PCA(mDS6, graph = FALSE)
#DS6_ev <- get_eigenvalue(DS6.pca)
#DS6_ev
#fviz_eig(DS6.pca, addlabels = TRUE, ylim = c(0, 50))
### --------------------------------------------------------------------------------------
##-----------END OF DS6 DATA ANALYSIS ----------------------------------------------------
##***************************************************************************************

#dataset - 7 Africa, South Asian (Yeruba,Nigeria,Kenya,Beb,ITU)
#DS7_YerBeb <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/yer/20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
DS7_YerBeb <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/yer/YerBeb_20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS7_YerBeb)), decreasing = TRUE)[1:5]    ## NO duplicates
dimension(DS7_YerBeb) ## 29323 X 495

## CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS7_YB <- t(extract.gt(DS7_YerBeb, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS7_YB)         ## 494 X 29323

##-----  PCA  -------------------------------------------------------------------------------------
pcaDS7_YB <- prcomp(mDS7_YB) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS7_YB) 
str(pcaDS7_YB)                          ## 29323 X 494

# another method to find the proportional covariance
DS7pov <- summary(pcaDS7_YB)$importance[2,] ## The same pov can be calculated through this formula also
DS7pov[1:25]                                 ## displays the first 25 pov values

## ---------------------------------------------------------------------------------

## DS7 :: PLOT - PCA coverage ---------------------------
plot(pcaDS7_YB, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, 
     main="DS7:: PCA variant coverage graph   \n in DS7 dataset, Yeruba, Kenya, Nigeria \n BEB, ITU population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))
## ~60 % of variants is covered by PC1 and ~37% covered by PC2

pca7.var     <- pcaDS7_YB$sdev^2
pca7.var.per <- round(pca7.var/sum(pca7.var)*100,1)
pca7.var[1:10]
dimension(pca7.var.per) ##494
pca7.var.per[1:200]

## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca7.var.per[1:20], main = "DS7 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca7.var.per[1:20], main = "DS7 :: Bar plot to show the coverage of \n each principal component, first 20,
        AFR , AMR , EUR dataset ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(20,8,"PC1 = ~10%, PC2 = ~6%", col = "blue", font = 2))

##-----------------------------------------------------------------------------------
#biplot(pcaDS7_YB, scale = 0)          ## plots biplot graph, takes time
##-----------------------------------------------------------------------------------
# DS7 :: Plot dataset 7 - variance explained for each principal component
plot(DS7pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS7:: PoV - Proportion of Variance Explained graph \n in DS7 dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2 )
#----------------------------------
# DS7 :: Plot - cumulative proportion of variance explained
DS7y = sum(DS7pov[1:20])
sum(DS7y)
plotds7 <- plot((cumsum(DS7pov[1:20])), xlab = "Principal Components ",
                ylab = "Cumulative PoV Explained  ",
                main="DS7:: CPV - Cumulative Proportion of Variance Explained \n in DS7 dataset for first 20 samples ", 
                col.main="blue", font.main=3, cex.main=1.2,
                ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)
##------------------------------------------------------------------------------------------
## create a combined file DS7 with GT data and population panel file--------------------
pcaDS7_YB$x <- cbind(as.data.frame(pcaDS7_YB$x), "sample" = rownames(pcaDS7_YB$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS7_YB$x[,1:10]
dimension(pcaDS7_YB$x)                 ## 494 X 495

DS7 <- merge(pcaDS7_YB$x[,c(1:6,495)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS7)               ## Output matrix 494 X 10 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS7)                 ## list
DS7[1:5,]
str(pcaDS7_YB)
pcaDS7_YB$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS7_YB$center                   ##provides mean of the variables
pcaDS7_YB$scale                    ## provides Std Dev of the variables
dimension(pcaDS7_YB$rotation)      ## 29323 X 495

#- DS7 :: PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## DS 7:: PLOT FOR POPULATION ; data set 7 = Yeruba, Kenya, Nigeria, BEB, ITU population -------------
library(ggplot2)
ggplot(DS7, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20: DS7 : Population distribution for Yeruba, Kenya, Nigeria, BEB, ITU
           super population from pure DS7 population dataset ")

#-- DS7---------- PLOTTING --------------------------------------------
##------------------------------------------------------------------------------------------
## PLOT FOR POPULATION ; data set 7 = Yeruba, Kenya, Nigeria, BEB, ITU populations-------------
library(ggplot2)
ggplot(DS7, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr20: DS7 : Population distribution for Yeruba, Kenya, Nigeria, BEB, ITU
           super population from pure DS7 population dataset ")


##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS7 :: Let us find top 10 measurements (genes) that contribute to PC1 
##DS7 :: loading score calc for PC1
pca7loading <- pcaDS7_YB$rotation[,1] #as PC1 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca7loading)   #29323
DS7genes <- abs(pca7loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS7genes[1:25] ## displays loading scores of first 25 genes
DS7gene_rank <- sort(DS7genes, decreasing = TRUE) ## sorts the genes based on score
DS7gene_rank[1:25]  ##displays the ranked genes, first 25
DS7_top_10 <- names(DS7gene_rank[1:10]) ## gets the names of top 10 genes and display
DS7_top_10 ## displays top 10 genes affecting the PC1 for variants in DS7
pcaDS7_YB$rotation[DS7_top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

## DS7 :: Let us find top 10 measurements (genes) that contribute to PC2
##DS7 :: loading score calc for PC2
pca7loading1 <- pcaDS7_YB$rotation[,2] #as PC2 accounts for second highest variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca7loading1)   #29323
DS7genes1 <- abs(pca7loading1) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS7genes1[1:25] ## displays loading scores of first 25 genes
DS7gene1_rank <- sort(DS7genes1, decreasing = TRUE) ## sorts the genes based on score
DS7gene1_rank[1:25]  ##displays the ranked genes, first 25
DS7_top1_10 <- names(DS7gene1_rank[1:10]) ## gets the names of top 10 genes and display
DS7_top1_10 ## displays top 10 genes affecting the PC1 for variants in DS7
pcaDS7_YB$rotation[DS7_top1_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##scree plot-----------------------------------------------------------
library("FactoMineR")
#DS7.pca <- PCA(mDS7, graph = FALSE)
#DS7_ev <- get_eigenvalue(DS7.pca)
#DS7_ev
#fviz_eig(DS7.pca, addlabels = TRUE, ylim = c(0, 50))
### --------------------------------------------------------------------------------------
##------END OF DS7 DATA ANALYSIS ---------------------------------------------------------
##****************************************************************************************
##---------------------------------------------------------------------------------------------
###*##--------------------------------------------------------------------------------------------
#------------ POPULATION dataset 9 :: pure dataset -------------------------------------------
##  DS9 :: dataset - 9 ACB, ESN, GBR, IBS (4 pop)  PURE dataset ---------------------------------
library(vcfR)
library(R.utils)
DS9_AfEu <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/1-afr-eur/pure/22.25000001-26500000.ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
sort(table(getID(DS9_AfEu)), decreasing = TRUE)[1:5]    ## NO duplicates
#DS9_AfEu <- vcfR::read.vcfR("D:/GENOMICS/1000genome/C22/1-afr-eur/cleanAfEu_20.2500001-3500000.ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", verbose = FALSE)
dimension(DS9_AfEu)                     ## 48784 X 394

# CREATE LIST TRANSPOSED ------------------------------
## parses GT portion of vcf file
mDS9 <- t(extract.gt(DS9_AfEu, element = 'GT', as.numeric = TRUE, IDtoRowNames = TRUE))  ##pass
dimension(mDS9)         ## 393 X 48784

##-----  PCA -------------------------------------------------------------------------------------
pcaDS9 <- prcomp(mDS9) # performs PCA on the data and returns results as an object of class prcomp
dimension(pcaDS9) 
str(pcaDS9)                          ## 48784 X 393

# another method to find the proportional covariance
DS9pov <- summary(pcaDS9)$importance[2,] ## The same pov can be calculated through this formula also
DS9pov[1:5]                                 ## displays the first 5 pov values
typeof(DS9pov)                              ## double
dimension(DS9pov)                           ## 393
##------------------------------------------------------------------------------------------
# Plot dataset 9 - variance explained for each principal component
## DS9 :: PLOT - PCA coverage ---------------------------
plot(pcaDS9, type = "b", xlab = "Principal Components  ",
     pch=20,  col="brown", lwd = 2, ylim=c(1,80),
     main="DS9:: PCA variant coverage graph   \n in DS9 pure dataset, ACB, ESN, GBR and IBS population ", 
     col.main="brown", font.main=3, cex.main=1.2)   ##plots a line graph with variance = (Sqr(SD))

## 4.4% of variants is covered by first 2 Principal components

plot(DS9pov[1:20], xlab = "Principal Components  ",
     ylab = "Proportion of Variance Explained ",
     ylim = c(0, 1), type = "b",pch=20, cex=1, col="red", lwd = 2,
     main="DS9:: PoV - Proportion of Variance Explained graph \n in pure AFR/EUR dataset for first 20 samples ", 
     col.main="red", font.main=3, cex.main=1.2)

# Plot - cumulative proportion of variance explained

DS9y = sum(DS9pov[1:20])
sum(DS9y)
plotds2 <- plot((cumsum(DS9pov[1:20])), xlab = "Principal Components ",
                ylab = "Cumulative PoV Explained  ",
                main="DS9:: CPV - Cumulative Proportion of Variance Explained \n in DS9 pure dataset for first 20 samples ", 
                col.main="blue", font.main=3, cex.main=1.2,
                ylim = c(0, 1), type = "b", pch=20, cex=1, col="blue", lty = 5)

## DS9 :: Let us find the BARPLOT of percent variation represented by each PCs
pca9.var     <- pcaDS9$sdev^2 ##Finds how much variation in the orig. data each PC accounts for!!
pca9.var.per <- round(pca9.var/sum(pca9.var)*100,1) ##Finds % of variation each PC accounts for!!
pca9.var[1:10]
dimension(pca9.var.per) ##393
pca9.var.per[1:200]
pca9.var.per[1:20]  ## decide the limit based on <1 criteria

## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
## Barplot for % of variations represented by each pc. 
## we do not consider the PCs whose % is <1, so consider first 20 only
barplot(pca9.var.per[1:20], main = "DS9 :: Bar plot to show the coverage of \n each principal component, first 20", 
        xlab = "Principal components ",
        ylab = "Percent variation ")  ##pass
## first run the graph without the text inside and then run the below code.
barplot(pca9.var.per[1:20], main = "DS9 :: Bar plot to show the coverage of \n each principal component, first 20,
        ACB, ESN, GBR and IBS population dataset ", 
        xlab = "Principal components ",
        ylab = "Percent variation ", text(20,7,"PC1 = ~9%, PC2 = ~3%", col = "blue", font = 2))

##-----------------------------------------------------------------------------------

## 17% of variations are covered by first 2 Principal components
##------------------------------------------------------------------------------------------
##------------------------------------------------------------------------------------------
## create a combined file DS1 with GT data and population panel file--------------------
pcaDS9$x <- cbind(as.data.frame(pcaDS9$x), "sample" = rownames(pcaDS9$x))  ## Adds a 'sample' column at end of Pca matrix
pcaDS9$x[,1:10]
dimension(pcaDS9$x)                 ## 393 X 394

DS9 <- merge(pcaDS9$x[,c(1:3,394)], Popdata, by = "sample")  #creates new data binding original vcf data to PC's 1 to 6
dimension(DS9)               ## Output matrix 393 X 7 column, 1-sample;2-3 = PCA;4,5,6 = Popdata values
typeof(DS9)                 ## list
DS9[1:5,]
str(pcaDS9)
pcaDS9$rotation                 ## Provides PC loadings, with each column of this matrix containing PC vector
pcaDS9$center                   ##provides mean of the variables
pcaDS9$scale                    ## provides Std Dev of the variables
dimension(pcaDS9$rotation)      ## 48784 X 393

#------------ PLOTTING dataset 3  --------------------------------------------
# dataset = DS9
##------------------------------------------------------------------------------------------
## DS9 :: PLOT FOR POPULATION ; data set 9 = ACB, ESN, GBR and IBS  population, populations-------------
library(ggplot2)
ggplot(DS9, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS9 - : Population distribution for ACB, ESN, GBR and IBS 
            population from pure dataset DS9 ")

## DS9 :: PLOT FOR POPULATION ; data set 9 = ACB, ESN, GBR and IBS  population, populations-------------
library(ggplot2)
ggplot(DS9, aes(PC1, PC2, col = super_pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS9 - : Population distribution for ACB, ESN, GBR and IBS 
            population from pure dataset DS9 ")

## Scree plot ------------------------------------
library("FactoMineR")
DS9.pca <- PCA(mDS9, graph = FALSE)
DS9_ev <- get_eigenvalue(DS9.pca)
DS9_ev
fviz_eig(DS9.pca, addlabels = TRUE, ylim = c(0, 50))


##------------------------------------------------------------------------------------------
## --------------------
## PLOT FOR POPULATION- pure data set 3.1(a)- ACB, ESN, GBR and IBS population ONLY -------------

library(ggplot2)
ggplot(DS9, aes(PC1, PC2, col = pop)) +
  geom_point(shape = 22)  +
  ggtitle("Chr22:: DS9 - : Population distribution for ACB, ESN, GBR and IBS 
            population from pure dataset DS9 ")

##---FINDING THE GENES THAT ARE RESPONSIBLE FOR MAXIMUM VARIANTIONS ----------------------------------
## DS9 :: Let us find which genes have the largest effect , finding top 10 genes 
##DS9 :: loading score calc for PC1
pca9loading <- pcaDS9$rotation[,1] #as PC1 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca9loading)   #48784
DS9genes <- abs(pca9loading) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS9genes[1:25] ## displays loading scores of all genes
DS9gene_rank <- sort(DS9genes, decreasing = TRUE) ## sorts the genes based on score
DS9gene_rank[1:25]  ##displays the ranked genes
DS9_top_10 <- names(DS9gene_rank[1:10]) ## gets the names of top 10 genes and display
DS9_top_10 ## displays top 10 genes affecting the PC1 for variants in DS9
pcaDS9$rotation[DS9_top_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left

##DS9 :: loading score calc for PC2
pca9loading1 <- pcaDS9$rotation[,2] #as PC2 accounts for max variation of data

# Genes have large -ve values that push samples to left and have large +ve values that push samples to right
dimension(pca9loading1)   #48784
DS9genes1 <- abs(pca9loading1) ##gets absolute value of the gene loading score as we are interested in both +/_ values
DS9genes1[1:25] ## displays loading scores of all genes
DS9gene1_rank <- sort(DS9genes1, decreasing = TRUE) ## sorts the genes based on score
DS9gene1_rank[1:25]  ##displays the ranked genes
DS9_top1_10 <- names(DS9gene1_rank[1:10]) ## gets the names of top 10 genes and display
DS9_top1_10 ## displays top 10 genes affecting the PC1 for variants in DS9
pcaDS9$rotation[DS9_top1_10,1] ## see which genes has + and which has _ loading score pushing
## the samples to right or left
##------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------
##------End of DS9 data analysis -----------------------------------------------------------
## ****************************************************************************************

