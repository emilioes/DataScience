#' ---
#' title: "Description of data from DisGeNET"
#' output: html_notebook
#' ---
#' 
#' # Introduction
#' 
#' [DisGeNET](https://www.disgenet.org/) is a discovery platform containing one of the largest publicly available collections of genes and variants associated to human diseases (Piñero et al., 2019; Piñero et al., 2016; Piñero et al., 2015). DisGeNET integrates data from expert curated repositories, GWAS catalogues, animal models and the scientific literature. DisGeNET data are homogeneously annotated with controlled vocabularies and community-driven ontologies. Additionally, several original metrics are provided to assist the prioritization of genotype–phenotype relationships.
#' 
#' The current version of DisGeNET (v6.0) contains 628,685 gene-disease associations (GDAs), between 17,549 genes and 24,166 diseases, disorders, traits, and clinical or abnormal human phenotypes, and 210,498 variant-disease associations (VDAs), between 117,337 variants and 10,358 diseases, traits, and phenotypes.
#' 
#' Here we study the dataset **"Curated gene-disease associations"**. The file contains Gene-Disease associations from UniProt, CGI, ClinGen, Genomics England Panel App, PsyGeNET, Orphanet, the HPO, and CTD (human data).
#' 
#' **The columns in the files are**:
#' 
#' * geneId 		-> NCBI Entrez Gene Identifier
#' * geneSymbol	-> Official Gene Symbol
#' * DSI		-> The Disease Specificity Index for the gene
#' * DPI		-> The Disease Pleiotropy Index for the gene
#' * diseaseId 	-> UMLS concept unique identifier
#' * diseaseName 	-> Name of the disease	
#' * diseaseType  	-> The DisGeNET disease type: disease, phenotype and group
#' * diseaseClass	-> The MeSH disease class(es)
#' * diseaseSemanticType	-> The UMLS Semantic Type(s) of the disease
#' * score		-> DisGENET score for the Gene-Disease association
#' * EI		-> The Evidence Index for the Gene-Disease association
#' * YearInitial	-> First time that the Gene-Disease association was reported
#' * YearFinal	-> Last time that the Gene-Disease association was reported
#' * NofPmids	-> Total number of publications reporting the Gene-Disease association
#' * NofSnps		-> Total number of SNPs associated to the Gene-Disease association
#' * source		-> Original source reporting the Gene-Disease association
#' 
#' # Inatall relevant packages
#' 
## -------------------------------------------------------------------------
if (!require("pROC")) install.packages("pROC")

#' 
#' 
#' # Loading data 
#' 
#' We first read the gzipped file. The option "stringsAsFactors" is always FALSE per default since a recent update. This is mostly useful, except when showing the summary of the data frame in the cell below. For all practical purposes it does not make a difference in the current analysis.
## -------------------------------------------------------------------------
# Loading the tsv (tab seaprated variable) file into a dataframe
df <- read.delim('curated_gene_disease_associations.tsv', stringsAsFactors = T)

#' 
#' We show the beginning of the file to check its content.
## -------------------------------------------------------------------------
head(df)

#' 
#' We can also look at **summary statistics** of the file. This contains information about quantiles for numeric variables, and top words for strings.
## -------------------------------------------------------------------------
summary(df)

#' 
#' 
#' 
#' # Temporal insights
#' 
#' There are two timestamps in the dataset:
#' 
#' * YearInitial	-> First time that the Gene-Disease association was reported
#' * YearFinal	-> Last time that the Gene-Disease association was reported
#' 
#' From these, we can compute things like:
#' 
#' * Growth of number of discoveries in time.
#' * How certain genes have been "forgotten"
#' 
#' ## Growth process
#' 
#' First it is always good to get an idea of **distributions**. A distribution corresponds to a count of the number (or probability, if normalised) of elements per bin.
#' 
## -------------------------------------------------------------------------
hist(df$YearInitial, col='gray', xlab='Initial year', breaks = 100)

#' 
#' 
#' What can be observed? There as been a surge in discoveries since the 1980s (beginning of molecular/cellular biology), a rapid increase with arrival of high throughput technics late 1990s, and a more recent decay (there is only a finite amount of possible associations!).
#' 
#' To quantify better this rise-and-fall curve, one can look at the cumulative function! It will show us the cumulative amount of associations found through time. In particular, it will show us if there was a change of rate of increase / decrease through time.
#' 
#' 
#' We first compute the number of observations per year
## -------------------------------------------------------------------------
tb <- table(df$YearInitial)
head(tb)

#' 
#' And then we can look at the cumulative sum:
## -------------------------------------------------------------------------
x <- as.numeric(names(tb))
y <- cumsum(as.numeric(tb))

# another way os to do a empirical cumulative distribution, which is normalised to 1
# y <- ecdf(df$YearInitial)(x)

#' 
#' 
## -------------------------------------------------------------------------
plot(x,y, pch=19, type='l', xlab='Year', ylab='Total number of associations', lwd=5)

#' 
#' To better understand if the rates of increase are changing, we need to use a log scale.
#' 
## -------------------------------------------------------------------------
plot(x,y, log='y', 
     pch=19, type='l', xlab='Year', ylab='Total number of associations', lwd=5)

#' 
#' The linearity of the curve tells us a story! The discovery rate of novel associations follows an exponental growth $y \sim e^{at}$. The rate $a$ gives us the rate of growth. For example to get the doubling rate (time to have twice the number of known associations) we compute
#' 
#' $$\frac{2y}{y} = e^{a(t_2-t_1)}$$
#' 
#' $$ t_2-t_1 =1/a * \log(2)$$
#' 
#' Let us fit the data with a linear fit!
#' The function *lm* generates a linear regression model:
#' 
#' $$y \sim x \Leftrightarrow y=ax+b$$
#' 
## -------------------------------------------------------------------------
# we are interested in the log here!
logy <- log(y)/log(2)

# we don't want to fit before 1945 
inds <- which(x>1945) 
# the lm function can be understood by looking at the manual through ?lm
fit <- lm(logy[inds]~x[inds])

#' 
#' Once the fit is created we can observe the statistics of the quality of fit (see ?summary.lm for help). Here, we are interested in
#' 
#' * *R-squared*, the ‘fraction of variance explained by the model’ --> the closer to 1, the better. Here 0.995 is very very good!
#' * the *coefficient estimate*: 0.125 is very well constrained (small amount of noise)
#' * the *precision of the estimate*: the t value or the p-value give a standardized insight into the quality of the test. If the t value is (in absolute value) above 2, then it is usually a sign that the association is significant.
#' 
## -------------------------------------------------------------------------
summary(fit)

#' 
#' Let us represent things a bit better. We can plot the fit on top of the other curve!
#' 
## -------------------------------------------------------------------------
plot(x,logy, type='l', cex=0.2, xlab='Year',ylab='Cumulative number (log2)')
abline(fit, lty=2, col='red', lwd=3)

#' 
#' 
#' The trend is observed to be a very good fit! So what does it give us in terms of doubling rate $t_2-t_1 =\frac{\log(2)}{a}$?
#' 
## -------------------------------------------------------------------------
a <- fit$coefficients[2]
# the fit already incliudes the term log(2) in our case
doubling_rate <- 1 / a
print(paste('The total number of known associations doubles every',
            signif(doubling_rate,2),
            'years!'))

#' 
#' ## Evolution of quantities
#' 
#' Now we are interested in the "Final" year an association has been observed. This gives some insights about which associations have been forgotten. Let us see first look at how this final year depends on the initial year.
#' 
## ---- fig.asp=1, fig.width=3----------------------------------------------
#fig.asp=1 makes a square, easier to visualise y=x
plot(df$YearInitial, df$YearFinal, 
     pch=19, # nicer looking
     cex=0.2, # many points so we reduce the size
     xlab='Initial year',
     ylab='Final year'
     )

#' 
#' Elements on the diagonal have only been observed once. The more "up" we are, the longer the observation has been observed. We want to understand the relative time an association has been observed. For this we compute a normalised duration:
#' 
#' 
## ---- warning=FALSE-------------------------------------------------------
dt <- df$YearFinal - df$YearInitial
dt_max <- max(df$YearInitial, df$YearFinal, na.rm=T) - df$YearInitial
survival <- dt / dt_max
hist(survival)

#' 
#' We see a lot of associations are at 0 (only observed once). We can filter them out to see better the rest of the distribution
#' 
## -------------------------------------------------------------------------
inds <- which(survival>0)
hist(survival[inds], xlab='Survival', main='Keeping only positive survivals')

#' 
#' Now we wonder what are the trends for survival? Did it change in time? Are the first associations found more "interesting" (i.e they have been surviving longer)?
## ---- fig.height=3--------------------------------------------------------
# this creates a list of survivals for each Year Initial
sp <- split(survival, df$YearInitial)
boxplot(sp, xlab='Year Initial', ylab='Survival',col='gray')

#' 
#' 
#' We can have a simpler way to look at this through mean and standard error:
#' 
## ---- warning=FALSE-------------------------------------------------------
# this library allows to plot error bars
require(gplots)

# we compute the mean and standard error for each element of the list
means <- sapply(sp, mean, na.rm=T)
stderrs <- sapply(sp, function(x) sd(x, na.rm=T) / sqrt(length(x)))

years <- as.numeric(names(means))

# barplot.split(df$YearInitial[inds], survival[inds])
# plot(means, border = NA, ylab='Survival', las =3)
plotCI(years, means, stderrs, 
       xlab='Year',
       ylab='Survival',
       gap=0, pch=19, type='o',cex=0.3, sfrac=0)

#' 
#' 
#' Let us try to distinguish two periods: before 1970 and after. We want to explore simple tests that these periods have different survival distributions.
#' 
## -------------------------------------------------------------------------
means1 <- means[which(years<1970)]
means2 <- means[which(years>=1970)]

boxplot(list('Before 1970'=means1,
             'After 1970'=means2),
        col=c('lightblue','gray'),
        ylab='Survival')

#' 
#' When we have two distributions like this, we want to test that they are different. 
#' 
#' To do this we can use a *non-parametric* method. By this we mean that this method looks at the **rank** of quantities, not their absolute values. In this context, we want to show that before 1970 survival is higher than after 1970. We convert survival to a rank (first largest, second largest etc) and we test that "Before 1970" has higher rank. 
#' 
#' This can be easily seen using a ROC curve. The ROC curve is used in classification schemes to test that one distribution is separated from another one. The more the ROC curves points towards the top-left corner the better.
#' 
## ---- fig.asp=1, fig.width=3, warning=FALSE-------------------------------
library(pROC)

responses <- c(means1, means2)
labels <- c(rep(1, length(means1)), rep(0, length(means2)))

r <- roc(labels, responses, direction="<", levels = c(0,1))
plot(r, cex.lab=1.5, cex.axis=1.5, lwd=4)

#' 
#' Here the ROC curve looks really good. We can quantify this in two ways. First, one can compute the area under the ROC curve (AUC), and if the AUC is above 0.5, this point towards distinct distirbutions.
#' 
## -------------------------------------------------------------------------
print(r)

#' 
#' Here an AUC of 0.98 is extremely good. However it is possible to have hjigh values by chance when we have low number of observations. To control for this, we compute a probability to observe such a high value by chance given the number of observations. To do this, we can use the *Mann Whitney U test* also call *wilcoxon test*
#' 
## -------------------------------------------------------------------------
wilcox.test(means1, means2, alternative = 'greater')

#' 
#' We see that the p-value that AUC is this big is ~2e-11, very significant as expected "by eye".
#' 
#' Tomorrow we will look at how to build networks from this dataset!
