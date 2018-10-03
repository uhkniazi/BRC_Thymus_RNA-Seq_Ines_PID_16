# File: 01_krt76Analysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: selection and modelling of the gene of interest
# Date: 28/09/2018

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 34) AND (MetaFile.comment like "%ines%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

dbDisconnect(db)

library(Biobase)
dfSample = pData(oExp)
mData = exprs(oExp)

library(org.Mm.eg.db)
head(rownames(mData))

# columns(org.Mm.eg.db)
# df = AnnotationDbi::select(org.Mm.eg.db, keys = 'Krt76', columns = 'ENSEMBL', keytype = 'SYMBOL')

#### present in this raw data
table(rownames(mData) %in% 'Krt76')

summary(mData['Krt76',])
plot(mData['Krt76',])

######## perform some EDA on the count matrix and normalisation
dim(mData)


########################################## single cell qc
library(scater)
oSce = SingleCellExperiment(assays = list(counts=mData), colData=pData(oExp))
dim(oSce)

## extract spikein and mitochondrial rows
iSpike = grep('^ERCC', rownames(oSce)); length(iSpike)

## calculate QC metrics 
oSce = calculateQCMetrics(oSce, feature_controls = list(ERCC=iSpike))
head(colnames(colData(oSce)))

## provide spike-in information
library(scran)
isSpike(oSce, 'ERCC') = iSpike

# Two common measures of cell quality are the library size and the number of expressed features in each library.
par(mfrow=c(1,2))
hist(oSce$total_counts/1e6, xlab="Library size in millions", main="QC Library Size", prob=T,
     ylab="")
hist(oSce$total_features, xlab="Number of expressed genes", main="QC Features with non-zero counts", 
     ylab="", prob=T)

# remove cells with log-library sizes that are more than 3 MADs 
# below the median log-library sizs
bDropLibsize = isOutlier(oSce$total_counts, nmads=3, type="lower", log=TRUE); table(bDropLibsize)
bDropFeature = isOutlier(oSce$total_features, nmads=3, type="lower", log=TRUE); table(bDropFeature)

# The quantity of spike-in RNA added to each cell/well should be constant, 
# which means that the proportion should increase upon loss of endogenous RNA in low-quality cells.
hist(oSce$pct_counts_feature_control, xlab="ERCC proportion (%)", 
     ylab="", main="QC Proportion of reads mapped to ERCC")

# cells with more endogenous RNA or that are assayed with protocols using less spike-in RNA will have lower spike-in proportions
bDropSpike = isOutlier(oSce$pct_counts_feature_control, nmads=3, type="higher"); table(bDropSpike)

# Removal of a substantial proportion of cells (> 10%) may be indicative of an overall issue with data quality.
oSce = oSce[,!(bDropFeature | bDropLibsize)]# | bDropSpike)]
data.frame(ByLibSize=sum(bDropLibsize), ByFeature=sum(bDropFeature),
           BySpike=sum(bDropSpike), Remaining=ncol(oSce))
# ByLibSize ByFeature BySpike Remaining
# 1         1         3      26       107
par(p.old)

dim(oSce)
dim(colData(oSce))

## normalise the data
oSce = computeSumFactors(oSce, sizes=c(10, 20, 25))
sizeFactors(oSce)

# we compute a separate set of size factors for the spike-in set. For each cell, 
# the spike-in-specific size factor is defined as the total count across all transcripts in the spike-in set.
# try both ways to see which method is better for normalization i.e. general.use=T and general.use=F
oSce.F = computeSpikeFactors(oSce, type='ERCC', general.use=F)
oSce.T = computeSpikeFactors(oSce, type='ERCC', general.use=T)

# The computed values are stored as an exprs matrix in addition to the other assay elements
# and are log-transformed
oSce.F = normalize(oSce.F)
oSce.T = normalize(oSce.T)
##########################################

## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

## drop the intestine cells
i = grep('Intestine', as.character(dfSample$characteristics_ch1.1))

oSce.F.bk = oSce.F
oSce.T.bk = oSce.T
oSce.F = oSce.F[,-i]
oSce.T = oSce.T[,-i]

oDiag.1 = CDiagnosticPlots(exprs(oSce.F), 'no spikein')
oDiag.2 = CDiagnosticPlots(exprs(oSce.T), 'spikein')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
dfSample = data.frame(colData(oSce.T))
# fBatch = dfSample$characteristics_ch1.4
# table(fBatch)
# fBatch = factor(fBatch)

fBatch = as.character(dfSample$characteristics_ch1.4)
fBatch = gsub('mouse age: ', '', fBatch)
fBatch = factor(fBatch, levels=c('E14.5', "E18.5", "6 days", "4-6 weeks"))
tx = gsub('selection marker: ', '', as.character(dfSample$characteristics_ch1.2))
tx = factor(tx)
st = gsub('strain: ', '', as.character(dfSample$characteristics_ch1))
st = factor(st)

dir.create('results')
pdf('results/matrixClustering.pdf')

## compare the 2 methods using various plots
par(mfrow=c(1,1))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)
plot.PCA(oDiag.1, fBatch, cex.main=1, csLabels = as.character(tx))
# 
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

table(fBatch)
fBatch = gsub('mouse age: ', '', fBatch)
tx = gsub('selection marker: ', '', dfSample$characteristics_ch1.2)
fBatch = factor(fBatch)

plot.PCA(oDiag.1, fBatch, csLabels = tx, cex.main=1)
plot.PCA(oDiag.1, factor(tx), csLabels = as.character(fBatch), cex.main=1)
# ## change parameters 
# l = CDiagnosticPlotsGetParameters(oDiag.1)
# l$PCA.jitter = F
# l$HC.jitter = F
# 
# oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
# oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
# plot.PCA(oDiag.1, fBatch)
# plot.PCA(oDiag.2, fBatch)
# plot.dendogram(oDiag.1, fBatch, labels_cex = 0.6)
# plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

##### fit a model and produce some plots
#levels(fBatch) = c('E14.5', "E18.5", "4-6 weeks", "6 days")
library(lattice)
fBatch = as.character(dfSample$characteristics_ch1.4)
fBatch = gsub('mouse age: ', '', fBatch)
fBatch = factor(fBatch, levels=c('E14.5', "E18.5", "6 days", "4-6 weeks"))
tx = gsub('selection marker: ', '', as.character(dfSample$characteristics_ch1.2))
tx = factor(tx)
st = gsub('strain: ', '', as.character(dfSample$characteristics_ch1))
st = factor(st)
dfData = data.frame(Krt76=exprs(oSce.F)['Krt76', ], Coef=fBatch, Coef.adj=st, Coef.adj2 = tx)

# dotplot(Krt76 ~ Coef, data=dfData, groups=Coef.adj, pch=c(1:4), col=1:4, key=list(space='top', columns=2, points=list(pch=1:4, col=1:4),
#                                                                                   text=list(levels(dfData$Coef.adj))),
#         main='Normalised Krt Expression on Log Scale')

dotplot(Krt76 ~ Coef | Coef.adj:Coef.adj2, data=dfData, type='p', par.strip.text=list(cex=0.7),
        main='Krt76 expression vs mouse age grouped by strain:selection marker ')

## drop some empty factor levels
dfData = dfData[dfData$Coef.adj2 %in% c('CD45-', 'CD45- Epcam+'),]
dfData = droplevels.data.frame(dfData)
dfData$Coef.adj = dfData$Coef.adj:dfData$Coef.adj2
dfData = droplevels.data.frame(dfData)


library(lme4)
fit.lme1 = lmer(Krt76 ~ 1  + (1 | Coef), data=dfData)
summary(fit.lme1)

fit.lme2 = lmer(Krt76 ~ 1  + (1 | Coef) + (1 | Coef.adj), data=dfData)
summary(fit.lme2)

anova(fit.lme1, fit.lme2)
## setup the stan model
detach('package:org.Mm.eg.db', unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='../tResponse2RandomEffectsNoFixed.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(dfData$Krt76+0.5), 2*sd(dfData$Krt76+0.5))


### set stan input data
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nclusters2=nlevels(dfData$Coef.adj),
                 NgroupMap1=as.numeric(dfData$Coef),
                 NgroupMap2=as.numeric(dfData$Coef.adj),
                 y=dfData$Krt76, 
                 gammaShape=l$shape, gammaRate=l$rate)

fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=4,
                    pars=c('betas', 'sigmaRan1', 'sigmaRan2', 
                           'sigmaPop', 'mu', 'nu',
                           'rGroupsJitter1', 'rGroupsJitter2'),
                    cores=4, control=list(adapt_delta=0.99, max_treedepth = 12))

print(fit.stan, c('betas', 'sigmaRan1', 'sigmaRan2', 'sigmaPop', 'rGroupsJitter1', 'rGroupsJitter2', 'nu'), digits=3)
traceplot(fit.stan, 'betas')
traceplot(fit.stan, 'sigmaRan1')

m = extract(fit.stan)$mu
plot(density(dfData$Krt76))
lines(density(fitted(fit.lme1)), col=1)
lines(density(fitted(fit.lme2)), col=2)
lines(density(colMeans(m)), col=3)
#apply(m, 1, function(x) lines(density(x), col=3))

## make some plots for results
library(lattice)
bwplot(Krt76 ~ Coef, data=dfData, panel=panel.violin, type='b',
       par.strip.text=list(cex=0.7), varwidth=F, main='Normalised expression for Krt76')

bwplot(Krt76 ~ Coef | Coef.adj, data=dfData, panel=panel.violin, type='b',
       par.strip.text=list(cex=0.7), varwidth=F, main='Normalised expression for Krt76')

bwplot(colMeans(m) ~ Coef | Coef.adj, data=dfData, panel=panel.violin, type='b',
       par.strip.text=list(cex=0.7), varwidth=F, main='Model Fitted Krt76')

bwplot(Krt76 ~ Coef, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
       par.strip.text=list(cex=0.7), varwidth=T, main='Normalised expression for Krt76')

xyplot(Krt76 ~ Coef, data=dfData, type='p')
dotplot(Krt76 ~ Coef, data=dfData, type='p')
dotplot(Krt76 ~ Coef, data=dfData, groups=Coef.adj, auto.key=list(space='top', columns=2))
dotplot(Krt76 ~ Coef, data=dfData, groups=Coef.adj, pch=c(1:4), col=1:4, key=list(space='top', columns=2, points=list(pch=1:4, col=1:4),
                                                                         text=list(levels(dfData$Coef.adj))),
        main='Normalised Krt Expression on Log Scale')

dotplot(Krt76 ~ Coef | Coef.adj, data=dfData, type='p', par.strip.text=list(cex=0.7))

table(dfData$Coef)
table(dfData$Coef.adj)
tapply(dfData$Krt76, dfData$Coef.adj, mean)
tapply(dfData$Krt76, dfData$Coef, mean)
# mCoef = extract(fit.stan)$rGroupsJitter1
# dim(mCoef)
# colnames(mCoef) = levels(dfData$Coef)
# 
# mCoef.adj = extract(fit.stan)$rGroupsJitter2
# dim(mCoef.adj)
# colnames(mCoef.adj) = levels(dfData$Coef.adj)

## get the intercept at population level
# iIntercept = as.numeric(extract(fit.stan)$betas)
# # add the intercept to each random effect variable, to get the full coefficient
# m = mCoef[,'4-6 weeks'] + mCoef.adj[,'CD45+ IL25+']
# m = m + iIntercept
# hist(m, xlab='', main='Fitted Value')

### create plot for coefficients
mModules = extract(fit.stan)$rGroupsJitter1
dim(mModules)
## get the intercept at population level
iIntercept = extract(fit.stan)$betas
## add the intercept to each random effect variable, to get the full coefficient
mModules = sweep(mModules, 1, iIntercept, '+')

## format data for plotting
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Coef)
d$mods = factor(d$mods, levels=c('E14.5', "E18.5", "6 days",  "4-6 weeks"))

dotplot(mods ~ m+s1+s2, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=0.5), main='Model Based Average for Krt76 ~ sample group', xlab='Average log Expression')
