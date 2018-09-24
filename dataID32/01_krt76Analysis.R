# File: 01_krt76Analysis.R
# Auth: umar.niazi@kcl.ac.uk
# DESC: selection and modelling of the gene of interest
# Date: 24/09/2018

## load the data
source('header.R')

## load the data
library(RMySQL)

db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)

# get the query
g_did
q = paste0('select MetaFile.* from MetaFile
           where (MetaFile.idData = 32) AND (MetaFile.comment like "%ines%")')
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

columns(org.Mm.eg.db)
df = AnnotationDbi::select(org.Mm.eg.db, keys = 'Krt76', columns = 'ENSEMBL', keytype = 'SYMBOL')

#### present in this raw data
table(rownames(mData) %in% df$ENSEMBL)

summary(mData[df$ENSEMBL,])
plot(mData[df$ENSEMBL,])

######## perform some EDA on the count matrix and normalisation
dim(mData)

# drop the genes where average across rows is less than 3
i = rowMeans(mData)
table( i < 3)
mData = mData[!(i< 3),]
dim(mData)

ivProb = apply(mData, 1, function(inData) {
  inData[is.na(inData) | !is.finite(inData)] = 0
  inData = as.logical(inData)
  lData = list('success'=sum(inData), fail=sum(!inData))
  return(mean(rbeta(1000, lData$success + 0.5, lData$fail + 0.5)))
})

table(ivProb < 0.8)
hist(ivProb)
mData = mData[!(ivProb < 0.8), ]
dim(mData)

library(DESeq2)
sf = estimateSizeFactorsForMatrix(mData)
mData.norm = sweep(mData, 2, sf, '/')

## load CDiagnostics and test
## compare the normalised and raw data
library(downloader)
url = 'https://raw.githubusercontent.com/uhkniazi/CDiagnosticPlots/master/CDiagnosticPlots.R'
download(url, 'CDiagnosticPlots.R')

# load the required packages
source('CDiagnosticPlots.R')
# delete the file after source
unlink('CDiagnosticPlots.R')

oDiag.1 = CDiagnosticPlots(log(mData.norm+1), 'Normalised')
oDiag.2 = CDiagnosticPlots(log(mData+1), 'Original')

# the batch variable we wish to colour by, 
# this can be any grouping/clustering in the data capture process
#### pre-Aire (MHCIIlo;RFP-), early-Aire (MHCIIhi;RFP-), late-Aire (MHCIIhi;RFP+), and post-Aire (MHCIIlo;RFP+). 
fBatch = gl(4, k = 4, labels = c('RFPn_hi', 'RFPn_low', 'RFPp_hi', 'RFPp_low'))
table(fBatch)

pdf('results/matrixClustering.pdf')

## compare the 2 methods using various plots
par(mfrow=c(1,2))
boxplot.median.summary(oDiag.1, fBatch, legend.pos = 'center', axis.label.cex = 0.7)
boxplot.median.summary(oDiag.2, fBatch, legend.pos = 'topright', axis.label.cex = 0.7)

plot.mean.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.mean.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.sigma.summary(oDiag.1, fBatch, axis.label.cex = 0.7)
plot.sigma.summary(oDiag.2, fBatch, axis.label.cex = 0.7)

plot.missing.summary(oDiag.1, fBatch, axis.label.cex = 0.7, cex.main=1)
plot.missing.summary(oDiag.2, fBatch, axis.label.cex = 0.7, cex.main=1)

plot.PCA(oDiag.1, fBatch, cex.main=1)
plot.PCA(oDiag.2, fBatch, cex.main=1)
# 
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.8, cex.main=0.7)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.8, cex.main=0.7)

## change parameters 
l = CDiagnosticPlotsGetParameters(oDiag.1)
l$PCA.jitter = F
l$HC.jitter = F

oDiag.1 = CDiagnosticPlotsSetParameters(oDiag.1, l)
oDiag.2 = CDiagnosticPlotsSetParameters(oDiag.2, l)
plot.PCA(oDiag.1, fBatch)
plot.PCA(oDiag.2, fBatch)
plot.dendogram(oDiag.1, fBatch, labels_cex = 0.6)
plot.dendogram(oDiag.2, fBatch, labels_cex = 0.7)

##### fit a model and produce some plots

dfData = data.frame(Krt76=round(mData.norm[df$ENSEMBL,],0), Coef=factor(fBatch, levels = c("RFPn_low", "RFPn_hi", "RFPp_hi", "RFPp_low")))

library(lme4)
fit.lme1 = glmer.nb(Krt76 ~ 1  + (1 | Coef), data=dfData)
summary(fit.lme1)

## setup the stan model
detach('package:org.Mm.eg.db', unload=T)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stanDso = rstan::stan_model(file='../nbResp1RandomEffect.stan')

## calculate hyperparameters for variance of coefficients
l = gammaShRaFromModeSD(sd(log(dfData$Krt76+0.5)), 2*sd(log(dfData$Krt76+0.5)))

### set stan input data
lStanData = list(Ntotal=nrow(dfData), Nclusters1=nlevels(dfData$Coef),
                 Nsizes = 1,
                 NgroupMap1=as.numeric(dfData$Coef),
                 NsizeMap=rep(1, times=nrow(dfData)),
                 Ncol=1,
                 y=dfData$Krt76, 
                 gammaShape=l$shape, gammaRate=l$rate,
                 intercept = mean(log(dfData$Krt76+0.5)), intercept_sd= sd(log(dfData$Krt76+0.5))*2)


fit.stan = sampling(stanDso, data=lStanData, iter=10000, chains=4,
                    pars=c('sigmaRan1', 'betas',
                           'iSize',
                           'rGroupsJitter1',
                           'mu'
                           ),
                    cores=4)#, init=initf)#, control=list(adapt_delta=0.99, max_treedepth = 15))


print(fit.stan, digits=3)
traceplot(fit.stan, c('betas'))
traceplot(fit.stan, 'sigmaRan1')

## make some plots for results
library(lattice)
bwplot(log(Krt76) ~ Coef, data=dfData, panel=panel.violin, type='b',
       par.strip.text=list(cex=0.7), varwidth=F, main='Normalised expression for Krt76')

bwplot(log(Krt76) ~ Coef, data=dfData, panel=function(x, y, ...) panel.bwplot(x, y, pch='|',...), type='b',
       par.strip.text=list(cex=0.7), varwidth=T, main='Normalised expression for Krt76')

### create plot for coefficients
mModules = extract(fit.stan)$rGroupsJitter1
dim(mModules)
## get the intercept at population level
iIntercept = extract(fit.stan)$betas[,1]
## add the intercept to each random effect variable, to get the full coefficient
mModules = sweep(mModules, 1, iIntercept, '+')

## format data for plotting
m = colMeans(mModules)
s = apply(mModules, 2, sd)*1.96
d = data.frame(m, s, s1=m+s, s2=m-s)
d$mods = levels(dfData$Coef)

dotplot(mods ~ m+s1+s2, data=d, panel=llines(d$s1, d$s2), cex=0.6, pch=20,
        par.strip.text=list(cex=0.5), main='Fitted Coefficients for Krt76 ~ sample group', xlab='log Expression')
