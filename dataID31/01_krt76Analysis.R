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
           where (MetaFile.idData = 31) AND (MetaFile.comment like "%ines%")')
dfSample = dbGetQuery(db, q)
dfSample
n = paste0(dfSample$location, dfSample$name)
load(n)

oExp = lExp$normalised
dfSample = pData(oExp)
mData = exprs(oExp)

library(org.Mm.eg.db)
head(rownames(mData))

columns(org.Mm.eg.db)
df = AnnotationDbi::select(org.Mm.eg.db, keys = 'Krt76', columns = 'ENSEMBL', keytype = 'SYMBOL')

#### not present in this normalised data
table(rownames(mData) %in% df$ENSEMBL)
