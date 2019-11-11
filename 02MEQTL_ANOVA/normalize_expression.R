"%&%" = function(a,b) paste(a,b,sep="")
library(tibble)
library(dplyr)
source('/home/rschubert1/data/GEU/rntransform.R')
source('/home/rschubert1/data/GEU/ztransform.R')
library(preprocessCore)
library(argparse)
library(methods)

parser <- ArgumentParser()
parser$add_argument("-e", "--expression", help="file path of the gene expression file")
parser$add_argument("-o", "--outputdir", help="output directory", type="character")
parser$add_argument("-t", "--tag", help="file tag for this run of samples")
args <- parser$parse_args()

expr = read.table(file=args$expression,header=T,check.names = F)
expmat = as.matrix(t(expr[,-1])) #need genes in cols, ids in rows
qn.expmat <- as.matrix(normalize.quantiles(expmat)) ##quantile normalize
qn<-as.matrix(t(qn.expmat))
rownames(qn) <- expr[,1]
colnames(qn) <- colnames(expr)[-1]
write.table(file = args$outputdir %&% "/" %&% args$tag %&% "QN.txt", x = qn, quote = F, sep = '\t')

rn.qn.expmat <- apply(qn.expmat,1,"rntransform") ##rank transform to normality & transposes##
rownames(rn.qn.expmat) <- expr[,1]
colnames(rn.qn.expmat) <- colnames(expr)[-1]
write.table(file = args$outputdir %&% "/" %&% args$tag %&% "RN.txt", x = rn.qn.expmat, quote = F, sep = '\t')

