#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
print(args[1])
print(args)

library(bereR)
library(bptools)

local.wd <- "~/gd/Harvard/Research/"
data.dir  <- "./data/Eclipse/"
setwd(local.wd)
setwd(data.dir)

copd.filename <- paste(local.wd,"bere_null/","null.networks",Sys.time(),".rds",sep="")


######################################################
##      Data Loading from ECLIPSE dataset          ###
##                                                 ###
######################################################
eclipse <- list()
eclipse$motif    <- read.table("ECLIPSE_Blood_Motif.txt",header=F)
eclipse$exp      <- read.table("ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
eclipse$ppi      <- read.table("OV_PPI.txt",header=F)
eclipse$clinical <- read.table("ECLIPSE_blood.txt",header=T,fill = TRUE, sep="\t",row.names=1)
eclipse$exp      <- eclipse$exp[,order(colnames(eclipse$exp))]  # Make sure expression and clinical is in same order
eclipse$clinical <- eclipse$clinical[colnames(eclipse$exp),]    # Make sure clinical only contains patients with expression data

# Specify the group partition
filter.vec.1 <- eclipse$clinical$Subject.type=="COPD Subjects"
filter.vec.2 <- eclipse$clinical$Subject.type=="Smoker Controls"


######################################################
##  Running observed networks                      ###
##                                                 ###
######################################################

res1 <- bere(eclipse$motif, eclipse$exp[,filter.vec.1],cpp=F)
res2 <- bere(eclipse$motif, eclipse$exp[,filter.vec.2],cpp=F)
tmObs <- transformation.matrix(res1, res2)

######################################################
##  Running null networks with improved algorithm  ###
##                 2/25/15    START                ###
######################################################

# Copy expression data for null network generation
null.exp <- eclipse$exp

#start time
strt<-Sys.time()
    
rownames(null.exp) <- rownames(null.exp)[sample(1:nrow(null.exp))]
res1 <- bere(eclipse$motif,null.exp[,filter.vec.1],cpp=F)
res2 <- bere(eclipse$motif,null.exp[,filter.vec.2],cpp=F)
null.networks <- list(res1,res2) 


print(Sys.time()-strt)

# Save the observed and null networks (as separate files)
saveRDS(eclipse,eclipse.filename)
saveRDS(null.networks,copd.filename)