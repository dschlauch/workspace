library(bereR)
library(pandaR)
library(bptools)
library(reshape2)
library(penalized)
library(Biobase)
library(org.Hs.eg.db)
library(foreach)
library(doParallel)
library(limma)
library(igraph)
library(ggrepel)


# Figure 2a-c -------------------------------------------------------------

load("~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/ECLIPSE_combined_runs/activeImage.RData")
outputDir <- "~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/"
analysisCode <- "ECLIPSE"
source("~/gd/Harvard/Research/R_workspace/process_TM.R")
supDTFIvsEXPdf <- cbind(plotDF, study="ECLIPSE")

load("~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/COPDGene_combined_runs/activeImage.RData")
outputDir <- "~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/"
analysisCode <- "COPDGene"
source("~/gd/Harvard/Research/R_workspace/process_TM.R")
supDTFIvsEXPdf <- rbind(supDTFIvsEXPdf,cbind(plotDF, study="COPDGene"))

load("~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/LGRC_combined_runs/activeImage.RData")
outputDir <- "~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/"
analysisCode <- "LGRC"
source("~/gd/Harvard/Research/R_workspace/process_TM.R")
supDTFIvsEXPdf <- rbind(supDTFIvsEXPdf,cbind(plotDF, study="LGRC"))


load("~/gd/Harvard/Research/TM_outputs/JASPAR2014/BERE/LTCOPD_combined_runs/activeImage.RData")
outputDir <- "~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/"
analysisCode <- "LTCOPD"
source("~/gd/Harvard/Research/R_workspace/process_TM.R")
supDTFIvsEXPdf <- rbind(supDTFIvsEXPdf,cbind(plotDF, study="LTCOPD"))


dTFI_LIMMA_all <- ggplot(data=supDTFIvsEXPdf, aes(x=limmanegLogPValues, y=negLogZPValues)) + facet_wrap( ~ study, ncol=2) + 
  geom_point(aes(col=logfoldchangeTF), size=7, alpha=.8) +
  geom_point(shape = 1,size = 7,colour = "black") +
  geom_text_repel(data=supDTFIvsEXPdf[supDTFIvsEXPdf$labels!="",], aes(limmanegLogPValues, negLogZPValues, label=labels), size = 10,box.padding = unit(1, 'lines')) + 
  ylab("Differential TF Involvement, -log(p-value)") + xlab("Differential Expression,  LIMMA -log(p-value)") + 
  labs(title=paste0("Differential Involvement vs Differential Expression \n(Smoker Controls to COPD Patients )")) +
  scale_colour_gradient2(limits=c(-.02,.02), oob = scales::squish, name="log(FC)", low = "blue", high = "yellow", mid="white") +
  theme_classic() + 
  theme(plot.title = element_text(size=60,hjust=.5), axis.text=element_text(size=25), axis.title=element_text(size=50), 
        strip.text.x = element_text(size=40),
        legend.text=element_text(size=30), legend.title=element_text(size=40, vjust=0), legend.key.size=unit(1,"in"))

cairo_pdf(file.path(outputDir,paste('dTFIvsLIMMAall.pdf', sep="")), width=24, height=24)
print(dTFI_LIMMA_all)
dev.off()

# Figure 3 a-b ------------------------------------------------------------

source('~/gd/Harvard/Research/R_workspace/consolidateResultTables.R')


# Combine plots into figures ----------------------------------------------

outputDir <- "~/gd/Harvard/Research/R_workspace/TM_manuscript/figures/"
setwd(outputDir)
system('pdflatex figure2.tex')
system('pdflatex figure3.tex')
system('pdfcrop figure2.pdf figure2.pdf')
system('pdfcrop figure3.pdf figure3.pdf')
system('pdfcrop "dTFI\ vs\ LIMMA\ ECLIPSE.pdf" figure4.pdf')


# Create pngs for google docs --------------------------------------------

system('convert -density 1200 figure2.pdf figure2.png')
system('convert -density 1200 figure3.pdf -resize 50% figure3.png')
system('convert -density 1200 figure4.pdf -resize 50% figure4.png')


# Supplement --------------------------------------------------------------

system('pdflatex figure2ECLIPSE.tex')
system('pdflatex figure2COPDGene.tex')
system('pdflatex figure2LGRC.tex')
system('pdflatex figure2LTCOPD.tex')
system('pdfcrop figure2ECLIPSE.pdf figure2ECLIPSE.pdf')
system('pdfcrop figure2COPDGene.pdf figure2COPDGene.pdf')
system('pdfcrop figure2LGRC.pdf figure2LGRC.pdf')
system('pdfcrop figure2LTCOPD.pdf figure2LTCOPD.pdf')

system('pdfcrop "dTFI\ vs\ LIMMA\ ECLIPSE.pdf" figure4ECLIPSE.pdf')
system('pdfcrop "dTFI\ vs\ LIMMA\ COPDGene.pdf" figure4COPDGene.pdf')
system('pdfcrop "dTFI\ vs\ LIMMA\ LGRC.pdf" figure4LGRC.pdf')
system('pdfcrop "dTFI\ vs\ LIMMA\ LTCOPD.pdf" figure4LTCOPD.pdf')

