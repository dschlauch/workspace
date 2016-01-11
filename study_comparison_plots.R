load("~/gd/Harvard/Research/TM_outputs/LTCOPD_bere_bare_92540/activeImage92540.RData")
ltcopd_TM <- transMatrices[[1]]
load("~/gd/Harvard/Research/TM_outputs/LGRC_bere_56432/activeImage56432.RData")
lgrc_TM <- transMatrices[[1]]
load("~/gd/Harvard/Research/TM_outputs/ECLIPSE_bere_bare_55557/activeImage55557.RData")
eclipse_TM <- transMatrices[[1]]
load("~/gd/Harvard/Research/TM_outputs/COPDGene_bere_70856/activeImage70856.RData")
COPDGene_TM <- transMatrices[[1]]

library(ggplot2)
qplot(c(ltcopd_TM),c(lgrc_TM)) + geom_smooth(method=lm)
cor(c(ltcopd_TM),c(lgrc_TM))
summary(lm(c(ltcopd_TM)~c(lgrc_TM)))

qplot(c(eclipse_TM),c(COPDGene_TM)) + geom_smooth(method=lm)
cor(c(eclipse_TM),c(COPDGene_TM))
summary(lm(c(eclipse_TM)~c(COPDGene_TM)))

qplot(c(lgrc_TM),c(COPDGene_TM)) + geom_smooth(method=lm)
cor(c(lgrc_TM),c(COPDGene_TM))
summary(lm(c(lgrc_TM)~c(COPDGene_TM)))
