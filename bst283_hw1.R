seqres <- read.table("~/gd/Harvard/BST283/BST283_Hw1/var_file_problem_1.txt", header=T)


e_i <- .0001
theta_t <- 6.3

seqres$t_af <- seqres$t_alt_count/(seqres$t_ref_count+seqres$t_alt_count)
seqres$n_af <- seqres$n_alt_count/(seqres$n_ref_count+seqres$n_alt_count)

seqres$LMFM <- (seqres$t_af*(e_i/3) + (1-seqres$t_af)*(1-e_i))^seqres$t_ref_count *
  (seqres$t_af*(1-e_i)+(1-seqres$t_af)*(e_i/3))^seqres$t_alt_count
seqres$LM0 <- (1-e_i)^seqres$t_ref_count * (1-seqres$t_af)*(e_i/3)^seqres$t_alt_count 
seqres$LODT <- log((seqres$LMFM)/(seqres$LM0),10)

seqres$Candidates <- seqres$LODT> theta_t
seqres$Candidates[is.na(seqres$Candidates)] <- T
sum(seqres$Candidates)


seqres$LMFM_N <- (seqres$n_af*(e_i/3) + (1-seqres$n_af)*(1-e_i))^seqres$n_ref_count *
  (seqres$n_af*(1-e_i)+(1-seqres$n_af)*(e_i/3))^seqres$n_alt_count
seqres$LM0_N <- (1-e_i)^seqres$n_ref_count * (1-seqres$n_af)*(e_i/3)^seqres$n_alt_count 
seqres$LODT_N <- log(seqres$LMFM_N/seqres$LM0_N,10)

seqres$GermlineMuts <- seqres$LODT_N>theta_t
seqres$GermlineMuts[is.na(seqres$GermlineMuts)] <- T
sum(seqres$GermlineMuts)

seqres$somatic <- seqres$Candidates&(!seqres$GermlineMuts)

somaticVariants <- seqres[seqres$somatic,]
somaticVariants[order(-somaticVariants$t_af),]
sum(seqres$somatic)

