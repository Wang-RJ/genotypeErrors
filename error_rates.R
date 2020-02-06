# Calculate error rates from owl monkey genotype-phase combinations

library(ggplot2)

# Start in directory with genotype-phase combination .txt files

child_ids <- as.character(c(39103, 39108, 39117, 39122))

owl_counts <- list()
for(id in child_ids) {
  files <- list.files()[grep("genophase", list.files())]
  owl_counts[[id]] <- read.table(files[grep(id, files)], header = TRUE, check.names = FALSE)
}

strip_phase <- function(gp_combos) {
  f_sims <- gp_combos[grep("F", names(gp_combos))]
  m_sims <- gp_combos[grep("M", names(gp_combos))]
  
  names(f_sims) <- gsub("F", "", names(f_sims))
  names(m_sims) <- gsub("M", "", names(m_sims))
  
  return(colSums(rbind(f_sims, m_sims)))
}

F_violation_genos <- c("01101", "01112", "01122", "21100", "21110", "21121", "02101", "02112", "02122", "20121", "20110", "20100", "12122", "12112", "12101", "10121", "10110", "10100")
M_violation_genos <- c("01100", "01110", "01121", "21101", "21112", "21122", "02100", "02110", "02121", "20122", "20112", "20101", "12121", "12110", "12100", "10122", "10112", "10101")

parse_violationsF <- function(row) { row[paste(F_violation_genos, "F", sep = "")] }
parse_violationsM <- function(row) { row[paste(M_violation_genos, "M", sep = "")] }

calc_Mmatrix <- function(gcombo_vector) {
  G <- gcombo_vector
  M <- matrix(c(
G["01100"],     0.5*(G["11101"]+G["01111"]), G["21101"] + G["01121"], 0,                     0,                           0,
0,              0.5*G["11112"],              G["21112"],              G["01110"],            0.5*G["01111"],              0,
0,              0.5*G["11122"],              G["21122"],              0,                     G["01121"],                  0,
0,              G["21101"],                  0,                       G["01100"],            0.5*G["11100"],              0,
0,              0.5*G["21111"],              G["21112"],              G["01110"],            0.5*G["11110"],              0,
0,              0,                           0,                       G["01121"]+G["21101"], 0.5*(G["21111"]+G["11121"]), G["21122"],
G["02100"],     0.5*G["02111"],              G["02121"],              0,                     0,                           0,
0,              0,                           0,                       G["02110"],            0.5*G["02111"],              0,
0,              0,                           0,                       0,                     G["02121"],                  0,
0,              0,                           0,                       G["20101"],            0.5*G["20111"],              G["20122"],
0,              0.5*G["20111"],              G["20112"],              0,                     0,                           0,
0,              G["20101"],                  0,                       0,                     0,                           0,    
0,              0,                           0,                       G["10122"],            G["12121"]+0.5*G["11122"],   0.5*G["12222"],
0,              0,                           0,                       G["12110"]+G["10112"], 0.5*(G["12111"]+G["11112"]), 0.5*G["12212"],
G["12100"],     0.5*G["12111"],              G["12121"],              G["10101"],            0.5*G["11101"],              0.5*G["12201"],
0.5*G["10021"], 0.5*G["11121"],              G["12121"],              G["10101"],            0.5*G["10111"],              G["10122"],
0.5*G["10010"], 0.5*(G["11110"]+G["10111"]), G["10112"]+G["12110"],   0,                     0,                           0,
0.5*G["10000"], G["10101"]+0.5*G["11100"],   G["12100"],              0,                     0,                           0
  ),
  ncol = 6, byrow = TRUE)
  
  return(M)
}

extract_oopVec <- function(gp_combos) {
  F_violations <- parse_violationsF(gp_combos)
  M_violations <- parse_violationsM(gp_combos)
  
  oop_vec <- as.vector(unlist(cbind(F_violations, M_violations)))
  return(oop_vec)
}

oop_matrices <- lapply(owl_counts, function(family_matrix) {
  oop_matrix <- t(apply(family_matrix, 1, extract_oopVec))
  return(oop_matrix)
})

optimize_errs <- function(gp_combos, W_vec) {
  stripped_genos <- strip_phase(gp_combos)
  
  M_a <- calc_Mmatrix(stripped_genos)
  M_b <- calc_Mmatrix(stripped_genos)[18:1,]
  
  M <- rbind(M_a, M_b)
  
  oop_vec <- extract_oopVec(gp_combos)
  
  return(optim(par = rep(1e-3, 6), function(x) { return(sum((oop_vec - M %*% x)^2)) }, method = "L-BFGS-B", lower = 0)$par)
}

errmat_list <- lapply(1:4, function(j) {
  return(t(sapply(1:4, function(i) { return(optimize_errs(owl_counts[[j]][i,], matrix(rep(1,144), ncol = 36)[i,])) })))
})

# Rebuild as data frame for ggplot
err_df <- lapply(errmat_list, function(mat) { return(data.frame(estimate = as.vector(mat),
                                                                err = as.vector(sapply(c("0>1", "1>0", "2>0", "0>2", "1>2", "2>1"), rep, 4)),
                                                                GQ = rep(c(0,20,40,60), 6))) })

err_df <- do.call(rbind, err_df)

ggplot(err_df, aes(x = factor(err, levels = c("0>1", "1>0", "2>0", "0>2", "1>2", "2>1")), y = estimate, fill = factor(GQ))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_dodge(0.4)) +
  scale_fill_brewer(type = "seq", palette = 16, direction = -1) +
  theme_bw()
