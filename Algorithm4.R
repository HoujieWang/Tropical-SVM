# Please input assignment in vector: c(i_P, j_P, i_Q, j_Q)
# Please input P and Q in matrix
algorithm4 = function(assignment, P, Q){
  options(warn = -1)
  i_P = assignment[1]; j_P = assignment[2]; i_Q = assignment[3]; j_Q = assignment[4] 
  p1 = P[, j_P] - P[, i_P]; q1 = Q[, j_P] - Q[, i_P]
  p_ind = order(p1, decreasing = T)[-which(sort(q1) - sort(p1, decreasing = T) < 0)]
  q_ind = order(q1)[-which(sort(q1) - sort(p1, decreasing = T) < 0)]
  if (length(p_ind) >= 10 & length(q_ind) >= 10){
    Psub = P[p_ind, ]; Qsub = Q[q_ind, ]
    Ptst_sub = sample(1: nrow(Psub), size = floor(nrow(Psub)/4)); Qtst_sub = sample(1: nrow(Qsub), size = floor(nrow(Qsub)/4))
    trn_P = Psub[-Ptst_sub, ]; trn_Q = Qsub[-Qtst_sub, ]; tst_P = Psub[Ptst_sub, ]; tst_Q = Qsub[Qtst_sub, ]
    
    z = (min(trn_Q[, j_P] - trn_Q[, i_P]) + min(trn_P[, i_P] - trn_P[, j_P]))/2
    omega = rep(0, ncol(P))
    if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, i_P]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) <
       (min(trn_Q[, j_P] - trn_Q[, i_P]) - min(trn_P[, i_P] - trn_P[, j_P]))/2){
      omega[c(i_P, j_P)] = c((min(trn_Q[, j_P] - trn_Q[, i_P]) - min(trn_P[, i_P] - trn_P[, j_P]))/2 + max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]), max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]))
      
    }
    if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, i_P]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) >=
       (min(trn_Q[, j_P] - trn_Q[, i_P]) - min(trn_P[, i_P] - trn_P[, j_P]))/2){
      omega[c(i_P, j_P)] = c(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, i_P]), max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, i_P]) - (min(trn_Q[, j_P] - trn_Q[, i_P]) - min(trn_P[, i_P] - trn_P[, j_P]))/2)
    }
    options(warn = 0)
    return(list("Train P" = trn_P, "Train Q" = trn_Q, "Test P" = tst_P, "Test Q" = tst_Q, "optimal margin" = z, "optimal normal vector" = omega, "assignment" = assignment))
  } else {return(0)}
}