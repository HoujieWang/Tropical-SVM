# Please input assignment in vector: c(i_P, j_P, i_Q, j_Q)
# Please input P and Q in matrix
algorithm3 = function(assignment, P, Q){
  options(warn = -1)
  i_P = assignment[1]; j_P = assignment[2]; i_Q = assignment[3]; j_Q = assignment[4] 
  p1 = P[, i_P] - P[, j_P]; p2 = P[, i_Q] - P[, j_Q]
  q1 = Q[, i_Q] - Q[, j_Q]; q2 = Q[, i_P] - Q[, j_P]
  p_ind1 = order(p1)[-which(sort(p1) - sort(q2, decreasing = T) < 0)]
  q_ind1 = order(q2, decreasing = T)[-which(sort(p1) - sort(q2, decreasing = T) < 0)]
  p_ind2 = order(p2, decreasing = T)[-which(sort(q1) - sort(p2, decreasing = T) < 0)]
  q_ind2 = order(q1)[-which(sort(q1) - sort(p2, decreasing = T) < 0)]
  p_ind = intersect(p_ind1, p_ind2); q_ind = intersect(q_ind1, q_ind2)
  if (length(p_ind) >= 10 & length(q_ind) >= 10){
    Psub = P[p_ind, ]; Qsub = Q[q_ind, ]
    Ptst_sub = sample(1: nrow(Psub), size = floor(nrow(Psub)/4)); Qtst_sub = sample(1: nrow(Qsub), size = floor(nrow(Qsub)/4))
    trn_P = Psub[-Ptst_sub, ]; trn_Q = Qsub[-Qtst_sub, ]; tst_P = Psub[Ptst_sub, ]; tst_Q = Qsub[Qtst_sub, ]
    z = min(c(min(trn_Q[, j_P] - trn_Q[, i_P]) + min(trn_P[, i_P] - trn_P[, j_P]), 
              min(trn_Q[, i_Q] - trn_Q[, j_P]) + min(trn_P[, j_P] - trn_P[, i_Q])))
    omega = rep(0, ncol(P))
    omega[j_P] = max(rbind(trn_P, trn_Q)[,-c(i_P, i_Q)]- rbind(trn_P, trn_Q)[, j_P])
    omega[c(i_P, i_Q)] = c(omega[j_P] + min(trn_Q[, j_P] - trn_Q[, i_P]), omega[j_Q] + min(trn_P[, j_Q] - trn_P[, i_Q]))
    options(warn = 0)
    return(list("Train P" = trn_P, "Train Q" = trn_Q, "Test P" = tst_P, "Test Q" = tst_Q, "optimal margin" = z, "optimal normal vector" = omega, "assignment" = assignment))
  } else {return(0)}
}