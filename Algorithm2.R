# Please input assignment in vector: c(i_P, j_P, i_Q, j_Q)
# Please input P and Q in matrix
algorithm2 = function(assignment, P, Q){
  options(warn = -1)
  i_P = assignment[1]; j_P = assignment[2]; i_Q = assignment[3]; j_Q = assignment[4] 
  p1 = P[, j_P] - P[, i_P]; p2 = P[, j_P] - P[, i_Q]
  q1 = Q[, j_Q] - Q[, i_Q]; q2 = Q[, j_P] - Q[, i_P]
  maxP = -which(sort(-q1) - sort(p1, decreasing = T) - sort(-p2, decreasing = T) < 0)
  maxQ = -which(sort(-q2) + sort(-q1) - sort(-p2, decreasing = T) < 0)
  p_ind1 = order(p1, decreasing = T)[maxP]; p_ind2 = order(-p2, decreasing = T)[maxP]
  q_ind1 = order(-q1)[maxP]
  q_ind2 = order(-q2)[maxQ]; q_ind3 = order(-q1)[maxQ]
  p_ind3 = order(-p2, decreasing = T)[maxQ]
  p_ind = Reduce(intersect, list(p_ind1, p_ind2, p_ind3)); q_ind = Reduce(intersect, list(q_ind1, q_ind2, q_ind3))
  if (length(p_ind) >= 10 & length(q_ind) >= 10){
    Psub = P[p_ind, ]; Qsub = Q[q_ind, ]
    Ptst_sub = sample(p_ind, size = floor(nrow(Psub)/4)); Qtst_sub = sample(q_ind, size = floor(nrow(Qsub)/4))
    Ptrn_sub = p_ind[!p_ind %in% Ptst_sub]; Qtrn_sub = q_ind[!q_ind %in% Qtst_sub]
    trn_P = P[Ptrn_sub, ]; trn_Q = Q[Qtrn_sub, ]; tst_P = P[Ptst_sub, ]; tst_Q = Q[Qtst_sub, ]
    A1 = min(-p1[Ptrn_sub]); B = min(p2[Ptrn_sub])
    A2 = min(-q2[Qtrn_sub]); C = min(-q1[Qtrn_sub])
    if ((A1 + B + C)/2 <= A2 + B + C){
      if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) <= (B + C - A1)/2){
        omega_j_P = max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]); omega_j_Q = omega_j_P + (B + C - A1)/2; omega_i_Q = omega_j_P + B
      } else {
        omega_j_Q = max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]); omega_j_P = omega_j_Q - (B + C - A1)/2; omega_i_Q = omega_j_P + B
      }
      z = (A1 + B + C)/2
    } 
    if (A2 + B + C <= (A1 + B + C)/2) {
      if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) <= -A2){
        omega_j_P = max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]); omega_j_Q = omega_j_P - A2; omega_i_Q = omega_j_P + B
      } else {
        omega_j_Q = max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]); omega_j_P = omega_j_Q + A2; omega_i_Q = omega_j_P + B
      }
      z = A2 + B + C
    }
    omega = rep(0, ncol(P))
    omega[assignment[-4]] = c(omega_j_Q, omega_j_P, omega_i_Q)
    options(warn = 0)
    return(list("Train P" = trn_P, "Train Q" = trn_Q, "Test P" = tst_P, "Test Q" = tst_Q, "optimal margin" = z, "optimal normal vector" = omega, "assignment" = assignment))
  } else {return(0)}
}