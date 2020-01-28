# Please input assignment in vector: c(i_P, j_P, i_Q, j_Q)
algorithm1 = function(assignment, P, Q){
  options(warn = -1)
  i_P = assignment[1]; j_P = assignment[2]; i_Q = assignment[3]; j_Q = assignment[4] 
  p1 = P[, i_P] - P[, j_P]; p2 = P[, j_P] - P[, j_Q]
  p3 = P[, i_Q] - P[, j_P]; p5 = P[, j_Q] - P[, i_Q]
  q1 = Q[, i_Q] - Q[, j_Q]; q2 = Q[, j_Q] - Q[, j_P]
  q5 = Q[, i_P] - Q[, j_P]
  p_ind1 = order(p2)[-which(sort(p2) + sort(q2) < 0)]
  q_ind1 = order(q2)[-which(sort(p2) + sort(q2) < 0)]
  p_ind2 = order(p1)[-which(sort(p1) - sort(q5, decreasing = T) < 0)]
  q_ind2 = order(q5, decreasing = T)[-which(sort(p1) - sort(q5, decreasing = T) < 0)]
  p_ind3 = order(p5)[-which(sort(p5) - sort(-q1, decreasing = T) < 0)]
  q_ind3 = order(-q1, decreasing = T)[-which(sort(p5) - sort(-q1, decreasing = T) < 0)]
  p_ind = Reduce(intersect, list(p_ind1, p_ind2, p_ind3)); q_ind = Reduce(intersect, list(q_ind1, q_ind2, q_ind3))
  if (length(p_ind) >= 10 & length(q_ind) >= 10){
    #Psub = P[p_ind, ]; Qsub = Q[q_ind, ]
    Ptst_sub = sample(p_ind, size = floor(nrow(Psub)/4)); Qtst_sub = sample(q_ind, size = floor(nrow(Qsub)/4))
    Ptrn_sub = p_ind[!p_ind %in% Ptst_sub]; Qtrn_sub = q_ind[!q_ind %in% Qtst_sub]
    trn_P = P[Ptrn_sub, ]; trn_Q = Q[Qtrn_sub, ]; tst_P = P[Ptst_sub, ]; tst_Q = Q[Qtst_sub, ]
    A = min(p1[Ptrn_sub]); B = min(trn_P[, j_P] - trn_P[, i_Q]); C = min(p2[Ptrn_sub])
    D = min(q1[Qtrn_sub]); E = min(trn_Q[, j_Q] - trn_Q[, i_P]); FF = min(q2[Qtrn_sub])
    if((A+B+D+E)/2 <= min(c(A+C+E, B+D+FF))){
      if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) <= (B+D-E-A)/2){
        omega_j_P = max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]); omega_j_Q = omega_j_P + (B+D-E-A)/2; omega_i_Q = omega_j_P + B; omega_i_P = omega_j_Q + E
      } else {
        omega_j_Q = max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]); omega_j_P = omega_j_Q - (B+D-E-A)/2; omega_i_Q = omega_j_P + B; omega_i_P = omega_j_Q + E
      }
      z = (A+B+D+E)/2
    } 
    if (A+C+E <= min((A+B+D+E)/2, B+D+FF)){
      if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) <= C){
        omega_j_P = max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]); omega_j_Q = omega_j_P + C; omega_i_P = omega_j_Q + E; omega_i_Q = omega_j_Q + A + C + E - D
      } else {
        omega_j_Q = max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]); omega_j_P = omega_j_Q - C; omega_i_P = omega_j_Q + E; omega_i_Q = omega_j_Q + A + C + E - D
      }
      z = A+C+E
    }
    if (B+D+FF <= min((A+B+D+E)/2, A+C+E)){
      if(max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]) - max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]) <= -FF){
        omega_j_P = max(trn_P[,-c(i_P,i_Q)]- trn_P[, j_P]); omega_j_Q = omega_j_P - FF; omega_i_P = omega_j_P + B + D + FF - A; omega_i_Q = omega_j_Q + B + FF
      } else {
        omega_j_Q = max(trn_Q[,-c(i_P,i_Q)]- trn_Q[, j_Q]); omega_j_P = omega_j_Q + FF; omega_i_P = omega_j_P + B + D + FF - A; omega_i_Q = omega_j_Q + B + FF
      }
      z = B+D+FF
    }
    omega = rep(0, ncol(P))
    omega[assignment] = c(omega_i_P, omega_j_P, omega_i_Q, omega_j_Q)
    options(warn = 0)
    return(list("Train P" = trn_P, "Train Q" = trn_Q, "Test P" = tst_P, "Test Q" = tst_Q, "optimal margin" = z, "optimal normal vector" = omega, "assignment" = assignment))
  } else {return(0)}
}