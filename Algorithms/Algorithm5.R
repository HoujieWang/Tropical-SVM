# omega is the optimal normal vector
# Please input assignment in vector: c(i_P, j_P, i_Q, j_Q)
algorithm5 = function(omega, tst_P, tst_Q, assignment){
  classification = apply(t(apply(rbind(tst_P, tst_Q), 1, function(x){x + omega})), 1, function(x){which.max(x)})
  return((sum(classification[1:nrow(tst_P)] == assignment[1])+sum(classification[-c(1:nrow(tst_P))] == assignment[3]))/(nrow(tst_P) + nrow(tst_Q)))
}