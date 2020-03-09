algorithm3 = function(assignment, data, tst_data, n, ntst, beta){
  k1 = assignment[1]; k2 = assignment[2]
  omega = rep(0, ncol(tst_data))
  f.obj = c(1, 0, 0, rep(-1, 4*n)*beta)
  f.conp = rbind(cbind(rep(1, n), matrix(rep(c(-1, 1), n), nrow = n, ncol = 2, byrow = T)), 
                 cbind(rep(0, n), matrix(rep(c(-1, 1), n), nrow = n, ncol = 2, byrow = T)))
  f.conq = rbind(cbind(rep(1, n), matrix(rep(c(1, -1), n), nrow = n, ncol = 2, byrow = T)), 
                 cbind(rep(0, n), matrix(rep(c(1, -1), n), nrow = n, ncol = 2, byrow = T)))
  f.con = cbind(rbind(f.conp, f.conq), diag(-1, nrow = 4*n, ncol = 4*n))
  f.rhs = c(rep(data[1: n, k1] - data[1: n, k2], 2), rep(data[-c(1: n), k2] - data[-c(1: n), k1], 2))
  f.dir=rep("<=",nrow(f.con))
  sol = lp("max", f.obj, f.con, f.dir, f.rhs)$solution
  omega[c(k1, k2)] = sol[2: 3]
  omega[-c(k1, k2)] = apply(data[, -c(k1, k2)], 2, function(x){min(c(min(data[, k2] - x)+omega[k2], min(data[, k1] - x)+omega[k1]))})
  classification = apply(t(apply(tst_data, 1, function(x){x+omega})), 1, which.max)
  return((sum(classification[1:ntst] == k1)+sum(classification[-c(1: ntst)] == k2))/length(classification))
}
