algorithm2 = function(assignment, data, tst_data, n, ntst, beta, method_ind){
  k = assignment[1]; iq = assignment[2]; jp = assignment[3]
  f.obj = c(1, 0, 0, 0, rep(-1, 6*n)*beta)
  f.conp = rbind(cbind(rep(1, n), rep(-1, n), rep(1, n), rep(0, n)), cbind(rep(0, n), rep(-1, n), rep(1, n), rep(0, n)), cbind(rep(0, n), rep(0, n), rep(-1, n), rep(1, n)))
  f.conq = rbind(cbind(rep(1, n), rep(1, n), rep(0, n), rep(-1, n)), cbind(rep(0, n), rep(1, n), rep(0, n), rep(-1, n)), cbind(rep(0, n), rep(-1, n), rep(1, n), rep(0, n)))
  f.con = cbind(rbind(f.conp, f.conq), diag(-1, nrow = 6*n, ncol = 6*n))
  f.rhs = c(rep(data[1: n, k] - data[1: n, jp], 2), data[1: n, jp] - data[1: n, iq], rep(data[-c(1: n), iq] - data[-c(1: n), k], 2), data[-c(1: n), k] - data[-c(1: n), jp])
  f.dir=rep("<=",nrow(f.con))
  omega = rep(0, ncol(data))
  sol = lp("max", f.obj, f.con, f.dir, f.rhs)$solution
  omega[c(k, jp, iq)] = sol[2: 4]
  omega[-c(k, jp, iq)] = apply(data[, -c(k, jp, iq)], 2, function(x){min(c(min(data[1: n, jp] - x[1:n])+omega[jp], min(data[-c(1: n), k] - x[-c(1:n)])+omega[k]))})
  classification = lapply(as.list(as.data.frame(apply(tst_data, 1, function(x){x+omega}))), function(x){which(x == max(x))})
  PQ_com = matrix(c(1, 0, 0,
                    0, 1, 0,
                    0, 0, 1,
                    1, 1, 0,
                    1, 0, 1,
                    0, 1, 1,
                    1, 1, 1,
                    0, 0, 0), ncol = 3, byrow = T)
  ind_matrix = combinations(8, 4)
  accuracy = c()
  for (i in method_ind){
    P = PQ_com[ind_matrix[i, ], ]; Q = PQ_com[-ind_matrix[i, ], ]
    accuracy = c(accuracy, sum(c(sapply(classification[1: ntst], function(x){
      v = c(k, jp, iq) %in% x;
      return(sum(colSums(t(P) == v) == ncol(P)))
    }), sapply(classification[-c(1: ntst)], function(x){
      v = c(k, jp, iq) %in% x;
      return(sum(colSums(t(Q) == v) == ncol(Q)))
    })))/length(classification))
  }
  return(accuracy)
}
