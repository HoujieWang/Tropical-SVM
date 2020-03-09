algorithm4 = function(assignment, data, tst_data, n, ntst, beta, method_ind){
  ip = assignment[1]; iq = assignment[3]; j = assignment[2]
  f.obj = c(1, 0, 0, 0, rep(-1, 6*n)*beta)
  f.conp = rbind(cbind(rep(1, n), rep(-1, n), rep(0, n), rep(1, n)), cbind(rep(0, n), rep(-1, n), rep(0, n), rep(1, n)), cbind(rep(0, n), rep(0, n), rep(1, n), rep(-1, n)))
  f.conq = rbind(cbind(rep(1, n), rep(0, n), rep(-1, n), rep(1, n)), cbind(rep(0, n), rep(0, n), rep(-1, n), rep(1, n)), cbind(rep(0, n), rep(1, n), rep(0, n), rep(-1, n)))
  f.con = cbind(rbind(f.conp, f.conq), diag(-1, nrow = 6*n, ncol = 6*n))
  f.rhs = c(rep(data[1: n, ip] - data[1: n, j], 2), data[1: n, j] - data[1: n, iq], rep(data[-c(1: n), iq] - data[-c(1: n), j], 2), data[-c(1: n), j] - data[-c(1: n), ip])
  f.dir=rep("<=",nrow(f.con))
  omega = rep(0, ncol(tst_data))
  omega[c(ip, iq, j)] = lp("max", f.obj, f.con, f.dir, f.rhs)$solution[2: 4]
  omega[-c(ip, iq, j)] = apply(data[, -c(ip, iq, j)], 2, function(x){min(data[, j] - x) + omega[j]})
  classification = apply(t(apply(tst_data, 1, function(x){x+omega})), 1, function(x){which(x == max(x))})
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
      v = c(ip, iq, j) %in% x
      return(sum(colSums(t(P) == v) == ncol(P)))
    }), sapply(classification[-c(1: ntst)], function(x){
      v = c(ip, iq, j) %in% x
      return(sum(colSums(t(Q) == v) == ncol(Q)))
    })))/length(classification))
  }
  return(accuracy)
}
