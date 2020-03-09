algorithm1 = function(assignment, data, tst_data, n, ntst, beta, method_ind){
  ip = assignment[1]; iq = assignment[2]; jp = assignment[3]; jq = assignment[4]
  f.obj = c(1, rep(0, 4), c(rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n), rep(-1, n))*beta)
  f.conp = rbind(cbind(rep(1, n), rep(-1, n), rep(1, n), rep(0, n), rep(0, n)), cbind(rep(0, n), rep(-1, n), rep(1, n), rep(0, n), rep(0, n)), 
                 cbind(rep(0, n), rep(0, n), rep(-1, n), rep(1, n), rep(0, n)), cbind(rep(0, n), rep(0, n), rep(-1, n), rep(0, n), rep(1, n)))
  f.conq = rbind(cbind(rep(1, n), rep(0, n), rep(0, n), rep(-1, n), rep(1, n)), cbind(rep(0, n), rep(0, n), rep(0, n), rep(-1, n), rep(1, n)), 
                 cbind(rep(0, n), rep(1, n), rep(0, n), rep(0, n), rep(-1, n)), cbind(rep(0, n), rep(0, n), rep(1, n), rep(0, n), rep(-1, n)))
  f.con = cbind(rbind(f.conp, f.conq), diag(-1, nrow = 8*n, ncol = 8*n))
  f.rhs = c(rep(data[1: n, ip] - data[1: n, jp], 2), data[1: n, jp] - data[1: n, iq], data[1: n, jp] - data[1: n, jq], 
            rep(data[-c(1: n), iq] - data[-c(1: n), jq], 2), data[-c(1: n), jq] - data[-c(1: n), ip], data[-c(1: n), jq] - data[-c(1: n), jp])
  f.dir=rep("<=",nrow(f.con))
  omega = rep(0, ncol(data))
  omega[c(ip, jp, iq, jq)] = lp("max", f.obj, f.con, f.dir, f.rhs)$solution[2: 5]
  omega[-c(ip, jp, iq, jq)] = apply(data[, -c(ip, jp, iq, jq)], 2, function(x){min(c(min(data[1: n, jp] - x[1:n])+omega[jp], min(data[-c(1: n), jq] - x[-c(1:n)])+omega[jq]))})
  classification = lapply(as.list(as.data.frame(apply(tst_data, 1, function(x){x+omega}))), function(x){which(x == max(x))})
  P_base = matrix(c(1, 0, 0, 0,
                    0, 1, 0, 0,
                    1, 1, 0, 0, 
                    1, 1, 1, 1), ncol = 4, byrow = T); 
  Q_base = matrix(c(0, 0, 1, 0,
                    0, 0, 0, 1,
                    0, 0, 1, 1,
                    0, 0, 0, 0), ncol = 4, byrow = T); 
  PQ_com = matrix(c(1, 0, 1, 0,
                    1, 0, 0, 1,
                    0, 1, 1, 0,
                    0, 1, 0, 1,
                    1, 1, 1, 0,
                    1, 1, 0, 1,
                    1, 0, 1, 1,
                    0, 1, 1, 1), ncol = 4, byrow = T)
  ind_matrix = combinations(8, 4)
  accuracy = c()
  for (l in method_ind){
    P = rbind(P_base, PQ_com[ind_matrix[l, ], ]); Q = rbind(Q_base, PQ_com[-ind_matrix[l, ], ])
    accuracy = c(accuracy, sum(c(sapply(classification[1: ntst], function(x){
      v = c(ip, jp, iq, jq) %in% x;
      return(sum(colSums(t(P) == v) == ncol(P)))
    }), sapply(classification[-c(1: ntst)], function(x){
      v = c(ip, jp, iq, jq) %in% x;
      return(sum(colSums(t(Q) == v) == ncol(Q)))
    })))/length(classification))
  }
  return(accuracy)
}
