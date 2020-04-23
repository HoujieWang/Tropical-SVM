graph_producer = function(ntst){
  library(ggplot2); library(e1071); library(lpSolve); library(ape); library(gtools); library(parallel)
  full_n = 100;
  d = 10;
  beta = 1
  n = full_n - ntst
  method_ind_1 = c(19, 31); method_ind_4 = c(31, 36); method_ind_2 = c(3, 31)
  # Algorithms
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
  
  # Load the data and assignments
  load(paste("data_",ntst, "%.RData", sep = ""))
  for (i in paste(paste("asgn_", 1: 4, "_", sep = ""), ntst, "%.RData", sep = "")){load(i)}
  
  accuracy3 = c()
  accuracy4 = c()
  accuracy2 = c()
  accuracy1 = c()
  accuracy0 = c()
  for (i in 1: 12){
    # print(i)
    a3 = assignmen3[[i]]
    a4 = assignmen4[[i]]
    a2 = assignmen2[[i]]
    a1 = assignmen1[[i]]
    temp3= c()
    temp4= c()
    temp2= c()
    temp1= c()
    temp0= c()
    for (j in 1: 10){
      # print(j)
      P_tst = data_matrix_list[[i]][[j]][[3]]; Q_tst = data_matrix_list[[i]][[j]][[4]]
      P_trn = data_matrix_list[[i]][[j]][[1]]; Q_trn = data_matrix_list[[i]][[j]][[2]]
      data = rbind(P_trn, Q_trn); tst_data = rbind(P_tst, Q_tst)
      if (i <= 8){
        temp3 = c(temp3, algorithm3(a3[j, ], data = data, tst_data = tst_data, n, ntst, beta))
        temp4 = c(temp4, algorithm4(a4[j, ], data = data, tst_data = tst_data, n, ntst, beta, method_ind_4[2]))
        temp2 = c(temp2, algorithm2(a2[j, ], data = data, tst_data = tst_data, n, ntst, beta, method_ind_2[1]))
        temp1 = c(temp1, algorithm1(a1[j, ], data = data, tst_data = tst_data, n, ntst, beta, method_ind_1[1]))
      }
      if (i > 8){
        temp3 = c(temp3, algorithm3(a3[j, ], data = data, tst_data = tst_data, n, ntst, beta))
        temp4 = c(temp4, algorithm4(a4[j, ], data = data, tst_data = tst_data, n, ntst, beta, method_ind_4[1]))
        temp2 = c(temp2, algorithm2(a2[j, ], data = data, tst_data = tst_data, n, ntst, beta, method_ind_2[2]))
        temp1 = c(temp1, algorithm1(a1[j, ], data = data, tst_data = tst_data, n, ntst, beta, method_ind_1[2]))
      }
      fit = svm(data, as.factor(c(rep("P", n), rep("Q", (nrow(data)-n)))))
      temp0 = c(temp0, sum(diag(table(predict(fit, tst_data), c(rep("P", ntst), rep("Q", (nrow(tst_data)-ntst))))))/nrow(tst_data))
    }
    accuracy3 = c(accuracy3, max(temp3))
    accuracy4 = c(accuracy4, max(temp4))
    accuracy2 = c(accuracy2, max(temp2))
    accuracy1 = c(accuracy1, max(temp1))
    accuracy0 = c(accuracy0, max(temp0))
  }
  result_table = cbind.data.frame("C" = c(seq(0.2,1,by=0.2), seq(1.2,6,by=1.2), 8, 10),
                                  "Algorithm 1" = accuracy1,
                                  "Algorithm 2" = accuracy2, 
                                  "Algorithm 3" = accuracy3,
                                  "Algorithm 4" = accuracy4, 
                                  "classical" = accuracy0)
  result_table = cbind.data.frame(c(seq(0.2,1,by=0.2), seq(1.2,6,by=1.2), 8, 10),
                                  matrix(unlist(lapply(apply(result_table[, -1], 2, function(x){lowess(result_table[, 1], x)}), function(x){return(x[[2]])})), ncol = 5))
  colnames(result_table) = c("C", "Algorithm 1", "Algorithm 2", "Algorithm 3", "Algorithm 4", "classical")
  C = rep(result_table[, 1], 5); accu = as.vector(matrix(as.matrix(result_table[, -1]), ncol = 1)); 
  type = rep(colnames(result_table)[-1], each = nrow(result_table))
  result_table = cbind.data.frame(C, accu, type)
  result_table[, 2][result_table[, 2] > 1] = 1
  attach(result_table)
  ggplot(data = result_table, aes(x = C, y = accu, group = type)) +
    geom_line(aes(color = type, size = type)) + geom_point(aes(shape = type, size = type)) +
    scale_color_manual(values =  c("royalblue3", "turquoise3","goldenrod3", "springgreen4", "grey25")) +
    scale_size_manual(values = rep(1.5, 5)) + scale_shape_manual(values = rep(17, 5)) +
    labs(x = "C", y = "accuracy") + ggtitle(paste("accuracy vs. C  (", ntst, "% test data)", sep = ""))  + scale_x_continuous(breaks = seq(0, 10,by = 2))
}

# Specify the percent of test data, please choose among 15, 20 and 25.
graph_producer(ntst = 15)


