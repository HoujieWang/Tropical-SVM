library(lpSolve)
library(ape)
library(Parallel)

# Number of leaves of each tree
N = 5

# Dimension of a tree after transformed into a vector
d = choose(N, 2)

# We first begin with a simluated data set at C=6
P  = matrix(unlist(lapply(read.tree("LittleTree6A.txt"), function(x){cophenetic.phylo(x)[upper.tri(cophenetic.phylo(x))]})), ncol = d, byrow = T)
Q  = matrix(unlist(lapply(read.tree("LittleTree6B.txt"), function(x){cophenetic.phylo(x)[upper.tri(cophenetic.phylo(x))]})), ncol = d, byrow = T)

############################################## All Possible Assignments for All Theorems ############################################################

# Theorem 4.8
all_i1i2j1j2 = expand.grid(1:d, 1:d, 1:d, 1:d)
all_i1i2j1j2 = all_i1i2j1j2[apply(all_i1i2j1j2, 1, function(x){return(length(unique(x)) == length(x))}), ]
all_i1i2j1j2 = as.list(as.data.frame(t(all_i1i2j1j2)))

# Theorem 4.10
all_ijl =  expand.grid(1:d, 1:d, 1:d)
all_ijl = all_ijl[-which(all_ijl[,1] == all_ijl[,2]), ]
all_ijl = all_ijl[-which(all_ijl[,1] == all_ijl[,3]), ]
all_ijl = all_ijl[-which(all_ijl[,2] == all_ijl[,3]), ]
all_ijl1 = lapply(as.list(as.data.frame(t(all_ijl))), function(x){return(c(x, x[1]))})

# Theorem 4.11

all_ij = expand.grid(1:d, 1:d)
all_ij = all_ij[-which(all_ij[, 1] == all_ij[, 2]), ]
all_ij = lapply(as.list(as.data.frame(t(all_ij))), function(x){return(c(x, rev(x)))})

# Theorem 4.12 

all_i1i2j = expand.grid(1:d, 1:d, 1:d)
all_i1i2j = all_i1i2j[-which(all_i1i2j[,1] == all_i1i2j[,2]), ]
all_i1i2j = all_i1i2j[-which(all_i1i2j[,1] == all_i1i2j[,3]), ]
all_i1i2j = all_i1i2j[-which(all_i1i2j[,2] == all_i1i2j[,3]), ]
all_i1i2j = lapply(as.list(as.data.frame(t(all_i1i2j))), function(x){return(c(x, x[2]))})

############################################################## Algorithms 1~5 #######################################################################

# Algorithm 1~4 output 0 if size of feasible data set for P and Q is smaller than 10
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
    Psub = P[p_ind, ]; Qsub = Q[q_ind, ]
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
algorithm5 = function(omega, tst_P, tst_Q, assignment){
  classification = apply(t(apply(rbind(tst_P, tst_Q), 1, function(x){x + omega})), 1, function(x){which.max(x)})
  return((sum(classification[1:nrow(tst_P)] == assignment[1])+sum(classification[-c(1:nrow(tst_P))] == assignment[3]))/(nrow(tst_P) + nrow(tst_Q)))
}

############################################################# Test for Accuracy #####################################################################

# Accuracy for Theorem 4.8
result_4.8 = mclapply(all_i1i2j1j2, function(x){algorithm1(x, P, Q)})
sapply(result_4.8[sapply(result_4.8, function(x){class(x) == "list"})], function(x){algorithm5(x[[6]], x[[3]], x[[4]], x[[7]])})

# Accuracy for Theorem 4.10
result_4.10 = mclapply(all_ijl1, function(x){algorithm2(x, P, Q)})
sapply(result_4.10[sapply(result_4.10, function(x){class(x) == "list"})], function(x){algorithm5(x[[6]], x[[3]], x[[4]], x[[7]])})

# Accuracy for Theorem 4.11
result_4.11 = mclapply(all_i1i2j, function(x){algorithm3(x, P, Q)})
sapply(result_4.11[sapply(result_4.11, function(x){class(x) == "list"})], function(x){algorithm5(x[[6]], x[[3]], x[[4]], x[[7]])})

# Accuracy for Theorem 4.12
result_4.12 = mclapply(all_ij, function(x){algorithm4(x, P, Q)})
sapply(result_4.12[sapply(result_4.12, function(x){class(x) == "list"})], function(x){algorithm5(x[[6]], x[[3]], x[[4]], x[[7]])})
