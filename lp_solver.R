lp_solver = function(data, assignment, np, nq, flatten){
  ip = assignment[1]; jp = assignment[2]; iq = assignment[3]; jq = assignment[4] 
  d = ncol(data)
  f.obj = c(1, rep(0, d))
  f.conp = array(0, dim = c(d, (d+1), np))
  f.conq = array(0, dim = c(d, (d+1), nq))
  f.rhsp = array(0, dim = c(d, 1, np))
  f.rhsq = array(0, dim = c(d, 1, nq))
  f.con = matrix(0, nrow = d*(np+nq), ncol = (d+1))
  f.rhs = rep(0, d*(np+nq))
  for (k in 1: np){
    f.conp[1, 1, k]  = 1
    f.conp[, , k][as.matrix(expand.grid(1:2, 1+c(jp, ip)))] = c(1,1,-1,-1)
    f.conp[-c(1:2), (1+jp), k] = rep(-1, (d-2))
    f.conp[, , k][cbind(c(3:d), 1+c(1:d)[-c(ip, jp)])] = 1
    f.con[(1+d*(k-1)):(d*k), ] = f.conp[, ,k]
    f.rhsp[, , k] = -f.conp[, -1,1] %*% data[k,]
    #f.rhsp[, , k][1] = 0
    f.rhs[(1+d*(k-1)):(d*k)] =  f.rhsp[, , k]
  }
  for (k in 1: nq){
    f.conq[1, 1, k]  = 1
    f.conq[, , k][as.matrix(expand.grid(1:2, 1+c(jq, iq)))] = c(1,1,-1,-1)
    f.conq[-c(1:2), (1+jq), k] = rep(-1, (d-2))
    f.conq[, , k][cbind(c(3:d), 1+c(1:d)[-c(iq, jq)])] = 1
    f.con[(np*d+1+d*(k-1)):(d*(k+np)), ] = f.conq[, ,k]
    f.rhsq[, , k] = -f.conq[, -1,1] %*% data[(np+k),]
    #f.rhsq[, , k][1] = 0
    f.rhs[(np*d+1+d*(k-1)):(d*(k+np))] =  f.rhsq[, , k]
  }
  if(flatten){f.rhs[f.rhs < 10 & f.rhs > -10] = 0}
  f.rhs = matrix(f.rhs)
  f.dir=rep("<=",(length(f.rhs)))
  return(lp("max", f.obj, f.con, f.dir, f.rhs))
}

