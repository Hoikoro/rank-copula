# maximum rank likelihood estimator (MRLE)
library(hgm)


#generating covariance matrix
cov.fun1 = function(th, d){ # AR(1) model (-1 < th < 1)
  S = matrix(0, d, d)
  for(i in 1:d){
    for(j in i:d){
      S[i,j] = S[j,i] = th^(j-i)
    }
  }
  S
}

emp.copula = function(R1, R2){ # empirical copula
  n = length(R1)
  if(n <= 1 || length(R2) != n) stop("aho")
  ec = matrix(0, n-1, n-1)
  s = 1:(n-1)
  for(t in 1:n){
    ec = ec + (1/n) * outer(s >= R1[t], s >= R2[t])
  }
  ec
}


# hgm.ncorthant(x, y) # x: covariance, y: mean, output: int_{t>0} phi_d(t|x,y) dt
loglike = function(R, th, cov.fun){
  # R: n*d matrix, R[i,j] is the rank of the i-th data in j-th variable
  n = nrow(R)
  d = ncol(R)
  S = solve(cov.fun(th, d))
  ll = (n-1)/2*log(det(S)) - d/2*log(n)

  BB = array(0, c(n-1, d, n-1, d))
  q = (1:(n-1))/n
  for(i in 1:d){
    for(j in 1:d){
      ec = emp.copula(R[,i], R[,j])
      BB[,i,,j] = S[i,j] * n * (ec - outer(q, q))
    }
  }

  BBmat = matrix(0, (n-1)*d, (n-1)*d)
  for(i in 1:d){
    for(j in 1:d){
      sub1 = (n-1)*(i-1) + (1:(n-1))
      sub2 = (n-1)*(j-1) + (1:(n-1))
      BBmat[sub1, sub2] = BB[,i,,j]
    }
  }
  x = solve(BBmat)
  x = (x + t(x)) / 2
  y = rep(0, (n-1)*d)
  nc = hgm.ncorthant(x, y) * det(x)^(1/2)

  ll + log(nc)
}

#d = 4      # dimension of copula
ranks = cbind(c(1,2,3,4),c(1,2,4,3),c(1,3,2,4),c(1,3,4,2),c(1,4,2,3),c(1,4,3,2),
               c(2,1,3,4),c(2,1,4,3),c(2,3,1,4),c(2,3,4,1),c(2,4,1,3),c(2,4,3,1),
               c(3,1,2,4),c(3,1,4,2),c(3,2,1,4),c(3,2,4,1),c(3,4,1,2),c(3,4,2,1),
               c(4,1,2,3),c(4,1,3,2),c(4,2,1,3),c(4,2,3,1),c(4,3,1,2),c(4,3,2,1)) #candidate for gamma(,i)

#ranks = cbind(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))

cov = cov.fun1                    # generating covariance matrix
d = 2                             # dimension of copula
n = dim(ranks)[1]                 # number of observation
rankcnt = dim(ranks)[2]           # equal to factorial(n)
grid = seq(-0.99,0.99,by=0.01)    # variable1

m=rankcnt                         # number of different rank functions

#by symmetry, first data can be fixed.

log_likelihood=rep(0,length(dataset))
for(i in 1:length(log_likelihood)){
    mat = cbind(ranks[,1],ranks[,i])
    print(mat)
    for(j in 1:length(grid)){
        param = grid[j]
        log_likelihood[j]=loglike(mat,param,cov)
    }
    pdf(paste("rank-loglike-", i, ".pdf", sep=""))
    plot(grid, log_likelihood, type="l", xlab=expression(theta), ylab="log-likelihood",main = paste(mat[,2]))
    #dev.copy2eps(file=paste("rank-loglike-", i, ".eps", sep=""))
    pos = which.max(log_likelihood)
    print(grid[pos])
    dev.off()
}