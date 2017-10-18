## Inputs ##
# Ymat: coefficient matrix of basis functions for functional Y's (N x K matrix)
# Xmat: matrix of scalar predictors (N x R matrix)
# penmat: penalty matrix (K x K matrix)
#		if no penalty, use 
#			penmat = matrix(0,nrow=dim(Ymat)[2],ncol=dim(Ymat)[2])
# penlamb: penalty lambda for each beta coefficients (a vector of length R)
#		if no penalty, use 
#			penlamb = rep(0,dim(Xmat)[2])

## Outputs ##
# Bhat: coefficient matrix for functional beta that takes same basis as Y's
#		(R x K matrix)
# pvals: test statistics size and p-values based on L2-test, PC-test, 
#		and Choi-test


FS_penreg <- function(Ymat, Xmat, penmat, penlamb) {

  N = dim(Ymat)[1]
  K = dim(Ymat)[2]
  R = dim(Xmat)[2]

  if (dim(Xmat)[1] != N) {
    stop("Check the dimension of Ymat and Xmat")
  }
  if (dim(penmat)[1] != K | dim(penmat)[2] != K) {
    stop("Check the dimension of penmat")
  }
  if (length(penlamb) != R) {
    stop("Check the length of penlamb")
  }

  Lamb <- diag(penlamb)

  # Penalized Estimation
  V = c(t(Ymat)%*%Xmat)
  A = ((t(Xmat)%*%Xmat) %x% diag(1,K)) + Lamb %x% penmat
  Ainv = chol2inv(chol(A))
  Bhat_v = Ainv%*%V
  Bhat_m_T = matrix(Bhat_v,nrow=K,ncol=R)

  Bhat = t(Bhat_m_T)	### Penalized LSE

  Eps_mat = Ymat - Xmat%*%Bhat
  Sig_mat = (1/(N - R)) * t(Eps_mat)%*%Eps_mat
  Big_Cov_mat = Ainv%*%((t(Xmat)%*%Xmat)%x%Sig_mat)%*%Ainv


  # Testing
  EV_PC<-.85
  EV_Choi<-.99
  EV_L2<-.999

  pvalmat <- data.frame(betanum=1:nrow(Bhat), Ti_L2 = NA, pval_L2 = NA, 
								Ti_PC = NA, pval_PC = NA,
								Ti_Choi = NA, pval_Choi = NA)
  require('CompQuadForm')
  for (i in 1:nrow(Bhat)) {
		
    covBi = Big_Cov_mat[(K*(i-1)+1):(K*i),(K*(i-1)+1):(K*i)]
    Eig = eigen(covBi, symmetric=T)
    eigvals_i <- Eig$values
    eigvecs_i <- Eig$vectors

    Scores<-Bhat%*%eigvecs_i

    k_PC<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_PC)
    k_Choi<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_Choi)
    k_L2<-1+sum(cumsum(eigvals_i)/sum(eigvals_i)<EV_L2)
  
    Ti_L2 = sum( Scores[i,1:k_L2]^2 )
    Ti_PC = sum( Scores[i,1:k_PC]^2/eigvals_i[1:k_PC] )
    Ti_Choi = sum( Scores[i,1:k_Choi]^2/sqrt(eigvals_i[1:k_Choi]) )
    pvalmat$Ti_L2[i] = Ti_L2
    pvalmat$Ti_PC[i] = Ti_PC
    pvalmat$Ti_Choi[i] = Ti_Choi
    pvalmat$pval_L2[i] = imhof(Ti_L2/eigvals_i[1], eigvals_i[1:k_L2]/eigvals_i[1], epsabs = 10^(-32), epsrel =10^(-32))$Qq
    pvalmat$pval_PC[i] = pchisq(Ti_PC, df = k_PC, lower.tail=FALSE)
    pvalmat$pval_Choi[i] = imhof(Ti_Choi/sqrt(eigvals_i[1]), sqrt(eigvals_i[1:k_Choi])/sqrt(eigvals_i[1]), epsabs = 10^(-32), epsrel =10^(-32))$Qq

  }

  return(list(Bhat = Bhat, pvals = pvalmat))
}
