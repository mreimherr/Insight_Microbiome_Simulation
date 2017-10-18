# F on Scalar Regression
FonS_reg<-function(Y.f,X,q=20,lamb_beta=1){
	EV_PC<-.95
	EV_CH<-.99
	EV_L2<-.999
	
	X_mod<-model.matrix(~X)
	p = dim(X_mod)[2]
	n = dim(X_mod)[1]
	Y.pc<-pca.fd(Y.f,nharm=q)
	v<-Y.pc$harmonics
	#Ymat<-Y.pc$scores
	Ymat<-inprod(Y.f,v)
	
	v2 = deriv.fd(v,2)
	U = inprod(v2,v2)
	
	Ik = diag(1,q, q)
	lamb_int = 0
	#lamb_beta = 10000
	Lamb = diag(c(lamb_int,lamb_beta))
	
	V = c(t(Ymat)%*%X_mod)
	A = ((t(X_mod)%*%X_mod) %x% Ik) + Lamb%x%U
	Ainv = solve(A)
	Bhat_v = Ainv%*%V
	Bhat_m_T = matrix(Bhat_v,nrow=q,ncol=p)
	Bhat_c = v$coefs%*%Bhat_m_T
	Bhat_f = fd(coef=Bhat_c,v$basis)
	
	#plot(Bhat_f)
	
	### Testing
	Eps_mat = Ymat - X_mod%*%t(Bhat_m_T)
	Sig_mat = (1/(n - p))*t(Eps_mat)%*%Eps_mat
	
	Big_Cov_mat = Ainv%*%((t(X_mod)%*%X_mod)%x%Sig_mat)%*%Ainv
	CV_x = Big_Cov_mat[(q+1):(q*p),(q+1):(q*p)]
	Bhat1_c = Bhat_m_T[,2]
	
	CV_E<-eigen(CV_x)
	CV_v<-CV_E$vectors
	CV_val<- CV_E$values
	k_PC<-1+sum(cumsum(CV_val)/sum(CV_val)<EV_PC)
	k_CH<-1+sum(cumsum(CV_val)/sum(CV_val)<EV_CH)
	k_L2<-1+sum(cumsum(CV_val)/sum(CV_val)<EV_L2)
	
	
	Scores<-CV_v%*%Bhat1_c
	
	T_PC = sum(Scores[1:k_PC]^2/CV_val[1:k_PC])
	T_CH = sum(Scores[1:k_CH]^2/sqrt(CV_val[1:k_CH]))
	T_L2 = sum(Scores[1:k_L2]^2)
	#pchisq(T,df=q,lower.tail=FALSE)
	pval_PC<- pchisq(T_PC,df=k_PC,lower.tail=FALSE)
	pval_CH<- imhof(T_CH,sqrt(CV_val[1:k_CH]))[[1]]
	pval_L2<- imhof(T_L2,CV_val[1:k_L2])[[1]]
	
	
	Cv_c<- v$coefs%*%CV_x%*%t(v$coefs)
	Cv_fun<-bifd(Cv_c,v$basis,v$basis)
	pts<-seq(v$basis$rangeval[1],v$basis$rangeval[2],length=50)
	Cv_mat = eval.bifd(pts,pts,Cv_fun)
	se_vec<-sqrt(diag(Cv_mat))
	
	se_fun<-Data2fd(pts,se_vec,Y.f$basis)

	pvals<-matrix(c(pval_PC,pval_CH,pval_L2),ncol=1)
	rownames(pvals) = c("PC","Choi","L2")
	return(list(beta_est = Bhat_f, se_fun = se_fun, pvals=pvals))
}