# load the necessary packages
# 
library(fda)
library (refund)

# FUNCTIONS:
# 
# funtion to indentify the correct value of lambda, using the GCV on the data
# or on the derivatives

choice_lambda_spline <- function(data, x_data,
                                 order = 4, lambda=(10^(seq(-5,5,by=0.5))),  n.breaks=100,
                                 GCV.derivatives = TRUE , plot.GCV = TRUE)
{
    NT=dim(data)[2]
    
    num_points_projection <- n.breaks
    
    
    GCV_d0 <- rep(NA, length(lambda))
    GCV_d1 <- rep(NA, length(lambda))
    
    m=order   # order of the  spline
    
    degree=m-1  # degree
    
    spline <- matrix(NA, nrow=dim(data)[1], ncol=NT)
    spline_der <- matrix(NA, nrow=dim(data)[1], ncol=NT)
    matrix_x <- x_data
    for (j in 1:length(lambda))
    {
        print(paste('iteration on lambda: ', j))
        breaks=seq(min(x_data, na.rm=TRUE),  max(x_data, na.rm=TRUE), 
                   length=num_points_projection) 
        
        basis=create.bspline.basis(breaks,norder=m)
        
        functionalPar = fdPar(fdobj=basis,Lfdobj=2,lambda=lambda[j])
        mod=smooth.basis(argvals=x_data, t(as.matrix(data)),functionalPar)
        
        GCV=mod$gcv
        GCV_d0[j] <- mean(GCV)
        spline = t(eval.fd(x_data, mod$fd, Lfdobj=0))
        spline_der = t(eval.fd(x_data, mod$fd, Lfdobj=1))
        SSE_d1 <- rep(0, dim(data)[1])
        GCV_d1_vett <- rep(0, dim(data)[1])
        for (t in 1:dim(data)[1])
        {
            # SSE is the L2 distance between the difference
            SSE_d1[t] <- sum((spline_der[t,(1:(NT-1))]-diff(as.vector(data[t,])))^2)*1 - 
                ((spline_der[t,(1:(NT-1))]-diff(as.vector(data[t,])))^2)[1]/2 - 
                ((spline_der[t,(1:(NT-1))]-diff(as.vector(data[t,])))^2)[NT-1]/2 #1 in the distance between 
            # two points (delta t). and then the Trapezi rule is applied
            GCV_d1_vett[t] <- NT*SSE_d1[t]/(NT-mod$df)^2
        }
        GCV_d1[j] <- mean(GCV_d1_vett)
        
    }
    
    
    if (plot.GCV)
    {
        par(mfrow=c(1,2))
        # plot of the GCV on the data
        plot(log10(lambda), GCV_d0, pch=19, type='b')
        
        # plot of the GCV on the data
        plot(log10(lambda), GCV_d1, pch=19, type='b')
    }
    if (GCV.derivatives){
        
        lambda_selected <- lambda[which.min(GCV_d1)]
        
    }else{
        lambda_selected <- lambda[which.min(GCV_d0)]
    }
    
    return(lambda_selected)
}


# function to define th spline, with a fixed lambda (lambda)
# it will retur the splines (and the derivatives of these splines)
# evaulated in the points of the grid

definition_spline <- function(data, x_data, x_out,
                              order = 4, lambda=0,  n.breaks=100)
{
    
    NT=length(x_out)
    
    num_points_projection <- n.breaks
    
    m=order   # order of the  spline
    
    degree=m-1  # degree
    
    spline <- matrix(NA, nrow=dim(data)[1], ncol=NT)
    spline_der <- matrix(NA, nrow=dim(data)[1], ncol=NT)
    matrix_x <- x_out
    
    breaks=seq(min(x_data),max(x_data),length=num_points_projection) 
    
    basis=create.bspline.basis(breaks,norder=m)
    
    functionalPar = fdPar(fdobj=basis, Lfdobj=2, lambda=lambda)
    mod=smooth.basis(argvals=x_data, t(as.matrix(data)), functionalPar)
    spline = t(eval.fd(x_out, mod$fd, Lfdobj=0))
    spline_der = t(eval.fd(x_out, mod$fd, Lfdobj=1))
    
    return(list(spline= spline, spline_der=spline_der, fd_obj=mod$fd, y2c = mod$y2cMap))
}


# define a smoothing for a dafa_frame 
# (different x grid values for different observations)



define_smoothing <- function(data_frame, col, lambda, breaks)
{
    
    norder=4
    
    basis <- create.bspline.basis(c(min(breaks), max(breaks)), norder= norder, breaks=breaks)
    
    start_element <- matrix(0, basis$nbasis, 1)
    Wfd0 <- fd(start_element, basis)
    
    growpar <- fdPar(Wfd0,2,lambda) # penalization on the 2rd derivative
    
    age_here <- data_frame$age[which(!is.na(data_frame$age))]
    y_here <- data_frame[which(!is.na(data_frame$age)), col]
    
    
    mod_w <- smooth.basis(age_here, y_here , growpar)
    
    result_w <- mod_w$fd
    
    spline = t(eval.fd(age_here, result_w, Lfd=0))
    
    SSE_d0 <- sum((spline-t(y_here))^2)
    
    gcv_2 <- mod_w$gcv # exactly the same as computed before!
    
    return(list(fd_smooth=result_w, gcv=gcv_2))
    
}