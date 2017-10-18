library(fda)

all_the_registrations <- function(funzioni, sex=NULL, show_plot = FALSE)
{
    if (is.null(sex))
    {
        sex <- 1:(dim(funzioni$coefs)[2])
    }
    
    rng.age <-  range(funzioni$basis$rangeval)
    Wnbasis <- 10
    Wbasis <- create.bspline.basis(rng.age, norder=4, nbasis=Wnbasis)
    Wfd0 <- fd(matrix(0, Wnbasis, 1), Wbasis)
    WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=100)
    
    ## 1) registration of the fuctions
    registration <- register.fd(mean.fd(funzioni), funzioni, WfdParobj, crit=1,  dbglev= 0)
  #  quartz()
    # pdf('Registration_BMI_FPCA.pdf', height=7, width = 20)
    if (show_plot)
    {
        par(mfrow=c(1,3))
        plot(funzioni, col=sex, main='Original Functions')
        plot(registration$regfd, col=sex, main='Registered Functions')
        plot(registration$warpfd, col=sex, main='Warping Functions')
   
    }
     # dev.off()
    
    #plotreg.fd(registration)
    
    ## 2) registration of the derivatives
    derivatives <- deriv(funzioni, Lfdobj=1)
    #plot(derivatives)
    
    registration_der <- register.fd(mean.fd(derivatives), derivatives, WfdParobj, crit=1,  dbglev= 0)
    
    #quartz()
    if(show_plot)
    {
        par(mfrow=c(1,3))
        plot(derivatives, col=sex, main='Original Derivatives')
        plot(registration_der$regfd,col=sex, main='Registered Derivatives')
        plot(registration_der$warpfd, col=sex, main='Warping Functions')
    
    }
    
    
    data_reg_warp_der <- register.newfd(funzioni, registration_der$warpfd, type = 'monotone')
    
  #  quartz()
    
    #pdf('Registration_derivatives_BMI_FPCA.pdf', height=14, width = 20)
    if(show_plot)
    {

        par(mfrow=c(2,3))
        plot(derivatives, col=sex, main='Original Derivatives')
        plot(registration_der$regfd,col=sex, main='Registered Derivatives')
        plot(registration_der$warpfd, col=sex, main='Warping Functions')
    
        plot(funzioni, col=sex, main='Original Functions')
        plot(data_reg_warp_der,col=sex, main='Registered Functions')
    }
    #dev.off()
}


all_the_registrations_derivatives <- function(funzioni, show_plot = FALSE)
{
    rng.age <-  range(funzioni$basis$rangeval)
    Wnbasis <- 10
    Wbasis <- create.bspline.basis(rng.age, norder=4, nbasis=Wnbasis)
    Wfd0 <- fd(matrix(0, Wnbasis, 1), Wbasis)
    WfdParobj <- fdPar(Wfd0, Lfdobj=2, lambda=100)
    
    
    ## 2) registration of the derivatives
    derivatives <- deriv(funzioni, Lfdobj=1)
    #plot(derivatives)
    
    registration_der <- register.fd(mean.fd(derivatives), derivatives, WfdParobj, crit=1, dbglev= 0)
    
    
    data_reg_warp_der <- register.newfd(funzioni, registration_der$warpfd)
    
    col = rainbow(dim(funzioni$coefs)[2])
    
  #  pdf('Registration_derivatives_BMI_FPCA_new.pdf', height=14, width = 20)
  if (show_plot)
  {
    par(mfrow=c(2,3))
    plot(derivatives, main='Original Derivatives', col= col)
    plot(registration_der$regfd, main='Registered Derivatives', col= col)
    plot(registration_der$warpfd, main='Warping Functions', col= col)
    
    plot(funzioni, main='Original Functions', col= col)
    plot(data_reg_warp_der, main='Registered Functions', col = col)
    
    }
 #   dev.off()
    #
    return(list(registered_fd = data_reg_warp_der, registered_der = registration_der$regfd, 
                warping = registration_der$warpfd))
}