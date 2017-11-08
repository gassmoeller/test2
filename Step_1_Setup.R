## Script to set up functions for cross-convolution shear wave splitting code
## Must run before sourcing Step_2_Data.R
## Last Edited: JAN 2017 - Derek Witt

######################################################################
######################################################################
######################################################################
# Set up manually the function rsac to open and read sac files in R
rsac<- function(files, endian = .Platform$endian)
{
  if(length(endian) == 1 & length(files) > 1)
    endian <- rep(endian, length(files))
  n <- length(files)
  data <- vector(mode = "list", length = n)
  for(i in 1:n)
  {
    file <- files[i]
    zz <- file(file, "rb")
    h1 <- readBin(con = zz, what = numeric(), n = 70, size = 4,
                  endian = endian[i])
    dim(h1) <- c(5, 14)
    h1 <- aperm(h1)
    # NA values:
    h1[h1 == -12345] <- NA
    h2 <- readBin(con = zz, what = integer(), n = 35, size = 4,
                  endian = endian[i])
    dim(h2) <- c(5, 7)
    h2 <- aperm(h2)
    # NA values:
    h2[h2 == -12345] <- NA
    h3 <- readBin(con = zz, what = logical(), n = 5, size = 4,
                  endian = endian[i])
    h4 <- readBin(con = zz, what = character(), n = 1, size = 4,
                  endian = endian[i])
    # Define header variables:
    dt <- h1[1, 1]
    depmin <- h1[1, 2]
    depmax <- h1[1, 3]
    scale <- h1[1, 4]
    odelta <- h1[1, 5]
    b <- h1[2, 1]
    e <- h1[2, 2]
    o <- h1[2, 3]
    a <- h1[2, 4]
    f <- h1[5,
    stla <- h1[7, 2]
    stlo <- h1[7, 3]
    stel <- h1[7, 4]
    stdp <- h1[7, 5]
    evla <- h1[8, 1]
    evlo <- h1[8, 2]
    evel <- h1[8, 3]
    evdp <- h1[8, 4]
    mag <- h1[8, 5]
    dist <- h1[11, 1]
    az <- h1[11, 2]
    baz <- h1[11, 3]
    gcarc <- h1[11, 4]
    cmpaz <- h1[12, 3]
    cmpinc <- h1[12, 4]
    nzyear <- h2[1, 1]
    nzjday <- h2[1, 2]
    nzhour <- h2[1, 3]
    nzmin <- h2[1, 4]
    nzsec <- h2[1, 5]
    nzmsec <- h2[2, 1]
    norid <- h2[2, 3]
    nevid <- h2[2, 4]
    N <- h2[2, 5]
    idep <- h2[4, 2]
    iztype <- h2[4, 3]
    leven <- h3[1]
    lpspol <- h3[2]
    kstnm <- substr(h4, 1, 8)
    kstnm <- sub("-12345", "      ", kstnm)
    kevnm <- substr(h4, 9, 24)
    kevnm <- sub("-12345", "      ", kevnm)
    khole <- substr(h4, 25, 32)
    khole <- sub("-12345", "      ", khole)
    ko <- substr(h4, 33, 40)
    ko <- sub("-12345", "      ", ko)
    ka <- substr(h4, 41, 48)
    ka <- sub("-12345", "      ", ka)
    kcmpnm <- substr(h4, 161, 168)
    kcmpnm <- sub("-12345", "      ", kcmpnm)
    knetwork <- substr(h4, 169, 176)
    knetwork <- sub("-12345", "      ", knetwork)
    kinst <- substr(h4, 185, 192)
    kinst <- sub("-12345", "      ", kinst)
    seek(con = zz, where = 632)
    x <- readBin(con = zz, what = numeric(), n = N,
                 size = 4, endian = endian[i])
    close(zz)
    data[[i]] <- list(amp = x, dt = dt, depmin = depmin, depmax = depmax,
                      scale = scale, odelta = odelta,
                      b = b, e = e, o = o, a = a, f = f,
                      stla = stla, stlo = stlo, stel = stel, stdp = stdp,
                      evla = evla, evlo = evlo, evel = evel, evdp = evdp,
                      mag = mag, dist = dist, az = az, baz = baz, gcarc = gcarc,
                      cmpaz = cmpaz, cmpinc = cmpinc,
                      nzyear = nzyear, nzjday = nzjday, nzhour = nzhour,
                      nzmin = nzmin, nzsec = nzsec,
                      nzmsec = nzmsec, norid = norid,
                      nevid = nevid, N = N,
                      units = idep, iztype = iztype,
                      leven = leven, lpspol = lpspol,
                      sta = kstnm, kevnm = kevnm, khole = khole,
                      ko = ko, ka = ka,
                      comp = kcmpnm, knetwork = knetwork, kinst = kinst)
  }
  class(data) <- "rsac"
  invisible(data)
}
## Code: Thompson, E. M. and Lees, J. M. \code{Rsac}
######################################################################
######################################################################
######################################################################
preprocess<-function(R,T,sigma=-9,difference=TRUE){
	# Center the data by subtracting off the mean
	if(difference){
		R<-diff(R,1,1)
		T<-diff(T,1,1)
	}# end if difference
	R<- R-mean(R)
	T<- T-mean(T)
	if(sigma<0){
		sigma_R<-sqrt(var(R))
		sigma_T<-sqrt(var(T))
		sigma<-(sigma_R+sigma_T)/2
	}# end if negative sigma
	R<-R/sigma 
	T<-T/sigma 
	out<-list(R=R,T=T,sigma=sigma)
	return(out)
}
######################################################################
######################################################################
######################################################################
unpreprocess<-function(R,T,Sigma=1,undifference=TRUE){
    if(undifference){
    	R <- cumsum(c(0,R))
    	T <- cumsum(c(0,T))
    }
    #R<- R+mean(R)
    #T<- T+mean(T)
    R <- Sigma*R
    T <- Sigma*T
    out<-list(R=R,T=T)
    return(out)
}
######################################################################
######################################################################
######################################################################
# Function to compute the Fourier frequencies.
Fourier_Freq<-function(n){
	Ln<- floor((n-1)/2) # this means take the integer part
	Un<-floor(n/2)
	Fn<- -Ln:Un
	omega<-2*pi*Fn/n
	return(omega)
}
######################################################################
######################################################################
######################################################################
# Periodogram/(2*pi) function.
Pgram<-function(y){
	n<-length(y)
	Ln<- floor((n-1)/2) # this means take the integer part
	Un<-floor(n/2)
	Fn<- -Ln:Un
	# dft in R does not divide by n; need to do that here.
	#
	dft<-fft(y)/sqrt(n*2*pi)
	#
	# fft orders these (0,2pi) instead of (-pi,pi), 
	# so (pi,2pi) needs to map to (-pi,0). 
	#
	tmp<-dft
	dft[1:Ln]<-tmp[(Un+2):n]
	dft[(Ln+1):(Ln+Un+1)]<-tmp[1:(Un+1)]
	#
	I_y<-Re(dft*Conj(dft))
	return(I_y)
}
######################################################################
######################################################################
######################################################################
fit_pspline<-function(I,omega,K=90,eps=0.005,penalty=1){
	# Need to delete frequency zero.
	zero<-(omega==0)
	I<-I[!zero]
	omega<-omega[!zero]
	knots<-log(eps)+(1:K)*(log(pi-eps)-log(eps))/K
	knots<-c(-rev(exp(knots)),exp(knots))
	ZZ<-outer(omega,knots,"-")
	ZZ<-ZZ*(ZZ>0)
	ZZ<-cbind(rep(1,length(omega)),omega,ZZ)
	fit<-solve(t(ZZ)%*%ZZ+diag(c(0,0,rep(penalty,2*K))))%*%t(ZZ)%*%cbind(log(I))
	return(fit)
}
######################################################################
######################################################################
######################################################################
eval_pspline<-function(fit,omega,K=90,eps=0.005,penalty=1){
	knots<-log(eps)+(1:K)*(log(pi-eps)-log(eps))/K
	knots<-c(-rev(exp(knots)),exp(knots))
	ZZ<-outer(omega,knots,"-")
	ZZ<-ZZ*(ZZ>0)
	ZZ<-cbind(rep(1,length(omega)),omega,ZZ)
	eval<-ZZ%*%cbind(c(fit))
	return(eval)
}
######################################################################
######################################################################
######################################################################
# Matrix of dft's.
DFT_Matrix<-function(y,Max_Lag){
    D<-c()
    N<-length(y)
      for(g in 0:Max_Lag){
		D<-cbind(fft(y[(g+1):(g+N-Max_Lag)]),D)
      }

    # Reorder D for consistency with periodogram.
    # fft orders these (0,2pi) instead of (-pi,pi), 
    # so (pi,2pi) needs to map to (-pi,0). 

    n<-N-Max_Lag
    Ln<- floor((n-1)/2) # this means take the integer part
    Un<-floor(n/2)
    tmp<-D
    D[1:Ln,]<-tmp[(Un+2):n,]
    D[(Ln+1):(Ln+Un+1),]<-tmp[1:(Un+1),]
    return(D)
}
######################################################################
######################################################################
######################################################################
# Plot a criterion surface
Surface_Plot<-function(DT,phi,criterion){
	G<-length(DT)
	h<-topo.colors(G)
	plot(range(DT),range(phi),type="n",xlab="Split Time",ylab="Fast Axis")
	r<-rank(criterion)
	delta_DT<-max(diff(DT,1,1))
	delta_phi<-max(diff(abs(phi),1,1)) # can be messed up depending on order of grid values
	for(i in 1:G){
		rect(DT[i],phi[i],DT[i]+delta_DT,phi[i]+delta_phi,col=h[r[i]],border=NA)
	}# end loop on i
	g_max<-(1:G)[r==G]
	points(DT[g_max],phi[g_max],pch=17,col="red",cex=1.2)
}
######################################################################
######################################################################
######################################################################
#Save the values of DT and phi in a list
Save_results<-function(DT,phi,criterion){
	G<-length(DT)
	r<-rank(criterion)
	g_max<-(1:G)[r==G]
	save<-list(DTm=DT[g_max],phim=phi[g_max])
	return(save)
}
######################################################################
######################################################################
######################################################################
# Create lookup table using Menke and Levin one-layer functional forms.
# Change manually the maximum and minimum dt & maximun and minimum fast axis:
Lookup_Table<-function(sqrt_G=200,scale=1,max_DT=4.0,min_DT=0,max_phi=180,min_phi=0,sample_interval=0.01,back_azimuth=baz){
	fast_axis<-rep((min_phi+((0:(sqrt_G-1))/sqrt_G)*(max_phi-min_phi)),sqrt_G)
	delta<-fast_axis-back_azimuth
	delta<-2*pi*delta/360 
	DT<-sort(rep((min_DT+((0:(sqrt_G-1))/sqrt_G)*(max_DT-min_DT)),sqrt_G))
	R1<-scale*cos(delta)^2
	R2<-scale*sin(delta)^2
	T1<-scale*cos(delta)*sin(delta)
	T2<--T1
	Lag<-(DT/sample_interval)%/%1 
	out<-list(fast_axis=fast_axis,DT=DT,R1=R1,R2=R2,T1=T1,T2=T2,Lag=Lag)
	return(out)
}
######################################################################
######################################################################
######################################################################
# Set up manually the function SurrogateData from WaveletComp package
SurrogateData <-
function(x, method = "white.noise", 
                          params=list(AR     = list(p=1),
                                      ARIMA  = list(p=1, q=1, include.mean=T, sd.fac=1, trim = F, trim.prop = 0.01)
#                                       ,
#                                       meboot = list(trim = 0.1, force.clt = F, expand.sd = T, fiv = 5)
                                     ) 
                          ){
                          
  if(method == "white.noise")  x.sur <- rnorm(length(x)) 
  if(method == "shuffle")      x.sur <- sample(x, length(x)) 
  if(method == "Fourier.rand") x.sur <- FourierRand(x) 
  
  if(method == "AR")           { 
 
     x.sur <- AR(x, params = params) 
     
  } 
  
#   if(method == "meboot")       { 
#   
#      trim      = params$meboot$trim
#      force.clt = params$meboot$force.clt
#      expand.sd = params$meboot$expand.sd
#      fiv       = params$meboot$fiv
#      
#      x.sur <- meboot(x, reps=2, trim = trim, force.clt = force.clt, expand.sd = expand.sd, fiv = fiv)$ensemble[,1]
#      
#   }
  
  if(method == "ARIMA")         {
  
     x.sur <- ARIMA(x, params = params)
 
  }
  
  return(invisible(x.sur))
}
## Code: Tian, H. and Cazelles, B., \code{WaveletCo}
######################################################################
######################################################################
######################################################################
# Set up manually the function FourierRand from WaveletComp package
#Needed for using SurrogateData with the Fourier randomization method
FourierRand <-
function(x){
  n <- length(x)
  z <- fft(x)

  if(n%%2 == 0){
    ph <- 2*pi*runif(n/2-1)
    ph <- c(0, ph, 0, -rev(ph))}

  if(n%%2 != 0){
    ph <- 2*pi*runif((n-1)/2)
    ph <- c(0, ph, -rev(ph))}

  ph <- complex(imaginary = ph)
  z <- z * exp(ph)
  x.sur <- Re(fft(z, inverse = TRUE)/n)
  
  return(invisible(x.sur))
}
## Code: Tian, H. and Cazelles, B., \code{WaveletCo}
######################################################################
######################################################################
######################################################################
## Requires several globally-defined variables: R1, R2, T1, T2, Outer_Real, Outer_Imag
## Evaluate log-like.
Whittle_fast<-function(DT,D_Rstar,D_Tstar,ln_f_R,ln_f_T,Omega,sample_interval=0.02){
    Max_Lag<-floor(4/sample_interval)
    N<-dim(D_Tstar)[1]
    SS<-(1:N)[Omega>0]
    lag<-DT/sample_interval
    Lag1<-floor(lag)
    ##
    ## Construct matrix with length(psi) rows and length(Omega) columns
    ## Vector is added to each column.
    rho2<-(R1^2+R2^2)+outer(2*R1*R2,cos(Lag1*Omega),"*")
    tau2<-(T1^2+T2^2)+outer(2*T1*T2,cos(Lag1*Omega),"*")
    ##
    ## Construct matrix of spectra for cross-convolved series.
    spectrum_cross1<-t(t(rho2)*c(exp(ln_f_T))+t(tau2)*c(exp(ln_f_R)))
    tmp_real<-Outer_Real+outer(R2,Re(D_Tstar[,Lag1+1]),"*")-outer(T2,Re(D_Rstar[,Lag1+1]),"*")
    tmp_imag<-Outer_Imag+outer(R2,Im(D_Tstar[,Lag1+1]),"*")-outer(T2,Im(D_Rstar[,Lag1+1]),"*")
    I_cross1<-(tmp_real^2+tmp_imag^2)/(2*pi*N)
    Lag2<-Lag1+1
    I_cross2<-0
    spectrum_cross2<-1
    like2<-rep(0,length(psi))
    if(lag>Lag1){
        tmp_real<-Outer_Real+outer(R2,Re(D_Tstar[,Lag2+1]),"*")-outer(T2,Re(D_Rstar[,Lag2+1]),"*")
        tmp_imag<-Outer_Imag+outer(R2,Im(D_Tstar[,Lag2+1]),"*")-outer(T2,Im(D_Rstar[,Lag2+1]),"*")
        I_cross2<-(tmp_real^2+tmp_imag^2)/(2*pi*N)
        rho2<-(R1^2+R2^2)+outer(2*R1*R2,cos(Lag2*Omega),"*")
        tau2<-(T1^2+T2^2)+outer(2*T1*T2,cos(Lag2*Omega),"*")
        ##
        ## Construct matrix of spectra for cross-convolved series.
        spectrum_cross2<-t(t(rho2)*c(exp(ln_f_T))+t(tau2)*c(exp(ln_f_R)))
        like2<-apply(-(I_cross2/spectrum_cross2)-log(spectrum_cross2),MAR=1,FUN="sum")
    }# end if lag
    ##
    like1<-apply(-(I_cross1/spectrum_cross1)-log(spectrum_cross1),MAR=1,FUN="sum")
    log_like<-(Lag2-lag)*like1+(lag-Lag1)*like2
    ##
    ##
    return(log_like)
}

######################################################################
######################################################################
######################################################################
estimate_S_star<-function(Rstar,Tstar,L=120,rho=c(cos(pi/9)^2,sin(pi/9)^2),tau=c(cos(pi/9)*sin(pi/9),-cos(pi/9)*sin(pi/9))){
  N<-length(Rstar) 
  Z_R<-matrix(0,N,N+L)
  Z_T<-matrix(0,N,N+L)

  tmp1<-cbind(matrix(0,N,L),diag(rep(rho[1],N)))
  tmp2<-cbind(diag(rep(rho[2],N)),matrix(0,N,L))
  Z_R<-tmp1+tmp2
  
  tmp1<-cbind(matrix(0,N,L),diag(rep(tau[1],N)))
  tmp2<-cbind(diag(rep(tau[2],N)),matrix(0,N,L))
  Z_T<-tmp1+tmp2
  
  YT<-cbind(Tstar)
  YR<-cbind(Rstar)
  Y<-rbind(YR,YT)
  Z<-rbind(Z_R,Z_T)
  S_star_hat<-lm(Y~-1+Z)$coef # estimate of the pre-processed signal, (1-B)S_t/sigma
  return(S_star_hat)
}
###
