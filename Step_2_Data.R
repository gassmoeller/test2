## Cross-convolution shear wave splitting code
## Implements new statistical assessment of splitting measurements via random-phase microseismic noise and bootstrapping
## Last Edited: SEPT 2017 - Derek Witt

################ First, source latest version of Step_1_Setup.R ####################
library(foreach) # foreach package installation is needed.
library(doParallel) # doParallel package installation is needed. It is a "parallel backend" for the foreach package using the %dopar% function. By default, doParallel uses multicore functionality on Unix-like systems and snow functionality on Windows.
num_cores<-4
registerDoParallel(cores=num_cores)# setup multi-core cluster

starttime<-Sys.time()

## 1ST.- Calculate the estimated values of {DT,phi} for a many events with known BAZ and sum likelihoods
########################################################################################################################
########################################################################################################################

## Set up station/event pairs to iterate through
pathname<-"C:/Users/ksabunis/Documents/R/"
setwd(pathname)
This is an improvement
stations<-read.table(paste(pathname,"stations.txt",sep=""))
S<-dim(stations)[1]
for (s in 1:S){

## Loop over stations
st<-stations[s,]

## The list(s) of events need to be made prior to running code
events<-read.table(paste(pathname, st, "_events.txt",sep=""))
E<-dim(events)[1]

setwd(paste(pathname,st,sep=""))

## Set up vectors of event-specific arrays
like_presum<-c()
bazvec<-c()
Ri_p<-c()
Ti_p<-c()
ln_f_Ri<-c()
ln_f_Ti<-c()
bazphi<-c()
bazdt<-c()

for (e in 1:E){
  
  ## Input event name, use consistent nomenclature for event naming  
  event_name<-c(paste(events[e,]))
  Rsacfile<-rsac(paste(events[e,],".BHR.sac",sep=""), endian = .Platform$endian)
  Tsacfile<-rsac(paste(events[e,],".BHT.sac",sep=""), endian = .Platform$endian)
  
  ## Get the back azimuth from the sac file
  bazvec<-cbind(bazvec,Rsacfile[[1]]$baz)
  
  ## Get the sampling interval from the sac file, calculate sqrt_G which is used to set up integer lag
  si<- round(as.numeric(Rsacfile[[1]]$dt),3)
  sqrt_G<-4/si
  
  ## Create lookup table of fast axis and split time parameters
  lu<-Lookup_Table(sqrt_G=sqrt_G, sample_interval=si,back_azimuth=bazvec[,e])
  
  ## Setup R and T sac files
  Rsacfile<-Rsacfile[[1]]$amp
  Tsacfile<-Tsacfile[[1]]$amp
  
  ## Cut off some initial and some final values, 
  ## then split into pre-SKS event (noise) and SKS event (signal plus noise)
  top<-length(Rsacfile)-1
  cut<-round(top/8) # may have to adjust "cut" intervals based on sample interval (8 is ideal for 50 Hz data, 10-16 for 20 Hz data)
  R<-Rsacfile[cut:(top-cut)] # set up the length manually to get the pre-SKS noise and the SKS signal 
  T<-Tsacfile[cut:(top-cut)] # set up the length manually to get the pre-SKS noise and the SKS signal 
  mid<-round((top-cut-cut)/2)
  
  ########################################################################################################################
  
  Ri<-R[1:mid] # Radial pre-SKS noise
  Ti<-T[1:mid] # Tangential pre-SKS noise
  
  outi<-preprocess(Ri,Ti,difference=TRUE)
  
  Ri_p<-cbind(Ri_p,outi$R)
  Ti_p<-cbind(Ti_p,outi$T)
  
  ni<-length(Ri_p[,e])
  omegai<-Fourier_Freq(ni)
  
  I_Ri<-Pgram(Ri_p[,e])
  fit_Ri<-fit_pspline(I_Ri,omegai,penalty=0.5)
  #omi<-omegai[omegai>0]
  #ln_f_Ri<-eval_pspline(fit_Ri,omi,penalty=0.5)
  
  I_Ti<-Pgram(Ti_p[,e])
  fit_Ti<-fit_pspline(I_Ti,omegai,penalty=0.5)
  #ln_f_Ti<-eval_pspline(fit_Ti,omi,penalty=0.5)
  
  ########################################################################################################################
  
  Rstari<-R[(mid+1):(top-cut-cut)] # Radial SKS signal
  Tstari<-T[(mid+1):(top-cut-cut)]  # Tangential SKS signal
  
  out_stari<-preprocess(Rstari,Tstari,sigma=outi$sigma,difference=TRUE)
  Rstari_p<-out_stari$R
  Tstari_p<-out_stari$T
  MLagi<-max(lu$Lag)
  Ni<-length(Rstari_p)-MLagi
  Omegai<-Fourier_Freq(Ni)
  
  ln_f_Ri<-cbind(ln_f_Ri,eval_pspline(fit_Ri,Omegai,penalty=0.5))
  ln_f_Ti<-cbind(ln_f_Ti,eval_pspline(fit_Ti,Omegai,penalty=0.5))
  
  ########################################################################################################################
  
  D_Rstari<-DFT_Matrix(Rstari_p,Max_Lag=MLagi)
  D_Tstari<-DFT_Matrix(Tstari_p,Max_Lag=MLagi)
  
  ########################################################################################################################
  
  ## Vectorized Whittle set-up
  phi<-sort(unique(lu$fast_axis))-bazvec[,e]
  psi<-2*pi*phi/360 # convert to radians
  R1<<-cos(psi)^2
  R2<<-sin(psi)^2
  T1<<-cos(psi)*sin(psi)
  T2<<--T1
  Delta_T<-sort(unique(lu$DT))
  Outer_Real<<-outer(R1,Re(D_Tstari[,1]),"*")-outer(T1,Re(D_Rstari[,1]),"*")
  Outer_Imag<<-outer(R1,Im(D_Tstari[,1]),"*")-outer(T1,Im(D_Rstari[,1]),"*")
  
  
  Gi<-length(unique(lu$DT))
  
  # Whittle Grid Search
  likevi<-c()
  # for(g in 1:Gi){
  #  likevi<-c(likevi,Whittle_fast(Delta_T[g],D_Rstar=D_Rstari,D_Tstar=D_Tstari,ln_f_R=ln_f_Ri[,e],ln_f_T=ln_f_Ti[,e],Omega=Omegai,sample_interval=si))
  # }
  
  likevi_par<-foreach(g=1:Gi) %dopar% {
    out<-Whittle_fast(Delta_T[g],D_Rstar=D_Rstari,D_Tstar=D_Tstari,ln_f_R=ln_f_Ri[,e],ln_f_T=ln_f_Ti[,e],Omega=Omegai,sample_interval=si)
  }
  likevi<-c(likevi,unlist(likevi_par))
  
  ## Smooth likelihood surface for roughness in the DT dimension
  x <- unique(lu$DT)   # Get unique DT values
  L <- matrix(likevi,sqrt_G,sqrt_G)  # Put likei into a matrix
  L_smoother <- apply(L,MAR=1,FUN=smooth.spline,x=x,spar=0.6) # Apply smooth.spline to DT values
  tmp<-c()   # Loop through all fast axis directions to smooth DT
  for (i in 1:sqrt_G){
    tmp <- c(tmp,L_smoother[[i]]$y)
  }
  smooth_likevi<-c(t(matrix(tmp,sqrt_G,sqrt_G))) # Surface comes out flipped so take the transpose...
  
  ## Create matrix of individual likelihoods to sum
  like_presum<-cbind(like_presum,c(smooth_likevi))
  
  ## save baz vs fast axis and split time for each event
  tmp<-Save_results(lu$DT,lu$fast_axis,smooth_likevi)
  bazphi<-rbind(bazphi,cbind(bazvec[,e],tmp$phim))
  bazdt<-rbind(bazdt,cbind(bazvec[,e],tmp$DTm))
} ## End of likelihood calculation 


## Sum likelihood surfaces to get final estimate
final_like<-apply(like_presum,MAR=1,FUN=sum)

## Plot and save results

pdf(paste(pathname,st,'_log_likelihood.pdf',sep=""))
Surface_Plot(lu$DT,lu$fast_axis,final_like)
title(main=paste("Log-Likelihood for ",st))
dev.off()

save<-Save_results(lu$DT,lu$fast_axis,final_like)
est_DT<-c(round(save$DTm,3)) # estimated value of DT with the code based on the cross-convolution method
est_phi<-c(round(save$phim,2)) # estimated value of phi with the code based on the cross-convolution method

list<-cbind(est_DT,est_phi)

## Save splitting parameters, fast axis vs. baz, and split time vs. baz
write.table(list, file=(paste(pathname,st,'_splitting_results.txt',sep="")), sep="\t", row.names=FALSE, col.names=FALSE) 
write.table(bazphi, file=(paste(pathname,st,'_baz_vs_fast_axis.txt',sep="")), sep="\t", row.names=FALSE, col.names=FALSE) 
write.table(bazdt, file=(paste(pathname,st,'_baz_vs_split_time.txt',sep="")), sep="\t", row.names=FALSE, col.names=FALSE) 

########################################################################################################################
########################################################################################################################
########################################################################################################################

## 2ND.- Create synthetic signal, radial, and transverse seismograms with the estimated values of {DT,phi} calculated on the 1st step 
########################################################################################################################
########################################################################################################################

# Set up earth model parameters
split_time<-est_DT
fast_axis<-est_phi
Lag<-floor(split_time/si) # integer lag based on sample interval

# Set up event-specific arrays that need to be filled after initial estimate
top_angle<-c()
Rho1<-c()
Rho2<-c()
Tau1<-c()
Tau2<-c()
Rsig<-c()
Tsig<-c()
Rstar_in<-c()
Tstar_in<-c()
St<-c()

# Loop through event list to populate arrays

for (e in 1:E){
  
  # Find angle between fast axis and baz in radians
  top_angle_degrees<-fast_axis-bazvec[,e] # angle between fast axis and baz in top layer
  top_angle<-cbind(top_angle,2*pi*top_angle_degrees/360)
  
  # Set up the amplitudes of impulses (from Levin and Menke)

  Rho1<-cbind(Rho1,cos(top_angle[,e])^2)
  Rho2<-cbind(Rho2,sin(top_angle[,e])^2)
  Tau1<-cbind(Tau1,cos(top_angle[,e])*sin(top_angle[,e]))
  Tau2<-cbind(Tau2,-Tau1[,e])
  
  # Set up R and T components for least-square regression for source signal estimation 
  Rsig<-cbind(Rsig,R[-1])
  Tsig<-cbind(Tsig,T[-1])
  
  # Estimate the true signal that produced the data to use it as the input wavelet
  true_signal<-estimate_S_star(Rsig[,e],Tsig[,e],L=Lag,rho=c(Rho1[,e],Rho2[,e]),tau=c(Tau1[,e],Tau2[,e]))
  St<-cbind(St,c(rep(0,cut+Lag),true_signal,rep(0,cut-Lag))) # Make the estimated signal the same length as the input seismograms + integer Lag, fill beginning and end with zeros
  
  # Convolve Menke and Levin coefficients with estimated source signal
  RDATA<-rep(0,top)
  TDATA<-RDATA
  
  for(i in (Lag+1):(top+Lag)){
    RDATA[i-Lag]=Rho1[,e]*St[i,e]+Rho2[,e]*St[i-Lag,e]
    TDATA[i-Lag]=Tau1[,e]*St[i,e]+Tau2[,e]*St[i-Lag,e]
  }  
  
  
  # Cut and scale the signal to the same length and same variance as the real signal (Rstari)
  Rboot<-RDATA[cut:(top-cut)]
  Tboot<-TDATA[cut:(top-cut)]
  
  Rstar_in<-cbind(Rstar_in,Rboot[(1+mid):(top-cut-cut)])
  Tstar_in<-cbind(Tstar_in,Tboot[(1+mid):(top-cut-cut)])
  
}


## 3RD.- Bootstrapping the code with the synthetic noise free seismogram created on the 2nd step and the noise file from 1st step
# Set up bootstrapping

B<-100 # Bootstrap iterations

########################################################################################################################
########################################################################################################################
# LOOP STARTS HERE
# Loop over the code to get the values of DT and phi for different repeated simulations of synthetic data foreach iterates over the variables in parallel

bootstart<-Sys.time()
  
  # Set up c vectors
  Rb_p<-c()
  R2b_p<-c()
  Tb_p<-c()
  T2b_p<-c()
  Rb_up<-c()
  R2b_up<-c()
  Tb_up<-c()
  T2b_up<-c()
  
  results<-c()
  indresults<-c()
  likev_presum<-c()
  final_likev<-c()
  
  # Results<-foreach(i=1:ncol(Rb_up), j=1:ncol(Tb_up), k=1:ncol(R2b_up), l=1:ncol(T2b_up)) %dopar% {
  
  for(b in 1:B){
  #Results<-foreach(b=1:B) %dopar% {  
    likev_presum<-c()
    final_likev<-c()
    for(e in 1:E){
      
    ## Generate random-phase noise surrogated from pre-SKS arrival noise 
    Rb_p<-c(SurrogateData(Ri_p[,e], method = "Fourier.rand"))
    R2b_p<-c(SurrogateData(Ri_p[,e], method = "Fourier.rand"))
    
    Tb_p<-c(SurrogateData(Ti_p[,e], method = "Fourier.rand"))
    T2b_p<-c(SurrogateData(Ti_p[,e], method = "Fourier.rand"))
    
    out<-unpreprocess(Rb_p,Tb_p,Sigma=outi$sigma,undifference=TRUE)
    Rb_up<-c(out$R)
    Tb_up<-c(out$T)
    
    out<-unpreprocess(R2b_p,T2b_p,Sigma=outi$sigma,undifference=TRUE)
    R2b_up<-c(out$R)
    T2b_up<-c(out$T)
    
    ########################################################################################################################
    
    out<-preprocess(Rb_up,Tb_up,difference=TRUE)
    R_out<-out$R
    T_out<-out$T
    n<-length(R_out)
    omega<-Fourier_Freq(n)
    I_R<-Pgram(R_out)
    fit_R<-fit_pspline(I_R,omega,penalty=0.5)
    om<-omega[omega>0]
    ln_f_R<-eval_pspline(fit_R,om,penalty=0.5)
    
    I_T<-Pgram(T_out)
    fit_T<-fit_pspline(I_T,omega,penalty=0.5)
    ln_f_T<-eval_pspline(fit_T,om,penalty=0.5)
    
    ########################################################################################################################
    
    Rstar<-Rstar_in[,e] + R2b_up
    Tstar<-Tstar_in[,e] + T2b_up
    
    out_star<-preprocess(Rstar,Tstar,sigma=out$sigma,difference=TRUE)
    Rstar_out<-out_star$R
    Tstar_out<-out_star$T
    MLag<-max(lu$Lag)
    N<-length(Rstar_out)-MLag
    Omega<-Fourier_Freq(N)
    
    ln_f_R<-eval_pspline(fit_R,Omega,penalty=0.5)
    ln_f_T<-eval_pspline(fit_T,Omega,penalty=0.5)
    
    ########################################################################################################################
    
    D_Rstar<-DFT_Matrix(Rstar_out,Max_Lag=MLag)
    D_Tstar<-DFT_Matrix(Tstar_out,Max_Lag=MLag)
    
    ########################################################################################################################
    ## Vectorized Whittle approximation set-up
    phi<-sort(unique(lu$fast_axis))-bazvec[,e]
    psi<-2*pi*phi/360 # convert to radians
    R1<<-cos(psi)^2
    R2<<-sin(psi)^2
    T1<<-cos(psi)*sin(psi)
    T2<<--T1
    Delta_T<-sort(unique(lu$DT))
    Outer_Real<<-outer(R1,Re(D_Tstar[,1]),"*")-outer(T1,Re(D_Rstar[,1]),"*")
    Outer_Imag<<-outer(R1,Im(D_Tstar[,1]),"*")-outer(T1,Im(D_Rstar[,1]),"*")
    
    ########################################################################################################################
    
    Gi<-length(unique(lu$DT))
    ## like<-rep(0,G)
    #SS<-(1:N)[Omega>0]
    likev<-c()
    # for(g in 1:Gi){
    #  likev<-c(likev,Whittle_fast(Delta_T[g],D_Rstar=D_Rstar,D_Tstar=D_Tstar,ln_f_R=ln_f_Ri[,e],ln_f_T=ln_f_Ti[,e],Omega=Omegai,sample_interval=si))
    # }
    likev_par<-foreach(g=1:Gi) %dopar% {
      out<-Whittle_fast(Delta_T[g],D_Rstar=D_Rstari,D_Tstar=D_Tstari,ln_f_R=ln_f_Ri[,e],ln_f_T=ln_f_Ti[,e],Omega=Omegai,sample_interval=si)
    }
    #
    likev<-c(likev,unlist(likev_par))
    
    ########################################################################################################################
    ## Smooth 3D likelihood surface in DT dimension for each bootstrap replicate
    
    x <- unique(lu$DT)   # Get unique DT values
    L_b <- matrix(likev,sqrt_G,sqrt_G)  # Put likei into a matrix
    L_smoother_b <- apply(L_b,MAR=1,FUN=smooth.spline,x=x,spar=0.6) # Apply smooth.spline to DT values
    
    tmp<-c()   # Loop through all fast axis directions to smooth DT
    for (i in 1:sqrt_G){
      tmp <- c(tmp,L_smoother_b[[i]]$y)
    }
    
    smooth_likev<-c(t(matrix(tmp,sqrt_G,sqrt_G))) # Surface comes out flipped so take the transpose...
    
    # Save every station-event bootstrapping result
    saveind<-Save_results(lu$DT,lu$fast_axis,smooth_likev)
    indresults<-rbind(indresults,c(saveind$DTm,saveind$phim))
    
    # Save pre-summed bootstrapping likelihoods
    likev_presum<-cbind(likev_presum,c(smooth_likev))
    
    } # end of event loop
    
    final_likev<-cbind(final_likev,apply(likev_presum,MAR=1,FUN=sum))
    
    # save<-Save_results(lu$DT,lu$fast_axis,final_likev[,b])
    saveboot<-Save_results(lu$DT,lu$fast_axis,final_likev)
    results<-rbind(results,c(saveboot$DTm,saveboot$phim))
    
    #Results<-matrix(data=Results,nrow=B,ncol=2)
    #Results[b,]<-c(results)
    
    ########################################################################################################################
    ########################################################################################################################
    ########################################################################################################################
    
  } # end of bootstrap loop 
  
 #} #dopar
  # Calculate stats for summed and unsummed likelihood realizations
  bootdtsd<-sd(results[,1])
  bootphisd<-sd(results[,2])
  indbootdtsd<-sd(indresults[,1])
  indbootphisd<-sd(indresults[,2])
  # bootdtmean<-mean(results[,1])
  # bootphimean<-mean(results[,2])
  # indbootdtmean<-mean(indresults[,1])
  # indbootphimean<-mean(indresults[,2])
  bootstats<-cbind(bootdtsd,bootphisd)
  indbootstats<-cbind(indbootdtsd,indbootphisd)
  stats<-cbind(est_DT,est_phi,indbootstats)
  
  #Results<-t(round(matrix(unlist(Results),2,B),2))
  write.table(results, file=(paste(pathname,st,'_bootstrapping_results.txt',sep="")), sep="\t", row.names=FALSE, col.names=FALSE) 
  write.table(indresults, file=(paste(pathname,st,'_stationevent_bootstrapping_results.txt',sep="")), sep="\t", row.names=FALSE, col.names=FALSE)
  write.table(stats, file=(paste(pathname,st,'_results_stats.txt',sep="")), sep="\t", row.names=FALSE, col.names=FALSE)
  boottime<-Sys.time()-bootstart
} #end of station loop (if being used)
totaltime<-Sys.time()-starttime
 
