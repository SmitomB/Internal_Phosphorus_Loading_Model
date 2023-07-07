# functions needed for Jordan Lake model
# Original author: Dario Del Giudice
# New author: Smitom Borah

# EDITS
## 1. 08-31-2022 07:00 am: edited the sysanal.plot.chains() to make it publishable in the paper

# calculate density of prior distribution
logpri  <- function(par,priors)
{
  logprior <- 0
  for ( i in 1:length(priors)) # for every parameters
  { 
    if ( priors[[i]][1] == "uniform" ) 
      logpdf <-dunif(par[i],as.numeric(priors[[i]][2]),as.numeric(priors[[i]][3]), log =T)
    if ( priors[[i]][1] == "truncnorm" ) 
    {
      logpdf <-log(dtruncnorm(par[i],as.numeric(priors[[i]][2]),as.numeric(priors[[i]][3]),as.numeric(priors[[i]][4]),as.numeric(priors[[i]][5])))
    }
    logprior <- logprior + logpdf
  }
  return(logprior)
}



# define objective function (only having par as argument, other "inputs" defined in previous functions & hardcoded)
logpost <- function(par,priors,...)
{
  if(any(par<0)){ logdens = -Inf}else{  
    logprior  <- logpri(par,priors)
    loglikeli <- loglik(par,...)
    logdens <- as.numeric(logprior+loglikeli)}
  if (rnorm(1)>2) {print(par);print(logdens)}
  return(logdens)
} 

# log-likelihood: evaluates model fit, returns log density
loglik <- function(par,mod_fun,yObs,lambda1=0,lambda2=0){
  yM <- mod_fun(CalPar=par) %>% pull(value)
  ll_trans <- sum(dnorm(x=sysanal.boxcox(yObs, lambda1,lambda2), 
                   mean =sysanal.boxcox(yM, lambda1,lambda2), 
                    sd = par["sy"], log=T)+
    log(sysanal.boxcox.deriv(yObs,lambda1,lambda2)))
  return(ll_trans)
  }

# BoxCox transformation 

sysanal.boxcox <- function(data,lambda1=1,lambda2=1)
{
  if ( lambda1 == 0 )
  {
    return(ifelse(data>-lambda2,log(data+lambda2),NA))
  }
  else
  {
    return(ifelse(data>=-lambda2,((data+lambda2)^lambda1 - 1)/lambda1,NA))
  }
}

sysanal.boxcox.deriv <- function(data,lambda1=1,lambda2=1)
{
  return(ifelse(data>-lambda2,(data+lambda2)^(lambda1 - 1),NA))
}


sysanal.boxcox.inv <- function(data,lambda1=1,lambda2=1)
{
   if ( lambda1 == 0 )
   {
      return(exp(data)-lambda2)
   }
   else
   {
      return(ifelse(lambda1*data>-1,(lambda1*data+1)^(1/lambda1)-lambda2,
                                    -lambda2))
   }
}

# calculates R2 as fraction of variance explained i.e. as N-S
# --------------------------------------

VarExp <- function(mod,obs)
{
  RSS = t(obs-mod)%*%(obs-mod); SYY = t(obs)%*%(obs)-length(obs)*mean(obs)^2 
  return(paste(signif(1-RSS/SYY,3)))
}

# function to plot markov chains:
# -------------------------------

sysanal.plot.chains <- function(postsamp,ncol=NA,mar=NA,
                                ylim=list(),
                                titles=list(),xlab="chain index",ylab=list())
{
   nvar <- ncol(postsamp)
   if ( is.na(ncol) ) nc <- floor(sqrt(nvar))
   nr <- ceiling(nvar/nc)
   marg <- mar
   if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,2) # c(bottom, left, top, right)
   par.def <- par(no.readonly=T)
   par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
   for ( i in 1:nvar )
   {
      name <- colnames(postsamp)[i]
      data <- postsamp[,i]
      ylim.i <- c(min(data),max(data))
      if ( ylim.i[1] == ylim.i[2] ) 
      {
         ylim.i[1] <- 0.5*ylim.i[1]
         ylim.i[2] <- 1.5*ylim.i[2]
      }
      if ( length(ylim) > 0 )
      {
         ind <- match(name,names(ylim))
         if ( !is.na(ind) )
         {
            ylim.i <- ylim[[ind]]
         }
      }
      title.i <- name
      if ( length(titles) > 0 )
      {
         ind <- match(name,names(titles))
         {
            if ( !is.na(ind) )
            {
               title.i <- titles[[ind]]
            }
         }
      }
      ylab.i <- name
      if ( length(ylab) > 0 )
      {
         ind <- match(name,names(ylab))
         if ( !is.na(ind) )
         {
            ylab.i <- ylab[[ind]]
         }
      }
      plot(data,ylim=ylim.i,type="l",main=title.i,xlab=xlab,ylab=ylab.i)
   }
   par(par.def)
}

#####################################################
# current model version to compile

lakemod_v0 <- odin::odin({
  #           Resuspension & tmp adjust    Wsh in           dwnstr in   upstr in     Settling& tmp adjust                                    outflow
  deriv(M1) <- (S1*Rs)*(ThetaR^(Tmp-20))+(Q1in*C1in)*(1-Sh)+(Q21*M2/V2)            -((M1*(Ks+Vs*A1/V1))*(ThetaV^(Tmp-20))+(Q12*M1/V1))
  deriv(M2) <- (S2*Rs)*(ThetaR^(Tmp-20))+(Q2in*C2in)*(1-Sh)+(Q32*M3/V3)+(Q12*M1/V1)-((M2*(Ks+Vs*A2/V2))*(ThetaV^(Tmp-20))+(Q23*M2/V2)+(Q21*M2/V2))
  deriv(M3) <- (S3*Rs)*(ThetaR^(Tmp-20))+(Q3in*C3in)*(1-Sh)+(Q43*M4/V4)+(Q23*M2/V2)-((M3*(Ks+Vs*A3/V3))*(ThetaV^(Tmp-20))+(Q34*M3/V3)+(Q32*M3/V3))
  deriv(M4) <- (S4*Rs)*(ThetaR^(Tmp-20))+(Q4in*C4in)*(1-Sh)+            (Q34*M3/V3)-((M4*(Ks+Vs*A4/V4))*(ThetaV^(Tmp-20))+(M4*Q4out/V4)+(Q43*M4/V4))
  #          Settle in&temp adjust          Burial   resuspension&temp adjust
  deriv(S1) <- (M1*(Ks+Vs*A1/V1))*(ThetaV^(Tmp-20))-((S1*Bs)+(S1*Rs)*(ThetaR^(Tmp-20)))
  deriv(S2) <- (M2*(Ks+Vs*A2/V2))*(ThetaV^(Tmp-20))-((S2*Bs)+(S2*Rs)*(ThetaR^(Tmp-20)))
  deriv(S3) <- (M3*(Ks+Vs*A3/V3))*(ThetaV^(Tmp-20))-((S3*Bs)+(S3*Rs)*(ThetaR^(Tmp-20)))
  deriv(S4) <- (M4*(Ks+Vs*A4/V4))*(ThetaV^(Tmp-20))-((S4*Bs)+(S4*Rs)*(ThetaR^(Tmp-20)))
  
  initial(M1) <- M1init
  initial(M2) <- M2init
  initial(M3) <- M3init
  initial(M4) <- M4init
  
  initial(S1) <- S1init*S1fac
  initial(S2) <- S2init*S1fac
  initial(S3) <- S3init*S1fac
  initial(S4) <- S4init*S1fac
  
  M1init = user()
  M2init = user()
  M3init = user()
  M4init = user()
  
  S1init = user()
  S2init = user()
  S3init = user()
  S4init = user()
  
  Vs = user()
  Ks = user()
  Bs = user()
  Rs = user()
  Sh = user()
  ThetaR = user()
  ThetaV = user()
  S1fac = user()
  # S2fac = user()
  # S3fac = user()
  # S4fac = user()
  timestep[] = user()
  V1 = interpolate(timestep,V1t,"constant")
  V2 = interpolate(timestep,V2t,"constant")
  V3 = interpolate(timestep,V3t,"constant")
  V4 = interpolate(timestep,V4t,"constant")
  Q1in = interpolate(timestep,Q1int,"constant")
  Q2in = interpolate(timestep,Q2int,"constant")
  Q3in = interpolate(timestep,Q3int,"constant")
  Q4in = interpolate(timestep,Q4int,"constant")
  Q12 = interpolate(timestep,Q12t,"constant")
  Q21 = interpolate(timestep,Q21t,"constant")
  Q23 = interpolate(timestep,Q23t,"constant")
  Q32 = interpolate(timestep,Q32t,"constant")
  Q34 = interpolate(timestep,Q34t,"constant")
  Q43 = interpolate(timestep,Q43t,"constant")
  Q4out = interpolate(timestep,Q4outt,"constant")
  A1 = interpolate(timestep,A1t,"constant")
  A2 = interpolate(timestep,A2t,"constant")
  A3 = interpolate(timestep,A3t,"constant")
  A4 = interpolate(timestep,A4t,"constant")
  C1in= interpolate(timestep,C1int,"constant")
  C2in= interpolate(timestep,C2int,"constant")
  C3in= interpolate(timestep,C3int,"constant")
  C4in= interpolate(timestep,C4int,"constant")
  Tmp = interpolate(timestep,Tmpt,"constant")
  
  V1t[] = user()
  V2t[] = user()
  V3t[] = user()
  V4t[] = user()
  Q1int[] = user()
  Q2int[] = user()
  Q3int[] = user()
  Q4int[] = user()
  Q12t[] = user()
  Q21t[] = user()
  Q23t[] = user()
  Q32t[] = user()
  Q34t[] = user()
  Q43t[] = user()
  Q4outt[] = user()
  A1t[] = user()
  A2t[] = user()
  A3t[] = user()
  A4t[] = user()
  C1int[] = user()
  C2int[] = user()
  C3int[] = user()
  C4int[] = user()
  Tmpt[] = user()
  
  dim(V1t) = user()
  dim(V2t) = user()
  dim(V3t) = user()
  dim(V4t) = user()
  dim(Q1int) = user()
  dim(Q2int) =user()
  dim(Q3int) = user()
  dim(Q4int) = user()
  dim(Q12t) =user()
  dim(Q21t) =user()
  dim(Q23t) =user()
  dim(Q32t) =user()
  dim(Q34t) = user()
  dim(Q43t) =user()
  dim(Q4outt) = user()
  dim(A1t) =user()
  dim(A2t) = user()
  dim(A3t) = user()
  dim(A4t) = user()
  dim(C1int)= user()
  dim(C2int)= user()
  dim(C3int)= user()
  dim(C4int)= user()
  dim(Tmpt) = user()
  dim(timestep)=user()
})    


#################################################################
# Create input list

creaInp <- function(inpdf,tlim=tmax)
{
  return(list(V1t = inpdf$V1[1:tlim],
                V2t = inpdf$V2[1:tlim],
                V3t = inpdf$V3[1:tlim],
                V4t = inpdf$V4[1:tlim],
                V1downt = inpdf$V1down[1:tlim],
                V2downt = inpdf$V2down[1:tlim],
                V3downt = inpdf$V3down[1:tlim],
                V4downt = inpdf$V4down[1:tlim],
                mot = inpdf$mo[1:tlim],
                Q1int = inpdf$Q1in[1:tlim],
                Q2int = inpdf$Q2in[1:tlim],
                Q3int = inpdf$Q3in[1:tlim],
                Q4int = inpdf$Q4in[1:tlim],
                Q12t = inpdf$Q12[1:tlim],
                Q21t = inpdf$Q21[1:tlim],
                Q23t = inpdf$Q23[1:tlim],
                Q32t = inpdf$Q32[1:tlim],
                Q34t = inpdf$Q34[1:tlim],
                Q43t = inpdf$Q43[1:tlim],
                Q4outt = inpdf$Q4out[1:tlim],
                A1t = inpdf$A1[1:tlim],
                A2t = inpdf$A2[1:tlim],
                A3t = inpdf$A3[1:tlim],
                A4t = inpdf$A4[1:tlim],
                C1int= inpdf$C1in[1:tlim],
                C2int= inpdf$C2in[1:tlim],
                C3int= inpdf$C3in[1:tlim],
                C4int= inpdf$C4in[1:tlim],
                Tmpt = inpdf$Temperature[1:tlim]))
}



######################################################
# processes model output, for predictions. Returns model results for every month

mod_wrap_pred <- function(inp_dat=inp_lis, #needs updgraded/extended inputs
                        mod_fun=lakemod_v0, 
                        CalPar=par_ini, 
                        FixInp_lis=Const_lis){
  
  CalPar_lis <- list(
    Vs= CalPar["Vs"]/vs_scal,Ks= CalPar["Ks"], Bs= CalPar["Bs"]/par_scal,Rs= CalPar["Rs"]/par_scal,
    Sh= CalPar["Sh"],
    ThetaR= CalPar["ThetaR"],ThetaV= CalPar["ThetaV"],
    S1fac= CalPar["S1fac"])

  
  OmnInp <- c(CalPar_lis,inp_dat,FixInp_lis)

  mod_cmpl <- suppressWarnings(mod_fun(user=OmnInp)) # DDG:key part missing
  y_c <- mod_cmpl$run(FixInp_lis$timestep)
  
  y_c <- data.frame(y_c)

# model would predict average water column conc 
y_c$C1 <- (y_c$M1)/inp_dat$V1t
y_c$C2 <- (y_c$M2)/inp_dat$V2t
y_c$C3 <- (y_c$M3)/inp_dat$V3t
y_c$C4 <- (y_c$M4)/inp_dat$V4t

#..........................

#get modifiers to convert avg conc to surf
CoAdj2 <- inp_dat$V2t/(inp_dat$V2t+inp_dat$V2downt*(-1+Rbs2[inp_dat$mot])) 
CoAdj3 <- inp_dat$V3t/(inp_dat$V3t+inp_dat$V3downt*(-1+Rbs3[inp_dat$mot])) 
CoAdj4 <- inp_dat$V4t/(inp_dat$V4t+inp_dat$V4downt*(-1+Rbs4[inp_dat$mot])) 


#obtain surface concentrations, useful to compare model and suface data
  y_c$C2sur=y_c$C2*CoAdj2
  y_c$C3sur=y_c$C3*CoAdj3
  y_c$C4sur=y_c$C4*CoAdj4

  #record month
  y_c$mo <- inp_dat$mot

  return(y_c %>% as_tibble())
}


####################################################
mon2seas <- function(mo){
    if (mo %in% c(12,1,2)) seas = "winter"
    if (mo %in% c(3:5)) seas = "spring"
    if (mo %in% c(6:8)) seas = "summer"
    if (mo %in% c(9:11)) seas = "autumn"
  return(seas)
}
