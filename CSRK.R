CSRK <- function(cpd,warn=FALSE){
   ## Volume and fugacity coefficient of Carnahan-Starling Redlich-Kwong
   ## EoS for fluids.
   ## Returns function f(P,T) for volume in J/bar and log(gamma)
   ##   f = gamma*P; input f(P,T) is T(K), P(bars)

   ## warn parameter TRUE enables warning if hard sphere density yields a solid

   crit <- switch(cpd,
   He  = c(  5.1953, 2.2746),
   Ne  = c( 44.49,  27.686),
   Ar  = c(150.687, 48.63),
   H2  = c( 41.2,   21.1),
   CO  = c(133.16,  34.98),
   N2  = c(126.20,  33.95),
   O2  = c(154.581, 50.430),
   CH4 = c(190.564, 45.922),
   CO2 = c(304.18,  73.80),
   NH3 = c(405.50, 112.77),
   H2O = c(647.3,  220.6),
   H2S = c(373.4,  89.63),
   SO2 = c(430.65, 78.84),
         c(NA,NA)
   )
   if (is.na(crit[1])) stop(sprintf('CSRK: %s unknown species',cpd),call.=FALSE)

   ab  <- switch(cpd,
   Ar =  list(a=c(-0.700190,0.722121,0),           b=0.153066),
   He =  list(a=c(-8.88469,0.191372,11.4976),      b=0.0676166),
   Ne =  list(a=c(0.552466,0.184632,-0.408789),    b=0.149340),
   H2 =  list(a=c(-2.17856,0,0),                   b=0.080398),
   CO =  list(a=c(-4.30923,1.29948,5.59317),       b=0.156382),
   N2 =  list(a=c(1.56616,0,-2.96891),             b=0.141462),
   O2 =  list(a=c(1.00862,0.113076,-0.328853),     b=0.156248),
   CH4 = list(a=c(-0.0798752,0.220581,0.474302),   b=0.133213),
   CO2 = list(a=c(0.389276,0.195472,0),            b=0.135557),
   NH3 = list(a=c(0.946982,0.0798373,-0.279543),   b=0.127271),
   H2O = list(a=c(1.08193,0.0296910,-0.337474),    b=0.117309),
   H2S = list(a=c(0.440691,0.0864720,0),           b=0.121762),
   SO2 = list(a=c(0.890050,0,0),                   b=0.154133)
   )
   ac <-                                                ## a scaling by Tc**(5/2)/Pc
      ab$a*crit[1]^(3/2)*c(crit[1],1,crit[1]^2)/crit[2] ## CORK a + 1/T term
   bc <-                                                ## b scaling by Tc/Pc
      ab$b*crit[1]/crit[2]  
   id <- cpd

   function(P,T,R=8.3144,PF=0.74048048,tol=5e-8,UR=FALSE){
      fxi <- function(xi){
         P - fRTob*xi/(1-xi)^3 * (1 + xi + xi^2 - xi^3) + faosT*xi^2/(1+4*xi)
      }
      dfxi <- function(xi){
         ## Derivative of fxi, dP/d(xi) for Newton-Raphson root finding
         -fobsrT/((4*xi+1)^2*(xi-1)^4) * (
            bRTsrT*(
               16*xi^6 - 56*xi^5 + 33*xi^4 + 92*xi^3 + 52*xi^2 + 12*xi + 1
            ) -
            8*av*xi*(2*xi^5 - 7*xi^4 + 8*xi^3 - 2*xi^2 - 2*xi + 1)
         )
      }
      av <-                      ## Redlich-Kwong a
         R^2*sum(ac*c(1,T,1/T))  ## CORK + 1/T recipe
      bv <-                      ## Redlich-Kwong b
         R*bc                    ## standard, also CORK
      RT <- R*T
      fRTob <- 4*RT/bv
      sqrT <- sqrt(T)
      Pbo4RT <- P*bv/(4*RT)
      faosT <- 16*av/(sqrT*bv^2)
      bRTsrT <- bv*RT*sqrT
      fobsrT <- 4/(sqrT*bv^2)
      xs <- c(0,0.2,0.4,0.5,0.6,0.8); xl <- length(xs)
      fx <- vapply(xs,fxi,pi); fs <- fx[1]*fx <= 0
      if(UR && any(fs)) {
         ## R root finder (sublinear convergence)
         xr <- sum(!fs)
         xi <- uniroot(fxi,xs[xr+c(0,1)],tol=tol)$root
      } else {
         ## Newton-Raphson (quadratic convergence)
         x <- 0.5; dfdx <- tol + 1; it <- 0; itmax <- 60
         while (abs(dfdx) >= tol && it < itmax){
            it <- it + 1
            dfdx <- -fxi(x)/dfxi(x)
            x <- min(1.15*PF,max(0,x + dfdx))
         }
         if (it >= itmax && abs(dfdx)/P > 1e-5)warning(
            sprintf("CSRK: Unconverged %s %.3g[%.0f,%.0f]",id,dfdx,P,T),
            call.=FALSE
         )
         xi <- x
      }
      if (xi > 0.58 && warn) warning(
         sprintf("CSRK: Metastable %s %.3f[%.0f,%.0f]",id,xi,P,T),
         call.=FALSE
      )
      vol <- bv/(4*xi)
      z <- vol*P/RT
      attr(vol,'gamma') <- z - 3 - log(z) - (xi^2 - 2)/(1 - xi)^2 -
         av/(bv*sqrT*RT) * log(4*xi + 1)
      vol
   }
}

# Test calculation with H2O at 15 GPa and 1500 K

H2O <- CSRK('H2O')

vol <- H2O(150000,1500) 

cat(sprintf("H2O volume at 15 GPa 1500 K: %.4f J/b\n",vol))
cat(sprintf("H2O fugacity at 15 GPa 1500 K: %.0f b\n",
   exp(attr(vol,'gamma') + log(150000))))
