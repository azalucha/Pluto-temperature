      PROGRAM YELLE_T

C---------------------------------------------------------------
C Calculate Pluto's radiative-conductive atmospheric
C temperature
C Author: Angela Zalucha
C Date: 22 July 2016
C Inputs: methane mixing ratio, surface temperature,
C surface pressure, surface albedo, season, latitude
C (see user inputs section)
C needs: get_qo.f, tridiag.f
C To run: ./yelle_T.f > output.txt (name out output file
C specified by user, otherwise will print to screen)
C Outputs: index, altitude (km), pressure (microbars), 
C temperature (k)
C About: calculates radiative-conductive temperature as
C shown in Yelle, R. V. and Lunine}, J. I.,
C "Evidence for a molecule heavier than methane in the 
C  atmosphere of {Pluto}",
C  Nature 339, 288-290 (1989), doi 10.1038/339288a0.
C See also Zalucha, et al., "An analysis of {P}luto occultation 
C light curves using an atmospheric radiative-conductive model",
C Icarus 211, 804-818 (2011), doi 10.1016/j.icarus.2010.08.018.
C---------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nlevs,k,niters,nmax,pertind
      PARAMETER(nlevs=140)
C      PARAMETER(nlevs=180) 

      DOUBLE PRECISION ko,alpha,wl,wn,piF,nlosch,a10,mun2
      DOUBLE PRECISION much4,molwgt,p10,bigGrav,bodyMass,gamma0
      DOUBLE PRECISION mi,bigR,poo,amu,z10,kb,g,c,h,cp,ts
      DOUBLE PRECISION rm(nlevs),rs,drmm(nlevs),drmp(nlevs)
      DOUBLE PRECISION temp(nlevs),zm(nlevs),scaleheight
      DOUBLE PRECISION p(nlevs),pint(nlevs),na(nlevs),n(nlevs)
      DOUBLE PRECISION q(nlevs),l(nlevs),epsilon(nlevs)
      DOUBLE PRECISION rnet(nlevs),xi(nlevs),bigG(nlevs)
      DOUBLE PRECISION drmb(nlevs),xiPrev(nlevs)
      DOUBLE PRECISION alpham(nlevs),lambdam(nlevs),flag
      DOUBLE PRECISION omegam(nlevs),qm(nlevs),deltaT,um(nlevs)
      DOUBLE PRECISION error,tempPrev(nlevs),bs,rho(nlevs)
      DOUBLE PRECISION bigC(nlevs),bigS(nlevs),littleG
      DOUBLE PRECISION moldiam,z10b,c10,epsln,eps
      DOUBLE PRECISION wll,wnl,nul,a10l,epslnl,epsl
      DOUBLE PRECISION fwhm,rp,a1,delT,temppert(nlevs)
      DOUBLE PRECISION dtdzhi,condtermhi,dtdzlo,condtermlo
      DOUBLE PRECISION condterm,ppert,napert,npert,rhopert
      DOUBLE PRECISION qpert,lpert,dfdz,tau,two
      DOUBLE PRECISION a10s,wls,wns,bss,epslns,epss,piFs,qs(nlevs)
      DOUBLE PRECISION c10l,c10s,p10l,p10s,lat,ls,piF0,piF1
      DOUBLE PRECISION albedo,axialTilt,lsp,eccentricity

C*************** User Inputs (in MKS)**************
C fraction of CH4 (dimensionless)
      gamma0=0.006
C surface pressure (kg/s^2/m=Pa)
C note: 1 microbar = 1 g/s^2/cm
C note: 1 microbar = 0.1 Pa
      poo=1.4
C surface temperature, remains fixed (K)
      ts=37.
C surface radius (km)
C Pluto
       rs=1190.
C latitude (degrees)
      lat=-32.
C ls (degrees)
      ls=226.17
C surface albedo (dimensionless)
      albedo=0.
C solar forcing at 3.3 um (W/m^2)
C Pluto
       piF0=13.23805d-3
C*************** END User Inputs*******************

       rs=rs*1000.

C 3.3 micron heating
C 7.6 micron cooling


C diffusion coefficient (J/m/s/K^(alpha+1))
      ko=5.63d-5
C exponent on temperature in diffusion term (dimensionless)
      alpha=1.12
C wavelength for heating (microns)
      wl=3.3
C convert to wavenumber for heating (m^-1)
      wn=1.d6/wl
C Loschmidt number (m^-3)
      nlosch=2.6868d25
C band strength (m)
      bs=30.*(wn/wl)/nlosch
C Einstein coefficient at 3.3 microns (s^-1)
      a10=4.24
C molecular weight of N2 (g/mol)
      mun2=28.
C molecular weight of CH4 (g/mol)
      much4=16.0426
C molecular weight (g/mol)
      molwgt=(mun2*much4)/(mun2+much4)
      mi=molwgt
C probability that a V-T transition occurs (3.3 microns) (dimensionless)
      p10=1.d-6
C probability that a V-T transition occurs (7.6	microns) (dimensionless)
      p10l=p10
C gravitational constant (m^3/kg/s^2)
      bigGrav=6.67d-11
C Body mass (kg)
      bodyMass=2.13975d22
C specific gass constant (J/K/kg)
      bigR=(8.314)/(mi*1.d-3)
C atmomic mass unit (kg)
      amu=1.66d-27
C collision rate (m^3/s)
      z10=1.5d-16
C Boltzmann Constant (J/K)
      kb=1.308658d-23
C ratio of statistical weights (dimensionless)
      g=3.
C speed of light (m/s)
      c=3.d8
C Planck's constant (J/s)
      h=6.626d-34
C specific heat at constant pressure (J/K/kg)
      cp=1043.
C gravitational acceration at the surface (m/s^2)
      littleG=0.7792
C axial tilt (obliquity, degrees)
      axialTilt=28.32
C orbital eccentricity
      eccentricity=0.0097
C ls of perihelion (degrees)
      lsp=37.9835
C wavelength for cooling (microns)
      wll=7.8
C convert to wavenumber (m^-1)
      wnl=1.d6/wll
C convert to frequency (1/s)
      nul=c*wnl
C Einstein coefficient at 7.6 microns (s^-1)
      a10l=2.56
C epsilon at 3.3 micron      
      c10=p10*z10*0.1/kb/100.
      epsln=c10/a10
      eps=epsln/(1.+epsln)
C epsilon at 7.6 micron
      c10l=p10l*z10*0.1/kb/100.
      epslnl=c10l/a10l
      epsl=epslnl/(1.+epslnl)

C timestep in numerical scheme (s)
      deltaT=1000000.
C tolerance for numerical convergence (K)
      error=0.0001
C initialize the number of iterations completed
      niters=0
C maximum number of iterations before stopping
      nmax=500000

C solar forcing at 3.3 um
      CALL GET_QO(lat,ls,piF0,albedo,axialTilt,
     &                            eccentricity,lsp,piF1)
      piF=piF1*(wl/wn)

C Set up the vertical grid
C rm is the radius in m
C drmm is delta r-minus from Zalucha et al. 2011
C drmp is delta r-plus from Zalucha et al. 2011
C drmb is delta r-bar from Zalucha et al. 2011
      rm(1)=rs
      rm(2)=rs+300.
      drmm(1)=0.
      drmm(2)=300.
      drmp(1)=300.

      DO k=2,nlevs-1
       drmp(k)=1.03*drmm(k)
       rm(k+1)=rm(k)+drmp(k)
       drmm(k+1)=drmp(k)
       drmb(k)=(rm(k+1)-rm(k-1))/2.
      ENDDO

C altitude
      DO k=1,nlevs
       zm(k)=rm(k)-rs
      ENDDO

C initialize temperatures.  Here level 1 is surface temperature,
C everywhere else is a constant
      temp(1)=ts
      DO k=2,nlevs
       temp(k)=80.
      ENDDO

      flag=1


C begin implicit time stepping
      DO WHILE(flag .GT. 0)
       DO k=1,nlevs
        xi(k)=temp(k)**(alpha+1.)
        tempPrev(k)=temp(k)
        xiPrev(k)=xi(k)
       ENDDO


C get pressure from integrating hydrostatic balance equation
C variable temperature, variable gravity with height
       pint(1)=0.
       DO k=2,nlevs
        pint(k)=pint(k-1)-bigGrav*bodyMass/rm(k)**2./bigR/
     &            temp(k)*drmm(k)
       ENDDO
       DO k=1,nlevs
        p(k)=poo*DEXP(pint(k))
C number density of atmosphere
        na(k)=p(k)/kb/temp(k)
C number density of CH4
        n(k)=gamma0*na(k)
C mass density of atmosphere
        rho(k)=mi*amu*na(k)
        epsilon(k)=eps
C heating rate (J/m^3/s) 
        q(k)=piF*epsilon(k)*n(k)*bs
C cooling rate (J/m^3/s)
        l(k)=g*epsl*n(k)*h*nul*a10l*
     &      DEXP(-h*nul/kb/temp(k))
C net heating rate
        rnet(k)=q(k)-l(k)
C net heating rate with 2.3 micron heating
C         rnet(k)=q(k)-l(k)+qs(k)
CCCCC
C leading coefficients in dxi/dr, see Zalucha et al. 2011
        bigG(k)=ko*temp(k)**alpha/cp/rho(k)
        bigC(k)=temp(k)**alpha*2.*ko/cp/rho(k)/rm(k)
        bigS(k)=temp(k)**alpha/cp/rho(k)*(1.+alpha)*rnet(k)
       ENDDO

C coefficients in xi, see Zalucha et al. 2011
       alpham(1)=0.
       lambdam(1)=1.
       omegam(1)=0.
       qm(1)=temp(1)**(alpha+1.)

       alpham(nlevs)=-1.
       lambdam(nlevs)=1.
       omegam(nlevs)=0.
       qm(nlevs)=0.

       DO k=2,nlevs-1
        alpham(k)=deltaT/drmb(k)*(-bigG(k)/drmm(k)+bigC(k)/2.)
        lambdam(k)=1.+2.*deltaT*bigG(k)/drmp(k)/drmm(k)
        omegam(k)=-deltaT/drmb(k)*(bigG(k)/drmp(k)+bigC(k)/2.)
        qm(k)=xi(k)+deltaT*bigS(k)
       ENDDO


C Call numerical recipes in fortran tridiagonal matrix solver
       CALL TRIDAG(alpham,lambdam,omegam,qm,um,nlevs)

       flag=0

C Convert back to temperature.  Check for convergence.
       DO k=1,nlevs
        temp(k)=um(k)**(1./(alpha+1.))
        IF(DABS(um(k)-xiPrev(k)) .GT. error) flag=flag+1
       ENDDO

       niters=niters+1
       IF(niters .GT. nmax) flag=0

      ENDDO

C write output 
      write(*,*) 'niters=',niters

      DO k=1,nlevs
       write(*,*) 'k=',k,', z=',(rm(k)-rs)/1000.,
     &            ' km, p=',p(k)*10.,' microbar, T=',temp(k),' K'
      ENDDO

      STOP
      END
