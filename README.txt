For: Fortran 77 or later

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
