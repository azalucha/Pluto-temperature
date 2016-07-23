      SUBROUTINE GET_QO(lat,ls,inpFlux,albedo,axialTilt,
     &                         eccentricity,lsp,qo)

C Subroutine to calculate solar forcing
C ==Inputs==
C lat: latitude (degrees)
C ls: "L sub s" (ecliptic longitude, degrees)
C inpFlux: mean solar constant (W/m^2)
C albedo: surface albedo
C axialTilt: obliquity in IAU convention where North Pole is 
C            in same hemisphere as Sun's (i.e. Pluto rotates "backwards")
C            (in degrees)
C lsp: ls of perihelion (degrees)
C ==Outputs==
C qo (output, W/m^2)

      IMPLICIT NONE

      DOUBLE PRECISION lat,ls,inpFlux,albedo,qo,axialTilt
      DOUBLE PRECISION solarDec,haossArg,hourAngleOfSunsetSunrise
      DOUBle PRECISION pi,lsp,eccentricity

      pi=3.14159265359

      lat=lat*2.*pi/360.
      ls=ls*2.*pi/360.
      axialTilt=axialTilt*2.*pi/360
      lsp=lsp*2.*pi/360.

      solarDec = DASIN(DSIN(axialTilt)*DSIN(ls))
      haossArg = -DTAN(lat)*DTAN(solarDec)
      IF(haossArg .LT. -1.) haossArg = -1.
      IF(haossArg .GT. 1.) haossArg = 1.
      hourAngleOfSunsetSunrise = DACOS(haossArg)
      qo = inpFlux/pi*(1. - albedo)*
     &       ((1. + eccentricity*DCOS(ls - lsp))
     &       /(1. - eccentricity**2.))**2.
     &       *(DSIN(lat)*DSIN(solarDec)*hourAngleOfSunsetSunrise
     &       + DCOS(lat)*DCOS(solarDec)*DSIN(hourAngleOfSunsetSunrise))

      RETURN
      END
