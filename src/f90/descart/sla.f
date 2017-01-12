      DOUBLE PRECISION FUNCTION sla_DBEAR (A1, B1, A2, B2)
*+
*     - - - - - -
*      D B E A R
*     - - - - - -
*
*  Bearing (position angle) of one point on a sphere relative to another
*  (double precision)
*
*  Given:
*     A1,B1    d    spherical coordinates of one point
*     A2,B2    d    spherical coordinates of the other point
*
*  (The spherical coordinates are RA,Dec, Long,Lat etc, in radians.)
*
*  The result is the bearing (position angle), in radians, of point
*  A2,B2 as seen from point A1,B1.  It is in the range +/- pi.  If
*  A2,B2 is due east of A1,B1 the bearing is +pi/2.  Zero is returned
*  if the two points are coincident.
*
*  P.T.Wallace   Starlink   23 March 1991
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT NONE

      DOUBLE PRECISION A1,B1,A2,B2

      DOUBLE PRECISION DA,X,Y


      DA=A2-A1
      Y=SIN(DA)*COS(B2)
      X=SIN(B2)*COS(B1)-COS(B2)*SIN(B1)*COS(DA)
      IF (X.NE.0D0.OR.Y.NE.0D0) THEN
         sla_DBEAR=ATAN2(Y,X)
      ELSE
         sla_DBEAR=0D0
      END IF

      END
      
      
      
      
      SUBROUTINE sla_DH2E (AZ, EL, PHI, HA, DEC)
*+
*     - - - - -
*      D E 2 H
*     - - - - -
*
*  Horizon to equatorial coordinates:  Az,El to HA,Dec
*
*  (double precision)
*
*  Given:
*     AZ      d     azimuth
*     EL      d     elevation
*     PHI     d     observatory latitude
*
*  Returned:
*     HA      d     hour angle
*     DEC     d     declination
*
*  Notes:
*
*  1)  All the arguments are angles in radians.
*
*  2)  The sign convention for azimuth is north zero, east +pi/2.
*
*  3)  HA is returned in the range +/-pi.  Declination is returned
*      in the range +/-pi/2.
*
*  4)  The latitude is (in principle) geodetic.  In critical
*      applications, corrections for polar motion should be applied.
*
*  5)  In some applications it will be important to specify the
*      correct type of elevation in order to produce the required
*      type of HA,Dec.  In particular, it may be important to
*      distinguish between the elevation as affected by refraction,
*      which will yield the "observed" HA,Dec, and the elevation
*      in vacuo, which will yield the "topocentric" HA,Dec.  If the
*      effects of diurnal aberration can be neglected, the
*      topocentric HA,Dec may be used as an approximation to the
*      "apparent" HA,Dec.
*
*  6)  No range checking of arguments is done.
*
*  7)  In applications which involve many such calculations, rather
*      than calling the present routine it will be more efficient to
*      use inline code, having previously computed fixed terms such
*      as sine and cosine of latitude.
*
*  P.T.Wallace   Starlink   21 February 1996
*
*  Copyright (C) 1996 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT NONE

      DOUBLE PRECISION AZ,EL,PHI,HA,DEC

      DOUBLE PRECISION SA,CA,SE,CE,SP,CP,X,Y,Z,R


*  Useful trig functions
      SA=SIN(AZ)
      CA=COS(AZ)
      SE=SIN(EL)
      CE=COS(EL)
      SP=SIN(PHI)
      CP=COS(PHI)

*  HA,Dec as x,y,z
      X=-CA*CE*SP+SE*CP
      Y=-SA*CE
      Z=CA*CE*CP+SE*SP

*  To HA,Dec
      R=SQRT(X*X+Y*Y)
      IF (R.EQ.0D0) THEN
         HA=0D0
      ELSE
         HA=ATAN2(Y,X)
      END IF
      DEC=ATAN2(Z,R)

      END
      
      
      
      
      
      SUBROUTINE sla_DE2H (HA, DEC, PHI, AZ, EL)
*+
*     - - - - -
*      D E 2 H
*     - - - - -
*
*  Equatorial to horizon coordinates:  HA,Dec to Az,El
*
*  (double precision)
*
*  Given:
*     HA      d     hour angle
*     DEC     d     declination
*     PHI     d     observatory latitude
*
*  Returned:
*     AZ      d     azimuth
*     EL      d     elevation
*
*  Notes:
*
*  1)  All the arguments are angles in radians.
*
*  2)  Azimuth is returned in the range 0-2pi;  north is zero,
*      and east is +pi/2.  Elevation is returned in the range
*      +/-pi/2.
*
*  3)  The latitude must be geodetic.  In critical applications,
*      corrections for polar motion should be applied.
*
*  4)  In some applications it will be important to specify the
*      correct type of hour angle and declination in order to
*      produce the required type of azimuth and elevation.  In
*      particular, it may be important to distinguish between
*      elevation as affected by refraction, which would
*      require the "observed" HA,Dec, and the elevation
*      in vacuo, which would require the "topocentric" HA,Dec.
*      If the effects of diurnal aberration can be neglected, the
*      "apparent" HA,Dec may be used instead of the topocentric
*      HA,Dec.
*
*  5)  No range checking of arguments is carried out.
*
*  6)  In applications which involve many such calculations, rather
*      than calling the present routine it will be more efficient to
*      use inline code, having previously computed fixed terms such
*      as sine and cosine of latitude, and (for tracking a star)
*      sine and cosine of declination.
*
*  P.T.Wallace   Starlink   9 July 1994
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT NONE

      DOUBLE PRECISION HA,DEC,PHI,AZ,EL

      DOUBLE PRECISION D2PI
      PARAMETER (D2PI=6.283185307179586476925286766559D0)

      DOUBLE PRECISION SH,CH,SD,CD,SP,CP,X,Y,Z,R,A


*  Useful trig functions
      SH=SIN(HA)
      CH=COS(HA)
      SD=SIN(DEC)
      CD=COS(DEC)
      SP=SIN(PHI)
      CP=COS(PHI)

*  Az,El as x,y,z
      X=-CH*CD*SP+SD*CP
      Y=-SH*CD
      Z=CH*CD*CP+SD*SP

*  To spherical
      R=SQRT(X*X+Y*Y)
      IF (R.EQ.0D0) THEN
         A=0D0
      ELSE
         A=ATAN2(Y,X)
      END IF
      IF (A.LT.0D0) A=A+D2PI
      AZ=A
      EL=ATAN2(Z,R)

      END


      DOUBLE PRECISION FUNCTION sla_PA (HA, DEC, PHI)
*+
*     - - -
*      P A
*     - - -
*
*  HA, Dec to Parallactic Angle (double precision)
*
*  Given:
*     HA     d     hour angle in radians (geocentric apparent)
*     DEC    d     declination in radians (geocentric apparent)
*     PHI    d     observatory latitude in radians (geodetic)
*
*  The result is in the range -pi to +pi
*
*  Notes:
*
*  1)  The parallactic angle at a point in the sky is the position
*      angle of the vertical, i.e. the angle between the direction to
*      the pole and to the zenith.  In precise applications care must
*      be taken only to use geocentric apparent HA,Dec and to consider
*      separately the effects of atmospheric refraction and telescope
*      mount errors.
*
*  2)  At the pole a zero result is returned.
*
*  P.T.Wallace   Starlink   16 August 1994
*
*  Copyright (C) 1995 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT NONE

      DOUBLE PRECISION HA,DEC,PHI

      DOUBLE PRECISION CP,SQSZ,CQSZ



      CP=COS(PHI)
      SQSZ=CP*SIN(HA)
      CQSZ=SIN(PHI)*COS(DEC)-CP*SIN(DEC)*COS(HA)
      IF (SQSZ.EQ.0D0.AND.CQSZ.EQ.0D0) CQSZ=1D0
      sla_PA=ATAN2(SQSZ,CQSZ)

      END

      SUBROUTINE sla_DS2TP (RA, DEC, RAZ, DECZ, XI, ETA, J)
*+
*     - - - - - -
*      D S 2 T P
*     - - - - - -
*
*  Projection of spherical coordinates onto tangent plane:
*  "gnomonic" projection - "standard coordinates" (double precision)
*
*  Given:
*     RA,DEC      dp   spherical coordinates of point to be projected
*     RAZ,DECZ    dp   spherical coordinates of tangent point
*
*  Returned:
*     XI,ETA      dp   rectangular coordinates on tangent plane
*     J           int  status:   0 = OK, star on tangent plane
*                                1 = error, star too far from axis
*                                2 = error, antistar on tangent plane
*                                3 = error, antistar too far from axis
*
*  P.T.Wallace   Starlink   18 July 1996
*
*  Copyright (C) 1996 Rutherford Appleton Laboratory
*  Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
*-

      IMPLICIT NONE

      DOUBLE PRECISION RA,DEC,RAZ,DECZ,XI,ETA
      INTEGER J

      DOUBLE PRECISION SDECZ,SDEC,CDECZ,CDEC,
     :                 RADIF,SRADIF,CRADIF,DENOM

      DOUBLE PRECISION TINY
      PARAMETER (TINY=1D-6)


*  Trig functions
      SDECZ=SIN(DECZ)
      SDEC=SIN(DEC)
      CDECZ=COS(DECZ)
      CDEC=COS(DEC)
      RADIF=RA-RAZ
      SRADIF=SIN(RADIF)
      CRADIF=COS(RADIF)

*  Reciprocal of star vector length to tangent plane
      DENOM=SDEC*SDECZ+CDEC*CDECZ*CRADIF

*  Handle vectors too far from axis
      IF (DENOM.GT.TINY) THEN
         J=0
      ELSE IF (DENOM.GE.0D0) THEN
         J=1
         DENOM=TINY
      ELSE IF (DENOM.GT.-TINY) THEN
         J=2
         DENOM=-TINY
      ELSE
         J=3
      END IF

*  Compute tangent plane coordinates (even in dubious cases)
      XI=CDEC*SRADIF/DENOM
      ETA=(SDEC*CDECZ-CDEC*SDECZ*CRADIF)/DENOM

      END
