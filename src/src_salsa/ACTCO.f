**************************************************************
C
        SUBROUTINE  ACTCO( ML1,MLT2,ML3,ML4,ML6,T,G1,G3,G4,G6,AW)
C
C
C  (the interface:)
C  input:     ML1,MLT2,ML3,ML4,ML6,T
C  output:   G1,G3,G4,G6,AW
C
C
        IMPLICIT DOUBLE PRECISION (A-Z)
C
C
C Debye-Hyckel constant for activity coefficient
        DOUBLE PRECISION AGAM
        PARAMETER (AGAM = 0.511 )
C
C MLi is molality of ion i
C indexies of ions: H = 1 NH4 = 3
C SO4 = 2 NO3 = 4 CL = 6  HSO4 = 8      ( 5,7 = nothing )
C
C
C stoichiometric ionic strength (no ML8)
C ( Z1=1.0, Z2=2.0, Z3=1.0, Z4=1.0, Z6=1.0, Z8=1.0 )
        I = ( ML1 + MLT2*4.0D0 + ML3 + ML4 + ML6)/2.0D0
C
C square root of ionic strength
        I12 = SQRT(I)
C
        TL = 298D0/T - 1.0D0
        TC = 1.0D0+DLOG(298D0/T)-298D0/T
C
C  the logarithms of activity coefficients of salts if alone at ionic strength I
C  from Jacobson's thesis
C  2H+/SO_4^2-  (6.0m)
        IF ( I .GT. 18.0D0 ) THEN
            IS = 18.0D0
            IS12 = SQRT(18.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F012    =  -2.316892E1  + 2.357379E-1*T  -  8.020047E-4*T*T
     J  + 9.096087E-7*T*T*T
        F112  =  (3.559700E2  -  3.656462E0*T  +   1.243911E-2*T*T
     J   -  1.409604E-5*T*T*T)/1.7320508D0
        F212  = (-6.199375E2  + 6.260075E0*T  -   2.108119E-2*T*T
     J  +  2.369106E-5*T*T*T)/3.0D0
        F312  = ( 3.107313E2  -  3.073566E0*T   +  1.018136E-2*T*T
     J   -  1.128845E-5*T*T*T)/5.19615D0
        F412  = ( -2.457443E1  + 2.147608E-1*T  -   6.330187E-4*T*T
     J   +   6.306867E-7*T*T*T)/9.0D0
        F512  = ( -6.413438E0  +  6.935867E-2*T  -  2.461235E-4*T*T
     J   +  2.877242E-7*T*T*T)/15.5885D0
C
        LG120 = (F012 + F112*IS12 + F212*IS + F312*IS*IS12 
     J  + F412*IS*IS +  F512*IS*IS*IS12) /2.30259
C
C
C  H+/HSO4 (6.0m)
        IF ( I .GT. 18.0D0 ) THEN
            IS = 18.0D0
            IS12 = SQRT(18.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  =  -1.225455E2  + 1.230707E0*T  -  4.119776E-3*T*T
     J    +  4.595935E-6*T*T*T
        F1  =  1.996164E3  -   1.997451E1*T  +  6.660462E-2*T*T
     J    -  7.403799E-5*T*T*T
        F2  =  -4.102617E3  + 4.080428E1*T  -  1.352508E-1*T*T
     J    + 1.494309E-4*T*T*T
        F3  =  2.665534E3  -  2.610425E1*T  +  8.525047E-2*T*T
     J  -   9.282615E-5*T*T*T
        F4  =  -5.421704E2  + 5.121146E0*T  -  1.611159E-2*T*T
     J   + 1.687197E-5*T*T*T
        F5  =  1.886410E1  -  1.415977E-1*T  +  3.214582E-4*T*T
     J   -  1.950665E-7*T*T*T
C
        LG180=(F0 + F1*IS12/1.7320508D0 + F2*IS/3.0D0
     J  + F3*IS*IS12/5.19615D0 + F4*IS*IS/9.0D0 
     J  + F5*IS*IS*IS12/15.5885D0 )/2.30259
C
C  the (approximate) HSO4/SO4  ratio
C
        RHSSO = 1.015D-2*TEMPCOR(8.85D0,25.14D0,T)*ML1*
     J       10**( 3*LG120 - 2*LG180 )
C
        ML2 = MLT2 /(1.0D0+RHSSO)
        ML8 = RHSSO*ML2
C
C ionic strength
        I = ( ML1 + ML2*4.0D0 + ML3 + ML4 + ML6 + ML8)/2.0D0
        I12 = SQRT(I)
C
C
C  2H+/SO_4^2-  (6.0m)
        IF ( I .GT. 18.0D0 ) THEN
            IS = 18.0D0
            IS12 = SQRT(18.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        LG120 = ( F012+F112*IS12+ F212*IS + F312*IS*IS12 
     J   + F412*IS*IS +  F512*IS*IS*IS12) /2.30259
C
C
C  H+/HSO4 (6.0m)
        IF ( I .GT. 6.0D0 ) THEN
            IS = 6.0D0
            IS12 = SQRT(6.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        LG180 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
C sum of molalities
        SUMML = ML1 + ML2 + ML3 + ML4 + ML6 + ML8
C
        IF (SUMML .EQ. 0.0 ) THEN
        AW = 1.0D0
        G1 = 1.0D0
        G3 = 1.0D0
        G4 = 1.0D0
        G6 = 1.0D0
        RETURN
        END IF
C
C
        FF1 =  2.25D0*F112*ML1*ML2  + F1*ML1*ML8
        FF2 =  2.25D0*F212*ML1*ML2  + F2*ML1*ML8
        FF3 =  2.25D0*F312*ML1*ML2  + F3*ML1*ML8
        FF4 =  2.25D0*F412*ML1*ML2  + F4*ML1*ML8
        FF5 =  2.25D0*F512*ML1*ML2  + F5*ML1*ML8
C
C
C HNO3 (28m)
        IF ( I .GT. 28.0D0 ) THEN
            IS = 28.0D0
            IS12 = SQRT(28.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  =  -2.388378D-2
        F1  =  -7.777787D-1  + 5.785894D-1*TL - 4.785171D-1*TC
        F2  =   5.950086D-1  -  9.860271D-1*TL + 6.521896D0*TC
        F3  =  -1.284278D-1  + 6.043012D-1*TL -  2.605544D0*TC
        F4  =  1.291734D-2   -  1.123169D-1*TL + 3.739984D-1*TC
        F5  =  -6.257155D-4  + 6.688134D-3*TL - 1.832646D-2*TC
C 
        LG140 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
        FF1 = FF1  + F1*ML1*ML4
        FF2 = FF2  + F2*ML1*ML4
        FF3 = FF3  + F3*ML1*ML4
        FF4 = FF4  + F4*ML1*ML4
        FF5 = FF5  + F5*ML1*ML4
C
C
C  HCl (16m)
        IF ( I .GT. 16.0D0 ) THEN
            IS = 16.0D0
            IS12 = SQRT(16.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  =  -1.998104E-2
        F1  =  -7.959068E-1 + 5.532198E-1*TL + 2.108728E0*TC
        F2  =   6.580198E-1 -  2.571126E-1*TL + 8.542292E-1*TC
        F3  =  -7.409941E-2 + 2.790048E-1*TL - 6.237459E-1*TC
        F4  =  1.345075E-2 - 4.691631E-2*TL + 1.935911E-1*TC
        F5  =  -2.248651E-3 + 2.382485E-3*TL - 2.037543E-2*TC
C
         LG160 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
        FF1 = FF1  + F1*ML1*ML6
        FF2 = FF2  + F2*ML1*ML6
        FF3 = FF3  + F3*ML1*ML6
        FF4 = FF4  + F4*ML1*ML6
        FF5 = FF5  + F5*ML1*ML6
C
C
C  (NH4)2SO4 (5.8m)
        IF ( I .GT. 17.4D0 ) THEN
            IS = 17.4D0
            IS12 = SQRT(17.4D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  = -2.163694D-2 + 2.297972D-1*TL
        F1=(-3.377941D0+4.255129D-1*TL+1.609902D-3*TC)/1.7320508D0
        F2=(3.118007D0  -  2.220594D0*TL + 4.437758D0*TC)/3.0D0
        F3=(-1.920544D0  + 2.607601D0*TL + 6.101756D-3*TC)/5.19615D0
        F4=(6.372975D-1 -  1.243384D0*TL  + 4.021805D-1*TC)/9.0D0
        F5=(-8.277292D-2 +2.102563D-1*TL + 4.375833D-4*TC)/15.5885D0
C
         LG320 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
        FF1 = FF1  + 2.25D0*F1*ML3*ML2
        FF2 = FF2  + 2.25D0*F2*ML3*ML2
        FF3 = FF3  + 2.25D0*F3*ML3*ML2
        FF4 = FF4  + 2.25D0*F4*ML3*ML2
        FF5 = FF5  + 2.25D0*F5*ML3*ML2
C
C
C NH4NO3 (25.9m)
        IF ( I .GT. 25.9D0 ) THEN
            IS = 25.9D0
            IS12 = SQRT(25.9D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  =  -1.044572D-2
        F1  =  -1.004940D0  + 4.362921D-1*TL + 2.611682D0*TC
        F2  =   4.674064D-1 - 1.455444D0*TL + 3.158677D0*TC
        F3  =  -1.750495D-1 + 6.282104D-1*TL - 2.005748D0*TC
        F4  =   3.253844D-2 - 1.123507D-1*TL + 4.113737D-1*TC
        F5  =  -2.276789D-3 + 7.438990D-3*TL - 2.820677D-2*TC
C
         LG340 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
        FF1 = FF1  + F1*ML3*ML4
        FF2 = FF2  + F2*ML3*ML4
        FF3 = FF3  + F3*ML3*ML4
        FF4 = FF4  + F4*ML3*ML4
        FF5 = FF5  + F5*ML3*ML4
C
C
C NH4CL (7.4m)
        IF ( I .GT. 7.4D0 ) THEN
            IS = 7.4D0
            IS12 = SQRT(7.4D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  =  -5.022484D-3
        F1  =  -1.037873D0  + 4.890513D-1*TL + 1.959107D0*TC
        F2  =   8.517483D-1 - 7.013315D-1*TL  +  9.894682D-1*TC
        F3  =  -4.225323D-1 + 4.682151D-1*TL - 1.024499D-1*TC
        F4  =   1.214996D-1 - 1.702461D-1*TL - 2.354376D-1*TC
        F5  =  -1.471525D-2 + 2.502643D-2*TL + 6.600384D-2*TC
C
         LG360 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
        FF1 = FF1  + F1*ML3*ML6
        FF2 = FF2  + F2*ML3*ML6
        FF3 = FF3  + F3*ML3*ML6
        FF4 = FF4  + F4*ML3*ML6
        FF5 = FF5  + F5*ML3*ML6
C
C
C NH4HSO4 (6.0m)
        IF ( I .GT. 6.0D0 ) THEN
            IS = 6.0D0
            IS12 = SQRT(6.0D0)
        ELSE
            IS = I
            IS12 = I12
        END IF
C
        F0  = -2.708121D-3
        F1  = -1.095646D0 + 4.890513D-1*TL + 1.959107D0*TC
        F2  =  1.042878D0 -  7.013315D-1*TL + 9.894682D-1*TC
        F3  = -6.289405D-1 + 4.682151D-1*TL - 1.024499D-1*TC
        F4  =  2.079631D-1 - 1.702461D-1*TL - 2.354376D-1*TC
        F5  = -2.776957D-2 + 2.502643D-2*TL + 6.600384D-2*TC
C
         LG380 = ( F0 + F1*IS12 + F2*IS + F3*IS*IS12 + F4*IS*IS 
     J        + F5*IS*IS*IS12 )/2.30259
C
        FF1 = FF1  + F1*ML3*ML8
        FF2 = FF2  + F2*ML3*ML8
        FF3 = FF3  + F3*ML3*ML8
        FF4 = FF4  + F4*ML3*ML8
        FF5 = FF5  + F5*ML3*ML8
C
C the final values
        FF1 =  2.0D0*FF1/(SUMML*I)
        FF2 =  2.0D0*FF2/(SUMML*I)
        FF3 =  2.0D0*FF3/(SUMML*I)
        FF4 =  2.0D0*FF4/(SUMML*I)
        FF5 =  2.0D0*FF5/(SUMML*I)
C
C
C  [Z+Z-]
        ZPLZMI=(ML1 + ML2*4.0D0 + ML3 + ML4 + ML6 + ML8)/SUMML
C
C [Z1Z2] = [Z3Z2] = 2
C [Z1Z4] = [Z1Z6] = [Z3Z4] = [Z3Z6] = [Z1Z8] = [Z3Z8] = 1
C
        ALPHA = -2.0D0 * (  4.5D0*ML1*ML2
     J   +   ML1*ML4 + ML1*ML6 +  ML1*ML8
     J   +   4.5D0*ML3*ML2 +  ML3*ML4
     J   +   ML3*ML6 + ML3*ML8 )    / (SUMML*I)
     J   +   ZPLZMI
C
C
C formula 6 in Bromleys article
        SIG = ( 3.0D0/( I12 * I12 * I12 ) ) *
     J (1.0D0 + I12 - 1.0D0/(1.0D0 + I12) - 2.0D0*DLOG(1.0D0 + I12))
C
C
C from the gibbs-duhem equation
C
        OSMCOE = 1.0D0
     J  -   2.30259D0*ALPHA*AGAM*SIG*I12/3.0D0
     J  +  FF1*I12/3.0D0 + FF2*I/2.0D0 + FF3*3.0D0*I*I12/5.0D0
     J  +  FF4*2.0D0*I*I/3.0D0 +  FF5*5.0D0*I*I*I12/7.0D0
C
C activity of water
        LNACT1 = (- 18.02D-3* OSMCOE * SUMML )
        IF (LNACT1 .LT. -1.0D99 ) LNACT1 = -1.0D99
        IF (LNACT1 .GT. 0.0D0 ) LNACT1 = 0.0D0
        AW = EXP( LNACT1 )
C
C
C calculation of activity coefficints of ions
C
C formula 27
C
        Y21 = 2.25D0*ML2/I              !  Z1 = 1  Z2 = 2
        Y41 =  ML4/I                    !  Z1 = Z4 = 1
        Y61 =  ML6/I                    !  Z1 = Z6 = 1
        Y81 =  ML8/I                    !  Z1 = Z8 = 1
        Y23 = 2.25D0*ML2/I              !  Z3 = 1  Z2 = 2
        Y43 =  ML4/I                    !  Z3 = Z4 = 1
        Y63 =  ML6/I                    !  Z3 = Z6 = 1
        Y83 =  ML8/I                    !  Z3 = Z8 = 1
C
C formula 28
C
        X14 = ML1/I             !  Z1 = Z4 = 1
        X34 = ML3/I             !  Z3 = Z4 = 1
        X16 = ML1/I             !  Z1 = Z6 = 1
        X36 = ML3/I             !  Z3 = Z6 = 1
C
C formula 25
C
        F1 = Y21*LG120 + Y41*LG140 + Y61*LG160 + Y81*LG180
     J  + AGAM * I12 * ( 2.0*Y21 + Y41 + Y61 + Y81 )/(1+I12)
C
        F3 = Y23*LG320 + Y43*LG340 + Y63*LG360 + Y83*LG380
     J  + AGAM * I12 * ( 2.0*Y23 + Y43 + Y63 + Y83 )/(1+I12)
C
C formula 26
C
        F4 = X14*LG140 + X34*LG340
     J    + AGAM * I12 * (  X14 + X34 )/(1+I12)
C
        F6 = X16*LG160 + X36*LG360
     J    + AGAM * I12 * (  X16 + X36 )/(1+I12)
C
C calculation of activity coefficients of single ions
C formula 16
C
        LG1 = ( - AGAM*I12/(1+I12)  + F1   )
        IF ( LG1 .LT. -99.0 ) LG1 = -99.0
        IF ( LG1 .GT.  99.0 ) LG1 =  99.0
        G1 = 10.0D0**( LG1  )
C
        LG3 = ( - AGAM*I12/(1+I12)  + F3   )
        IF ( LG3 .LT. -99.0 ) LG3 = -99.0
        IF ( LG3 .GT.  99.0 ) LG3 =  99.0
        G3 = 10.0D0**( LG3 )
C
        LG4 = ( - AGAM*I12/(1+I12)  + F4   )
        IF ( LG4 .LT. -99.0 ) LG4 = -99.0
        IF ( LG4 .GT.  99.0 ) LG4 =  99.0
        G4 = 10.0D0**( LG4 )
C
        LG6 = ( - AGAM*I12/(1+I12)  + F6   )
        IF ( LG6 .LT. -99.0 ) LG6 = -99.0
        IF ( LG6 .GT.  99.0 ) LG6 =  99.0
        G6 = 10.0D0**( LG6 )
C
C
        RETURN
        END
C
*************************************************************
       DOUBLE PRECISION FUNCTION TEMPCOR(A,B,T)
       IMPLICIT DOUBLE PRECISION (A-Z)
C
C this function calculates the temperature corrections for
C equilibrium constants at
C temperatures different from 298 K
C according to the Van't Hoff equation
C
C inputs
C    A  is delta(H0)/(RT0)
C    B   is delta(CP0/R)
C    T   is the actual temperature in kelvins
C
      T0=298D0
      TR=T0/T
      COR=A*(TR-1.0D0)+B*(1.0D0+DLOG(TR)-TR)
      TEMPCOR=DEXP(COR)
C
C
      RETURN
      END
C                                                                        
