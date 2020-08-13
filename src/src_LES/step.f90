!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE step

  USE mo_submctl, ONLY : spec
  USE util, ONLY : getMassIndex, calc_correlation 
  
  IMPLICIT NONE
  
  INTEGER :: istpfl = 1
  REAL    :: timmax = 18000.
  LOGICAL :: corflg = .FALSE.
  
  REAL    :: frqhis =  9000.
  REAL    :: frqanl =  3600.
  REAL    :: radfrq =  0.
  
  REAL    :: time   =  0.
  REAL    :: strtim =  0.0    ! In decimal days, 0.5 mid-day
  LOGICAL :: outflg = .TRUE.
  
CONTAINS
   !
   ! ----------------------------------------------------------------------
   ! Subroutine model:  This is the main driver for the model's time
   ! integration.  It calls the routine tstep, which steps through the
   ! physical processes active on a time-step and updates variables.  It
   ! then checks to see whether or not different output options are
   ! satisfied.
   SUBROUTINE stepper

      USE mpi_interface, ONLY : myid, double_scalar_par_max

      USE grid, ONLY : dtl, dzt, zt, zm, nzp, dn0, u0, v0, a_up, a_vp, a_wp, &
                       a_uc, a_vc, a_wc, write_hist, write_anal, close_anal, dtlt,  &
                       dtlv, dtlong, nzp, nyp, nxp, level,                          &
                       ! For mass budged
                       a_rp, a_rc, a_srp, a_dn

      USE stat, ONLY : sflg, savg_intvl, ssam_intvl, write_ps, close_stat, mcflg, acc_massbudged,  &
                       write_massbudged
      USE thrm, ONLY : thermo

      LOGICAL, PARAMETER :: StopOnCFLViolation = .FALSE.
      REAL, PARAMETER :: cfl_upper = 0.50, cfl_lower = 0.30

      REAL         :: t1,t2,tplsdt,begtime
      REAL(kind=8) :: cflmax,gcflmax
      INTEGER      :: istp, iret
      LOGICAL :: cflflg
      !
      ! Timestep loop for program
      !
      begtime = time
      istp = 0

      CALL cpu_time(t1)

      DO WHILE (time + 0.1*dtl < timmax)

         istp = istp+1
         tplsdt = time + dtl + 0.1*dtl
         sflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dtl  &
            .OR. tplsdt >= timmax  .OR. tplsdt < 2.*dtl)

         CALL t_step(cflflg,cflmax)

         time = time + dtl

         CALL double_scalar_par_max(cflmax,gcflmax)
         cflmax = gcflmax

         IF (cflmax > cfl_upper .OR. cflmax < cfl_lower) THEN
            CALL tstep_reset(nzp,nxp,nyp,a_up,a_vp,a_wp,a_uc,a_vc,a_wc,     &
                             dtl,dtlong,cflmax,cfl_upper,cfl_lower)
            dtlv = 2.*dtl
            dtlt = dtl
         END IF

         !
         ! output control
         !
         IF (mod(tplsdt,savg_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_ps(nzp,dn0,u0,v0,zm,zt,time)

         IF ((mod(tplsdt,frqhis) < dtl .OR. time >= timmax) .AND. outflg)   &
            CALL write_hist(2, time)
         IF (mod(tplsdt,savg_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_hist(1, time)

         IF ((mod(tplsdt,frqanl) < dtl .OR. time >= timmax) .AND. outflg) THEN
            CALL thermo(level)
            CALL write_anal(time)
         END IF

         IF (cflflg) THEN
            cflflg = .FALSE.
            IF (StopOnCFLViolation) CALL write_hist(-1,time)
         END IF

         IF(myid == 0) THEN
            CALL cpu_time(t2) ! t2-t1 is the actual time from output
            IF (mod(istp,istpfl) == 0 ) THEN
               PRINT "('   Timestep # ',i6," //     &
                  "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3)",     &
                  istp, time, t2-t1
               CALL cpu_time(t1)
            END IF
         END IF

      END DO

      IF (mcflg) THEN
         !
         ! Juha:
         ! Get the final statistics of atmospheric water for mass budged
         CALL acc_massbudged(nzp,nxp,nyp,1,dtlt,dzt,a_dn,    &
                             rv=a_rp,rc=a_rc,prc=a_srp)

         CALL write_massbudged

      END IF ! mcflg

      CALL write_hist(1, time)
      iret = close_anal()
      iret = close_stat()

   END SUBROUTINE stepper
   !
   !----------------------------------------------------------------------
   ! Subroutine tstep_reset: Called to adjust current velocity and reset
   ! timestep based on cfl limits
   !
   SUBROUTINE tstep_reset(n1,n2,n3,up,vp,wp,uc,vc,wc,dtl,dtmx,cfl,c1,c2)

      INTEGER, INTENT (in)      :: n1,n2,n3
      REAL, INTENT (in)         :: up(n1,n2,n3),vp(n1,n2,n3),wp(n1,n2,n3),dtmx,c1,c2
      REAL(kind=8), INTENT (IN) :: cfl
      REAL, INTENT (inout)      :: uc(n1,n2,n3),vc(n1,n2,n3),wc(n1,n2,n3),dtl

      INTEGER :: i,j,k
      REAL    :: cbar, dtl_old

      cbar = (c1+c2)*0.5
      dtl_old = dtl

      IF (cfl > c1) dtl = min(dtmx,dtl*cbar/c1)
      IF (cfl < c2) dtl = min(dtmx,dtl*cbar/c2)

      DO j = 1, n3
         DO i = 1, n2
            DO k = 1, n1
               uc(k,i,j) = up(k,i,j) + (uc(k,i,j)-up(k,i,j))*dtl/dtl_old
               vc(k,i,j) = vp(k,i,j) + (vc(k,i,j)-vp(k,i,j))*dtl/dtl_old
               wc(k,i,j) = wp(k,i,j) + (wc(k,i,j)-wp(k,i,j))*dtl/dtl_old
            END DO
         END DO
      END DO

   END SUBROUTINE tstep_reset

   ! 
   !----------------------------------------------------------------------
   ! Subroutine set_LES_runtime: Set the status of process switches e.g.
   ! if they have a defined spinup time etc.
   !
   SUBROUTINE set_LES_runtime(time)
     USE mcrp, ONLY : sed_aero,   &
                      sed_cloud,  &
                      sed_precp,  &
                      sed_ice,    &
                      sed_snow,   &
                      bulk_autoc
     USE grid, ONLY : level
     IMPLICIT NONE

     REAL, INTENT(in) :: time

     IF ( sed_aero%switch .AND. time > sed_aero%delay ) sed_aero%state = .TRUE.
     IF ( sed_cloud%switch .AND. time > sed_cloud%delay ) sed_cloud%state = .TRUE.
     IF ( sed_precp%switch .AND. time > sed_precp%delay ) sed_precp%state = .TRUE.
     IF ( sed_ice%switch .AND. time > sed_ice%delay ) sed_ice%state = .TRUE.
     IF ( sed_snow%switch .AND. time > sed_snow%delay ) sed_snow%state = .TRUE.
     IF ( bulk_autoc%switch .AND. time > bulk_autoc%delay ) bulk_autoc%state = .TRUE.
     IF (level < 5) THEN
        sed_ice%state = .FALSE.
        sed_snow%state = .FALSE.
     END IF

   END SUBROUTINE set_LES_runtime

   !
   !----------------------------------------------------------------------
   ! Subroutine t_step: Called by driver to timestep through the LES
   ! routines.  Within many subroutines, data is accumulated during
   ! the course of a timestep for the purposes of statistical analysis.
   !
   SUBROUTINE t_step(cflflg,cflmax)

      USE grid, ONLY : level, dtl, dtlt,                                         &
                       ! Added parameters for interfacing with SALSA
                       nxp, nyp, nzp, a_press, a_temp, a_rsl,                             &
                       a_rc, a_wp, a_rp, a_rt, a_rh,                                      &
                       a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                       a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                       a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                       a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                       a_gaerop, a_gaerot, a_dn,  a_nactd,  a_vactd,            &
                       a_rsi

      USE stat, ONLY : sflg, statistics
      USE sgsm, ONLY : diffuse
      USE srfc, ONLY : surface
      USE thrm, ONLY : thermo
      USE mcrp, ONLY : micro
      USE prss, ONLY : poisson
      USE advf, ONLY : fadvect, newdroplet
      USE advl, ONLY : ladvect
      USE forc, ONLY : forcings
      USE util, ONLY : maskactiv !Juha: Included for SALSA

      USE mo_salsa_driver, ONLY : run_SALSA

      USE constrain_SALSA, ONLY : tend_constrain, SALSA_diagnostics

      LOGICAL, INTENT (out)      :: cflflg
      REAL(KIND=8), INTENT (out) :: cflmax

      LOGICAL :: zactmask(nzp,nxp,nyp)
      REAL    :: zwp(nzp,nxp,nyp)  !! FOR SINGLE-COLUMN RUNS

      INTEGER :: n4
      
      CALL set_LES_runtime(time)

      zwp = 0.5  ! single column run vertical velocity

      cflflg = .FALSE.

      ! Reset ALL tendencies here.
      !----------------------------------------------------------------
      ! "Scalar" timestep
      CALL tend0(.FALSE.)

      ! Put the newly activated to zero
      IF (level >= 4) THEN
         a_vactd = 0.
         a_nactd = 0.
      END IF

      CALL surface()

      CALL diffuse

      CALL sponge(0)

      IF (level >= 1) THEN

         CALL thermo(level)

         CALL forcings(time,strtim)

         IF (level >= 4) THEN

            n4 = spec%getNSpec() ! Aerosol components + water

            CALL tend_constrain(n4)
            CALL update_sclrs
            CALL tend0(.TRUE.)

            IF ( nxp == 5 .AND. nyp == 5 ) THEN
               ! 1D -runs
               CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,zwp,a_dn,  &
                              a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                              a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                              a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                              a_nicep,   a_nicet,   a_micep,   a_micet,    &
                              a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                              a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                              dtlt, time, level, .FALSE.)
            ELSE
               !! for 2D or 3D runs
               CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,a_wp,a_dn,  &
                              a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                              a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                              a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                              a_nicep,   a_nicet,   a_micep,   a_micet,    &
                              a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                              a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                              dtlt, time, level, .FALSE.)
            END IF !nxp==5 and nyp == 5
            CALL tend_constrain(n4)
         END IF

      END IF ! level

      CALL update_sclrs

      !-------------------------------------------
      ! "Deposition" timestep
      ! -- Reset only scalar tendencies
      CALL tend0(.TRUE.)

      ! Dont perform sedimentation or level 3 autoconversion during spinup
      CALL micro(level)

      IF (level >= 4) CALL tend_constrain(n4)
      CALL update_sclrs

      !-------------------------------------------
      ! "Advection" timestep
      ! -- Reset only scalar tendencies
      CALL tend0(.TRUE.)

      ! Mask for cloud base activation
      IF (level >= 4) CALL maskactiv(zactmask,nxp,nyp,nzp,2,a_rh,rc=a_rc,w=a_wp)
      ! Get tendencies from cloud base activation
      IF (level >= 4) CALL newdroplet(zactmask)

      CALL fadvect

      IF (level >= 4)  &
         CALL tend_constrain(n4)

      CALL update_sclrs

      CALL thermo(level)

      IF (level >= 4)  THEN
         CALL SALSA_diagnostics
         CALL thermo(level)
      END IF

      CALL corlos

      CALL ladvect

      CALL buoyancy

      CALL sponge(1)

      CALL poisson

      CALL cfl (cflflg, cflmax)

      CALL thermo(level)

      IF (level >= 4)  THEN
         CALL SALSA_diagnostics
         call thermo(level)
      ENDIF

      IF (sflg) THEN
         CALL statistics (time+dtl)
      END IF

   END SUBROUTINE t_step
   !
   !----------------------------------------------------------------------
   ! Subroutine tend0: sets all tendency arrays to zero
   !
   SUBROUTINE tend0(sclonly)

      USE grid, ONLY : a_ut, a_vt, a_wt, nscl, a_st, newsclr

      LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

      INTEGER :: n

      IF( .NOT. sclonly) THEN
         a_ut = 0.; a_vt = 0.; a_wt = 0.
      END IF
      DO n = 1, nscl
         CALL newsclr(n)
         a_st = 0.
      END DO

   END SUBROUTINE tend0

   !
   !----------------------------------------------------------------------
   ! Subroutine cfl: Driver for calling CFL computation subroutine
   !
   SUBROUTINE cfl(cflflg,cflmax)

      USE grid, ONLY : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dtlt
      USE stat, ONLY : fill_scalar

      LOGICAL, INTENT(out) :: cflflg
      REAL(KIND=8), INTENT (out)   :: cflmax
      REAL, PARAMETER :: cflnum = 0.95

      cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dtlt)

      cflflg = (cflmax > cflnum)
      IF (cflflg) PRINT *, 'Warning CFL Violation :', cflmax
      CALL fill_scalar(1,REAL(cflmax))

   END SUBROUTINE cfl
   !
   !----------------------------------------------------------------------
   ! Subroutine cfll: Checks CFL criteria, brings down the model if the
   ! maximum thershold is exceeded
   !
   REAL(KIND=8) FUNCTION cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

      INTEGER, INTENT (in) :: n1, n2, n3
      REAL, DIMENSION (n1,n2,n3), INTENT (in) :: u, v, w
      REAL, INTENT (in)    :: dxi,dyi,dzt(n1),dtlt

      INTEGER :: i, j, k
      cfll = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               cfll = max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                      abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
            END DO
         END DO
      END DO

   END FUNCTION cfll
   !
   !----------------------------------------------------------------------
   ! Subroutine update_sclrs:  Updates scalars by applying tendency and
   ! boundary conditions
   !
   SUBROUTINE update_sclrs

      USE grid, ONLY : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
                       dtlt, newsclr, isgstyp
      USE sgsm, ONLY : tkeinit
      USE util, ONLY : sclrset

      INTEGER :: n

      DO n = 1, nscl
         CALL newsclr(n)
         CALL update(nzp,nxp,nyp,a_sp,a_st,dtlt)
         CALL sclrset('mixd',nzp,nxp,nyp,a_sp,dzt)
      END DO

      IF (isgstyp == 2) THEN
         CALL tkeinit(nxyzp,a_qp)
      END IF

   END SUBROUTINE update_sclrs
   !
   ! ----------------------------------------------------------------------
   ! Subroutine update:
   !
   SUBROUTINE update(n1,n2,n3,a,fa,dt)

      INTEGER, INTENT(in)   :: n1, n2, n3
      REAL, INTENT (in)     :: fa(n1,n2,n3),dt
      REAL, INTENT (inout)  :: a(n1,n2,n3)
      INTEGER :: i, j, k

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               a(k,i,j) = a(k,i,j) + fa(k,i,j)*dt
            END DO
         END DO
      END DO

   END SUBROUTINE update
   !
   ! ----------------------------------------------------------------------
   ! Subroutine buoyancy:
   !
   SUBROUTINE buoyancy
     
     USE grid, ONLY : a_uc, a_vc, a_wc, a_wt, a_rv, a_rc, a_theta, &
          a_rp, a_srp, a_ri, a_srs, nxp, nyp, nzp, dzm, th00, level, pi1
     USE stat, ONLY : sflg, comp_tke
     USE util, ONLY : ae1mm
     USE thrm, ONLY : update_pi1
     
     REAL :: awtbar(nzp), a_tmp1(nzp,nxp,nyp), rv(nzp,nxp,nyp), rc(nzp,nxp,nyp)
     
     IF (level < 4) THEN
        rv = a_rv ! Water vapor
        rc = a_rp - a_rv ! Total condensate (cloud + precipitation)
     ELSE IF (level >= 4) THEN
        rv = a_rp ! Water vapor
        rc = a_rc + a_srp + a_ri + a_srs ! Total condensed water (aerosol+cloud+precipitation+ice+snow)
     END IF
     call boyanc(nzp,nxp,nyp,a_wt,a_theta,rv,th00,a_tmp1,rc)
     
     CALL ae1mm(nzp,nxp,nyp,a_wt,awtbar)
     CALL update_pi1(nzp,awtbar,pi1)
     
     IF (sflg)  CALL comp_tke(nzp,nxp,nyp,dzm,th00,a_uc,a_vc,a_wc,a_tmp1)
     
   END SUBROUTINE buoyancy
   !
   ! ----------------------------------------------------------------------
   ! Subroutine boyanc:
   !
   SUBROUTINE boyanc(n1,n2,n3,wt,th,rv,th00,scr,rc)

      USE defs, ONLY : g, ep2

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: th00,th(n1,n2,n3),  &
                             rv(n1,n2,n3)  ! water vapor
                                      
      REAL, INTENT(in)    :: rc(n1,n2,n3)  ! Total condensed water (aerosol, cloud, rain, ice and snow) mixing ratio

      REAL, INTENT(inout) :: wt(n1,n2,n3)
      REAL, INTENT(out)   :: scr(n1,n2,n3)

      INTEGER :: k, i, j
      REAL    :: gover2

      gover2 = 0.5*g

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rv(k,i,j))-th00)/th00-rc(k,i,j))
          end do

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)+scr(k,i,j)+scr(k+1,i,j)
          end do
       end do
    end do

   END SUBROUTINE boyanc
   !
   ! ----------------------------------------------------------------------
   ! Subroutine corlos:  This is the coriolis driver, its purpose is to
   ! from the coriolis accelerations for u and v and add them into the
   ! accumulated tendency arrays of ut and vt.
   !
   SUBROUTINE corlos

      USE defs, ONLY : omega
      USE grid, ONLY : a_uc, a_vc, a_ut, a_vt, nxp, nyp, nzp, u0, v0, cntlat

      LOGICAL, SAVE :: initialized = .FALSE.
      REAL, SAVE    :: fcor

      INTEGER :: i, j, k

      IF (corflg) THEN
         IF (.NOT. initialized) fcor = 2.*omega*sin(cntlat*0.01745329)
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = 2, nzp
                  a_ut(k,i,j) = a_ut(k,i,j) - fcor*(v0(k)-0.25*                   &
                                (a_vc(k,i,j)+a_vc(k,i+1,j)+a_vc(k,i,j-1)+a_vc(k,i+1,j-1)))
                  a_vt(k,i,j) = a_vt(k,i,j) + fcor*(u0(k)-0.25*                   &
                                (a_uc(k,i,j)+a_uc(k,i-1,j)+a_uc(k,i,j+1)+a_uc(k,i-1,j+1)))
               END DO
            END DO
         END DO
         initialized = .TRUE.
      END IF

   END SUBROUTINE corlos
   !
   ! ----------------------------------------------------------------------
   ! Subroutine sponge: does the rayleigh friction for the momentum terms,
   ! and newtonian damping of thermal term the damping is accumulated with the
   ! other tendencies
   !
   SUBROUTINE sponge (isponge)

      USE grid, ONLY : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
                       nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00

      INTEGER, INTENT (in) :: isponge

      INTEGER :: i, j, k, kk

      IF (maxval(spng_tfct) > epsilon(1.) .AND. nfpt > 1) THEN
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = nzp-nfpt, nzp-1
                  kk = k+1-(nzp-nfpt)
                  IF (isponge == 0) THEN
                     a_tt(k,i,j) = a_tt(k,i,j) - spng_tfct(kk)*                   &
                                   (a_tp(k,i,j)-th0(k)+th00)
                  ELSE
                     a_ut(k,i,j) = a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-u0(k))
                     a_vt(k,i,j) = a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-v0(k))
                     a_wt(k,i,j) = a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                  END IF
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE sponge


   ! NOTE: SALSA_DIAGNOSTICS AND OTHER SIMILAR SUBROUTINES ARE NOW FOUND IN CONSTRAIN_SALSA.F90



END MODULE step
