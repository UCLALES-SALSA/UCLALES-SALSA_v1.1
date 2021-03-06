
!****************************************************************
!*                                                              *
!*   MODULE MO_SALSA_DYNAMICS                               *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************

MODULE mo_salsa_dynamics
      USE classSection, ONLY : Section
      USE mo_salsa_types, ONLY : allSALSA
   IMPLICIT NONE


CONTAINS

   ! this calculated for empty bins too!!!
   ! fxm: test well, esp. self-coagulation (but other bits too!)
   ! AL_note: Diagnostic variables of cond and nucl mass
   !********************************************************************
   !
   ! Subroutine COAGULATION(kproma,kbdim,klev, &
   !       pnaero,pvols,pdwet, &
   !       pcore, ptstep)
   !
   !********************************************************************
   !
   ! Purpose:
   ! --------
   ! Calculates particle loss and change in size distribution
   !  due to (Brownian) coagulation
   !
   !
   ! Method:
   ! -------
   ! Semi-implicit, non-iterative method:
   !  Volume concentrations of the smaller colliding particles
   !  added to the bin of the larger colliding particles.
   !  Start from first bin and use the updated number and volume
   !  for calculation of following bins. NB! Our bin numbering
   !  does not follow particle size in regime 2.
   !
   !Schematic for bin numbers in different regimes:
   !             1                            2
   !    +-------------------------------------------+
   !  a | 1 | 2 | 3 || 4 | 5 | 6 | 7 |  8 |  9 | 10||
   !  b |           ||11 |12 |13 |14 | 15 | 16 | 17||
   !    +-------------------------------------------+
   !
   ! Exact coagulation coefficients for each pressure level
   !  are calculated in subroutine SET_COAGC (in mo_salsa_init)
   !  which is called once at the beginning of the simulation
   !  from model driver. In subroutine COAGULATION, these exact
   !  coefficients are scaled according to current particle wet size
   !  (linear scaling).
   !
   ! Juha: Now also considers coagulation between hydrometeors,
   !       and hydrometeors and aerosols.
   !
   !       Since the bins are organized in terms of the dry size of
   !       of the condensation nucleus, while coagulation kernell is
   !       calculated with the actual hydrometeor size, some assumptions
   !       are laid out:
   !                 1. Cloud droplets from each size bin are lost by
   !                    coagulation with other cloud droplets that have
   !                    larger condensation nucleus.
   !
   !                 2. Cloud droplets from each size bin are lost by
   !                    coagulation with all drizzle bins, regardless of
   !                    the nucleus size in the latter (collection of cloud
   !                    droplets by rain).
   !
   !                 3. Coagulation between drizzle bins acts like 1.
   !
   !       ISSUES:
   !           Process selection should be made smarter - now just lots of ifs
   !           inside loops. Bad.
   !
   !
   ! Interface:
   ! ----------
   ! Called from main aerosol model
   !
   !
   ! Coded by:
   ! ---------
   ! Hannele Korhonen (FMI) 2005
   ! Harri Kokkola (FMI) 2006
   ! Tommi Bergman (FMI) 2012
   ! Matti Niskanen(FMI) 2012
   ! Anton Laakso  (FMI) 2013
   ! Juha Tonttila (FMI) 2014
   !
   !---------------------------------------------------------------------


   SUBROUTINE coagulation(kproma,kbdim,klev,   &
                          ptstep,ptemp,ppres   )

     USE mo_salsa_types, ONLY : aero, cloud, precp, ice, snow
      USE mo_submctl, ONLY: ntotal,nbins,ncld,nprc,nice,nsnw,spec,   &
                            lscgaa, lscgcc, lscgpp, lscgii, lscgss,  &
                            lscgca, lscgpa, lscgia, lscgsa,          &
                            lscgpc, lscgic, lscgsc,                  &
                            lscgip, lscgsp,                          &
                            lscgsi

      USE mo_salsa_coagulation_kernels

      USE mo_salsa_coagulation_processes
      USE classSection, ONLY : Section
      USE mo_salsa_types, ONLY : allSALSA

      IMPLICIT NONE


      !-- Input and output variables -------------
      INTEGER, INTENT(IN) ::        &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical klev

      REAL, INTENT(IN) ::         &
         ptstep,                    & ! time step [s]
         ptemp(kbdim,klev),         &
         ppres(kbdim,klev)
      !-- Local variables ------------------------

      INTEGER :: nspec

      LOGICAL :: any_aero, any_cloud, any_precp, any_ice, any_snow
      
      REAL :: zccaa(kbdim,klev,nbins,nbins),    & ! updated coagulation coefficients [m3/s]
              zcccc(kbdim,klev,ncld,ncld),      & ! - '' - for collision-coalescence between cloud droplets [m3/s]
              zccca(kbdim,klev,nbins,ncld),     & ! - '' - for cloud collection of aerosols [m3/s]
              zccpc(kbdim,klev,ncld,nprc),      & ! - '' - for collection of cloud droplets by precip [m3/s]
              zccpa(kbdim,klev,nbins,nprc),     & ! - '' - for collection of aerosols by precip
              zccpp(kbdim,klev,nprc,nprc),      & ! - '' - for collision-coalescence between precip particles 
              zccia(kbdim,klev,nbins,nice),     & ! - '' - for collection of aerosols by ice 
              zccic(kbdim,klev,ncld,nice),      & ! - '' - for collection of cloud particles droplets by ice 
              zccii(kbdim,klev,nice,nice),      & ! - '' - for aggregation between ice 
              zccip(kbdim,klev,nprc,nice),      & ! - '' - for collection of precip by ice
              zccsa(kbdim,klev,nbins,nsnw),     & ! - '' - for collection of aerosols by snow 
              zccsc(kbdim,klev,ncld,nsnw),      & ! - '' - for collection of cloud droples by snow 
              zccsi(kbdim,klev,nice,nsnw),      & ! - '' - for collection of ice by snow 
              zccsp(kbdim,klev,nprc,nsnw),      & ! - '' - for collection of precip by snow 
              zccss(kbdim,klev,nsnw,nsnw)         ! - '' - for aggregation between snow

      INTEGER :: ii,jj,kk

      !-----------------------------------------------------------------------------
      !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
      !      CoagSink ~ Dp in continuum regime, thus we calculate
      !      'effective' number concentration of coarse particles

      zccaa(:,:,:,:) = 0.; zcccc(:,:,:,:) = 0.; zccca(:,:,:,:) = 0.; zccpc(:,:,:,:) = 0.
      zccpa(:,:,:,:) = 0.; zccpp(:,:,:,:) = 0.; zccia(:,:,:,:) = 0.; zccic(:,:,:,:) = 0.
      zccii(:,:,:,:) = 0.; zccip(:,:,:,:) = 0.; zccsa(:,:,:,:) = 0.; zccsc(:,:,:,:) = 0.
      zccsi(:,:,:,:) = 0.; zccsp(:,:,:,:) = 0.; zccss(:,:,:,:) = 0.
            
      !-- 2) Updating coagulation coefficients -------------------------------------
      nspec = spec%getNSpec()


      CALL update_coagulation_kernels(kbdim,klev,ppres,ptemp,    &
                                      zccaa, zcccc, zccca, zccpc, zccpa,  &
                                      zccpp, zccia, zccic, zccii, zccip,  &
                                      zccsa, zccsc, zccsi, zccsp, zccss)

      any_aero = ANY( aero(:,:,:)%numc > aero(:,:,:)%nlim ) .AND. &
                 ANY( [lscgaa,lscgca,lscgpa,lscgia,lscgsa] )
      any_cloud = ANY( cloud(:,:,:)%numc > cloud(:,:,:)%nlim ) .AND. &
                  ANY( [lscgcc,lscgca,lscgpc,lscgic,lscgsc] ) 
      any_precp = ANY( precp(:,:,:)%numc > precp(:,:,:)%nlim ) .AND. &
                  ANY( [lscgpp,lscgpa,lscgpc,lscgip,lscgsp])
      any_ice = ANY( ice(:,:,:)%numc > ice(:,:,:)%nlim ) .AND. &
                ANY( [lscgii,lscgia,lscgic,lscgip,lscgsi] )
      any_snow = ANY( snow(:,:,:)%numc > snow(:,:,:)%nlim ) .AND. &
                 ANY( [lscgss,lscgsa,lscgsc,lscgsp,lscgsi] )
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !-- 3) New particle and volume concentrations after coagulation -------------
      IF (any_precp) &
           CALL coag_precp(kbdim,klev,nspec,ptstep,zccpp,zccpa,zccpc,zccip,zccsp)

      IF (any_aero) &
           CALL coag_aero(kbdim,klev,nspec,ptstep,zccaa,zccca,zccpa,zccia,zccsa)

      IF (any_cloud) &
           CALL coag_cloud(kbdim,klev,nspec,ptstep,zcccc,zccca,zccpc,zccic,zccsc)

   END SUBROUTINE coagulation


   ! fxm: calculated for empty bins too
   ! fxm: same diffusion coefficients and mean free paths used for sulphuric acid
   !      and organic vapours (average values? 'real' values for each?)
   !********************************************************************
   !
   ! Subroutine CONDENSATION(kproma, kbdim,  klev,        &
   !                         pnaero, pvols,  pdwet, plwc, &
   !                         pcsa,   pcocnv, pcocsv,      &
   !                         ptemp,  ppres,  ptstep)
   !
   !********************************************************************
   !
   ! Purpose:
   ! --------
   ! Calculates the increase in particle volume and
   !  decrease in gas phase concentrations due to condensation
   !  of sulphuric acid and two organic compounds (non-volatile
   !  and semivolatile)
   !
   !
   ! Method:
   ! -------
   ! Regime 3 particles only act as a sink for condensing vapours
   !  while their size and composition does not change.
   ! Exception: Soluble fraction of regime 3c particles can change
   !  and thus they can be moved to regime 3b
   !
   ! New gas and aerosol phase concentrations calculated according
   !  to Jacobson (1997): Numerical techniques to solve
   !  condensational and dissolutional growth equations
   !  when growth is coupled to reversible reactions,
   !  Aerosol Sci. Tech., 27, pp 491-498.
   !
   ! fxm: one should really couple with vapour production and loss terms as well
   !      should nucleation be coupled here as well????
   !
   ! Juha: Now does the condensation of water vapour on hydrometeors as well,
   !       + the condensation of semivolatile aerosol species on hydromets.
   !       Modified for the new aerosol datatype. LWC is obtained from %volc(8)
   !
   !
   ! Interface:
   ! ----------
   ! Called from main aerosol model
   !
   !
   ! Coded by:
   ! ---------
   ! Hannele Korhonen (FMI) 2005
   ! Harri Kokkola (FMI) 2006
   ! Juha Tonttila (FMI) 2014
   !
   !---------------------------------------------------------------
   !
   ! Following parameterization has been used:
   ! ------------------------------------------
   !
   ! Molecular diffusion coefficient of condensing vapour [m2/s]
   !  (Reid et al. (1987): Properties of gases and liquids,
   !   McGraw-Hill, New York.)
   !
   ! D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / &
   !  {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
   !
   ! M_air = 28.965 : molar mass of air [g/mol]
   ! d_air = 19.70  : diffusion volume of air
   ! M_h2so4 = 98.08  : molar mass of h2so4 [g/mol]
   ! d_h2so4 = 51.96  : diffusion volume of h2so4
   !
   !---------------------------------------------------------------

   SUBROUTINE condensation(kproma,  kbdim,  klev,   krow,       &
                            level,  paero,   pcloud, pprecp,            &
                            pice,    psnow,                     &
                            pcsa, pcocnv,  pcocsv, pchno3, 	&
			    pcnh3, prv,prs, prsi,ptemp,  ppres,  	&
			    ptstep,  ppbl)

     USE mo_salsa_nucleation
     USE mo_submctl, ONLY :      &
          ntotal,                    &
          lscndgas,                  & 
          lscndh2oae, lscndh2ocl, lscndh2oic, & ! Condensation to aerosols, clouds and ice particles
          nsnucl,                  &                     ! nucleation
          fn2b,nbins ,                     &
          ncld,nprc,                  &
          nice,nsnw
     IMPLICIT NONE
     
     !-- Input and output variables ----------
     INTEGER, INTENT(IN) ::      &
          kproma,                    & ! number of horiz. grid kproma
          kbdim,                     & ! dimension for arrays
          klev,                      & ! number of vertical klev
          krow
     INTEGER, INTENT(in) :: level
     
      REAL, INTENT(IN) ::         &
           ptemp(kbdim,klev),         & ! ambient temperature [K]
           ppres(kbdim,klev),         & ! ambient pressure [Pa]
           ptstep,                    & ! timestep [s]
           prs(kbdim,klev),           & ! Water vapor saturation mixing ratio
           prsi(kbdim,klev)              ! Saturation mixing ratio    [kg/m3]
      
      INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level
      
      REAL, INTENT(INOUT) ::     &
           prv(kbdim,klev),          & ! Water vapor mixing ratio [kg/kg]
           pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
           pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
           pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
           pchno3(kbdim,klev),       & ! nitric acid concentration [#/m3]
           pcnh3(kbdim,klev)          ! ammonia concentration [#/m3]
      
      TYPE(section), INTENT(inout) :: &
           pcloud(kbdim,klev,ncld),     & ! Hydrometeor properties
           paero(kbdim,klev,fn2b),      & ! Aerosol properties
           pprecp(kbdim,klev,nprc),     & ! rain properties
           pice(kbdim,klev,nice),       & ! ice properties
           psnow(kbdim,klev,nsnw)        ! snowing properties
      
      REAL :: zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
           zxsa(kbdim,klev),           & ! ratio of sulphuric acid and organic vapor in 3nm particles
           zxocnv(kbdim,klev),         &
           zrh(kbdim,klev),zsubtime,ztstep2
      INTEGER :: ll
      
      zxocnv = 0.
      zxsa = 0.
      zj3n3 = 0.
      zrh(1:kbdim,:) = prv(1:kbdim,:)/prs(1:kbdim,:)
      
      !------------------------------------------------------------------------------

      ! Nucleation
      IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,  &
           ptemp,  zrh,    ppres,  &
           pcsa,   pcocnv, ptstep, zj3n3,  &
           zxsa,   zxocnv, ppbl            )
      
      ! Condensation of H2SO4 and organic vapors
      IF (lscndgas) CALL condgas(kproma,  kbdim,  klev,    krow,      &
           pcsa, pcocnv, pcocsv,     &
           zxsa, ptemp,  ppres, ptstep )
      
      ! Condensation of water vapour, HNO3 and NH3
      
      ztstep2=1e-4
      zsubtime = 0.
      DO ll = 1,15
         
         zsubtime = zsubtime + ztstep2
         if(zsubtime>=ptstep) cycle
         IF (lscndh2ocl .OR. lscndh2oae .OR. lscndh2oic) &
              CALL partitioning(kproma,kbdim,klev,krow,ppres,ptemp,paero,pcloud,    &
              pprecp, pchno3,pcnh3,prv,prs,ztstep2)
         ztstep2=min(0.2,ztstep2*10,ptstep- zsubtime)
         
      ENDDO
    END SUBROUTINE condensation

!
! ----------------------------------------------------------------------------------------------------------
!
    SUBROUTINE partitioning(kproma,kbdim,klev,krow,ppres,ptemp,aero,cloud,    &
         precp,pghno3,pgnh3,prv,prs,ptstep)
      
      USE classSection
      USE mo_submctl, ONLY : spec,           &
           nbins, ncld, nprc,   &
           in1a, fn2b,          &
           rhowa, mwa, mair,    &
           surfw0, mvno, mvnh,  &
           mvwa, boltz, rg,     &
           massacc,      &
           rhono, mno,          &
           rhonh, mnh,          &
           rhosu, msu,          &
           avog, pi,pi6,            &
           pstand,              &
           nlim, prlim, d_sa,ntotal, &
           alv
      
      
      USE aerosol_thermodynamics, ONLY : inorganic_pdfite
      
      IMPLICIT NONE
      
      ! Equation numbers refer to those in Jacobson: Fundamentals of Atmospheric Modelling, Second Edition (2005)
      
      INTEGER, INTENT(in)  :: kproma,kbdim,klev,krow
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev)
      REAL, INTENT(in) :: prs(kbdim,klev)
      
      
      TYPE(section), INTENT(inout) :: &
           cloud(kbdim,klev,ncld),     & ! Hydrometeor properties
           aero(kbdim,klev,fn2b),      & ! Aerosol properties
           precp(kbdim,klev,nprc)
      REAL, INTENT(inout) :: prv(kbdim,klev),      &
           pghno3(kbdim,klev),   &
           pgnh3(kbdim,klev)
      REAL :: zph(kbdim,klev,nbins+ncld+nprc)
      
      REAL :: zkelwaae(nbins),  zkelwacd(ncld),  zkelwapd(nprc)      ! Kelvin effects for water
      REAL :: zkelno3ae(nbins), zkelno3cd(ncld), zkelno3pd(nprc)     ! Kelvin effects for HNO3
      REAL :: zkelnh3ae(nbins), zkelnh3cd(ncld), zkelnh3pd(nprc)     ! Kelvin effects for NH3
      
      REAL :: zcwacae(nbins),   zcwanae(nbins),   & ! Current, intermediate and new water in aerosols
           zcno3cae(nbins),  zcno3nae(nbins),  & !  -  ''  - HNO3
           zcnh3cae(nbins),  zcnh3nae(nbins),  & !  -  ''  - NH3
           
           zcwaccd(ncld),    zcwancd(ncld),   & ! -  ''  - water in cloud droplets
           zcno3ccd(ncld),   zcno3ncd(ncld),  & ! -  ''  - HNO3 in cloud droplets
           zcnh3ccd(ncld),   zcnh3ncd(ncld),  & ! -  ''  - NH3 
           
           zcwacpd(nprc),    zcwanpd(nprc),   & ! -  ''  - water in precipitation
           zcno3cpd(nprc),   zcno3npd(nprc),  & ! -  ''  - HNO3 in precipitation
           zcnh3cpd(nprc),   zcnh3npd(nprc)     ! -  ''  - NH3
      
      REAL :: zcwatot,zcnh3tot,zcno3tot
      
      REAL :: zcwac, zcwan                                         ! Current, intermediate and new water gas concentration
      REAL :: zcno3c, zcno3n                                       ! -  ''  - HNO3
      REAL :: zcnh3c, zcnh3n                                       ! -  ''  - NH3
      
      REAL :: zcgnh3eqae(nbins), zcgno3eqae(nbins), zcgwaeqae(nbins), & ! Equilibrium gas concentrations 
           zcgnh3eqcd(ncld), zcgno3eqcd(ncld),   zcgwaeqcd(ncld), &
           zcgnh3eqpd(nprc), zcgno3eqpd(nprc),   zcgwaeqpd(nprc)
      
      REAL :: zmtwaae(nbins),  zmtwacd(ncld),  zmtwapd(nprc)  ! Mass transfer coefficients for H2O
      REAL :: zmtno3ae(nbins), zmtno3cd(ncld), zmtno3pd(nprc) ! Mass transfer coefficients for HNO3
      REAL :: zmtnh3ae(nbins), zmtnh3cd(ncld), zmtnh3pd(nprc) ! Mass transfer coefficients for NH3
      
      REAL :: zsathno3ae(nbins), zsathno3cd(ncld), zsathno3pd(nprc)
      REAL :: zsatnh3ae(nbins), zsatnh3cd(ncld), zsatnh3pd(nprc)
      
      REAL :: zbeta  ! transition correction factor
      REAL :: zdfvap ! Diffusion coefficient for vapors
      
      REAL :: zaw(nbins) ! water activity
      
      REAL :: zHp_ae(nbins,3),    & ! H' (Eq (17.99)) for aerosol, 1 = NO3-, 2=NH4+, 3=H2O
           zHp_cd(ncld,3),     & ! H' (Eq (17.99)) for clouds
           zHp_pd(nprc,3),     & ! H' (Eq (17.99)) for precipitation
           zexpterm_ae(nbins), & ! exponent term in (17.104) for aerosol bins
           zexpterm_cd(ncld), & ! exponent term in (17.104) for cloud bins
           zexpterm_pd(nprc)    ! exponent term in (17.104) for precipitation bins
      
      
      REAL :: zrhoair,            & ! air density [kg m-3]
           zthcond,            & ! thermal conductivity of air
           zdfh2o,zdfno3,zdfnh3                ! diffusion coefficient of water
      
      REAL :: adt ! timestep
      REAL :: adtc2(2,nbins)
      REAL :: telp,ttot ! Elapsed time
      REAL :: zsum1, zsum2, zhlp3 ! temporary variables 
      REAL :: zDeff_ae(nbins), zDeff_cd(ncld), zDeff_pd(nprc), & ! effective diffusion coefficient
           zDp_ae(nbins),   zDp_cd(ncld),   zDp_pd(nprc), &
           dwet_ae(nbins),   dwet_cd(ncld),   dwet_pd(nprc)
      INTEGER :: nstr,ph_switch= 0
      INTEGER :: ii,jj,kk,cc,tt,new,water, nspec

      REAL :: zbetaae(fn2b),zbetacd(ncld),zbetapd(nprc),fmax,fmin
      REAL :: zknae(fn2b),zkncd(ncld),zknpd(nprc),zknud(fn2b),h2oactae(fn2b),h2oactcd(ncld),h2oactpd(nprc)
      REAL :: zmfp_h2o, zmfp_no3, zmfp_nh3
      REAL :: zmfp,rhlimx, alpha_h2o, alpha_no3, alpha_nh3
      INTEGER :: iso4,ioc,ibc,idu,iss,iwa,ino,inh

      iso4 = spec%getIndex("SO4",notFoundValue = 0)
      ioc = spec%getIndex("OC",notFoundValue = 0)
      ibc = spec%getIndex("BC",notFoundValue = 0)
      idu = spec%getIndex("DU",notFoundValue = 0)
      iss = spec%getIndex("SS",notFoundValue = 0)
      ino = spec%getIndex("NO",notFoundValue = 0)
      inh = spec%getIndex("NH",notFoundValue = 0)
      iwa = spec%getIndex("H2O")
      nspec = spec%getNSpec()

      nstr = 1
      ! initialize
      adt = ptstep
    
      DO jj = 1,klev

         DO ii = 1,kproma
            
            ! density of air
            zrhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))
            zkelno3pd = 1.
            zkelno3cd = 1.
            zkelno3ae = 1.
            zkelnh3pd = 1.
            zkelnh3cd = 1.
            zkelnh3ae = 1.
            zsathno3ae = 1.
            zsathno3cd = 1.
            zsathno3pd = 1.
            zsatnh3ae = 1.
            zsatnh3cd = 1.
            zsatnh3pd = 1.
            zexpterm_ae = 0.
            zexpterm_cd = 0.
            zexpterm_pd = 0.
            
            ! Kelvin effects
            ! Wet diameter
            
            DO cc = 1, nbins
               if (aero(ii,jj,cc)%numc >nlim) then
                  dwet_ae(cc) = ( SUM(aero(ii,jj,cc)%volc(1:nspec))/aero(ii,jj,cc)%numc/pi6 )**(1./3.)
                  zkelwaae(cc) = min(2.,exp( 4.*surfw0*mvwa /  &
                       (boltz*ptemp(ii,jj)*dwet_ae(cc)) ) )
                  zkelno3ae(cc) = min(2.,exp( 4.*surfw0*mvno /  &
                       (boltz*ptemp(ii,jj)*dwet_ae(cc)) ) )
                  zkelnh3ae(cc) = min(2.,exp( 4.*surfw0*mvnh /  &
                       (boltz*ptemp(ii,jj)*dwet_ae(cc)) ) )
               endif
            enddo
                        
            ! Wet diameter
            DO cc = 1, ncld
               if (cloud(ii,jj,cc)%numc >nlim) then 
                  dwet_cd(cc) = ( SUM(cloud(ii,jj,cc)%volc(1:nspec))/cloud(ii,jj,cc)%numc/pi6 )**(1./3.)
                  zkelwacd(cc) = min(2.,exp( 4.*surfw0*mvwa /  & 
                       (boltz*ptemp(ii,jj)*dwet_cd(cc)) ) )
                  zkelno3cd(cc) = min(2.,exp( 4.*surfw0*mvno /  & 
                       (boltz*ptemp(ii,jj)*dwet_cd(cc)) ) )
                  zkelnh3cd(cc) = min(2.,exp( 4.*surfw0*mvnh /  &
                       (boltz*ptemp(ii,jj)*dwet_cd(cc)) ) )
               endif
            enddo
            
            ! Wet diameter
            DO cc = 1, nprc
               if (precp(ii,jj,cc)%numc >prlim) then 
                  dwet_pd(cc) = ( SUM(precp(ii,jj,cc)%volc(1:nspec))/precp(ii,jj,cc)%numc/pi6 )**(1./3.)
                  zkelwapd(cc) = min(2.,exp( 4.*surfw0*mvwa /  &
                       (boltz*ptemp(ii,jj)*dwet_pd(cc)) ) )
                  zkelno3pd(cc) = min(2.,exp( 4.*surfw0*mvno /  &
                       (boltz*ptemp(ii,jj)*dwet_pd(cc)) ) )
                  zkelnh3pd(cc) = min(2.,exp( 4.*surfw0*mvnh /  &
                       (boltz*ptemp(ii,jj)*dwet_pd(cc)) ) )
               endif
            enddo
            ! Current gas concentrations
            zcno3c = pghno3(ii,jj)/avog
            zcnh3c = pgnh3(ii,jj)/avog
            zcwac  = prv(ii,jj)/mwa*zrhoair
            
            ! Current particle concentrations
            if (ino>0) then
               zcno3cae(1:nbins) = aero(ii,jj,1:nbins)%volc(ino)*rhono/mno
               zcno3ccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(ino)*rhono/mno
               zcno3cpd(1:nprc) = precp(ii,jj,1:nprc)%volc(ino)*rhono/mno
            endif
            if (inh>0) then
               zcnh3cae(1:nbins) = aero(ii,jj,1:nbins)%volc(inh)*rhonh/mnh
               zcnh3ccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(inh)*rhonh/mnh
               zcnh3cpd(1:nprc) = precp(ii,jj,1:nprc)%volc(inh)*rhonh/mnh
            endif
            if (iwa>0) then
               zcwacae(1:nbins) = aero(ii,jj,1:nbins)%volc(iwa)*rhowa/mwa
               zcwaccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(iwa)*rhowa/mwa
               zcwacpd(1:nprc) = precp(ii,jj,1:nprc)%volc(iwa)*rhowa/mwa
            endif
            zcwatot=zcwac + sum(zcwacae(1:nbins)) + sum(zcwaccd(1:ncld)) + sum(zcwacpd(1:nprc)) ! total amount of water
            zcnh3tot=zcnh3c+sum(zcnh3cae(1:nbins)) + sum(zcnh3ccd(1:ncld)) + sum(zcnh3cpd(1:nprc))
            zcno3tot=zcno3c+sum(zcno3cae(1:nbins)) + sum(zcno3ccd(1:ncld)) + sum(zcno3cpd(1:nprc))

            ! Mass transfer coefficients 
            zmtwaae = 0.; zmtwaae = 0.
            zmtno3ae = 0.; zmtnh3ae = 0.
            zmtno3cd = 0.; zmtnh3cd = 0.
            zmtno3pd = 0.; zmtnh3pd = 0.
            
            ! Get the equilibrium concentrations before condensation of water
            ! aerosols
            zcgno3eqae=0.
            zcgnh3eqae=0.
            zcgwaeqae=0.
            zcgno3eqcd=0.
            zcgnh3eqcd=0.
            zcgwaeqcd=0.
            zcgno3eqpd=0.
            zcgnh3eqpd=0.
            zcgwaeqpd=0.
            h2oactae=0.
            h2oactcd=0.
            h2oactpd=0.
            ! aerosol droplets
            CALL thermoequil(aero(ii,jj,:),nbins,nlim,ptemp(ii,jj),ppres(ii,jj), zcgno3eqae, zcgnh3eqae, 	&
                 zcgwaeqae,h2oactae,zpH(ii,jj,1:nbins),ph_switch,zkelwaae,prv(ii,jj)/prs(ii,jj),ptstep)
            ! cloud droplets
            
            CALL thermoequil(cloud(ii,jj,:),ncld,nlim,ptemp(ii,jj),ppres(ii,jj), zcgno3eqcd, zcgnh3eqcd, 	&
		 zcgwaeqcd,h2oactcd,zpH(ii,jj,nbins+1:nbins+ncld),ph_switch,zkelwacd,prv(ii,jj)/prs(ii,jj),ptstep)
            ! precipitation
            CALL thermoequil(precp(ii,jj,:),nprc,prlim,ptemp(ii,jj),ppres(ii,jj), zcgno3eqpd, zcgnh3eqpd, 	&
		 zcgwaeqpd,h2oactpd,zpH(ii,jj,nbins+ncld+1:nbins+ncld+nprc),ph_switch,zkelwapd,prv(ii,jj)/prs(ii,jj),ptstep)
            
            ! 1) Condensation / evaporation of water
            zdfh2o = DIFFC(ptemp(ii,jj),ppres(ii,jj),1) 
            zdfno3 = DIFFC(ptemp(ii,jj),ppres(ii,jj),2) 
            zdfnh3 = DIFFC(ptemp(ii,jj),ppres(ii,jj),4) 
            
            ! write(*,*) zdfh2o,zdfno3,zdfnh3
            ! pause
            zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air 
            
            ! Initialization of variables
            zsum1 = 0.
            zsum2 = 0.
            zcwan = 0.
            zcwanae = 0.
            zcwancd = 0.
            zcwanpd = 0.
            rhlimx = 1.
            !-- 2) Transition regime correction factor for particles ---------
            !  
            !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
            !  Topics in current aerosol research, Pergamon.  
            !
            !  Size of condensing molecule considered only for 
            !  nucleation mode (3 - 20 nm) 
            !
            
            zmfp   = 3.*zdfno3*sqrt(pi*msu/(8.*rg*ptemp(ii,jj)))                      ! mean free path [m]
            zmfp_h2o = 1.e-2*(mair*1.e3/(pi*avog*zrhoair*1.e-3*(3.11e-8)**2))*	&
                 SQRT(mair/(mair+mwa))
            zmfp_no3 = 1.e-2*(mair*1.e3/(pi*avog*zrhoair*1.e-3*(3.07e-8)**2))*	&
                 SQRT(mair/(mair+mno))
            zmfp_nh3 = 1.e-2*(mair*1.e3/(pi*avog*zrhoair*1.e-3*(4.32e-8)**2))*	&
                 SQRT(mair/(mair+mnh))
            alpha_h2o = 1.!mass accomodation coeficient	Davidovits, P., et al. 2004
            alpha_no3 = 1. !jacobson
            alpha_nh3 = 1. !jacobson
            
            ! 1) Condensation / evaporation of water
            ! aerosol bins
            DO cc = nstr, nbins
               IF (aero(ii,jj,cc)%numc > nlim) THEN 
                  zcgwaeqae(cc)  = h2oactae(cc)*MAX(prs(ii,jj),prv(ii,jj)/1.01)*zrhoair/mwa
                  zknae(cc) = zmfp_h2o/(0.5*dwet_ae(cc))			!16.20
                  zbetaae(cc) = 1. + zknae(cc)*(( 1.33 + (0.71/zknae(cc)))/( 1. + (1./zknae(cc))) + 	&
                       4.*(1-alpha_h2o)/(3.*alpha_h2o) )			
                  zbetaae(cc) = 1./zbetaae(cc)	
                  zDp_ae(cc) = zdfh2o*zbetaae(cc)
                  zDeff_ae(cc) = zDp_ae(cc)/                                              & ! (16.55)
                       (mwa*zDp_ae(cc)*alv*zkelwaae(cc)*zcgwaeqae(cc)/                     &
                       (zthcond*ptemp(ii,jj))*(alv*mwa/(rg*ptemp(ii,jj)) - 1.) + 1.)**.5
                  zmtwaae(cc) = aero(ii,jj,cc)%numc*2.*pi*dwet_ae(cc)*        & ! (16.64)
                       zDeff_ae(cc)
                  zHp_ae(cc,3) = 1.e0                                                     ! initialization
                  IF(zcgwaeqae(cc) > 0. .and. zcwacae(cc) > 0.) zHp_ae(cc,3) = zcwacae(cc)/zcgwaeqae(cc)         ! (17.99)
                  
                  zhlp3 = max(-200.,-adt*zkelwaae(cc)*zmtwaae(cc)/zHp_ae(cc,3))           ! prevent underflow problem on some compilers
                  zexpterm_ae(cc) = exp(zhlp3)                                               ! exponent term in Eq (17.104)
                  IF(prv(ii,jj)/prs(ii,jj) < rhlimx) THEN                     ! APD
                     zsum1 = zsum1 + zcwacae(cc)*(1.-zexpterm_ae(cc))                     ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + zHp_ae(cc,3)/zkelwaae(cc)*(1.-zexpterm_ae(cc))       ! sum term in Eq (17.104) denominator
                  ELSE                                                                       ! APC
                     zsum1 = zsum1 + adt*zmtwaae(cc)*zkelwaae(cc)*zcgwaeqae(cc)              ! sum term in Eq (16.71) numerator
                     zsum2 = zsum2 + adt*zmtwaae(cc)                                         ! sum term in Eq (16.71) denominator
                  END IF
                  
               END IF
            END DO
            
            ! cloud bins
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > nlim) THEN
                  zcgwaeqcd(cc)  = h2oactcd(cc)*MAX(prs(ii,jj),prv(ii,jj)/1.01)*zrhoair/mwa
                  zkncd(cc) = zmfp_h2o/(0.5*dwet_cd(cc))		
                  zbetacd(cc) = 1. + zkncd(cc)*(( 1.33 + (0.71/zkncd(cc)))/( 1. + (1./zkncd(cc))) & 
                       + 4.*(1-alpha_h2o)/(3.*alpha_h2o) )			
                  
                  
                  zbetacd(cc) = 1./zbetacd(cc)	
                  zDp_cd(cc) = zdfh2o*zbetacd(cc)
                  zDeff_cd(cc) = zDp_cd(cc) /                                              & ! (16.55)
                       (mwa*zDp_cd(cc)*alv*zkelwacd(cc)*zcgwaeqcd(cc)/                     &
                       (zthcond*ptemp(ii,jj))*(alv*mwa/(rg*ptemp(ii,jj)) - 1.) + 1.)**.5
                  zmtwacd(cc) =  cloud(ii,jj,cc)%numc*2.*pi*dwet_cd(cc)*     & ! (16.64)
                       zDeff_cd(cc)
                  zHp_cd(cc,3) = 1.e0                                                     ! initialization
                  IF(zcgwaeqcd(cc) > 0. .and. zcwaccd(cc) > 0.) zHp_cd(cc,3) = zcwaccd(cc)/zcgwaeqcd(cc)         ! (17.99)
                  zhlp3 = max(-200.,-adt*zkelwacd(cc)*zmtwacd(cc)/zHp_cd(cc,3))
                  zexpterm_cd(cc) = exp(zhlp3)                                               ! exponent term in Eq (17.104)
                  IF(prv(ii,jj)/prs(ii,jj) < rhlimx) THEN                         ! APD
                     zsum1 = zsum1 + zcwaccd(cc)*(1.-zexpterm_cd(cc))                     ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + zHp_cd(cc,3)/zkelwacd(cc)*(1.-zexpterm_cd(cc))       ! sum term in Eq (17.104) denominator
                  ELSE                                                                       ! APC
                     zsum1 = zsum1 + adt*zmtwacd(cc)*zkelwacd(cc)*zcgwaeqcd(cc)              ! sum term in Eq (16.71) numerator
                     zsum2 = zsum2 + adt*zmtwacd(cc)                                         ! sum term in Eq (16.71) denominator
                  END IF
               END IF
            END DO
            
            ! precipitation bins
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > prlim) THEN
                  zcgwaeqpd(cc)  = h2oactpd(cc)*MAX(prs(ii,jj),prv(ii,jj)/1.01)*zrhoair/mwa
                  zknpd(cc) = zmfp_h2o/(0.5*dwet_pd(cc))			!16.20
                  zbetapd(cc) = 1. + zknpd(cc)*(( 1.33 + (0.71/zknpd(cc)))/( 1. + (1./zknpd(cc))) + 	&
                       4.*(1-alpha_h2o)/(3.*alpha_h2o) )			
                  zbetapd(cc) = 1./zbetapd(cc)	
                  zDp_pd(cc) = zdfh2o*zbetapd(cc)
                  zDeff_pd(cc) = zDp_pd(cc) /                                              & ! (16.55)
                       (mwa*zDp_pd(cc)*alv*zkelwapd(cc)*zcgwaeqpd(cc)/                     &
                       (zthcond*ptemp(ii,jj))*(alv*mwa/(rg*ptemp(ii,jj)) - 1.) + 1.)
                  zmtwapd(cc) =  precp(ii,jj,cc)%numc*2.*pi*dwet_pd(cc)*     & ! (16.64)
                       zDeff_pd(cc)
                  zHp_pd(cc,3) = 1.e0                                                     ! initialization
                  IF(zcgwaeqpd(cc) > 0. .and. zcwacpd(cc) > 0.) zHp_pd(cc,3) = zcwacpd(cc)/zcgwaeqpd(cc)         ! (17.99)
                  zhlp3 = max(-200.,-adt*zkelwapd(cc)*zmtwapd(cc)/zHp_pd(cc,3))
                  zexpterm_pd(cc) = exp(zhlp3)                                               ! exponent term in Eq (17.104)
                  IF(prv(ii,jj)/prs(ii,jj) < rhlimx) THEN                                     ! APD
                     zsum1 = zsum1 + zcwacpd(cc)*(1.-zexpterm_pd(cc))                     ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + zHp_pd(cc,3)/zkelwapd(cc)*(1.-zexpterm_pd(cc))       ! sum term in Eq (17.104) denominator
                  ELSE                                                                       ! APC
                     zsum1 = zsum1 + adt*zmtwapd(cc)*zkelwapd(cc)*zcgwaeqpd(cc)              ! sum term in Eq (16.71) numerator
                     zsum2 = zsum2 + adt*zmtwapd(cc)                                         ! sum term in Eq (16.71) denominator
                  END IF
               END IF
            END DO
            fmax = 0.3
            fmin = -0.2
            ! update the gas phase concentration [mol/m3] of water
            zcwan = MIN(zcwatot,(zcwac + zsum1)/(1. + zsum2)) ! Eq (17.104)
            ! update the particle phase concentration of water in each bin
            !aerosol bins
            DO cc = 1, nbins
               IF (aero(ii,jj,cc)%numc > nlim) THEN
                  IF(prv(ii,jj)/prs(ii,jj) < rhlimx .and. (zcwan-zkelwaae(cc)*zcgwaeqae(cc)) > 0.) THEN                                   ! APD
                     zcwanae(cc) = zHp_ae(cc,3)*zcwan/zkelwaae(cc) +                      & ! (17.102)
                          (zcwacae(cc) - zHp_ae(cc,3)*zcwan/zkelwaae(cc))*zexpterm_ae(cc)
                  ELSE                                                                     ! APC
                     zcwanae(cc) = zcwacae(cc) + min(max(adt*zmtwaae(cc)*(zcwan-zkelwaae(cc)*zcgwaeqae(cc)), &
                          fmin*zcwacae(cc)),fmax*zcwacae(cc))  				!16.69
                  END IF
                  
                  aero(ii,jj,cc)%volc(iwa) = zcwanae(cc)*mwa/rhowa          ! convert to volume concentration
               END IF
            END DO
            !cloud bins
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > nlim) THEN
                  IF(prv(ii,jj)/prs(ii,jj) < rhlimx .and. (zcwan-zkelwacd(cc)*zcgwaeqcd(cc)) > 0.) THEN                                   ! APD
                     zcwancd(cc) = zHp_cd(cc,3)*zcwan/zkelwacd(cc) +                     & ! (17.102)
                          (zcwaccd(cc) - zHp_cd(cc,3)*zcwan/zkelwacd(cc))*zexpterm_cd(cc)
                  ELSE                                                                     ! APC
                     zcwancd(cc) = zcwaccd(cc) + min(max(adt*zmtwacd(cc)*(zcwan-zkelwacd(cc)*zcgwaeqcd(cc)), &
                          fmin*zcwaccd(cc)),fmax*zcwaccd(cc))  				!16.69
                  END IF
                  cloud(ii,jj,cc)%volc(iwa) = zcwancd(cc)*mwa/rhowa          ! convert to volume concentration
               END IF
            END DO
            
            !precipitation bins
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > prlim) THEN
                  IF(prv(ii,jj)/prs(ii,jj) < rhlimx .and. (zcwan-zkelwapd(cc)*zcgwaeqpd(cc)) > 0.) THEN                                    ! APD
                     zcwanpd(cc) = zHp_pd(cc,3)*zcwan/zkelwapd(cc) +                      & ! (17.102)
                          (zcwacpd(cc) - zHp_pd(cc,3)*zcwan/zkelwapd(cc))*zexpterm_pd(cc)
                  ELSE
                     zcwanpd(cc) = zcwacpd(cc) + min(max(adt*zmtwapd(cc)*(zcwan-zkelwapd(cc)*zcgwaeqpd(cc)), &
                          fmin*zcwacpd(cc)),fmax*zcwacpd(cc))  				!16.69
                  END IF
                  precp(ii,jj,cc)%volc(iwa) = zcwanpd(cc)*mwa/rhowa                          ! convert to volume concentration
               END IF
            END DO
            
            zcwan = zcwatot - (sum(zcwanae)+sum(zcwancd)+sum(zcwanpd))
            ! END IF
            
            
            ! 1) Condensation / evaporation of HNO3 
            ! Initialization of variables
            zsum1 = 0.
            zsum2 = 0.
            zcno3nae = 0.
            zcno3ncd = 0.
            zcno3npd = 0.
            ! aerosol bins
            DO cc = nstr, nbins
               IF (aero(ii,jj,cc)%numc > nlim .and. ino>0.) THEN
                  zknae(cc) = zmfp_no3/(0.5*dwet_ae(cc))			!16.20
                  zbetaae(cc) = 1. + zknae(cc)*(( 1.33 + (0.71/zknae(cc)))/( 1. + (1./zknae(cc))) + 	&
                       4.*(1-alpha_no3)/(3.*alpha_no3) )			
                  zbetaae(cc) = 1./zbetaae(cc)	
                  zHp_ae(cc,1) = 1.e0                                                ! initialization
                  IF(zcgno3eqae(cc) > 0. .and. zcno3cae(cc) > 0.) &
                       zHp_ae(cc,1) = zcno3cae(cc)/zcgno3eqae(cc) ! (17.99)
                  zmtno3ae(cc) = 2.*pi*dwet_ae(cc) *  &                     ! (16.64)
                       zdfno3*aero(ii,jj,cc)%numc*zbetaae(cc)
                  zhlp3 = max(-200.,-adt*zkelno3ae(cc)*zmtno3ae(cc)/zHp_ae(cc,1))
                  zexpterm_ae(cc) = exp(zhlp3)                                          ! exponent term in Eq (17.104)
                  zsum1 = zsum1 + zcno3cae(cc)*(1.-zexpterm_ae(cc))                  ! sum term in Eq (17.104) numerator
                  zsum2 = zsum2 + zHp_ae(cc,1)/zkelno3ae(cc)*(1.-zexpterm_ae(cc))    ! sum term in Eq (17.104) denominator
               END IF
            END DO
            
            ! cloud bins
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > nlim .and. ino>0.) THEN
                  zkncd(cc) = zmfp_no3/(0.5*dwet_cd(cc))			!16.20
                  zbetacd(cc) = 1. + zkncd(cc)*(( 1.33 + (0.71/zkncd(cc)))/( 1. + (1./zkncd(cc))) + 	&
                       4.*(1-alpha_no3)/(3.*alpha_no3) )			
                  zbetacd(cc) = 1./zbetacd(cc)	
                  zHp_cd(cc,1) = 1.e0                                                ! initialization
                  IF(zcgno3eqcd(cc) > 0. .and. zcno3ccd(cc) > 0.) zHp_cd(cc,1) = zcno3ccd(cc)/zcgno3eqcd(cc) ! (17.99)
                  zmtno3cd(cc) = 2.*pi*dwet_cd(cc) *  &
                       zdfno3*cloud(ii,jj,cc)%numc*zbetacd(cc)
                  zhlp3 = max(-200.,-adt*zkelno3cd(cc)*zmtno3cd(cc)/zHp_cd(cc,1))
                  zexpterm_cd(cc) = exp(zhlp3)                                          ! exponent term in Eq (17.104)
                  zsum1 = zsum1 + zcno3ccd(cc)*(1.-zexpterm_cd(cc))                  ! sum term in Eq (17.104) numerator
                  zsum2 = zsum2 + zHp_cd(cc,1)/zkelno3cd(cc)*(1.-zexpterm_cd(cc))    ! sum term in Eq (17.104) denominator
               END IF
            END DO
          
            ! precipitation bins
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > prlim .and. ino>0.) THEN
                  zknpd(cc) = zmfp_no3/(0.5*dwet_pd(cc))			!16.20
                  zbetapd(cc) = 1. + zknpd(cc)*(( 1.33 + (0.71/zknpd(cc)))/( 1. + (1./zknpd(cc))) + 	&
                       4.*(1-alpha_no3)/(3.*alpha_no3) )			
                  zbetapd(cc) = 1./zbetapd(cc)
                  zHp_pd(cc,1) = 1.e0                                                ! initialization
                  IF(zcgno3eqpd(cc) > 0. .and. zcno3cpd(cc) > 0.) zHp_pd(cc,1) = zcno3cpd(cc)/zcgno3eqpd(cc) ! (17.99)
                  zmtno3pd(cc) = 2.*pi*dwet_pd(cc) *  &
                       zdfno3*precp(ii,jj,cc)%numc*zbetapd(cc)
                  zhlp3 = max(-200.,-adt*zkelno3pd(cc)*zmtno3pd(cc)/zHp_pd(cc,1))           
                  zexpterm_pd(cc) = exp(zhlp3)                                          ! exponent term in Eq (17.104)
                  zsum1 = zsum1 + zcno3cpd(cc)*(1.-zexpterm_pd(cc))                  ! sum term in Eq (17.104) numerator
                  zsum2 = zsum2 + zHp_pd(cc,1)/zkelno3pd(cc)*(1.-zexpterm_pd(cc))    ! sum term in Eq (17.104) denominator
               END IF
            END DO
            
            ! update the gas phase concentration [mol/m3] of HNO3
            zcno3n = (zcno3c + zsum1)/(1.  + zsum2) ! Eq (17.104)
            ! update the particle phase concentration of NO3- in each bin
            !aerosol bins
            DO cc = nstr, nbins
               IF (aero(ii,jj,cc)%numc > nlim .and. ino>0.) THEN
                  zcno3nae(cc) = max(zcnh3cae(cc)-aero(ii,jj,cc)%volc(1)*rhosu/msu*1.99,&
                       zHp_ae(cc,1)*zcno3n/zkelno3ae(cc) +                & ! (17.102)
                       (zcno3cae(cc) - zHp_ae(cc,1)*zcno3n/zkelno3ae(cc))*zexpterm_ae(cc))
                  aero(ii,jj,cc)%volc(ino) = zcno3nae(cc)*mno/rhono          ! convert to volume concentration
               END IF
            END DO
            
            !cloud bins
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > nlim .and. ino>0.) THEN
                  zcno3ncd(cc) = zHp_cd(cc,1)*zcno3n/zkelno3cd(cc) +                & ! (17.102)
                       (zcno3ccd(cc) - zHp_cd(cc,1)*zcno3n/zkelno3cd(cc))*zexpterm_cd(cc)
                  cloud(ii,jj,cc)%volc(ino) = zcno3ncd(cc)*mno/rhono         ! convert to volume concentration
               END IF
            END DO
            
            !precipitation bins
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > prlim .and. ino>0.) THEN
                  zcno3npd(cc) = zHp_pd(cc,1)*zcno3n/zkelno3pd(cc) +                & ! (17.102)
                       (zcno3cpd(cc) - zHp_pd(cc,1)*zcno3n/zkelno3pd(cc))*zexpterm_pd(cc)
                  precp(ii,jj,cc)%volc(ino) = zcno3npd(cc)*mno/rhono                   ! convert to volume concentration
               END IF
            END DO
            !CALL THERMODYNAMICS AGAIN
            zcgno3eqae=0.
            zcgnh3eqae=0.
            zcgwaeqae=0.
            zcgno3eqcd=0.
            zcgnh3eqcd=0.
            zcgwaeqcd=0.
            zcgno3eqpd=0.
            zcgnh3eqpd=0.
            zcgwaeqpd=0.
            ! aerosol droplets
            CALL thermoequil(aero(ii,jj,:),nbins,nlim,ptemp(ii,jj),ppres(ii,jj), zcgno3eqae, zcgnh3eqae, 	&
                 zcgwaeqae,h2oactae,zpH(ii,jj,1:nbins),ph_switch,zkelwaae,prv(ii,jj)/prs(ii,jj),ptstep)
            ! cloud droplets
            CALL thermoequil(cloud(ii,jj,:),ncld,nlim,ptemp(ii,jj),ppres(ii,jj), zcgno3eqcd, zcgnh3eqcd, 	&
                 zcgwaeqcd,h2oactcd,zpH(ii,jj,nbins+1:nbins+ncld),ph_switch,zkelwacd,prv(ii,jj)/prs(ii,jj),ptstep)
            ! precipitation
            CALL thermoequil(precp(ii,jj,:),nprc,prlim,ptemp(ii,jj),ppres(ii,jj), zcgno3eqpd, zcgnh3eqpd, 	&
                 zcgwaeqpd,h2oactpd,zpH(ii,jj,nbins+ncld+1:nbins+ncld+nprc),ph_switch,zkelwapd,prv(ii,jj)/prs(ii,jj),ptstep)
            
            ! 2) Condensation / evaporation of NH3
            ! Initialization of variables
            zsum1 = 0.
            zsum2 = 0.
            zcnh3nae = 0.
            zcnh3ncd = 0.
            zcnh3npd = 0.
            ! aerosol bins
            DO cc = nstr, nbins
               IF (aero(ii,jj,cc)%numc > nlim .and. inh>0.) THEN
                  zknae(cc) = zmfp_nh3/(0.5*dwet_ae(cc))			!16.20
                  zbetaae(cc) = 1. + zknae(cc)*(( 1.33 + (0.71/zknae(cc)))/( 1. + (1./zknae(cc))) + 	&
                       4.*(1-alpha_nh3)/(3.*alpha_nh3) )			
                  zbetaae(cc) = 1./zbetaae(cc)	
                  zHp_ae(cc,2) = 1.e0                                ! initialize
                  IF(zcgnh3eqae(cc) > 0. .and. zcnh3cae(cc) > 0.) &
                       zHp_ae(cc,2) = zcnh3cae(cc)/zcgnh3eqae(cc)   ! (17.99)
                  zmtnh3ae(cc) = 2.*pi*dwet_ae(cc) *  &
                       zdfnh3*aero(ii,jj,cc)%numc*zbetaae(cc)
                  zhlp3 = max(-200.,-adt*zkelnh3ae(cc)*zmtnh3ae(cc)/zHp_ae(cc,2))
                  zexpterm_ae(cc)=exp(zhlp3)                                              ! exponent term in Eq (17.104)
                  zsum1 = zsum1 + zcnh3cae(cc)*(1.-zexpterm_ae(cc))                    ! sum term in Eq (17.104) numerator
                zsum2 = zsum2 + zHp_ae(cc,2)/zkelnh3ae(cc)*(1.-zexpterm_ae(cc))      ! sum term in Eq (17.104) denominator
             END IF
          END DO
          
          !cloud bins
          DO cc = 1, ncld
             IF (cloud(ii,jj,cc)%numc > nlim .and. inh>0.) THEN
                zkncd(cc) = zmfp_nh3/(0.5*dwet_cd(cc))			!16.20
                zbetacd(cc) = 1. + zkncd(cc)*(( 1.33 + (0.71/zkncd(cc)))/( 1. + (1./zkncd(cc))) + 	&
                     4.*(1-alpha_nh3)/(3.*alpha_nh3) )			
                zbetacd(cc) = 1./zbetacd(cc)	
                zHp_cd(cc,2) = 1.e0                                                  ! initialize
                IF(zcgnh3eqcd(cc) > 0. .and. zcnh3ccd(cc) > 0.) zHp_cd(cc,2) = zcnh3ccd(cc)/zcgnh3eqcd(cc)   ! (17.99)
                zmtnh3cd(cc) = 2.*pi*dwet_cd(cc) *  &
                     zdfnh3*cloud(ii,jj,cc)%numc*zbetacd(cc)
                zhlp3 = max(-200.,-adt*zkelnh3cd(cc)*zmtnh3cd(cc)/zHp_cd(cc,2))
                zexpterm_cd(cc)=exp(zhlp3)                                              ! exponent term in Eq (17.104)
                zsum1 = zsum1 + zcnh3ccd(cc)*(1.-zexpterm_cd(cc))                    ! sum term in Eq (17.104) numerator
                zsum2 = zsum2 + zHp_cd(cc,2)/zkelnh3cd(cc)*(1.-zexpterm_cd(cc))      ! sum term in Eq (17.104) denominator
             END IF
          END DO
          
          ! precipitation bins
          DO cc = 1, nprc
             IF (precp(ii,jj,cc)%numc > prlim .and. inh>0.) THEN
                zknpd(cc) = zmfp_nh3/(0.5*dwet_pd(cc))			!16.20
                zbetapd(cc) = 1. + zknpd(cc)*(( 1.33 + (0.71/zknpd(cc)))/( 1. + (1./zknpd(cc))) + 	&
                     4.*(1-alpha_nh3)/(3.*alpha_nh3) )			
                zbetapd(cc) = 1./zbetapd(cc)	
                zHp_pd(cc,2) = 1.e0                                                  ! initialize
                IF(zcgnh3eqpd(cc) > 0. .and. zcnh3cpd(cc) > 0.) zHp_pd(cc,2) = zcnh3cpd(cc)/zcgnh3eqpd(cc)   ! (17.99)
                zmtnh3pd(cc) = 2.*pi*dwet_pd(cc) *  &
                     zdfnh3*precp(ii,jj,cc)%numc*zbetapd(cc)
                zhlp3 = max(-200.,-adt*zkelnh3pd(cc)*zmtnh3pd(cc)/zHp_pd(cc,2))
                zexpterm_pd(cc)=exp(zhlp3)                                              ! exponent term in Eq (17.104)
                zsum1 = zsum1 + zcnh3cpd(cc)*(1.-zexpterm_pd(cc))                    ! sum term in Eq (17.104) numerator
                zsum2 = zsum2 + zHp_pd(cc,2)/zkelnh3pd(cc)*(1.-zexpterm_pd(cc))      ! sum term in Eq (17.104) denominator
             END IF
          END DO
          
          ! update the gas phase concentration [mol/m3] of NH3
          zcnh3n = (zcnh3c + zsum1)/(1.  + zsum2)	! (17.104)
          ! update the particle phase concentration of NH3 in each bin
          
          !aerosol bins
          DO cc = nstr, nbins
             IF (aero(ii,jj,cc)%numc > nlim .and. inh>0.) THEN
                
                zcnh3nae(cc) = max(min(zHp_ae(cc,2)/zkelnh3ae(cc)*zcnh3n +               & ! (17.102)
                     (zcnh3cae(cc) - zHp_ae(cc,2)/zkelnh3ae(cc)*zcnh3n)*zexpterm_ae(cc), &
                     zcno3nae(cc)+(2.-1.e-5)*aero(ii,jj,cc)%volc(1)*rhosu/msu) &
                     ,0.95*zcnh3cae(cc))
                
                aero(ii,jj,cc)%volc(inh) = zcnh3nae(cc)*mnh/rhonh          ! convert to volume concentration
             END IF
          END DO
          
          ! cloud bins
          DO cc = 1, ncld
             IF (cloud(ii,jj,cc)%numc > nlim .and. inh>0.) THEN
                zcnh3ncd(cc) = max(min(zHp_cd(cc,2)/zkelnh3cd(cc)*zcnh3n +               & ! (17.102)
                     (zcnh3ccd(cc) - zHp_cd(cc,2)/zkelnh3cd(cc)*zcnh3n)*zexpterm_cd(cc), &
                     zcno3ncd(cc)+(2.-1.e-5)*cloud(ii,jj,cc)%volc(1)*rhosu/msu) &
                     ,0.95*zcnh3ccd(cc))

                cloud(ii,jj,cc)%volc(inh) = zcnh3ncd(cc)*mnh/rhonh           ! convert to volume concentration
             END IF
          END DO
          
          ! precipitation bins
          DO cc = 1, nprc
            IF (precp(ii,jj,cc)%numc > prlim .and. inh>0.) THEN
                zcnh3npd(cc) = max(min(zHp_pd(cc,2)/zkelnh3pd(cc)*zcnh3n +               & ! (17.102)
                     (zcnh3cpd(cc) - zHp_pd(cc,2)/zkelnh3pd(cc)*zcnh3n)*zexpterm_pd(cc), &
                     zcno3npd(cc)+(2.-1.e-5)*precp(ii,jj,cc)%volc(1)*rhosu/msu) &
                     ,0.95*zcnh3cpd(cc))

                precp(ii,jj,cc)%volc(inh) = zcnh3npd(cc)*mnh/rhonh           ! convert to volume concentration
             END IF
          END DO

!Make sure that mass balance holds in condensation
          
          zcnh3n = zcnh3tot - &
               (sum(aero(ii,jj,1:nbins)%volc(inh))+  &
               sum(cloud(ii,jj,1:ncld)%volc(inh))+   & 
               sum( precp(ii,jj,1:nprc)%volc(inh)))*rhonh/mnh
         
          zcno3n = zcno3tot -  &
               (sum(aero(ii,jj,1:nbins)%volc(ino))+  &
               sum(cloud(ii,jj,1:ncld)%volc(ino))+   & 
               sum( precp(ii,jj,1:nprc)%volc(ino)))*rhono/mno

          zcwatot=zcwac + sum(zcwacae(1:nbins)) + sum(zcwaccd(1:ncld)) + sum(zcwacpd(1:nprc))
          prv(ii,jj) = zcwan*mwa/zrhoair !      "            water concentration to kg/kg
          pghno3(ii,jj) = zcno3n*avog       ! convert gas phase concentration to #/m3
          pgnh3(ii,jj) = zcnh3n*avog        !      "

       END DO
    END DO
  
  END SUBROUTINE partitioning

 SUBROUTINE thermoequil(ppart,nb,nlim,ptemp,ppress, chno3g, cnh3g, ch2og, h2oact,acidity,ph_switch,zkelvin,prh,ptstep)


    USE mo_salsa_types, ONLY : section
    USE mo_submctl, ONLY : rhosu,msu,    &
                           rhoss,mss,    &
                           rhono,mno,    &
                           rhonh,mnh,    &
                           rhowa,mwa,    &
                           rg, avog, spec


    USE aerosol_thermodynamics, ONLY : inorganic_pdfite

    IMPLICIT NONE

!<testing
    INTEGER, PARAMETER :: NCmax=3, NAmax=5, NNmax=1
    REAL :: MOLAL(-NAmax:NCmax) 
!>
    LOGICAL  :: detailed_thermo
    REAL :: zions(7)                ! mol/m3

    REAL :: zwatertotal,        &   ! Total water in particles (mol/m3) ???
               chcl,                &   ! dummy variable for HCl concentration (not in use)
               zgammas(7)               ! Activity coefficients

    INTEGER :: cc
    INTEGER, INTENT(in) :: nb
    REAL,INTENT(in) :: nlim
    REAL :: pmols(nb,7), paw(nb),zkelvin(nb),prh
    REAL :: c_ions(7),SO4_rat(nb),acidity(nb)
    REAL, INTENT(in)  :: ptemp,ppress ,ptstep
    ! equilibrium gas phase concentrations over a flat surface
    REAL, INTENT(out) :: chno3g(nb), &
                             cnh3g(nb),  &
                             ch2og(nb) , h2oact(nb) 
    TYPE(section), INTENT(in) :: ppart(nb)
    
    REAL :: zKr,               & ! Equilibrium constants (see 
                zKeq,              & ! Jacobson (1999),  Atmos Environ 33, 3635 - 3649), Table 3
                dx,                & ! Change in ion concentration
                zcwl,              & ! Liquid water mol/m3-air
                zhlp,esl,rsl,       H_eff_no3          !

    REAL, PARAMETER ::ztemp0   = 298.15      ! Reference temperature (K)    

                      ! FROM/TO ISOROPIA
         REAL :: zwi_S(5)          ! Total species concentrations in moles/m**3 air
         REAL :: zcntrl_s(2)       ! nug for different types of problems solved
                                       ! different state of aerosols (deliquescent or
                                       ! metastable)
         REAL :: zwt_s(5)          ! ?
         REAL :: zaerliq_S(12)     ! Aqueous-phase concentration array in moles/m**3air
         REAL :: zaersld_S(9)      ! Solid-phase concentration array in moles/m**3 air  
         REAL :: zother_s(6)       ! Solution information array
         REAL :: zgas_s(3)         ! Gas-phase concentration array in moles/m**3 air
         REAL :: Kwa               
         CHARACTER*15 scasi_s          ! Returns the subcase which the input corresponds to
         REAL :: GAMA(6), gamma_no3, gamma_nh3, atm = 101325., H_eff_nh3
	 INTEGER :: thermo_model=3 !1 = ISOROPIA, 2 = PDFiTE, !3 = AIM !0 = No thermodynamic model
	 INTEGER :: ph_switch
      INTEGER :: iso4,ioc,ibc,idu,iss,iwa,ino,inh

      iso4 = spec%getIndex("SO4",notFoundValue = 0)
      ioc = spec%getIndex("OC",notFoundValue = 0)
      ibc = spec%getIndex("BC",notFoundValue = 0)
      idu = spec%getIndex("DU",notFoundValue = 0)
      iss = spec%getIndex("SS",notFoundValue = 0)
      ino = spec%getIndex("NO",notFoundValue = 0)
      inh = spec%getIndex("NH",notFoundValue = 0)
      iwa = spec%getIndex("H2O",notFoundValue = 0)

    ! choose how thermodynamical equilibrium is calculated
    ! .true. for detailed PD-Fite calculations of activity coefficients
    ! .false. for fast calculation, assuming ideality for all species

    	detailed_thermo = .true.
                   
     
    ! If there is no HNO3, PD-Fite is not required
    IF (ino /= 0) then 
        IF(sum(ppart(:)%volc(ino)*rhono/mno) == 0. .and. thermo_model==2) detailed_thermo = .false.
    ENDIF



    ! initialize
    pmols = 0.
    h2oact  = 0.
	ch2og = 0.
	chno3g = 0.
	cnh3g = 0.
	acidity = 0.
    DO cc = 1, nb
       IF(ppart(cc)%numc > nlim .and. ppart(cc)%volc(iwa)>0.) THEN
          
          !Calculate initial concentrations of individual ions
	if (iss>0) then
          zhlp = 2.*ppart(cc)%volc(iso4)*rhosu/msu + ppart(cc)%volc(iss)*rhoss/mss + &
               ppart(cc)%volc(ino)*rhono/mno - ppart(cc)%volc(iss)*rhoss/mss -  &
               ppart(cc)%volc(inh)*rhonh/mnh
	else
		if (ino > 0 .and. inh > 0) then
          		zhlp = 2.*ppart(cc)%volc(iso4)*rhosu/msu  + &
               			ppart(cc)%volc(ino)*rhono/mno  -  &
               			ppart(cc)%volc(inh)*rhonh/mnh
		elseif (ino == 0 .and. inh > 0) then
		          	zhlp = 2.*ppart(cc)%volc(iso4)*rhosu/msu    -  &
               			ppart(cc)%volc(inh)*rhonh/mnh
		elseif (ino == 0 .and. inh == 0) then
		          	zhlp = 2.*ppart(cc)%volc(iso4)*rhosu/msu   
		endif
	endif
          
          zhlp = MAX(zhlp, 1.e-30)
          c_ions= 0.
          c_ions(1) = zhlp                              ! H+
          c_ions(2) = ppart(cc)%volc(iso4)*rhosu/msu       ! SO42-
          c_ions(3) = 0.                             ! HSO4
          IF (ino > 0) c_ions(4) = ppart(cc)%volc(ino)*rhono/mno       ! NO3
          IF (inh > 0) c_ions(5) = MIN(ppart(cc)%volc(inh)*rhonh/mnh,c_ions(2)*2.+c_ions(4)-zhlp)       ! NH4
          if (iss>0) c_ions(6) = ppart(cc)%volc(iss)*rhoss/mss       ! Na
          if (iss>0) c_ions(7) = ppart(cc)%volc(iss)*rhoss/mss       ! Cl
 	
	  !getting the activity of water
	  MOLAL     = 0.
	  MOLAL(-2) = c_ions(2) ! SO42-
	  MOLAL(-3) = c_ions(4) ! NO3-
	  MOLAL(2)  = MIN((2.-1.e-10)*MOLAL(-2)+MOLAL(-3),c_ions(5)) ! NH4+ 
	  MOLAL(1)=   MOLAL(-3)+2.*MOLAL(-2)-MOLAL(2)
	  MOLAL = MOLAL/(ppart(cc)%volc(iwa)*rhowa) 

          zcwl = ppart(cc)%volc(iwa)*rhowa/mwa
          ! Reaction HSO4-    = H+ + SO42-               
          !          D        = A  + B
          

          IF (detailed_thermo) THEN

              zcwl = ppart(cc)%volc(iwa)*rhowa/mwa
              paw(cc) = 0.
              
              if(zcwl/(zcwl+c_ions(2)) > 0.00001) THEN
		   select case(thermo_model)
		   case(1) 
               		! WI moles/m3, aerosol phase, NA, SO4, HN4, NO3, CL
			zwi_s=0
                	if (iss>0) zwi_s(1)= ppart(cc)%volc(iss)*rhoss/mss                   !Na
                	zwi_s(2)= ppart(cc)%volc(iso4)*rhosu/msu               !SO4
                	IF (inh > 0) zwi_s(3)= ppart(cc)%volc(inh)*rhonh/mnh                   !NH4
                	IF (ino > 0) zwi_s(4)= ppart(cc)%volc(ino)*rhono/mno                   !NO3
                	if (iss>0) zwi_s(5)= ppart(cc)%volc(iss)*rhoss/mss                   !Cl

                	zcntrl_s(1)=1 ! reverse problem solved with WI only aerosol phase 
                	zcntrl_s(2)=1  ! The aerosol can have both solid+liquid phases (deliquescent) 

			!We are not licensed to publish ISORROPIA, hence we leaving out all the components of it from this distribution.
			!If interested in using ISORROPIA with this code, download the ISORROPIA source codes and 
			!place them in src/src_salsa/ directory.
			!and uncomment the call to ISORROPIA given below 
			!and also this SRCSF77 = ACTCO.f #isocom.f isofwd.f isorev.f thermo.f in the src/Makefile
			!CALL ISOROPIA(zwi_s, min(prh,1.), ptemp,  zcntrl_s,   &
			!zwt_s, zgas_s, zaerliq_s, zaersld_s, scasi_s, zother_s)

                        !ISORROPIA

			chno3g(cc) = zgas_s(2)
			cnh3g(cc)  = zgas_s(1)

			if (zcwl/(zcwl+sum(c_ions)) > 0.995) then
				!chno3g(cc) = chno3g(cc)
				H_eff_no3 = 3.2e6*exp(8700*(1./ptemp-1./298.))/(c_ions(1)/ppart(cc)%volc(iwa)*1.e-3)
				chno3g(cc) = (zwi_s(4)/ppart(cc)%volc(iwa)*1.e-3)/H_eff_no3*atm/(rg*ptemp)
                                Kwa=1.0e-14

				H_eff_nh3 = 62.0*1.7e-5*exp(450*(1./ptemp-1./298.))*(c_ions(1)/ppart(cc)%volc(iwa)*1.e-3)/Kwa
				cnh3g(cc) = (zwi_s(3)/ppart(cc)%volc(iwa)*1.e-3)/H_eff_nh3*atm/(rg*ptemp)
			endif

		   case(2)
			zions = 0.
          	        zions(1) = zhlp                        ! H+
             	        IF (inh > 0) zions(2) = ppart(cc)%volc(inh)*rhonh/mnh ! NH4
             	        if (iss>0) zions(3) = ppart(cc)%volc(iss)*rhoss/mss ! Na
             	        zions(4) = ppart(cc)%volc(iso4)*rhosu/msu ! SO4
             	        zions(5) = 0. ! HSO4
             	        IF (ino > 0) zions(ino) = ppart(cc)%volc(ino)*rhono/mno ! NO3
             	        if (iss>0) zions(inh) = ppart(cc)%volc(iss)*rhoss/mss ! Cl

             	        zcwl = ppart(cc)%volc(iwa)*rhowa/mwa

             	        ! calculate thermodynamical equilibrium in the particle phase
             	       CALL inorganic_pdfite(.9,ptemp,zions,zcwl,chno3g(cc),chcl,cnh3g(cc),zgammas,pmols(cc,:))

             	       ch2og(cc) =ch2og(cc)/(rg*ptemp)
             	       chno3g(cc)=chno3g(cc)/(rg*ptemp)
             	       cnh3g(cc) =cnh3g(cc)/(rg*ptemp)
		   case(3)

       
	                CALL THERMO(ptemp,MOLAL,ch2og(cc),chno3g(cc),cnh3g(cc),paw(cc))
             	        ch2og(cc) =ch2og(cc)/(rg*ptemp)
             	        chno3g(cc)=chno3g(cc)/(rg*ptemp)
             	        cnh3g(cc) =cnh3g(cc)/(rg*ptemp)
		   case(4)
			!do nothing
		   end select
	      end if

	      h2oact(cc) = zcwl/(zcwl+sum(c_ions))
          ELSE
             ! Equilibrium coefficient (mol/kg)
             zKeq = 1.02e-2*EXP(8.84*(ztemp0/ptemp-1.0)+25.15*(1.-ztemp0/ptemp+LOG(ztemp0/ptemp))) ! Table B7 
          
             ! Convert zKeq to zKr (molm-3) , Equation from Table 3 in Jacobson (1999) 
             !                                Atmos Environ 33, 3635 - 3649                                      
             zKr   = zKeq*(zcwl*mwa)
             ! Eq (17.13)
             dx = (-c_ions(2)-c_ions(1)-zKr &                                   ! Eq (17), Jacobson (1999)
               +SQRT((c_ions(2)+c_ions(1)+zKr)**2. & 
               -4.0*(c_ions(1)*c_ions(2)-c_ions(3)*zKr)))/2.0
              cnh3g(cc) = 0.
              chno3g(cc) = 0.

              !Change concentrations
              c_ions(1) = c_ions(1) + dx
              c_ions(2) = c_ions(2) + dx
              c_ions(3) = c_ions(3) - dx
              ch2og(cc) = zcwl/(zcwl+sum(c_ions))*satvaph2o(ptemp)/(rg*ptemp)
             
              ! vapor pressure of HNO3 at the droplet surface:
              zKeq=2.5e6*EXP(29.17*(ztemp0/ptemp-1.0)+16.83*(1.-ztemp0/ptemp+LOG(ztemp0/ptemp)))/101325. ! Table B.7 
              zKr = (zKeq*(zcwl*mwa)**2*rg*ptemp)                                                ! Table 3 in Jacobson (1999) 
              chno3g(cc) = 0.
              IF(c_ions(4) > 0.) chno3g(cc) = c_ions(4)*c_ions(1)/zKr                                   !    "
              ! vapor pressure of NH3 at the droplet surface:
              zKeq=2.58e17*EXP(64.02*(ztemp0/ptemp-1.0)+11.44*(1.-ztemp0/ptemp+LOG(ztemp0/ptemp)))/101325.**2 ! Table B.7 
              zKr = (zKeq*(zcwl*mwa*rg*ptemp)**2)                              ! Table 3 in Jacobson (1999) 
              cnh3g(cc) = 0.
              IF(chno3g(cc) > 0.) cnh3g(cc) = c_ions(4)*c_ions(1)/(chno3g(cc)*zKr)                      ! 
	      h2oact(cc) = zcwl/(zcwl+sum(c_ions))
	      if (ph_switch == 1) then
		CALL THERMO(ptemp,MOLAL,ch2og(cc),chno3g(cc),cnh3g(cc),paw(cc))
	           acidity(cc) = -log10(MOLAL(1))
	      endif
          END IF
       END IF

    END DO

  END SUBROUTINE thermoequil


   SUBROUTINE condgas(kproma,  kbdim,  klev,    krow,      &
                      pcsa,    pcocnv, pcocsv,  zxsa,      &
                      ptemp,   ppres,  ptstep              )

     USE mo_salsa_types, ONLY : aero, cloud, precp, ice, snow, allSALSA
      USE mo_submctl, ONLY :      &
         pi,                        &
         in1a, in2a,                & ! size bin indices
         fn2b,                      &
         ncld,                      &
         nprc,                      &
         nice,                      &
         nsnw,                      &
         ntotal,                    &
         nlim,                      &
         prlim,                     &
         boltz,                     & ! Boltzmann constant [J/K]
         rg,                        & ! molar gas constant [J/(mol K)]
         pstand,                    & ! standard pressure [Pa]
         mvsu, mvoc,                & ! molecular volumes of sulphate and OC [m3]
         spec,                      &
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         massacc,                   & ! mass accomodation coefficients in each bin
         n3                           ! number of molecules in one 3 nm particle [1]

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::      &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow

      REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep                       ! timestep [s]

      REAL, INTENT(INOUT) ::     &
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         zxsa(kbdim,klev)            ! ratio of sulphuric acid and organic vapor in 3nm particles

      !-- Local variables ----------------------
      INTEGER :: ii, jj    ! loop indices

      REAL ::                        &
         zvisc,                      & ! viscosity of air [kg/(m s)]
         zdfvap,                     & ! air diffusion coefficient [m2/s]
         zmfp,                       & ! mean free path of condensing vapour [m]
         zcs_tot,                    & ! total condensation sink [1/s] (gases)
         zcs_ocsv,                   & ! condensation sink for semivolatile organics [1/s]
         zcs_su,                     & ! condensation sink for sulfate [1/s]
         zcs_ocnv,                   & ! condensation sink for nonvolatile organics [1/s]
                                       ! vapour concentration after time step [#/m3]
         zcvap_new1,                 & ! sulphuric acid
         zcvap_new2,                 & ! nonvolatile organics
         zcvap_new3,                 & ! semivolatile organics
                                       ! change in vapour concentration [#/m3]
         zdvap1,                     & ! sulphuric acid
         zdvap2,                     & ! nonvolatile organics
         zdvap3,                     & ! semivolatile organics

         zdfpart(in1a+1),            & ! particle diffusion coefficient

         zknud(fn2b),                & ! particle Knudsen number
         zknca(ncld),                & ! Knudsen number for cloud droplets and aerosol vapours
         zknpa(nprc),                & ! Knudsen number for rain drops and aerosol vapours
         zknia(nice),                & ! Knudsen number for ice particles and aerosol vapours
         zknsa(nsnw),                & ! Knudsen number for snow flakes and aerosol vapours

         zbeta(fn2b),                & ! transitional correction factor for aerosols
         zbetaca(ncld),              & ! - '' - for condensing aerosol vapours on clouds (is this needed?)
         zbetapa(nprc),              & ! - '' - for condensing aerosol vapours on rain drops
         zbetaia(nice),              & ! - '' - for condensing aerosol vapours on ice (is this needed?)
         zbetasa(nsnw),              & ! - '' - for condensing aerosol vapours on snow flakes

         zcolrate(fn2b),             & ! collision rate of molecules to particles [1/s]
         zcolrate_ocnv(fn2b),        & ! collision rate of organic molecules to particles [1/s]
         zcolrateca(ncld),           & ! Collision rate of aerosol vapour molecules to cloud drops
         zcolratepa(nprc),           & ! Collision rate of gases to rain drops
         zcolrateia(nice),           & ! Collision rate of aerosol vapour molecules to ice particles
         zcolratesa(nsnw),           & ! Collision rate of gases to rain drops

         zdvolsa(fn2b),              & ! change of sulphate volume in each bin [fxm]
         zdvoloc(fn2b),              & !    - " - organics

         zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zn_vs_c,                    & ! ratio of nucleation of all mass transfer in the smallest bin
         zxocnv(kbdim,klev)

      INTEGER :: ioc, iso4


      ioc = spec%getIndex("OC",notFoundValue = 0)
      iso4 = spec%getIndex("SO4",notFoundValue = 0)

      zj3n3 = 0.
      zxocnv = 0.

      zdvolsa = 0.
      zn_vs_c = 0.
      DO jj = 1, klev
         DO ii = 1, kbdim

            zdvoloc = 0.

            !-- 1) Properties of air and condensing gases --------------------
            zvisc  = (7.44523e-3*SQRT(ptemp(ii,jj)**3))/(5093.*(ptemp(ii,jj)+110.4))      ! viscosity of air [kg/(m s)]
            zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]
            zmfp   = 3.*zdfvap*sqrt(pi*spec%msu/(8.*rg*ptemp(ii,jj)))                      ! mean free path [m]

            !-- 2) Transition regime correction factor for particles ---------
            !
            !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
            !  Topics in current aerosol research, Pergamon.
            !
            !  Size of condensing molecule considered only for
            !  nucleation mode (3 - 20 nm)
            !

            !-- particle Knudsen numbers
            zknud(in1a:in1a+1) = 2.*zmfp/(aero(ii,jj,in1a:in1a+1)%dwet+d_sa)              ! Gases on aerosols
            zknud(in1a+2:fn2b) = 2.*zmfp/aero(ii,jj,in1a+2:fn2b)%dwet

            zknca(1:ncld) = 2.*zmfp/cloud(ii,jj,1:ncld)%dwet          ! Knudsen number for gases on cloud drplets

            zknpa(1:nprc) = 2.*zmfp/precp(ii,jj,1:nprc)%dwet          ! Knudsen number for gases on rain drops

            zknia(1:nice) = 2.*zmfp/ice(ii,jj,1:nice)%dwet          ! Knudsen number for gases on ice particles

            zknsa(1:nsnw) = 2.*zmfp/snow(ii,jj,1:nsnw)%dwet          ! Knudsen number for gases on snow flakes

            !-- transitional correction factor
            zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &     ! Aerosol + gas
               (3.*massacc)*(zknud+zknud**2))

            zbetaca = 1. + zknca*( 1.33 + (0.71/zknca) )/( 1. + (1./zknca) ) ! Hydrometeor + gas
            zbetaca = 1./zbetaca

            zbetapa = 1. + zknpa*( 1.33 + (0.71/zknpa) )/( 1. + (1./zknpa) ) ! Rain drop + gas
            zbetapa = 1./zbetapa

            zbetaia = 1. + zknia*( 1.33 + (0.71/zknia) )/( 1. + (1./zknia) ) ! ice + gas
            zbetaia = 1./zbetaia

            zbetasa = 1. + zknsa*( 1.33 + (0.71/zknsa) )/( 1. + (1./zknsa) ) ! Rain drop + gas
            zbetasa = 1./zbetasa
            !-- 3) Collision rate of molecules to particles -------------------
            !
            !  Particle diffusion coefficient considered only for
            !  nucleation mode (3 - 20 nm)
            !

            !-- particle diffusion coefficient [m2/s]
            zdfpart = boltz*ptemp(ii,jj)*zbeta(in1a:in1a+1)/ &
                      (3.*pi*zvisc*aero(ii,jj,in1a:in1a+1)%dwet)

            !-- collision rate (gases on aerosols) [1/s]
            zcolrate = 0.
            zcolrate(in1a:in1a+1) = MERGE( 2.*pi*(aero(ii,jj,in1a:in1a+1)%dwet+d_sa)*    &
                                          (zdfvap+zdfpart)*zbeta(in1a:in1a+1)*           &
                                          aero(ii,jj,in1a:in1a+1)%numc,                 &
                                          0.,                                            &
                                          aero(ii,jj,in1a:in1a+1)%numc > nlim        )

            zcolrate(in1a+2:fn2b) = MERGE( 2.*pi*aero(ii,jj,in1a+2:fn2b)%dwet*zdfvap*       &
                                           zbeta(in1a+2:fn2b)*aero(ii,jj,in1a+2:fn2b)%numc, &
                                           0.,                                               &
                                           aero(ii,jj,in1a+2:fn2b)%numc > nlim        )

            !-- gases on hydrometeors
            zcolrateca = 0.
            zcolrateca(1:ncld) = MERGE( 2.*pi*cloud(ii,jj,1:ncld)%dwet*zdfvap*         &
                                        zbetaca(1:ncld)*cloud(ii,jj,1:ncld)%numc,      &
                                        0.,                                             &
                                        cloud(ii,jj,1:ncld)%numc > nlim           )

            ! Gases on rain drops
            zcolratepa = 0.
            zcolratepa(1:nprc) = MERGE( 2.*pi*precp(ii,jj,1:nprc)%dwet*zdfvap*    &
                                        zbetapa(1:nprc)*precp(ii,jj,1:nprc)%numc, &
                                        0.,                                        &
                                        precp(ii,jj,1:nprc)%numc > prlim       )
            !-- gases on ice particles
            zcolrateia = 0.
            zcolrateia(1:nice) = MERGE( 2.*pi*ice(ii,jj,1:nice)%dwet*zdfvap*      &
                                        zbetaia(1:nice)*ice(ii,jj,1:nice)%numc,   &
                                        0.,                                        &
                                        ice(ii,jj,1:ncld)%numc > prlim         )

            ! Gases on snow flakes
            zcolratesa = 0.
            zcolratesa(1:nsnw) = MERGE( 2.*pi*snow(ii,jj,1:nsnw)%dwet*zdfvap*     &
                                        zbetasa(1:nsnw)*snow(ii,jj,1:nsnw)%numc,  &
                                        0.,                                        &
                                        snow(ii,jj,1:nsnw)%numc > prlim        )

            !-- 4) Condensation sink [1/s] -------------------------------------

            zcs_tot = sum(zcolrate) + sum(zcolrateca) + sum(zcolratepa)+ sum(zcolrateia) + sum(zcolratesa)  ! total sink

            !-- 5) Changes in gas-phase concentrations and particle volume -----
            !
            !--- 5.1) Organic vapours ------------------------

            !---- 5.1.1) Non-volatile organic compound: condenses onto all bins
            IF(pcocnv(ii,jj) > 1.e-10 .AND. zcs_tot > 1.e-30 .AND. ioc > 0) THEN

               zn_vs_c = 0.

               IF(zj3n3(ii,jj,2) > 1.) zn_vs_c = (zj3n3(ii,jj,2))/(zj3n3(ii,jj,2) + &
                                                 pcocnv(ii,jj) * zcolrate(in1a))

               !   collision rate in the smallest bin, including nucleation and condensation
               !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
               !   equation (16.73)
               zcolrate_ocnv = zcolrate
               zcolrate_ocnv(in1a) = zcolrate_ocnv(in1a) + zj3n3(ii,jj,2)/pcocnv(ii,jj)

               zcs_ocnv = zcs_tot + zj3n3(ii,jj,2)/pcocnv(ii,jj)                   ! total sink for organic vapor

               zcvap_new2 = pcocnv(ii,jj)/(1.+ptstep*zcs_ocnv)                  ! new gas phase concentration [#/m3]
               zdvap2 = pcocnv(ii,jj) - zcvap_new2                                 ! change in gas concentration [#/m3]
               pcocnv(ii,jj) = zcvap_new2                                          ! updating vapour concentration [#/m3]

               zdvoloc = zcolrate_ocnv(in1a:fn2b)/zcs_ocnv*mvoc*zdvap2             ! volume change of particles
                                                                                   !  [m3(OC)/m3(air)]

               aero(ii,jj,in1a:fn2b)%volc(ioc) = aero(ii,jj,in1a:fn2b)%volc(ioc) + & !-- change of volume
                                                 zdvoloc                             !   due to condensation in 1a-2b

               ! Condensation on hydromets
               cloud(ii,jj,1:ncld)%volc(ioc) = cloud(ii,jj,1:ncld)%volc(ioc) +  &
                                              zcolrateca(1:ncld)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on rain drops
               precp(ii,jj,1:nprc)%volc(ioc) = precp(ii,jj,1:nprc)%volc(ioc) +  &
                                              zcolratepa(1:nprc)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on ice particles
               ice(ii,jj,1:nice)%volc(ioc) = ice(ii,jj,1:nice)%volc(ioc) +  &
                                            zcolrateia(1:nice)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on snow
               snow(ii,jj,1:nsnw)%volc(ioc) = snow(ii,jj,1:nsnw)%volc(ioc) +  &
                                             zcolratesa(1:nsnw)/zcs_ocnv*mvoc*zdvap2

               !-- Change of number concentration in the smallest bin caused by nucleation
               !   Jacobson (2005), equation (16.75)
               ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
               ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
               IF (zxocnv(ii,jj) > 0.) THEN
                  aero(ii,jj,in1a)%numc = aero(ii,jj,in1a)%numc + &
                     zn_vs_c * zdvoloc(in1a)/mvoc/(n3*zxocnv(ii,jj))
               END IF

            END IF


            !---- 5.1.2) Semivolatile organic compound: regimes 1, 2 and 3
            zcs_ocsv = sum(zcolrate(in2a:fn2b)) +  &       ! sink for semivolatile organics
                       sum(zcolrateca(1:ncld))  +  &       ! ... including condensation on cloud droplets
                       sum(zcolratepa(1:nprc))  +  &       ! and rain drops
                       sum(zcolrateia(1:nice))  +  &       ! and ice particles
                       sum(zcolratesa(1:nsnw))             ! and snow particles

            IF(pcocsv(ii,jj) > 1.e-10 .AND. zcs_ocsv > 1.e-30 .AND. ioc > 0) THEN


               zcvap_new3 = pcocsv(ii,jj)/(1.+ptstep*zcs_ocsv)   ! new gas phase concentration [#/m3]
               zdvap3 = pcocsv(ii,jj) - zcvap_new3                  ! change in gas concentration [#/m3]
               pcocsv(ii,jj) = zcvap_new3                           ! updating gas concentration [#/m3]

               zdvoloc(in2a:fn2b) = zdvoloc(in2a:fn2b) + &                     ! volume change of particles
                                    zcolrate(in2a:fn2b)/zcs_ocsv*mvoc*zdvap3   !  [m3(OC)/m3(air)]

               aero(ii,jj,in1a:fn2b)%volc(ioc) = &                   !-- change of volume due
                  aero(ii,jj,in1a:fn2b)%volc(ioc) + zdvoloc        !   due to condensation in 1a-2b

               ! Condensation on hydromets
               cloud(ii,jj,1:ncld)%volc(ioc) = cloud(ii,jj,1:ncld)%volc(ioc)  +  &
                                              zcolrateca(1:ncld)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on rain drops
               precp(ii,jj,1:nprc)%volc(ioc) = precp(ii,jj,1:nprc)%volc(ioc)  +  &
                                              zcolratepa(1:nprc)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on ice particles
               ice(ii,jj,1:nice)%volc(ioc) = ice(ii,jj,1:nice)%volc(ioc)  +  &
                                            zcolrateia(1:nice)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on snow particles
               snow(ii,jj,1:nsnw)%volc(ioc) = snow(ii,jj,1:nsnw)%volc(ioc)  +  &
                                             zcolratesa(1:nprc)/zcs_ocsv*mvoc*zdvap3

            END IF


            ! ---- 5.2) Sulphate -------------------------------------------
            IF(pcsa(ii,jj) > 1.e-10 .AND. zcs_tot > 1.e-30 .AND. iso4 > 0) THEN

               !-- Ratio of mass transfer between nucleation and condensation

               zn_vs_c = 0.

               IF(zj3n3(ii,jj,1) > 1.) zn_vs_c = (zj3n3(ii,jj,1)) / &
                                                 (zj3n3(ii,jj,1) +  &
                                                 pcsa(ii,jj) * zcolrate(in1a))

               !   collision rate in the smallest bin, including nucleation and condensation
               !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
               !   equation (16.73)
               zcolrate(in1a) = zcolrate(in1a) + zj3n3(ii,jj,1) / pcsa(ii,jj)

               zcs_su = zcs_tot + zj3n3(ii,jj,1) / pcsa(ii,jj)      ! total sink for sulfate

               !--- Sulphuric acid -------------------------
               !
               zcvap_new1 = pcsa(ii,jj) /(1.+ptstep*zcs_su)         ! new gas phase concentration [#/m3]
               zdvap1 = pcsa(ii,jj) - zcvap_new1                    ! change in gas concentration [#/m3]
               pcsa(ii,jj) = zcvap_new1                             ! updating vapour concentration [#/m3]

               zdvolsa = zcolrate(in1a:fn2b)/zcs_su*mvsu*zdvap1     ! volume change of particles
               ! [m3(SO4)/m3(air)] by condensation

               !-- Change of volume concentration of sulphate in aerosol [fxm]
               aero(ii,jj,in1a:fn2b)%volc(iso4) = aero(ii,jj,in1a:fn2b)%volc(iso4) + zdvolsa

               !-- Clouds
               cloud(ii,jj,1:ncld)%volc(iso4) = cloud(ii,jj,1:ncld)%volc(iso4)  +  &
                                              zcolrateca(1:ncld)/zcs_su*mvsu*zdvap1

               ! Rain drops
               precp(ii,jj,1:nprc)%volc(iso4) = precp(ii,jj,1:nprc)%volc(iso4)  +  &
                                              zcolratepa(1:nprc)/zcs_su*mvsu*zdvap1

               !-- Ice clouds
               ice(ii,jj,1:nice)%volc(iso4) = ice(ii,jj,1:nice)%volc(iso4)  +  &
                                            zcolrateia(1:nice)/zcs_su*mvsu*zdvap1

               ! Snow particles
               snow(ii,jj,1:nsnw)%volc(iso4) = snow(ii,jj,1:nsnw)%volc(iso4)  +  &
                                             zcolratesa(1:nsnw)/zcs_su*mvsu*zdvap1

               !-- Change of number concentration in the smallest bin caused by nucleation
               !   Jacobson (2005), equation (16.75)
               IF (zxsa(ii,jj) > 0.) THEN
                  aero(ii,jj,in1a)%numc = aero(ii,jj,in1a)%numc +          &
                                           zn_vs_c * zdvolsa(in1a)/mvsu/(n3*zxsa(ii,jj))
               END IF

            END IF

         END DO ! kbdim

      END DO ! klev

   END SUBROUTINE condgas

   !
   ! ----------------------------------------------------------------------------------------------------------
   !

   SUBROUTINE gpparth2o(kproma, kbdim,  klev,  krow,  &
                        level,  ptemp,  ppres, prs,   &
                        prsi,   prv,    ptstep        )
    
     USE mo_salsa_types, ONLY : aero, cloud, precp, ice, snow
     USE mo_submctl, ONLY : nbins, ncld, nprc,    &
          nice, nsnw, ntotal,            &
          spec,                           &
          mair,                         &
          surfw0, surfi0, rg,           &
          pi, pi6, prlim, nlim,      &
          massacc, avog,  &
          in1a, in2a,  &
          fn2b,            &
          lscndh2oae, lscndh2ocl, lscndh2oic, &
          alv, als 
     USE mo_particle_external_properties, ONLY : calcDiamSALSA
     USE mo_salsa_properties, ONLY : equilibration
     IMPLICIT NONE

      INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
      INTEGER, INTENT(in) :: level
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)

      REAL, INTENT(inout) :: prv(kbdim,klev)

      REAL :: zkelvin(nbins), zkelvincd(ncld), zkelvinpd(nprc), &  ! Kelvin effects
              zkelvinid(nice), zkelvinsd(nsnw)                      ! Kelvin effects ice'n'snow
      REAL :: zcwsurfae(nbins), zcwsurfcd(ncld), zcwsurfpd(nprc), & ! Surface mole concentrations
              zcwsurfid(nice), zcwsurfsd(nsnw)                 ! surface mole concentrations ice'n'snow
      REAL :: zmtae(nbins), zmtcd(ncld), zmtpd(nprc),      & ! Mass transfer coefficients
              zmtid(nice), zmtsd(nsnw)
      REAL :: zwsatae(nbins), zwsatcd(ncld), zwsatpd(nprc), &  ! Water saturation ratios above
              zwsatid(nice), zwsatsd(nsnw)                      ! ice'n'snow
      REAL :: zcwtot                                        ! Total water mole concentration
      REAL :: zcwc, zcwn, zcwint                            ! Current and new water vapour mole concentrations
      REAL :: zcwcae(nbins), zcwnae(nbins), zcwintae(nbins) ! Current and new water mole concentrations in aerosols
      REAL :: zcwccd(ncld), zcwncd(ncld), zcwintcd(ncld)    !     -  ''  -     in cloud droplets
      REAL :: zcwcpd(nprc), zcwnpd(nprc), zcwintpd(nprc)    !     -  ''  -     in rain drops
      REAL :: zcwcid(nice), zcwnid(nice), zcwintid(nice)    !     -  ''  -     in ice particles
      REAL :: zcwcsd(nsnw), zcwnsd(nsnw), zcwintsd(nsnw)    !     -  ''  -     in snow particles
      REAL :: zdfh2o, zthcond,rhoair
      REAL :: zbeta,zknud,zmfph2o
      REAL :: zact, zhlp1,zhlp2,zhlp3
      REAL :: adt,ttot
      REAL :: dwet, cap
      REAL :: zrh(kbdim,klev)
      REAL :: dwice(nice), dwsnow(nsnw)

      REAL :: zaelwc1(kbdim,klev), zaelwc2(kbdim,klev)

      INTEGER :: nstr
      INTEGER :: ii,jj,cc
      INTEGER :: counter

      INTEGER :: iwa,nspec

      zrh(:,:) = prv(:,:)/prs(:,:)
      
      iwa = spec%getIndex("H2O")
      nspec = spec%getNSpec()

      ! Calculate the condensation only for 2a/2b aerosol bins
      nstr = in2a

      ! Save the current aerosol water content
      zaelwc1(:,:) = SUM(aero(:,:,in1a:fn2b)%volc(iwa),DIM=3)*spec%rhowa

      ! For 1a bins do the equilibrium calculation
      CALL equilibration(kproma,kbdim,klev,      &
                         zrh,ptemp,.FALSE. )

      ! If RH < 98 % OR dynamic condensation for aerosols switched off, do equilibrium for all bins
      IF (zrh(1,1) < 0.98 .OR. .NOT. lscndh2oae)  CALL equilibration(kproma,kbdim,klev,      &
                                                                     zrh,ptemp,.TRUE. )

      ! The new aerosol water content after equilibrium calculation
      zaelwc2(:,:) = SUM(aero(:,:,in1a:fn2b)%volc(iwa),DIM=3)*spec%rhowa

      prv(:,:) = prv(:,:) - ( zaelwc2(:,:) - zaelwc1(:,:) )/(ppres(:,:)*mair/(rg*ptemp(:,:)))

      DO jj = 1, klev
         DO ii = 1, kbdim
			! Necessary?
            IF ( .NOT. ( &
                (ANY(cloud(ii,jj,:)%numc > nlim) .OR. ANY(precp(ii,jj,:)%numc > prlim) .AND. lscndh2ocl) .OR. &
                (ANY(ice(ii,jj,:)%numc > prlim) .OR. ANY(snow(ii,jj,:)%numc > prlim) .AND. lscndh2oic) .OR. &
                (ANY(aero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) &
                ) ) CYCLE

            rhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))

            zdfh2o = ( 5./(16.*avog*rhoair*1.e-3*(3.11e-8)**2) ) * &
               SQRT( rg*1e7*ptemp(ii,jj)*mair*1.e3*(spec%mwa+mair)*1.e3/( 2.*pi*spec%mwa*1.e3 ) )
            zdfh2o = zdfh2o*1.e-4

            zmfph2o = 3.*zdfh2o*sqrt(pi*spec%mwa/(8.*rg*ptemp(ii,jj)))
            zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air

            ! -- Water vapour (Follows the analytical predictor method by Jacobson 2005)
            zkelvinpd = 1.; zkelvincd = 1.; zkelvin = 1.; zkelvinid = 1.; zkelvinsd = 1.

            zcwc = 0.; zcwint = 0.; zcwn = 0.
            zcwcae = 0.; zcwccd = 0.; zcwcpd = 0.; zcwcid = 0.; zcwcsd = 0.;
            zcwintae = 0.; zcwintcd = 0.; zcwintpd = 0.; zcwintid = 0.; zcwintsd = 0.
            zcwnae = 0.; zcwncd = 0.; zcwnpd = 0.; zcwnid = 0.; zcwnsd = 0.
            zwsatae = 0.; zwsatcd = 0.; zwsatpd = 0.; zwsatid = 0.; zwsatsd = 0.
            
            zmtpd(:) = 0.
            zcwsurfpd(:) = 0.
            zmtcd(:) = 0.
            zcwsurfcd(:) = 0.
            zmtid(:) = 0.
            zcwsurfid(:) = 0.
            zmtsd(:) = 0.
            zcwsurfsd(:) = 0.
            zmtae(:) = 0.
            zcwsurfae(:) = 0.

            ! Cloud droplets --------------------------------------------------------------------------------
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > nlim .AND. lscndh2ocl) THEN

                  ! Wet diameter
                  dwet = ( SUM(cloud(ii,jj,cc)%volc(1:nspec))/cloud(ii,jj,cc)%numc/pi6 )**(1./3.)

                  ! Activity + Kelvin effect
                  zact = acth2o(cloud(ii,jj,cc))
                  zkelvincd(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentration over flat surface
                  zcwsurfcd(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatcd(cc) = zact*zkelvincd(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.)*(zknud+zknud**2))

                  ! Mass transfer according to Jacobson
                  zhlp1 = cloud(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatcd(cc)*zcwsurfcd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtcd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Rain drops --------------------------------------------------------------------------------
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > prlim .AND. lscndh2ocl) THEN
                  ! Wet diameter
                  dwet = ( SUM(precp(ii,jj,cc)%volc(1:nspec))/precp(ii,jj,cc)%numc/pi6 )**(1./3.)

                  ! Activity + Kelvin effect
                  zact = acth2o(precp(ii,jj,cc))
                  zkelvinpd(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentrations over flat surface
                  zcwsurfpd(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatpd(cc) = zact*zkelvinpd(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.)*(zknud+zknud**2))

                  ! Mass transfer according to Jacobson
                  zhlp1 = precp(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatpd(cc)*zcwsurfpd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtpd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Ice particles --------------------------------------------------------------------------------
            ! Dimension
            CALL calcDiamSALSA(nice,ice(ii,jj,:),dwice)
            DO cc = 1, nice
               IF (ice(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN
                  dwet=dwice(cc)
                     
                  ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                  cap=0.5*dwet
                     
                  ! Activity + Kelvin effect - edit when needed
                  !   Can be calculated just like for sperical homogenous particle or just ignored,
                  !   because these are not known for solid, irregular and non-homogenous particles.
                  !   Ice may not be that far from a sphere, but most particles are large and at least
                  !   growing particles are covered by a layer of pure ice.
                  zact = 1.0 
                  zkelvinid(cc) = exp( 4.*surfi0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )
                  
                  ! Saturation mole concentration over flat surface
                  zcwsurfid(cc) = prsi(ii,jj)*rhoair/spec%mwa
                  
                  ! Equilibrium saturation ratio
                  zwsatid(cc) = zact*zkelvinid(cc)
                  
                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                       (3.)*(zknud+zknud**2))
                  
                  ! Mass transfer according to Jacobson
                  zhlp1 = ice(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*als*zwsatid(cc)*zcwsurfid(cc)/(zthcond*ptemp(ii,jj)) 
                  zhlp3 = ( (als*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.
                  
                  zmtid(cc) = zhlp1/( zhlp2*zhlp3 + 1. )
                  
               END IF
            END DO
            
            ! Snow particles --------------------------------------------------------------------------------
            ! Dimension
            CALL CalcDiamSALSA(nsnw,snow(ii,jj,:),dwsnow)
            DO cc = 1, nsnw
               IF (snow(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN
                  dwet=dwsnow(cc)
                  
                  ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                  cap=0.5*dwet
                  
                  ! Activity + Kelvin effect
                  !   Can be calculated just like for sperical homogenous particle or just ignored,
                  !   because these are not known for solid, irregular and non-homogenous particles.
                  !   Especially snow is typically highly irregular (e.g. dendrite).
                  zact = 1.0 
                  zkelvinsd(cc) = exp( 4.*surfi0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )
                  
                  ! Saturation mole concentrations over flat surface
                  zcwsurfsd(cc) = prsi(ii,jj)*rhoair/spec%mwa
                  
                  ! Equilibrium saturation ratio
                  zwsatsd(cc) = zact*zkelvinsd(cc)
                  
                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                       (3.)*(zknud+zknud**2))
                  
                  ! Mass transfer according to Jacobson
                  zhlp1 = snow(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*als*zwsatsd(cc)*zcwsurfsd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (als*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.
                  
                  zmtsd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )
                  
               END IF
            END DO
            
            ! -- Aerosols: ------------------------------------------------------------------------------------
            DO cc = 1, nbins
               IF (aero(ii,jj,cc)%numc > nlim .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) THEN
                  ! Wet diameter
                  dwet = ( SUM(aero(ii,jj,cc)%volc(1:nspec))/aero(ii,jj,cc)%numc/pi6 )**(1./3.)

                  ! Water activity + Kelvin effect
                  zact = acth2o(aero(ii,jj,cc))
                  zkelvin(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentration over flat surface
                  ! Limit the supersaturation to max 1.01 for the mass transfer
                  ! EXPERIMENTAL
                  zcwsurfae(cc) = MAX(prs(ii,jj),prv(ii,jj)/1.01)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatae(cc) = zact*zkelvin(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.*massacc(cc))*(zknud+zknud**2))

                  ! Mass transfer
                  zhlp1 = aero(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatae(cc)*zcwsurfae(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtae(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Current mole concentrations
            zcwc = prv(ii,jj)*rhoair/spec%mwa
            zcwcae(1:nbins) = aero(ii,jj,1:nbins)%volc(iwa)*spec%rhowa/spec%mwa
            zcwccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(iwa)*spec%rhowa/spec%mwa
            zcwcpd(1:nprc) = precp(ii,jj,1:nprc)%volc(iwa)*spec%rhowa/spec%mwa
            zcwcid(1:nice) = ice(ii,jj,1:nice)%volc(iwa)*spec%rhoic/spec%mwa
            zcwcsd(1:nsnw) = snow(ii,jj,1:nsnw)%volc(iwa)*spec%rhosn/spec%mwa


            zcwtot = zcwc + SUM(zcwcae) + &
                     SUM(zcwccd) + &
                     SUM(zcwcpd) + &
                     SUM(zcwcid) + &
                     SUM(zcwcsd)
            ttot = 0.

            zcwintae = zcwcae; zcwintcd = zcwccd; zcwintpd = zcwcpd; zcwintid = zcwcid; zcwintsd = zcwcsd

            ! Substepping loop
            ! ---------------------------------
            zcwint = 0.
            counter = 0
            DO WHILE (ttot < ptstep)

               adt = 2.e-2
               ! New vapor concentration
               zhlp1 = zcwc + adt * ( SUM(zmtae(nstr:nbins)*zwsatae(nstr:nbins)*zcwsurfae(nstr:nbins))  + &
                                      SUM(zmtcd(1:ncld)*zwsatcd(1:ncld)*zcwsurfcd(1:ncld))              + &
                                      SUM(zmtpd(1:nprc)*zwsatpd(1:nprc)*zcwsurfpd(1:nprc))              + &
                                      SUM(zmtid(1:nice)*zwsatid(1:nice)*zcwsurfid(1:nice))              + &
                                      SUM(zmtsd(1:nsnw)*zwsatsd(1:nsnw)*zcwsurfsd(1:nsnw))              )

               zhlp2 = 1. + adt * ( SUM(zmtae(nstr:nbins)) + SUM(zmtcd(1:ncld)) + SUM(zmtpd(1:nprc)) &
                                  + SUM(zmtid(1:nice)) + SUM(zmtsd(1:nsnw)) )
               zcwint = zhlp1/zhlp2
               zcwint = MIN(zcwint,zcwtot)

               IF ( ANY(aero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98 ) THEN
                  DO cc = nstr, nbins
                     zcwintae(cc) = zcwcae(cc) + min(max(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae(cc)), &
                                                     -0.02*zcwcae(cc)),0.05*zcwcae(cc))
                     zwsatae(cc) = acth2o(aero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
                  END DO
               END IF
               IF ( ANY(cloud(ii,jj,:)%numc > nlim) ) THEN
                  DO cc = 1, ncld
                     zcwintcd(cc) = zcwccd(cc) + min(max(adt*zmtcd(cc)*(zcwint - zwsatcd(cc)*zcwsurfcd(cc)), &
                                                     -0.02*zcwccd(cc)),0.05*zcwccd(cc))
                     zwsatcd(cc) = acth2o(cloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
                  END DO
               END IF
               IF ( ANY(precp(ii,jj,:)%numc > prlim) ) THEN
                  DO cc = 1, nprc
                     zcwintpd(cc) = zcwcpd(cc) + min(max(adt*zmtpd(cc)*(zcwint - zwsatpd(cc)*zcwsurfpd(cc)), &
                                                     -0.02*zcwcpd(cc)),0.05*zcwcpd(cc))
                     zwsatpd(cc) = acth2o(precp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
                  END DO
               END IF
               IF (ANY(ice(ii,jj,:)%numc > prlim) ) THEN
                  DO cc = 1, nice
                     zcwintid(cc) = zcwcid(cc) + min(max(adt*zmtid(cc)*(zcwint - zwsatid(cc)*zcwsurfid(cc)), &
                          -0.02*zcwcid(cc)),0.05*zcwcid(cc))
                     zwsatid(cc) = zkelvinid(cc)
                  END DO
               END IF
               IF (ANY(snow(ii,jj,:)%numc > prlim) ) THEN
                  DO cc = 1, nsnw
                     zcwintsd(cc) = zcwcsd(cc) + min(max(adt*zmtsd(cc)*(zcwint - zwsatsd(cc)*zcwsurfsd(cc)),&
                          -0.02*zcwcsd(cc)),0.05*zcwcsd(cc))
                     zwsatsd(cc) = zkelvinsd(cc)
                  END DO
               END IF

               zcwintae(nstr:nbins) = MAX(zcwintae(nstr:nbins),0.)
               zcwintcd(1:ncld) = MAX(zcwintcd(1:ncld),0.)
               zcwintpd(1:nprc) = MAX(zcwintpd(1:nprc),0.)
               zcwintid(1:nice) = MAX(zcwintid(1:nice),0.)
               zcwintsd(1:nsnw) = MAX(zcwintsd(1:nsnw),0.)

               ! Update vapor concentration for consistency
               zcwint = zcwtot - SUM(zcwintae(1:nbins)) - &
                        SUM(zcwintcd(1:ncld))  - &
                        SUM(zcwintpd(1:nprc))  - &
                        SUM(zcwintid(1:nice))  - &
                        SUM(zcwintsd(1:nsnw))

               ! Update "old" values for next cycle
               zcwcae = zcwintae; zcwccd = zcwintcd; zcwcpd = zcwintpd; zcwcid =zcwintid; zcwcsd = zcwintsd;
               zcwc = zcwint

               ttot = ttot + adt

               counter = counter + 1
               IF (counter > 100) WRITE(*,*) "CONDENSATION SUBSTEPING: SOMETHING'S WRONG"

            END DO ! ADT

            zcwn = zcwint
            zcwnae = zcwintae
            zcwncd = zcwintcd
            zcwnpd = zcwintpd
            zcwnid = zcwintid
            zcwnsd = zcwintsd

            prv(ii,jj) = zcwn*spec%mwa/rhoair

            aero(ii,jj,1:nbins)%volc(iwa) = max(0.,zcwnae(1:nbins)*spec%mwa/spec%rhowa)
            cloud(ii,jj,1:ncld)%volc(iwa) = max(0.,zcwncd(1:ncld)*spec%mwa/spec%rhowa)
            precp(ii,jj,1:nprc)%volc(iwa) = max(0.,zcwnpd(1:nprc)*spec%mwa/spec%rhowa)
            ice(ii,jj,1:nice)%volc(iwa) = max(0.,zcwnid(1:nice)*spec%mwa/spec%rhoic)
            snow(ii,jj,1:nsnw)%volc(iwa) = max(0.,zcwnsd(1:nsnw)*spec%mwa/spec%rhosn)

         END DO !kproma

      END DO ! klev

   END SUBROUTINE gpparth2o
 
      FUNCTION DIFFC(Temp,Press,I)
	!
	! This function is used to calculate the diffusion coefficient from the
	! vapor in the gas as a function of temperature and pressure after expans.
	! I is species number.
	!
      Real :: Temp, Press
	!     
      Real :: MYY(6), DN(6), DIFFC
	INTEGER :: I,J
      DATA ( DN(J),MYY(J), J=1,6 ) /	&
         22.0748D-6, 1.6658D+0, 12.98D-6, 1.75D+0, 9.380307D-6, &
         1.75D+0, 20.0D-6, 1.75D+0, 14.76518D-6, 1.75D+0, 0.0, 0.0 /
	!   
      DIFFC = DN(I) * ( (TEMP/273.15)**MYY(I) ) * 101325. / Press
	!     
      RETURN
      END

  FUNCTION pH2O0(T)
    IMPLICIT NONE
!
!     Below 2oC an expression for pH2O0 based on Goff-Gratch is used
!     but modified so that the heat capacities of (supercooled) pure
!     water are those implied by an extrapolation of Hill's equation
!     of state. Above 2oC, an empirical fit to pH2O0 given in the CRC 
!     Handbook is used, valid to 100oC. At the 'join', the fit gives 
!     6.96289E-3 atm, and the lower temperature expression yields 
!     6.96157E-3 atm, a difference of 0.019%.
!
      REAL, intent(in) :: T
      REAL :: LNK1,LNK2, Agas, Bgas, Cgas, Dgas, Egas, delHr, R, T1, T2
      REAL :: Aliq,Bliq,Cliq,Dliq,Eliq,delA,delB,delC,delD,delE,DUM,pH2O0
      
      DATA LNK1, delHr/1.809429005D0, 45.077441D3/
      DATA R/8.3144D0/
      DATA Agas,Bgas,Cgas,Dgas,Egas/33.269811D0,0.00113261D0, &
                                   -1.09982D-5,               &
                                    3.573575D-8,0.D0/
!     Hill-based liquid phase heat capacities:
      DATA Aliq,Bliq,Cliq,Dliq,Eliq/295.1612D0,-1.540498D0,   &
                                    2.7023D-3,                &
                                    0.0D0,0.0D0/
!
!
!
      IF(T .LE. 275.15D0) THEN
!     ..below 2 degrees C, used the Goff-Gratch based expression:
!
        delA=Agas - Aliq
        delB=Bgas - Bliq
        delC=Cgas - Cliq
        delD=Dgas - Dliq
        delE=Egas - Eliq
!
        T2=T
        T1=273.15D0
        DUM=delA/R*LOG(T2/T1) + delB/(2*R)*(T2-T1)                     &
           + delC/(6*R)*(T2**2-T1**2)                                  &
           + delD/(12*R)*(T2**3-T1**3)                                 &
           + delE/(2*R)*(1.D0/T2**2-1.D0/T1**2)                        &
           + (1.D0/R)*(-delHr + delA*T1 + delB/2*T1**2 + delC/3*T1**3  &
                     + delD/4*T1**4 - delE/T1)*(1.D0/T2 - 1.D0/T1)
!
        LNK2 = DUM + LNK1
        pH2O0 = EXP(LNK2) * 0.000986923D0
!                           ^- converts to atm.
!
      ELSE
!     ..else use empirical fit to CRC tabulated vapour pressures:
!
        DUM=23.54872D0 - 6459.987931D0/T - 0.022791752D0*T             &
           + 1.6290826D-5*T**2
        pH2O0 = EXP(DUM) 
      ENDIF
!
      pH2O0 = pH2O0 *101325.
      END



!
! ---------------------------------------------------------------------
! This function calculates the liquid saturation vapor mixing ratio as
! a function of temperature and pressure
!
  subroutine rslf(p,t,esl,rsl)

  real, intent (in) :: p, t
  real, parameter :: c0=0.6105851e+03, c1=0.4440316e+02,    &
                     c2=0.1430341e+01, c3=0.2641412e-01,    &
                     c4=0.2995057e-03, c5=0.2031998e-05,    &
                     c6=0.6936113e-08, c7=0.2564861e-11,    &
                     c8=-.3704404e-13

  real ::   x
  real, intent (out) ::  esl, rsl
  x=min(max(-80.,t-273.16),50.)

  ! esl=612.2*exp(17.67*x/(t-29.65))
  esl=c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
  rsl=.622*esl/(p-esl)

  end subroutine rslf

   !-------------------------------------------------------
   REAL FUNCTION acth2o(ppart,pcw)

      USE classSection, ONLY : Section
      USE mo_submctl, ONLY : eps, spec
      IMPLICIT NONE

      TYPE(Section), INTENT(in) :: ppart
      REAL, INTENT(in), OPTIONAL  :: pcw

      REAL :: zns, znw
      INTEGER :: ndry, iwa, nn
      
      ndry = spec%getNSpec(type="dry")
      iwa = spec%getIndex("H2O")
      
      ! This is only relevant for solution particles so use rholiq
      zns = 0.
      DO nn = 1,ndry  ! Leaves out water and non-soluble species (zero dissociation factor)
         zns = zns + spec%diss(nn)*ppart%volc(nn)*spec%rholiq(nn)/spec%MM(nn)
      END DO 

      IF (present(pcw)) THEN
         znw = pcw
      ELSE
         znw = ppart%volc(iwa)*spec%rholiq(iwa)/spec%MM(iwa)
      END IF

      ! Assume activity coefficient of 1 for water...
      acth2o = MAX(0.1,znw/max(eps,(znw+zns)))

   END FUNCTION acth2o

   ! ------------------------------------------------------------------

   FUNCTION satvaph2o(ptemp) RESULT(psat)
      !-----------------------------------------------------------------
      ! Saturation vapour pressure of water vapour
      ! This is a local function for the subroutine *cloud_condensation*
      !
      ! J. Tonttila, FMI, 03/2014
      !-----------------------------------------------------------------
    
      IMPLICIT NONE

      REAL, INTENT(in) :: ptemp

      REAL, PARAMETER ::       &
         za0 = 6.107799961,    &
         za1 = 4.436518521e-1, &
         za2 = 1.428945805e-2, &
         za3 = 2.650648471e-4, &
         za4 = 3.031240396e-6, &
         za5 = 2.034080948e-8, &
         za6 = 6.136820929e-11

      REAL :: zt

      REAL :: psat

      zt = ptemp - 273.15

      psat = za0 + za1*zt + za2*zt**2 + za3*zt**3 +   &
             za4*zt**4 + za5*zt**5 + za6*zt**6

      ! To Pascals
      psat = psat*100.

   END FUNCTION satvaph2o
 
END MODULE mo_salsa_dynamics
