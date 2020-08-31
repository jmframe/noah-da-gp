program run_timestep
 use noahmp_veg_parameters
 use noahmp_globals
 use type_decs

! -- Variable Declarations ----------------------------------------------
 implicit none

 ! output parameters
 logical, parameter :: write_model_output = .true.
 logical, parameter :: write_ensemble_output = .false.
 logical, parameter :: calc_obj_fun = .true.

 ! simulation parameters, 
 !jmframe added sig_sm & Nlag as read in values, not set in code.
 integer :: sig_sm, Nlag
 integer :: Ne, Nt
 logical :: perturb, data_assim

 ! file I/O
 integer, parameter :: fid = 15
 character(1000)    :: fname, fdir
 character(1)       :: es1
 character(2)       :: es2
 character(3)       :: es3

 ! internal indexes
 integer :: s, t, e, ee, l, d, lags
 real    :: dummy
 integer :: vegtyp
 real    :: rz
 logical :: fixed, fexists, fixensemble = .true.
 integer, allocatable, dimension(:,:) :: date
 real, allocatable, dimension(:) :: time

 ! model data
 type(forcing_data), allocatable, dimension(:,:) :: forcing
 type(state_data),   allocatable, dimension(:,:) :: state, background
 type(state_data) :: state_tmp
 type(setup_data) :: setup
 type(setup_data),   allocatable, dimension(:)   :: setup_tmp
 type(output_data),  allocatable, dimension(:,:) :: output

 ! obs data
 integer, parameter                  :: Dz = 1
 real, allocatable, dimension(:)   :: obs
 real, dimension(Dz)                 :: zcov
 real, allocatable, dimension(:,:,:) :: X
 real, allocatable, dimension(:,:)   :: Y
 real, allocatable, dimension(:)     :: Z, R

 !CRAIG
 real :: noise
 logical :: isOneStep
 !JMFRAME
 logical :: perturb_one_step = .false.
 real(kind=8), dimension(3) :: lagged

 ! Forcing perturbations
 real, allocatable, dimension(:,:) :: P_forcing
 integer :: info
 real :: rn, rat, rnp, diff
 real, dimension(3,3) :: R_forcing
 real, dimension(3) :: S_forcing
 real :: random_normal

 ! mean output data
 real ::  mean_sfctmp,mean_sfcspd,mean_sfcprs,mean_q2,mean_lwrad, &
 mean_swrad,mean_prcprate,mean_smc1,mean_smc2,mean_smc3,mean_smc4,&
 mean_wood,mean_lfmass,mean_stmass,mean_rtmass,mean_qe,mean_qh,mean_nee      

! --- Set Up Run --------------------------------------------------------
! setup simulation

 call sim_init(setup,state_tmp,perturb,data_assim,Nt,Ne,sig_sm,Nlag,isOneStep)
 allocate(state(Nt,Ne))
 allocate(background(Nt,Ne))
 allocate(setup_tmp(Ne))
 allocate(output(Nt,Ne))
 do e = 1,Ne
   call sim_init(setup_tmp(e),state(1,e),perturb,data_assim,Nt,Ne,sig_sm,Nlag,isOneStep)
   do t = 2,Nt
     allocate(state(t,e)%stc(-setup%nsnow+1:setup%nsoil))
     allocate(state(t,e)%zsnso(-setup%nsnow+1:setup%nsoil))
     allocate(state(t,e)%tsno(setup%nsnow))
     allocate(state(t,e)%snice(setup%nsnow))
     allocate(state(t,e)%snliq(setup%nsnow))
     allocate(state(t,e)%sh2o(setup%nsoil))
     allocate(state(t,e)%smc(setup%nsoil))
     allocate(background(t,e)%stc(-setup%nsnow+1:setup%nsoil))
     allocate(background(t,e)%zsnso(-setup%nsnow+1:setup%nsoil))
     allocate(background(t,e)%tsno(setup%nsnow))
     allocate(background(t,e)%snice(setup%nsnow))
     allocate(background(t,e)%snliq(setup%nsnow))
     allocate(background(t,e)%sh2o(setup%nsoil))
     allocate(background(t,e)%smc(setup%nsoil))
   enddo
 enddo

! --- Load Observations -------------------------------------------------
 ! We need to load the observation for one step state update and data assimilation
 ! But there may be other instances where we'd want the observation loaded as well
 ! If so, will remove later.
 if ((data_assim).or.(isOneStep)) then
   ! jmframe July 2019, moved this up before the read-in from the observations file.
   ! It is kind of strange that it was below in the first place. Not sure why... 
   allocate(obs(Nt))
   obs = -9999.
   ! Take out of the data assimilation logic, because we use it for other things. 
   open(fid,file='obs.txt')
     do t = 1,Nt
       ! jmframe: Make sure matches observation file. Added 1 dummy.
       !Year Day hour   ?     ?     Soil Moisture? Soil Moisture? Precip Rate? 
       !2002 1   0.000 -0.92 -26.23 0.4086666870   0.4216666794   0.0000559161
       read(fid,*) dummy,dummy,dummy,dummy,dummy,obs(t),dummy,dummy
     enddo ! time
   close(fid)
 endif
 ! We only need the covariance for data assimilation, not for one step.
 if((data_assim).or.(perturb)) then
   fname = 'obs_cov.txt'
   open(fid,file=trim(fname))
    read(fid,*) zcov
   close(fid)
 endif

 !CRAIG add noise to observation
 if (perturb_one_step) then
  do t=1,Nt
   obs(t) = obs(t) + sqrt(noise)*random_normal()
  enddo
  open(fid,file='obs.txt')
  do t = 1,Nt
   write(fid,*) obs(t)
  enddo
  close(fid)
 endif

! forcing from file
 allocate(forcing(Nt,Ne))
 allocate(date(Nt,2))
 allocate(time(Nt))
 fdir = './forcing'

 !! Maintain Ensemble number 1 without any perturbations to the forcing data
 e = 1
 ! There was logic to name the ensembles with equal of digits, but I only print 9.
 fname = trim(fdir)//'.txt'
 open(fid,file=trim(fname))
 do t = 1,Nt
   ! the humidity here is kg/kg, not % and not relative humidity.
   read(fid,*) date(t,:),time(t),forcing(t,e)%sfcspd,dummy,   &
               forcing(t,e)%sfctmp,forcing(t,e)%q2,                 &
               forcing(t,e)%sfcprs,forcing(t,e)%swrad,              &
               forcing(t,e)%lwrad,forcing(t,e)%prcprate
 enddo ! times
 close(fid)

 if (Ne > 1) then
   !!!!!!!!!!!!!!!!!! Reference for Forcing Data perturbations !!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2011WR011420
   ! Perturbations on solar radiation and precipitation were 
   ! multiplicative lognormally distributed
   ! with a mean of 1 and standard deviations of 0.3 and 0.5,
   ! respectively; perturbations on temperature were additive
   ! Gaussian with a mean of 0 and unit variance, and the same
   ! daily perturbation was applied to daily maximum and minimum temperature.
  
   ! R = covariance matrix: always positive semidefinite and usually positive definite 
   ! https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2011WR011420 Page 7
   !                      Radiation    | Temperature     |  Precipitation
   R_forcing = reshape((/1.0, -0.8, 0.4, -0.8, 1.0, -0.32, 0.4, -0.32, 1.0/),shape(R_forcing))
   !R_forcing = reshape((/1.0, -0.4, 0.2, -0.4, 1.0, -0.151, 0.2, -0.151, 1.0/),shape(R_forcing))
   ! Cholesky decomposition
   ! DPOTRF VARIANT: top-looking block version of the algorithm, calling Level 3 BLAS.
    call dpotrf('U', 3, R_forcing, 3, info)
  
   ! Make consistent units.
   ! radiation, temp, and precip
   ! W/m2     ,C or K, kg/m2/s or mm/s 
   ! https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2017WR020991
   !     Shortwave radiation: 0.3
   !     Longwave Raduation: 50 W/m2
   !     Precipitation: 0.5
   ! Convert standard deviations according to paper. 
   ! 0.3 Mj/m2/d = 0.3 * 11.574 W/m2
   ! 0.5 mm/day = 0.5 * ( 1 mm/s * 60 s/min * 30 min/timestep) 
   S_forcing = (/0.3, 1.0, 0.35/)
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Ensemble loop to add perturbations to the foring data !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do e = 2,Ne
  
    ! Generate random numbers for the log normal transformed perts
    allocate(P_forcing(Nt,3))
    do t = 1, Nt
     P_forcing(t,1) = random_normal() !Radiation
     P_forcing(t,2) = random_normal() !Temperature
     P_forcing(t,3) = random_normal() !Precipitation
    enddo ! end TIME LOOP
  
   ! Multiply the perturbation and the covariance matrix.
    P_forcing = matmul(P_forcing, R_forcing)
  
    ! Scale perturbations based on standard deviations
    do t = 1, Nt
     P_forcing(t,1) = exp(S_forcing(1) * P_forcing(t,1)) !Radiation
     P_forcing(t,2) = S_forcing(2) * P_forcing(t,2)      !Temperature
     P_forcing(t,3) = exp(S_forcing(3) * P_forcing(t,3)) !Precipitation
    enddo ! end TIME LOOP
  
    ! Adding perturbations to the forcing data at each time step for the ensembles.
    do t = 1,Nt
  
     ! PRECIPITATION RATE
     !jmframe: I couldn't figure out why, but a few perturbed precipitation values
     !were very low magnitude negative values. They shouldn't be, based on the
     !math, but I added in this max() logic just in case.
     forcing(t,e)%prcprate = forcing(t,1)%prcprate * P_forcing(t,3)
     if (forcing(t,e)%prcprate.lt.0.0) then
      forcing(t,e)%prcprate = max(forcing(t,1)%prcprate, 0.0)
     endif
  
     ! Perturb the SURFACE TEMPERATURE
     forcing(t,e)%sfctmp = forcing(t,1)%sfctmp + P_forcing(t,2)
  
     ! RADIATION, combine the SW + LW before perturbing, then split back up.
     rn = forcing(t,1)%lwrad + forcing(t,1)%swrad
     rat = forcing(t,1)%lwrad/rn
     rnp = rn * P_forcing(t,1)
     diff = rnp - rn
     forcing(t,e)%lwrad = forcing(t,1)%lwrad + diff * rat
     forcing(t,e)%swrad = forcing(t,1)%swrad + diff * (1-rat)
     if (forcing(t,e)%lwrad.lt.0.0.or.forcing(t,e)%swrad.lt.0.0) then
      forcing(t,e)%lwrad = max(forcing(t,1)%lwrad, 0.0)
      forcing(t,e)%swrad = max(forcing(t,1)%swrad, 0.0)
     endif
  
    enddo ! end TIME LOOP
    deallocate(P_forcing)
  
    ! Set the forcing's that are not perturbed
    ! This seperate loop is uncessessary, but it will keep things organized
    do t = 1,Nt
     forcing(t,e)%sfcspd = forcing(t,1)%sfcspd
     forcing(t,e)%sfcprs = forcing(t,1)%sfcprs
     forcing(t,e)%q2 = forcing(t,1)%q2
    enddo
  
   enddo ! END Ensemble loop, Done perturbing forcing data
 endif ! Logic for running through ensembles greater than 1

! prescribed shade fraction
 forcing%shdfac = -9999.
 inquire(file='shdfac.txt',exist=fexists)
 if ((setup%dveg.eq.1).and.(fexists)) then
   do e = 1,Ne
    fname = 'shdfac.txt'
    open(fid,file=trim(fname))
      do t = 1,Nt
        read(fid,*) dummy,dummy,dummy,dummy,forcing(t,e)%shdfac
      enddo ! times
    close(fid)
   enddo ! ensembles
 endif
 
 ! parameters
 ! This standalone version does not use the parameter tables.
 ! The one place where this may cause a problem is on...
 ! line 8866 of module_sf_noahmplsm.f90 where carbon partitioning 
 ! to the leaf is different for elbforest than for other vegetation types. 
 ! We have set a constant vegetation type so that isurban, 
 ! iswater, issnow, and isbaren are not triggered.

 do e = 1,Ne
  call redprm(setup%nsoil,setup%tbot,vegtyp)
  setup%vegtyp = vegtyp ! this should !not! be a parameter 
 enddo

! --- Run the Model -----------------------------------------------------

 ! initial timestep
 t = 1
 do e = 1,Ne
   call driver(t,setup,forcing(t,e),state(t,e),output(t,e))
   if (data_assim) then
     background(t,e) = state(t,e)
   endif 
 enddo

 ! time loop
 lagged(1)  = state(1,1)%smc(1) !store gpr prediction before cut
 lagged(2)  = state(1,1)%smc(1) !store gpr prediction before cut
 lagged(3)  = state(1,1)%smc(1) !store gpr prediction before cut
 ! time loop
 do t = 2,Nt

   ! Ensemble loop
   do e = 1,Ne
     !store gpr lagged values before cut
     lagged(1) = lagged(2)
     lagged(2) = lagged(3)
     lagged(3)  = state(t-1,e)%smc(1)

     ! timestep
     state(t,e) = state(t-1,e)
     !JF:
     ! I finally figured out the trouble with the jumps in the Noah results when running the "onestep" 
     ! state update simulations. There are two states that need to be updated with the soil moisture observation.
     ! The SMC (liquid water) and SH2O (Liquid + Ice water). So below I am updating the liquid water directly
     ! with the observation, and the sH2O state with the observation plus the previous ice content. 
     if ((isOneStep).and.(obs(t-1)>0)) then
      ! Soil moisture states 2-4 should change in proportion to SMC1
!      state(t,e)%smc(2) = state(t-1,e)%smc(2) + (obs(t-1) - state(t,e)%smc(1)) * (state(t,e)%smc(2)/state(t,e)%smc(1)) 
!      state(t,e)%smc(3) = state(t-1,e)%smc(3) + (obs(t-1) - state(t,e)%smc(1)) * (state(t,e)%smc(3)/state(t,e)%smc(1)) 
!      state(t,e)%smc(4) = state(t-1,e)%smc(4) + (obs(t-1) - state(t,e)%smc(1)) * (state(t,e)%smc(4)/state(t,e)%smc(1)) 
      state(t,e)%smc(1) = obs(t-1)
      state(t,e)%sh2o(1) = obs(t-1) + (state(t-1,e)%sh2o(1) - state(t-1,e)%smc(1))

     ! add random perturbation
     !jmframe: why perturb the states of ensemble 1? added logic to prevent.
     if ((perturb).and.(e>1)) then 
     !if (perturb) then
       call perturb_state(setup,state(t,e),sig_sm)
     endif

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!! CRAIG moved this "CUT" up here after the one st`ep update !!!!!
     !!!!!! Still done below in Data Assimilation                 !!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (isOneStep) then
      if (state(t,e)%lfmass.le.50/SLA)     &
                       state(t,e)%lfmass = 50/SLA+0.01
      if (state(t,e)%lfmass.ge.5000/SLA)   &
                       state(t,e)%lfmass = 5000/SLA
      if (state(t,e)%stmass.le.0.05/0.003) &
                       state(t,e)%stmass = 0.05/0.003+0.01
      if (state(t,e)%rtmass.le.5)          &
                       state(t,e)%rtmass = 5.01 
      state(t,e)%lai = MAX(state(t,e)%lfmass*SLA/1000,0.05)
      state(t,e)%sai = MAX(state(t,e)%stmass*0.003,0.05)
      do d = 1,setup%nsoil
        if (state(t,e)%smc(d).gt.smcmax)   &
                       state(t,e)%smc(d) = smcmax
        if (state(t,e)%smc(d).lt.smcdry)     &
                       state(t,e)%smc(d) = smcdry
        if (state(t,e)%sh2o(d).gt.smcmax)  &
                       state(t,e)%sh2o(d) = smcmax
        if (state(t,e)%sh2o(d).lt.smcdry)    &
                       state(t,e)%sh2o(d) = smcdry
      enddo ! soil dimension
     endif

      if (t>2732) then
       print*, "time step", t
       print*, "obs(t-1)", obs(t-1)
       print*, "(state(t-1,e)%smc(1)", (state(t-1,e)%smc(1))
       print*, "(state(t-1,e)%sh2o(1)", (state(t-1,e)%sh2o(1))
       print*, "obs(t)", obs(t)
       print*, "forcing(t,1)%lwrad", forcing(t,1)%lwrad
       print*, "forcing(t,1)%swrad", forcing(t,1)%swrad
       print*, "forcing(t,1)%sfctmp", forcing(t,1)%sfctmp
       print*, "forcing(t,1)%sfcspd", forcing(t,1)%sfcspd
       print*, "forcing(t,1)%sfcprs", forcing(t,1)%sfcprs
       print*, "forcing(t,1)%q2",     forcing(t,1)%q2 
       print*, "---------- States Before driver"
       print*, "state(t,e)%smc(1)", state(t,e)%smc(1)
       print*, "state(t,e)%smc(2)", state(t,e)%smc(2)
       print*, "state(t,e)%smc(3)", state(t,e)%smc(3)
       print*, "state(t,e)%smc(4)", state(t,e)%smc(4)
       print*, "state(t,e)%sh2o(1)", state(t,e)%sh2o(1)
       print*, "state(t,e)%sh2o(2)", state(t,e)%sh2o(2)
       print*, "state(t,e)%sh2o(3)", state(t,e)%sh2o(3)
       print*, "state(t,e)%sh2o(4)", state(t,e)%sh2o(4)
      endif
     endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!! run model at timestep !!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call driver(t,setup,forcing(t,e),state(t,e),output(t,e))
      if (t>2732) then
       print*, "------- States After driver"
       print*, "state(t,e)%smc(1)", state(t,e)%smc(1)
       print*, "state(t,e)%smc(2)", state(t,e)%smc(2)
       print*, "state(t,e)%smc(3)", state(t,e)%smc(3)
       print*, "state(t,e)%smc(4)", state(t,e)%smc(4)
       print*, "state(t,e)%sh2o(1)", state(t,e)%sh2o(1)
       print*, "state(t,e)%sh2o(2)", state(t,e)%sh2o(2)
       print*, "state(t,e)%sh2o(3)", state(t,e)%sh2o(3)
       print*, "state(t,e)%sh2o(4)", state(t,e)%sh2o(4)
       print*, "OUTPUT VALUES:"
       print*, "output(t,1)%qe", output(t,1)%qe
       print*, "output(t,1)%qh", output(t,1)%qh
       print*, "output(t,1)%nee",output(t,1)%nee 
       print*, "---------------------------------------"
      endif

     ! Save background State
     if (data_assim) then
       background(t,e) = state(t,e)
     endif 

     ! error check in dynamic veg state 
     if (isnan(state(t,e)%stmass).or.isnan(state(t,e)%lfmass)) then
       if (fixensemble == .false.) then
        print*, 'WARNING:'
        print*, 'fixing at least one ensemble members'
        fixensemble = .true.
       endif
       fixed = .false.
       do ee = e-1,1,-1
         if ((.not.(fixed)).and.(.not.(isnan(state(t,ee)%stmass))) &
              .and.(.not.(isnan(state(t,ee)%lfmass)))) then
           state(t,e) = state(t,ee)
           output(t,e) = output(t,ee)
           background(t,e) = background(t,ee)
           fixed = .true.
         endif
       enddo
       if (.not.(fixed)) stop 'Error in Noah-MP veg state'
     endif

   enddo ! ensemble
  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !!!!!!!!! CALL DATA ASSIMILATION !!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ((data_assim).and.(obs(t).gt.0)) then
     lags = min(Nlag,t-1)
 
     ! soil moisture state updating
     ! X( lags, Soil Moisture Layers, No. Ensembles)
     allocate(X(lags,4,Ne))
     do l = 1,lags
       do e = 1,Ne
         X(l,:,e) = state(t-l+1,e)%smc(1:4)
       enddo
     enddo

     ! Y(DZ,Ne) etc... 
     ! Fix later on. take out the hard coded DZ values
     !jmframe: might as well fix now. Only dealing with one (SM1)
     allocate(Y(Dz,Ne))
     allocate(Z(Dz))
     allocate(R(Dz))
     do e = 1,Ne
       Y(Dz,e) = state(t,e)%smc(1)
      enddo
     Z = obs(t)
     R = zcov(Dz) ! read in as observation covariance from obs_cov.txt 

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!  Ensemle Kalman Smother !!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!     subroutine enks(X,Y,Zbar,Zsig,Nt,Ne,Dx,Dz)   !!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call enks(X,Y,Z,R,lags,Ne,4,Dz) 
     deallocate(Y)
     deallocate(Z)
     deallocate(R)

     ! This updates the model states based on the X result from enks() 
     do l = 1,lags
        do e = 1,Ne
         !Updating all 4 soil moisture states based on the lag.
         do d = 1,4
           state(t-l+1,e)%sh2o(d) = &
               state(t-l+1,e)%sh2o(d)+X(l,d,e)-state(t-l+1,e)%smc(d)
             state(t-l+1,e)%smc(d) = X(l,d,e)
           if (state(t-l+1,e)%smc(d).gt.smcmax)  &
                               state(t-l+1,e)%smc(d) = smcmax
           if (state(t-l+1,e)%smc(d).lt.smcdry)    &
                               state(t-l+1,e)%smc(d) = smcdry
           if (state(t-l+1,e)%sh2o(d).gt.smcmax) &
                                 state(t-l+1,e)%sh2o(d) = smcmax
           if (state(t-l+1,e)%sh2o(d).lt.smcdry)   &
                               state(t-l+1,e)%sh2o(d) = smcdry
         enddo ! sm dimension
       enddo ! ensemble
     enddo ! lag
      
     deallocate(X)
   endif ! DA + SM obs present 

   if (isOneStep == .false.) then
    do e = 1,Ne 
      if (state(t,e)%lfmass.le.50/SLA)     &
                       state(t,e)%lfmass = 50/SLA+0.01
      if (state(t,e)%lfmass.ge.5000/SLA)   &
                       state(t,e)%lfmass = 5000/SLA
      if (state(t,e)%stmass.le.0.05/0.003) &
                       state(t,e)%stmass = 0.05/0.003+0.01
      if (state(t,e)%rtmass.le.5)          &
                       state(t,e)%rtmass = 5.01 
      state(t,e)%lai = MAX(state(t,e)%lfmass*SLA/1000,0.05)
      state(t,e)%sai = MAX(state(t,e)%stmass*0.003,0.05)
      do d = 1,setup%nsoil
        if (state(t,e)%smc(d).gt.smcmax)   &
                       state(t,e)%smc(d) = smcmax
        if (state(t,e)%smc(d).lt.smcdry)     &
                       state(t,e)%smc(d) = smcdry
        if (state(t,e)%sh2o(d).gt.smcmax)  &
                       state(t,e)%sh2o(d) = smcmax
        if (state(t,e)%sh2o(d).lt.smcdry)    &
                       state(t,e)%sh2o(d) = smcdry
      enddo ! soil dimension
    enddo ! ensemble for bounds checking
   endif
 
 enddo ! time loop

! ------- Write Output --------------------------------------------
 if (write_model_output) then

  ! Write the main output file reguardless of the other options.
  fname = 'output.out'
  open(fid,file=trim(fname),status='replace')
  do t = 1,Nt
   write(fid,'( i7,i5,f7.3,                                          &
                f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                f17.6,f17.6,f17.6,f17.6,f17.6,f17.6)')               & 
    date(t,:),time(t), &
    forcing(t,1)%sfctmp, forcing(t,1)%sfcspd, forcing(t,1)%sfcprs,   &
    forcing(t,1)%q2, forcing(t,1)%lwrad, forcing(t,1)%swrad,         &
    forcing(t,1)%prcprate, state(t,1)%smc(1:4),                      &
    state(t,1)%rtmass, state(t,1)%wood,                              &
    state(t,1)%lfmass, state(t,1)%stmass,                            &
    output(t,1)%qe, output(t,1)%qh, output(t,1)%nee
  enddo
  close(fid)

  ! Now go through the DA/Perturbation options.
  if (perturb) then
   if (data_assim) then
    fname = 'enks_mean.out'
   else
    fname = 'open_mean.out'
   endif ! If data_assim
   open(fid,file=trim(fname),status='replace')

   do t = 1,Nt 
    mean_sfctmp   = sum(forcing(t,:)%sfctmp)   /Ne
    mean_sfcspd   = sum(forcing(t,:)%sfcspd)   /Ne
    mean_sfcprs   = sum(forcing(t,:)%sfcprs)   /Ne
    mean_q2       = sum(forcing(t,:)%q2)       /Ne
    mean_lwrad    = sum(forcing(t,:)%lwrad)    /Ne
    mean_swrad    = sum(forcing(t,:)%swrad)    /Ne
    mean_prcprate = sum(forcing(t,:)%prcprate) /Ne

    mean_smc1     = 0
    mean_smc2     = 0
    mean_smc3     = 0
    mean_smc4     = 0
    do e = 1,Ne
     mean_smc1     = mean_smc1 + state(t,e)%smc(1)
     mean_smc2     = mean_smc1 + state(t,e)%smc(2)
     mean_smc3     = mean_smc1 + state(t,e)%smc(3)
     mean_smc4     = mean_smc1 + state(t,e)%smc(4)
    enddo
    mean_smc1     = mean_smc1/Ne
    mean_smc2     = mean_smc2/Ne
    mean_smc3     = mean_smc3/Ne
    mean_smc4     = mean_smc4/Ne

    mean_wood     = sum(state(t,:)%wood)       /Ne
    mean_lfmass   = sum(state(t,:)%lfmass)     /Ne
    mean_stmass   = sum(state(t,:)%stmass)     /Ne
    mean_rtmass   = sum(state(t,:)%rtmass)     /Ne

    mean_qe       = sum(output(t,:)%qe)        /Ne
    mean_qh       = sum(output(t,:)%qh)        /Ne
    mean_nee      = sum(output(t,:)%nee)       /Ne
    write(fid,'(                                 &
                 f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                 f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                 f17.6,f17.6,f17.6,f17.6,f17.6,f17.6)')               & 
     mean_sfctmp,mean_sfcspd,mean_sfcprs,mean_q2,mean_lwrad,          &
     mean_swrad,mean_prcprate,                                        &
     mean_smc1,mean_smc2,mean_smc3,mean_smc4,                         &
     mean_rtmass,mean_wood,mean_lfmass,mean_stmass,                   &
     mean_qe,mean_qh,mean_nee
   enddo
   close(fid)
  endif ! If Perturb

  if (data_assim) then
   fname = 'back_mean.out'
   open(fid,file=trim(fname),status='replace')
  
   do t = 1,Nt 
    mean_sfctmp   = sum(forcing(t,:)%sfctmp)   /Ne
    mean_sfcspd   = sum(forcing(t,:)%sfcspd)   /Ne
    mean_sfcprs   = sum(forcing(t,:)%sfcprs)   /Ne
    mean_q2       = sum(forcing(t,:)%q2)       /Ne
    mean_lwrad    = sum(forcing(t,:)%lwrad)    /Ne
    mean_swrad    = sum(forcing(t,:)%swrad)    /Ne
    mean_prcprate = sum(forcing(t,:)%prcprate) /Ne
    
    mean_smc1     = 0
    mean_smc2     = 0
    mean_smc3     = 0
    mean_smc4     = 0
    do e = 1,Ne
     mean_smc1     = mean_smc1 + background(t,e)%smc(1)
     mean_smc2     = mean_smc1 + background(t,e)%smc(2)
     mean_smc3     = mean_smc1 + background(t,e)%smc(3)
     mean_smc4     = mean_smc1 + background(t,e)%smc(4)
    enddo
    mean_smc1     = mean_smc1/Ne
    mean_smc2     = mean_smc2/Ne
    mean_smc3     = mean_smc3/Ne
    mean_smc4     = mean_smc4/Ne

    mean_wood     = sum(background(t,:)%wood)       /Ne
    mean_lfmass   = sum(background(t,:)%lfmass)     /Ne
    mean_stmass   = sum(background(t,:)%stmass)     /Ne
    mean_rtmass   = sum(background(t,:)%rtmass)     /Ne

    mean_qe       = sum(output(t,:)%qe)        /Ne
    mean_qh       = sum(output(t,:)%qh)        /Ne
    mean_nee      = sum(output(t,:)%nee)       /Ne
    write(fid,'(                                &
                 f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                 f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                 f17.6,f17.6,f17.6,f17.6,f17.6,f17.6)')               & 
     mean_sfctmp,mean_sfcspd,mean_sfcprs,mean_q2,mean_lwrad,          &
     mean_swrad,mean_prcprate,                                        &
     mean_smc1,mean_smc2,mean_smc3,mean_smc4,                         &
     mean_rtmass,mean_wood,mean_lfmass,mean_stmass,                   &
     mean_qe,mean_qh,mean_nee
   enddo
   close(fid)
  endif ! If data_assim
 
  !ensemble output files
  if (write_ensemble_output) then
   do e = 1,min(9,Ne)
    
    if (perturb) then
     if (data_assim) then
      fname = 'enks_'
     else
      fname = 'open_'
     endif
     if (e.lt.10) then
      write(es1,'(i1)') e
      fname = trim(fname)//es1
     elseif (e.lt.100) then
      write(es2,'(i2)') e
      fname = trim(fname)//es2
     elseif (e.lt.1000) then
      write(es3,'(i3)') e
      fname = trim(fname)//es3
     endif 
     fname = trim(fname)//'.out'
    endif
    open(fid,file=trim(fname),status='replace')
  
    do t = 1,Nt
      write(fid,'( i7,i5,f7.3,                                          &
                   f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                   f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                   f17.6,f17.6,f17.6,f17.6,f17.6,f17.6)')               & 
       date(t,:),time(t), &
       forcing(t,e)%sfctmp, forcing(t,e)%sfcspd, forcing(t,e)%sfcprs,   &
       forcing(t,e)%q2, forcing(t,e)%lwrad, forcing(t,e)%swrad,         &
       forcing(t,e)%prcprate, state(t,e)%smc(1:4),                      &
       state(t,e)%rtmass, state(t,e)%wood,                              &
       state(t,e)%lfmass, state(t,e)%stmass,                            &
       output(t,e)%qe, output(t,e)%qh, output(t,e)%nee
    enddo
    close(fid)
  
    if (data_assim) then
     fname = 'back_'
     if (e.lt.10) then
      write(es1,'(i1)') e
      fname = trim(fname)//es1
     elseif (e.lt.100) then
      write(es2,'(i2)') e
      fname = trim(fname)//es2
     elseif (e.lt.1000) then
      write(es3,'(i3)') e
      fname = trim(fname)//es3
     endif 
     fname = trim(fname)//'.out'
     open(fid,file=trim(fname),status='replace')
  
     do t = 1,Nt
      write(fid,'( i7,i5,f7.3,                                          &
                   f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                   f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6,f17.6      & 
                   f17.6,f17.6,f17.6,f17.6,f17.6,f17.6)')               & 
       date(t,:),time(t), &
       forcing(t,e)%sfctmp, forcing(t,e)%sfcspd, forcing(t,e)%sfcprs,   &
       forcing(t,e)%q2, forcing(t,e)%lwrad, forcing(t,e)%swrad,         &
       forcing(t,e)%prcprate, state(t,e)%smc(1:4),                      &
       state(t,e)%rtmass, state(t,e)%wood,                              &
       state(t,e)%lfmass, state(t,e)%stmass,                            &
       output(t,e)%qe, output(t,e)%qh, output(t,e)%nee
     enddo
     close(fid)
    endif
   enddo ! open ensemble output files
  endif ! if write ensemble output
 endif ! if write model output

 if (calc_obj_fun) then
  call rmse_obj_fun(date,time,state(:,1),output(:,1),Nt)
 endif
 print*, 'The program has finished'
! -----------------------------------------------------------------------
end program

! ----------------------------------------------------------------
! ----------------------------------------------------------------
subroutine redprm(nsoil,tbot,vegtyp)
 use noahmp_globals
 use noahmp_veg_parameters
 implicit none

 ! inputs 
 integer, intent(in)    :: nsoil

 ! Locals
 integer, parameter :: fid = 134
 real    :: refdk
 real    :: refkdt
 real    :: frzk
 real    :: frzfact
 real    :: tbot
 character(1000) :: fname
 character(1) :: es1
 character(2) :: es2
 integer :: vegtyp

! ----Read in Paramter File----------------------------------------------

 fname = 'cal_parms.txt'
 open(fid,file=trim(fname),action='read')
  read(fid,*) Z0MVT
  read(fid,*) HVT
  read(fid,*) HVB
  read(fid,*) LTOVRC
  read(fid,*) DILEFW
  read(fid,*) RMF25
  read(fid,*) SLA
  read(fid,*) VCMX25
  read(fid,*) QE25
  read(fid,*) BEXP
  read(fid,*) DKSAT
! jmframe June 2019: taking out SMCDRY from the calibration,
! Since it is hard coded in for data assimilation
!  read(fid,*) SMCDRY
  read(fid,*) SMCMAX
  read(fid,*) SMCREF
  read(fid,*) SMCWLT
 close(fid)

 ! turn conductivity into real-valued space
 DKSAT = 10**DKSAT

 ! jmframe June 2019: Moved fractional soil moisture calcs below... 
 ! the rest of the parm inputs.

 ! ensure that canopy top is not lower than canopy bottom
 !HVB = 0.1 + HVB*(HVT-0.1)
 HVB = HVB*HVT

 fname = 'parms.txt'
 open(fid,file=trim(fname),action='read')

 ! veg parms
  read(fid,*) CH2OP
  read(fid,*) DLEAF
  !read(fid,*) Z0MVT
  !read(fid,*) HVT
  !read(fid,*) HVB
  read(fid,*) RC
  read(fid,*) RHOL(1)
  read(fid,*) RHOL(2)
  read(fid,*) RHOS(1)
  read(fid,*) RHOS(2)
  read(fid,*) TAUL(1)
  read(fid,*) TAUL(2)
  read(fid,*) TAUS(1)
  read(fid,*) TAUS(2)
  read(fid,*) XL
  read(fid,*) CWPVT
  read(fid,*) C3PSN
  read(fid,*) KC25
  read(fid,*) AKC
  read(fid,*) KO25
  read(fid,*) AKO
  read(fid,*) AVCMX
  !read(fid,*) LTOVRC
  read(fid,*) DILEFC
  !read(fid,*) DILEFW
  !read(fid,*) RMF25
  !read(fid,*) SLA
  read(fid,*) FRAGR
  read(fid,*) TMIN
  !read(fid,*) VCMX25
  read(fid,*) TDLEF
  read(fid,*) BP
  read(fid,*) MP
  !read(fid,*) QE25
  read(fid,*) RMS25
  read(fid,*) RMR25
  read(fid,*) ARM
  read(fid,*) FOLNMX
  read(fid,*) WDPOOL
  read(fid,*) WRRAT
  read(fid,*) MRP
  SAIM = 0.
  LAIM = 0.
  !read(fid,*) SAIM
  !read(fid,*) LAIM
  read(fid,*) SLAREA
  !read(fid,*) EPS
  read(fid,*) VEGTYP

 ! gen parms
  read(fid,*) csoil
  !read(fid,*) bexp
  !read(fid,*) dksat
  read(fid,*) dwsat
  read(fid,*) f1
  read(fid,*) psisat
  read(fid,*) quartz
  !jmframe: Not calibrating. Hardwired in as 0.02
  ! Added into parms.txt file as line 40
  read(fid,*) smcdry
  !read(fid,*) smcmax
  !read(fid,*) smcref
  !read(fid,*) smcwlt

  read(fid,*) zbot      ! 55
  read(fid,*) czil
  read(fid,*) frzk
  read(fid,*) refdk
  read(fid,*) refkdt
  read(fid,*) slope     ! 60
  read(fid,*) topt      ! 61
  read(fid,*) rgl       ! 62
  read(fid,*) rsmax     ! 63
  read(fid,*) rsmin     ! 64
  read(fid,*) hs        ! 65
  read(fid,*) nroot     
 close(fid)

 open(fid,file='tbot.txt',action='read')
  read(fid,*) tbot ! 67
 close(fid)

 open(fid,file='time_parms.txt',action='read')
  read(fid,*) LAIM   
  read(fid,*) SAIM    
!  read(fid,*) EPS
 close(fid)

 ! jmframe June 2019: moved from above, to take SMCDRY out of calibration.
 ! expand soil parameters as fractional intervals
 smcmax = smcdry + smcmax
 smcref = smcdry + (smcmax-smcdry)*smcref
 smcwlt = smcdry + (smcref-smcdry)*smcwlt

 ! some basic manipulations
 kdt = refkdt * dksat / refdk
 frzfact = (smcmax / smcref) * (0.412 / 0.468)
 frzx = frzk * frzfact

 ! error check on rooting layers
 if (nroot.gt.nsoil) nroot = nsoil

end subroutine redprm

subroutine sim_init(setup,state,perturb,data_assim,Ntimes,Nens,sig_sm,Nlag,isOneStep)
 use type_decs

 integer, parameter :: fid = 14
 type(state_data)   :: state
 type(setup_data)   :: setup
 logical, intent(out) :: perturb,data_assim,isOneStep
 integer, intent(out) :: Nens
 integer, intent(out) :: Ntimes
 !jmframe added sig_sm & Nlag here, instead of hard coded. July 2019
 integer, intent(out) :: sig_sm
 integer, intent(out) :: Nlag
 integer :: da_flag
 logical :: fexists

! simulation setup
 open(fid,file='init.txt',action='read')
  read(fid,*) setup%nsoil     
  read(fid,*) setup%nsnow     
 
! allocate dimensions
  allocate(state%stc(-setup%nsnow+1:setup%nsoil))
  allocate(state%zsnso(-setup%nsnow+1:setup%nsoil))
  allocate(state%tsno(setup%nsnow))
  allocate(state%snice(setup%nsnow))
  allocate(state%snliq(setup%nsnow))
  allocate(state%sh2o(setup%nsoil))
  allocate(state%smc(setup%nsoil))
  allocate(setup%sldpth(setup%nsoil))

! setup parameters
  read(fid,*) setup%zlvl      
  read(fid,*) setup%dt        
  read(fid,*) setup%opt_crs 
  read(fid,*) setup%opt_btr 
  read(fid,*) setup%opt_run 
  read(fid,*) setup%opt_sfc 
  read(fid,*) setup%opt_frz 
  read(fid,*) setup%opt_inf 
  read(fid,*) setup%opt_rad 
  read(fid,*) setup%opt_alb 
  read(fid,*) setup%opt_snf 
  read(fid,*) setup%opt_tbot 
  read(fid,*) setup%opt_stc 
  read(fid,*) setup%dveg     
  read(fid,*) setup%sldpth    

! initial state
  read(fid,*) state%stc(1:setup%nsoil) 
  read(fid,*) state%snowh   
  read(fid,*) state%sneqv   
  read(fid,*) state%canliq  
  read(fid,*) state%rtmass  
  read(fid,*) state%albold  
  read(fid,*) state%lai     
  read(fid,*) state%tv      
  read(fid,*) state%tg      
 close(fid)

 open(fid,file='startdate.txt')
  read(fid,*) setup%startdate
 close(fid)

 open(fid,file='da_flag.txt')
  read(fid,*) da_flag
 close(fid)
 if (da_flag.gt.0) then
  perturb = .true.
  data_assim = .true.
  isOneStep = .false.
  Nens = da_flag
 elseif (da_flag.lt.-1) then
  perturb = .true.
  data_assim = .false.
  isOneStep = .false.
  Nens = -da_flag
 elseif (da_flag.eq.0) then
  perturb = .false.
  data_assim = .false.
  isOneStep = .false.
  Nens = 1
elseif (da_flag.eq.-1) then
  perturb = .false.
  data_assim = .false.
  isOneStep = .true.
  Nens = 1
 else
  stop 9813
 endif

 open(fid,file='num_times.txt')
  read(fid,*) ntimes
 close(fid)

 ! values needed for data assimilation
 if (data_assim) then
  open(fid,file='sig_sm.txt')
   read(fid,*) sig_sm
  close(fid)
 
  open(fid,file='Nlag.txt')
   read(fid,*) Nlag
  close(fid)
 endif !data_assim

 open(fid,file='lat_lon.txt')
  read(fid,*) setup%latitude
  read(fid,*) setup%longitude
 close(fid)

 open(fid,file='plant_init.txt')
  read(fid,*) state%rtmass
  read(fid,*) state%wood
  read(fid,*) state%lfmass
  read(fid,*) state%stmass
 close(fid)

 open(fid,file='soil_init.txt')
  read(fid,*) state%smc(1)
  read(fid,*) state%smc(2)
  read(fid,*) state%smc(3)
  read(fid,*) state%smc(4)
  state%sh2o = state%smc
 close(fid)

 inquire(file='shdfac.txt',exist=fexists)
 if (fexists) then
   open(fid,file='shdfac.txt')
     read(fid,*) setup%shdfac_monthly
   close(fid)
 else
  setup%shdfac_monthly = (/0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98/)
 endif

end subroutine sim_init

subroutine dealoc(setup,state)
 use type_decs

 type(state_data)   :: state
 type(setup_data)   :: setup

 deallocate(state%stc)
 deallocate(state%zsnso)
 deallocate(state%tsno)
 deallocate(state%snice)
 deallocate(state%snliq)
 deallocate(state%sh2o)
 deallocate(state%smc)
 deallocate(setup%sldpth)

end subroutine dealoc

subroutine perturb_state(setup,state,sig_sm)
 use type_decs
 use noahmp_veg_parameters
 use noahmp_globals
 implicit none

 type(state_data), intent(inout) :: state
 type(setup_data), intent(in)    :: setup
 integer, intent(in)    :: sig_sm
 ! jmframe: changed to read in value.
 !real, parameter :: sig_sm = 0.005
! real, parameter :: sig_veg = 0.01
 integer, parameter :: N = 1
 real :: eta
 real :: random_normal
 integer :: d

 ! soil moisture
 do d = 1,setup%nsoil
  eta = random_normal()
  state%smc(d) = state%smc(d)+eta*sig_sm
  state%sh2o(d) = state%sh2o(d)+eta*sig_sm
 enddo

 return
end subroutine perturb_state

function random_normal() result(fn_val)
 implicit none
 real     :: half = 0.5
 real     :: fn_val
 real     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,    &
             r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
 integer :: i, n, clock

 do
   call random_number(u)
   call random_number(v)
   v = 1.7156 * (v - half)

   x = u - s
   y = ABS(v) - t
   q = x**2 + y*(a*y - b*x)

   ! Accept P if inside inner ellipse
   if (q < r1) EXIT
   ! Reject P if outside outer ellipse
   if (q > r2) CYCLE
   ! Reject P if outside acceptance region
   if (v**2 < -4.0*LOG(u)*u**2) EXIT
 enddo

 ! Return ratio of P's coordinates as the normal deviate
 fn_val = v/u

 return
end function random_normal
