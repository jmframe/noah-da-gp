program run_timestep

  use noahmp_veg_parameters
  use noahmp_globals
  use type_decs
 
 ! -- Variable Declarations ----------------------------------------------
  implicit none
 
  ! for output (formerly, these were parameters)
  logical :: write_model_output
  logical :: write_ensemble_output
  logical :: calc_obj_fun
  logical :: fixensemble
 
  ! SY: for obs data (formerly, this was parameter)
  integer                             :: Dz
 
  ! simulation parameters !SY: & data assimilation parameter also
  !jmframe added sig_sm & Nlag as read in values, not set in code.
  real    :: sig_sm, sig_veg
  integer :: Nlag, EQs
!  integer, allocatable, dimension(:) :: firstObs ! SY: REMEMBER to 
      !uncomment this if states (smc etc) are available as observations
  !Similarly, Soni added Threshold_on_CC
  real    :: Threshold_on_CC
  integer :: initCycle 
  integer :: Ne, Nt ! # of ensemble members, and # of timesteps
  logical :: if_state_perturb, data_assim, doEQ
 
  ! file I/O
  integer, parameter :: fid = 15
  character(1000)    :: fname, fdir, ofile
  character(1)       :: es1
  character(2)       :: es2
  character(3)       :: es3
 
  ! internal indexes
  integer :: t, e, ee, l, d, lags, isoil, tt, zz ! SY: removed s since 
                           !not used, added tt and zz
  real    :: dummy
  integer :: vegtyp
  logical :: fixed, fexists 
  integer, allocatable, dimension(:, :) :: date
  integer, allocatable, dimension(:, :) :: time
 
  ! model data
  type(forcing_data), allocatable, dimension(:,:) :: forcing
  type(state_data),   allocatable, dimension(:,:) :: state, background
  type(setup_data) :: setup
  type(output_data),  allocatable, dimension(:,:) :: output
 
  ! obs data
  real, allocatable, dimension(:,:)   :: obs
  real, allocatable, dimension(:)   :: NEE, Qle, Qh

  ! data assimilation-related
  real, allocatable, dimension(:)     :: zcov
  real, allocatable, dimension(:,:,:) :: X
  real, allocatable, dimension(:,:)   :: Y
  real, allocatable, dimension(:)     :: Z, R
  logical :: IfAssimAtThisTimestep 
  integer :: NumValidObsAtThisTimeStep, WhichValidObsAtThisTimeStep
 
  !CRAIG
  logical :: isOneStep
  !JMFRAME
  real(kind=8), dimension(3) :: lagged
 
  ! Forcing perturbations
  real, allocatable, dimension(:,:) :: P_forcing
  integer :: info
  real :: rn, rat, rnp, diff
  real, dimension(3,3) :: R_forcing
  real, dimension(3) :: S_forcing
 
  ! mean output data
  real ::  mean_sfctmp,mean_sfcspd,mean_sfcprs,mean_q2,mean_lwrad, &
  mean_swrad,mean_prcprate,mean_smc1,mean_smc2,mean_smc3,mean_smc4,&
  mean_wood,mean_lfmass,mean_stmass,mean_rtmass,mean_qe,mean_qh,mean_nee      

  !for very basic info about setup
  integer :: da_flag

  !function
  real :: random_normal
 
 ! --- Set Up Run --------------------------------------------------------
 ! setup simulation 
  inquire(file='da_flag.txt',exist=fexists)
  if (fexists) then
    open(fid, file = 'da_flag.txt')
    read(fid, *) da_flag
    close(fid)
    if (da_flag .gt. 0) then
      if_state_perturb = .true.
      data_assim = .true.
      isOneStep = .false.
      Ne = da_flag
      ofile = 'output.da'
    elseif (da_flag .eq. 0) then
      if_state_perturb = .false.
      data_assim = .false.
      isOneStep = .false.
      Ne = 1
      ofile = 'output.noah'
    elseif (da_flag .eq. -1) then
      if_state_perturb = .false.
      data_assim = .false.
      isOneStep = .true.
      Ne = 1
      ofile = 'output.onestep'
    elseif (da_flag .lt. -1) then
      if_state_perturb = .true.
      data_assim = .false.
      isOneStep = .false. 
      Ne = -da_flag
      ofile = 'output.perturbed'
    else
      stop 9813 ! SY: Actually, not required if above 4 cases exist
    endif
  else
    if_state_perturb = .false.
    data_assim = .false.
    isOneStep = .false.
    Ne = 1
    ofile = 'output.out'
  endif
 
  open(fid, file = 'num_times.txt')
  read(fid, *) Nt
  close(fid)
 
  allocate(state(Nt,Ne))
  allocate(background(Nt,Ne))
  allocate(output(Nt,Ne))

  call sim_init(setup, state(1, 1), data_assim, sig_sm, sig_veg, Nlag, Threshold_on_CC, doEQ, EQs, Dz, write_model_output, write_ensemble_output, calc_obj_fun, fixensemble)

  do e = 2, Ne

    ! allocate state component dimensions
    allocate( state(1, e)%stc( -setup%nsnow+1: setup%nsoil ) )
    allocate( state(1, e)%zsnso( -setup%nsnow+1: setup%nsoil ) )
    allocate( state(1, e)%tsno( setup%nsnow ) )
    allocate( state(1, e)%snice( setup%nsnow ) )
    allocate( state(1, e)%snliq( setup%nsnow ) )
    allocate( state(1, e)%sh2o( setup%nsoil ) )
    allocate( state(1, e)%smc( setup%nsoil ) )

    ! initial states
    do isoil = 1, setup%nsoil
      state(1, e)%stc(isoil) = state(1, 1)%stc(isoil)
    enddo
    state(1, e)%snowh = state(1, 1)%snowh
    state(1, e)%sneqv = state(1, 1)%sneqv
    state(1, e)%canliq = state(1, 1)%canliq
    state(1, e)%rtmass = state(1, 1)%rtmass !SY: NOTE that this is 
        ! repeated in initial plant-related carbon states section below
    state(1, e)%albold = state(1, 1)%albold
    state(1, e)%lai = state(1, 1)%lai
    state(1, e)%tv = state(1, 1)%tv
    state(1, e)%tg = state(1, 1)%tg

    ! initial plant-related carbon states
    !state(1, e)%rtmass = state(1, 1)%rtmass !SY: NOTE that I'm 
        !commenting this out for now since this has already been done 
        ! in initial states section above
    state(1, e)%wood = state(1, 1)%wood
    state(1, e)%lfmass = state(1, 1)%lfmass
    state(1, e)%stmass = state(1, 1)%stmass

    ! initial soil-related moisture states
    do isoil = 1, 4
      state(1, e)%smc(isoil) = state(1, 1)%smc(isoil)
    enddo
    state(1, e)%sh2o = state(1, 1)%sh2o

  enddo ! do e = 2, Ne

  do e = 1,Ne
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
  enddo ! do e = 1,Ne

  if (Dz .ne. 3) then
    stop 'Num_ObsVariables.txt has a different number than next section below!'
  endif

 ! --- Load Observations -------------------------------------------------
  ! We need to load the observation for one step state update and data assimilation
  ! But there may be other instances where we'd want the observation loaded as well
  ! If so, will remove later.
  if ((data_assim).or.(isOneStep).or.(doEQ).or.(calc_obj_fun)) then

    ! jmframe July 2019, moved this up before the read-in from the observations file.
    ! It is kind of strange that it was below in the first place. Not sure why...
    allocate(obs(Nt, Dz))
    allocate(NEE(Nt))
!    allocate(GPP(Nt)) ! SY: Jonathan changed to saying that this isn't used in NoahMP
    allocate(Qle(Nt))
    allocate(Qh(Nt))
    obs = -9999.
    ! Take out of the data assimilation logic, because we use it for other things.
    open(fid,file='obs.txt')
    do t = 1,Nt
      ! jmframe: PLUMBER-2 flux data Ordered from text file make from NetCDF data.
      ! year,month,day,hour,minute,NEE,GPP,Qle,Qh
      ! source: /discover/nobackup/jframe/data/plumber-2-flux-txt/
      read(fid,*) dummy, dummy, dummy, dummy, dummy, NEE(t), dummy, Qle(t), Qh(t)
      obs(t,1) = NEE(t)
!      obs(t,2) = GPP(t) ! SY: Jonathan says that this isn't used in NoahMP 
      obs(t,2) = Qle(t)
      obs(t,3) = Qh(t)
    enddo ! time
    close(fid)

    ! SY: REMEMBER: beginning of section to uncomment if states (smc etc) are 
          !available as observations
!    allocate(firstObs(Dz))
!    do zz = 1,Dz
!      do t = 1,Nt
!        if ( obs(t, zz) .gt. 0 ) then
!          firstObs(zz) = t
!          exit
!        endif
!      enddo
!    enddo
    ! SY: REMEMBER: ending of section to uncomment if states (smc etc) are 
          !available as observations

  endif ! if ((data_assim).or.(isOneStep).or.(doEQ).or.(calc_obj_fun)) then
  ! We only need the covariance for data assimilation, not for one step.
 
  if( (data_assim) .or. (if_state_perturb) ) then
    allocate(zcov(Dz))
    fname = 'obs_cov.txt'
    open(fid,file=trim(fname))
    read(fid,*) zcov
    close(fid)
  endif

 ! forcing from file
  allocate(forcing(Nt,Ne))
  allocate(date(Nt,3))
  allocate(time(Nt,2))
  fdir = './forcing'
 
 !! Maintain Ensemble number 1 without any perturbations to the forcing data
  e = 1
  ! There was logic to name the ensembles with equal of digits, but I only print 9.
 
  ! plumber-r met value order:
  ! year,month,day,hour,minute,Tair,SWdown,LWdown,VPD,Qair,Psurf,Precip,Wind,RH,CO2air,LAI_alternative,LAI,IGBP_veg_long
  fname = trim(fdir)//'.txt'
  open(fid,file=trim(fname))
  do t = 1,Nt
    ! the humidity here is kg/kg, not % and not relative humidity.
    read(fid,*) date(t,:),time(t,:), &
                forcing(t,e)%sfctmp,forcing(t,e)%swrad,forcing(t,e)%lwrad,dummy, &
                forcing(t,e)%q2,forcing(t,e)%sfcprs,forcing(t,e)%prcprate,forcing(t,e)%sfcspd, &
                dummy, dummy, dummy, dummy, dummy
    ! SY: Begin code portion that checks for any non-negative precip etc to separate this 
           !from the capping-to-0 issue during forcing ensembles creation further below 
    if (forcing(t,e)%prcprate .lt. 0.0) then
      print*, 'The forcing file precip value is a negative value of ',forcing(t,e)%prcprate,' at time ',t
      stop 'Forcing file contains negative value/s for precip!!!'
    endif
    if (forcing(t,e)%swrad .lt. 0.0) then
      print*, 'The forcing file swrad value is a negative value of ',forcing(t,e)%swrad,' at time ',t
      stop 'Forcing file contains negative value/s for swrad!!!'
    endif
    if (forcing(t,e)%lwrad .lt. 0.0) then
      print*, 'The forcing file lwrad value is a negative value of ',forcing(t,e)%lwrad,' at time ',t
      stop 'Forcing file contains negative value/s for lwrad!!!'
    endif
    ! SY: End code portion that checks for any non-negative precip etc to separate this 
           !from the capping-to-0 issue during forcing ensembles creation further below 
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
    
        ! RADIATION, combine the SW + LW before perturbing, then split back up.
        forcing(t,e)%lwrad = forcing(t,1)%lwrad * P_forcing(t,1)
        forcing(t,e)%swrad = forcing(t,1)%swrad * P_forcing(t,1)
        !SY: following should be because of a real-valued approximation-to-zero, 
        !    if at all, similar to forcing(t,e)%prcprate below 
        if (forcing(t,e)%lwrad.lt.0.0) then
          forcing(t,e)%lwrad = max(forcing(t,e)%lwrad, 0.0)
        endif
        if (forcing(t,e)%swrad.lt.0.0) then
          forcing(t,e)%swrad = max(forcing(t,e)%swrad, 0.0)
        endif
    
        ! Perturb the SURFACE TEMPERATURE
        forcing(t,e)%sfctmp = forcing(t,1)%sfctmp + P_forcing(t,2)
     
        ! PRECIPITATION RATE
        !jmframe: I couldn't figure out why, but a few perturbed precipitation values
        !were very low magnitude negative values. They shouldn't be, based on the
        !math, but I added in this max() logic just in case.
        !SY: this should be because of a real-valued approximation-to-zero, if at all 
        forcing(t,e)%prcprate = forcing(t,1)%prcprate * P_forcing(t,3)
        if (forcing(t,e)%prcprate.lt.0.0) then
          forcing(t,e)%prcprate = max(forcing(t,e)%prcprate, 0.0)
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
  inquire(file='shdfacAtTimesteps.txt',exist=fexists)
  if ((setup%dveg.eq.1).and.(fexists)) then

    e = 1
    fname = 'shdfacAtTimesteps.txt'
    open(fid,file=trim(fname))
    do t = 1,Nt
      read(fid,*) dummy,dummy,dummy,dummy,forcing(t,e)%shdfac
    enddo ! times
    close(fid)

    do e = 2,Ne

      do t = 1,Nt
        forcing(t,e)%shdfac = forcing(t,1)%shdfac  
      enddo ! do t = 1,Nt

    enddo ! do e = 2,Ne

  endif ! if ((setup%dveg.eq.1).and.(fexists)) then
  
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

  ! equilibriation+run loop
  do initCycle = 1,EQs+1 ! EQs number of equilibriation loops, and 1 run loop

    ! time loop
    do t = 1,Nt

      if ( ( data_assim ) .and. (initCycle .eq. EQs+1) ) then 
        if ( ( t .eq. 1 ) .and. ( initCycle .eq. 1 ) ) then 
          lags = min(Nlag,t)
        else
          lags = Nlag
        endif   
      endif   
 
      ! Ensemble loop
      do e = 1,Ne

        if ( ( t .eq. 1 ) .and. ( initCycle .eq. 1 ) ) then 
        ! initial timestep and initCycle

          call driver(t, setup, forcing(t,e), state(t,e), output(t,e), initCycle)
          if ( ( data_assim ) .and. (initCycle .eq. EQs+1) ) then 
            background(t,e) = state(t,e)
          endif 

          !SY: NOTE: Next 3 lines should  actually be initial conditions set by model, 
          !         not the states at the end of this timestep, but that might be a 
          !         little complicated to do (driver needs extra arguments to bring 
          !         those out here). Plus, this is only for this one timestep, so 
          !         maybe that's why such an implementation is skipped for now
          lagged(1)  = state(1,e)%smc(1) !store gpr prediction before cut
          lagged(2)  = state(1,e)%smc(1) !store gpr prediction before cut
          lagged(3)  = state(1,e)%smc(1) !store gpr prediction before cut

        else

          !store gpr lagged values before cut
          lagged(1) = lagged(2)
          lagged(2) = lagged(3)
          if ( ( t .eq. 1 ) .and. ( initCycle .gt. 1 ) ) then 
            lagged(3)  = state(Nt,e)%smc(1)
          else
            lagged(3)  = state(t-1,e)%smc(1)
          endif 
          
          ! value from last timestep
          if ( ( t .eq. 1 ) .and. ( initCycle .gt. 1 ) ) then 
            state(t,e) = state(Nt,e)
          else
            state(t,e) = state(t-1,e)
          endif 
          
          ! JMFRAME here we are not going to use the observation for the onestep,
          ! We are going to use the output from Data Assimilation
     !     if ((isOneStep).and.(obs(t-1)>0)) then
     !      ! Soil moisture states 2-4 should change in proportion to SMC1
     !      state(t,e)%smc(1) = obs(t-1)
     !      state(t,e)%sh2o(1) = obs(t-1) + (state(t-1,e)%sh2o(1) - state(t-1,e)%smc(1))
!            ! SY: Begin adding the cut lines of code here. They should have been here
                ! also inside this isOneStep decision and before the perturb_state below
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            !!!!!! CRAIG moved this "CUT" up here after the one step update  !!!!!
!            !!!!!! Still done below in Data Assimilation                 !!!!!!!!!
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            if (state(t,e)%lfmass .le. 50/SLA)     &
!                             state(t,e)%lfmass = 50/SLA + 0.01
!            if (state(t,e)%lfmass .gt. 5000/SLA)   &
!                             state(t,e)%lfmass = 5000/SLA
!            if (state(t,e)%stmass .le. 0.05/0.003) &
!                             state(t,e)%stmass = 0.05/0.003 + 0.01
!            if (state(t,e)%rtmass.le.5)          &
!                             state(t,e)%rtmass = 5.01 
!            state(t,e)%lai = MAX(state(t,e)%lfmass*SLA/1000,0.05)
!            state(t,e)%sai = MAX(state(t,e)%stmass*0.003,0.05)
!            do d = 1,setup%nsoil
!              if (state(t,e)%smc(d) .gt. smcmax)   &
!                             state(t,e)%smc(d) = smcmax
!              if (state(t,e)%smc(d) .lt. smcdry)     &
!                             state(t,e)%smc(d) = smcdry
!              if (state(t,e)%sh2o(d) .gt. smcmax)  &
!                             state(t,e)%sh2o(d) = smcmax
!              if (state(t,e)%sh2o(d) .lt. smcdry)    &
!                             state(t,e)%sh2o(d) = smcdry
!            enddo ! soil dimension
          ! SY: End adding the cut lines of code here. They should have been here 
             ! also inside this isOneStep decision and before the perturb_state below
     !     endif
     
          ! add random perturbation
          !jmframe: why perturb the states of ensemble 1? added logic to prevent.
          !SY: another reason to not perturb the states of ensemble 1 might be
          ! to have no error in dynamic veg states of ensemble 1 so that other 
          ! ensembles close to ensemble 1 can get assigned the valid dynamic 
          ! veg states of ensemble 1 in the "error check in dynamic veg state"
          ! section further below
          if ( if_state_perturb .and. ( e .gt. 1) ) then 
          !if (perturb) then

            call perturb_state(setup, state(t,e), sig_sm, sig_veg)

            !SY: Begin a cut here because of the perturb_state just above
            if (state(t,e)%lfmass .le. 50/SLA)     &
                             state(t,e)%lfmass = 50/SLA + 0.01
            if (state(t,e)%lfmass .gt. 5000/SLA)   &
                             state(t,e)%lfmass = 5000/SLA
            if (state(t,e)%stmass .le. 0.05/0.003) &
                             state(t,e)%stmass = 0.05/0.003 + 0.01
            if (state(t,e)%rtmass.le.5)          &
                             state(t,e)%rtmass = 5.01 
            state(t,e)%lai = MAX(state(t,e)%lfmass*SLA/1000,0.05)
            state(t,e)%sai = MAX(state(t,e)%stmass*0.003,0.05)
            do d = 1,setup%nsoil
              if (state(t,e)%smc(d) .gt. smcmax)   &
                             state(t,e)%smc(d) = smcmax
              if (state(t,e)%smc(d) .lt. smcdry)     &
                             state(t,e)%smc(d) = smcdry
              if (state(t,e)%sh2o(d) .gt. smcmax)  &
                             state(t,e)%sh2o(d) = smcmax
              if (state(t,e)%sh2o(d) .lt. smcdry)    &
                             state(t,e)%sh2o(d) = smcdry
            enddo ! soil dimension
            !SY: End  a cut here because of the perturb_state just above

          endif ! if ( if_state_perturb .and. ( e .gt. 1) ) then 
     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!! run model at timestep !!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call driver(t,setup,forcing(t,e),state(t,e),output(t,e), initCycle)
          
          ! Save background State
          if ( ( data_assim ) .and. (initCycle .eq. EQs+1) ) then 
            background(t,e) = state(t,e)
          endif 
     
        endif ! end of if ( ( t .eq. 1 ) .and. ( initCycle .eq. 1 ) ) then 

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
            endif ! if ((.not.(fixed)).and. ....
          enddo ! do ee = e-1,1,-1
          if (.not.(fixed)) stop 'Error in Noah-MP veg state'
        endif
   
      enddo ! ensemble

      ! SY: REMEMBER: beginning of section to uncomment  if states (smc etc) are 
             !available as observations
!      if ( ( t .eq. 1 ) .and. ( initCycle .eq. 1 ) ) then 
!      ! initial timestep and initCycle
!
!        do zz = 1,Dz
!          if (firstObs(zz) .eq. t) then      
!            if (zz .eq. 1) then
!              do e = 1,Ne
!                do d = 1,setup%nsoil
!                  BLAHBLAH state(t,e)%smc(d) = obs(t, zz) ! Copied from what Jonathan had; are all layers supposed to be set to obs, or only top layer?
!                enddo 
!                state(t,e)%sh2o = state(t,e)%smc
!              enddo 
!            elseif (zz .eq. 2) then
!              do e = 1,Ne
!                BLAHBLAH
!              enddo 
!            .
!            .
!            . 
!            endif ! if (zz .eq. 1) then
!          endif ! if (firstObs(zz) .eq. t) then      
!        enddo ! do zz = 1,Dz
!   
!      endif ! end of if ( ( t .eq. 1 ) .and. ( initCycle .eq. 1 ) ) then 
      ! SY: REMEMBER: ending of section to uncomment  if states (smc etc) are 
             !available as observations
     
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      !!!!!!!!! CALL DATA ASSIMILATION !!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( ( data_assim ) .and. (initCycle .eq. EQs+1) ) then 

        IfAssimAtThisTimestep = .false.
        NumValidObsAtThisTimeStep = 0
        do zz = 1,Dz
          if (obs(t,zz) .gt. 0) then
            if (.not. IfAssimAtThisTimestep) then
              IfAssimAtThisTimestep = .true.
            endif
            NumValidObsAtThisTimeStep = NumValidObsAtThisTimeStep + 1
          endif !  if (obs(t,zz) .gt. 0) then
        enddo !  do zz = 1,Dz
 
        if (IfAssimAtThisTimestep) then

          allocate(Y(NumValidObsAtThisTimeStep,Ne))
          allocate(Z(NumValidObsAtThisTimeStep))
          allocate(R(NumValidObsAtThisTimeStep))
      
          WhichValidObsAtThisTimeStep = 0
          do zz = 1,Dz
            if (obs(t,zz) .gt. 0) then
              WhichValidObsAtThisTimeStep = WhichValidObsAtThisTimeStep + 1    
              do e = 1,Ne
                if (zz .eq. 1) then
                  Y(WhichValidObsAtThisTimeStep,e) = output(t,e)%nee
                elseif (zz .eq. 2) then
                  Y(WhichValidObsAtThisTimeStep,e) = output(t,e)%qe
                elseif (zz .eq. 3) then
                  Y(WhichValidObsAtThisTimeStep,e) = output(t,e)%qh
                endif ! if (zz .eq. 1) then
              enddo
              Z(WhichValidObsAtThisTimeStep) = obs(t,zz)
              R(WhichValidObsAtThisTimeStep) = zcov(zz) ! read in as 
                     ! observation covariance from obs_cov.txt 
            endif !  if (obs(t,zz) .gt. 0) then
          enddo !  do zz = 1,Dz

          ! soil moisture state updating
          ! X( lags, Soil Moisture Layers, No. Ensembles)
          allocate(X(lags,14,Ne))
          do l = 1,lags
            do e = 1,Ne
              X(l,1:4,e) = state(t-l+1,e)%smc(1:4)
              X(l,5:8,e) = state(t-l+1,e)%stc(1:4)
              X(l,9,e) = state(t-l+1,e)%lfmass
              X(l,10,e) = state(t-l+1,e)%stmass
              X(l,11,e) = state(t-l+1,e)%rtmass
              X(l,12,e) = state(t-l+1,e)%wood
              X(l,13,e) = state(t-l+1,e)%lai
              X(l,14,e) = state(t-l+1,e)%sai
            enddo
          enddo

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!  Ensemle Kalman Smoother !!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!     subroutine enks(X,Y,Zbar,Zsig,Nt,Ne,Dx,Dz)   !!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          call enks(X, Y, Z, R, lags, Ne, 14, &
                    NumValidObsAtThisTimeStep, Threshold_on_CC) 
          deallocate(Y)
          deallocate(Z)
          deallocate(R)
   
          ! This section updates the model states based on the X result from 
          !   enks(), including performing cuts 
          do l = 1,lags
            do e = 1,Ne

              !Updating each of the 4 moisture, water and temperature states 
              !  of soil based on the lag. -- SY
              do d = 1,4 ! Alternately, do d = 1,setup%nsoil

                state(t-l+1,e)%sh2o(d) = &
                    state(t-l+1,e)%sh2o(d)+X(l,d,e)-state(t-l+1,e)%smc(d)
                if (state(t-l+1,e)%sh2o(d).gt.smcmax) &
                                      state(t-l+1,e)%sh2o(d) = smcmax
                if (state(t-l+1,e)%sh2o(d).lt.smcdry)   &
                                    state(t-l+1,e)%sh2o(d) = smcdry

                state(t-l+1,e)%smc(d) = X(l,d,e)
                if (state(t-l+1,e)%smc(d).gt.smcmax)  &
                                    state(t-l+1,e)%smc(d) = smcmax
                if (state(t-l+1,e)%smc(d).lt.smcdry)    &
                                    state(t-l+1,e)%smc(d) = smcdry

                state(t-l+1,e)%stc(d) = X(l,d+4,e)

              enddo ! soil layers dimension

              state(t-l+1,e)%lfmass = X(l,9,e)
              if (state(t-l+1,e)%lfmass .le. 50/SLA)     &
                               state(t-l+1,e)%lfmass = 50/SLA + 0.01
              if (state(t-l+1,e)%lfmass .gt. 5000/SLA)   &
                               state(t-l+1,e)%lfmass = 5000/SLA

              state(t-l+1,e)%stmass = X(l,10,e)
              if (state(t-l+1,e)%stmass .le. 0.05/0.003) &
                               state(t-l+1,e)%stmass = 0.05/0.003 + 0.01

              state(t-l+1,e)%rtmass = X(l,11,e)
              if (state(t-l+1,e)%rtmass.le.5)          &
                               state(t-l+1,e)%rtmass = 5.01

              state(t-l+1,e)%wood = X(l,12,e)

              state(t-l+1,e)%lai = X(l,13,e)
              state(t-l+1,e)%lai = MAX(state(t-l+1,e)%lfmass*SLA/1000, &
                                       0.05)

              state(t-l+1,e)%sai = X(l,14,e)
              state(t-l+1,e)%sai = MAX(state(t-l+1,e)%stmass*0.003,0.05)

            enddo ! ensemble
          enddo ! lag
           
          deallocate(X)

        endif ! if (IfAssimAtThisTimestep) then

      endif ! if ( ( data_assim ) .and (initCycle .eq. EQs+1) ) then 
   
      ! Eqlibration process (i.e., update time 1 with the value at time end.)
      if (t.eq.Nt) then

!        if (initCycle .le. EQs+1) then
!
!          do e = 1,Ne
!            call printstatevalues(t, e, state, Nt, Ne) 
!          enddo ! do e = 1,Ne
!
!        endif !  if (initCycle .le. EQs+1) then

        if (initCycle .eq. EQs) then
 
          ! write the plant and soil states:
          open(fid,file='plant_init.txt')
            write(fid,*) state(Nt,1)%rtmass
            write(fid,*) state(Nt,1)%wood
            write(fid,*) state(Nt,1)%lfmass
            write(fid,*) state(Nt,1)%stmass
          close(fid)

          open(fid,file='soil_init.txt')
            write(fid,*) state(Nt,1)%smc(1)
            write(fid,*) state(Nt,1)%smc(2)
            write(fid,*) state(Nt,1)%smc(3)
            write(fid,*) state(Nt,1)%smc(4)
          close(fid)

        endif !  if (initCycle .eq. EQs) then
        
        if ((calc_obj_fun).and.(initCycle.le.EQs+1)) then 
          call rmse_obj_fun(date,time,state(:,1),output(:,1),Nt)
          if (initCycle.le.EQs) then
            print*, "Initiliztion cycle", initCycle
          else
            print*, "Run cycle"
          endif
        endif

      endif ! if (t.eq.Nt) then
   
    enddo ! time loop: do t = 1,Nt

  enddo ! equilibriation+run loop: do initCycle = 1,EQs+1

! ------- Write Output --------------------------------------------
  if (write_model_output) then
 
    ! Write the main output file reguardless of the other options.
    fname = ofile
    open(fid,file=trim(fname),status='replace')
    do t = 1,Nt
      write(fid,'(i7,i7,i7,i7,i7,                                    & 
                  f17.6,f17.6,f17.6                                    & 
                  f17.6,f17.6                                          & 
                  f17.6,f17.6                                          & 
                  f17.6,f17.6,f17.6,f17.6                              & 
                  f17.6,f17.6,f17.6,f17.6                              & 
                  f17.6,f17.6                                          & 
                  f17.6,f17.6                                          & 
                  f17.6,f17.6,f17.6)')                                 & 
      date(t,:),time(t,:), &
      forcing(t,1)%sfctmp, forcing(t,1)%sfcspd, forcing(t,1)%sfcprs,   &
      forcing(t,1)%q2, forcing(t,1)%lwrad,                             &
      forcing(t,1)%swrad, forcing(t,1)%prcprate,                       & 
      state(t,1)%smc(1:4),                                             &
      state(t,1)%sh2o(1:4),                                            &
      state(t,1)%rtmass, state(t,1)%wood,                              &
      state(t,1)%lfmass, state(t,1)%stmass,                            &
      output(t,1)%qe, output(t,1)%qh, output(t,1)%nee
    enddo ! do t = 1,Nt
    close(fid)
  
    ! Now go through the DA/Perturbation options.
    if (if_state_perturb) then

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
                     f17.6,f17.6,f17.6,f17.6,f17.6, &
                     f17.6,f17.6, &
                     f17.6,f17.6,f17.6,f17.6, &
                     f17.6,f17.6,f17.6,f17.6, &
                     f17.6,f17.6,f17.6)')               & 
          mean_sfctmp,mean_sfcspd,mean_sfcprs,mean_q2,mean_lwrad, &
          mean_swrad,mean_prcprate,                               &
          mean_smc1,mean_smc2,mean_smc3,mean_smc4,                &
          mean_rtmass,mean_wood,mean_lfmass,mean_stmass,          &
          mean_qe,mean_qh,mean_nee

      enddo ! do t = 1,Nt 

      close(fid)

    endif ! if (if_state_perturb) then
  
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
        write(fid,'(                                 &
                     f17.6,f17.6,f17.6,f17.6,f17.6, &
                     f17.6,f17.6, &
                     f17.6,f17.6,f17.6,f17.6, &
                     f17.6,f17.6,f17.6,f17.6, &
                     f17.6,f17.6,f17.6)')               & 
          mean_sfctmp,mean_sfcspd,mean_sfcprs,mean_q2,mean_lwrad, &
          mean_swrad,mean_prcprate,                               &
          mean_smc1,mean_smc2,mean_smc3,mean_smc4,                &
          mean_rtmass,mean_wood,mean_lfmass,mean_stmass,          &
          mean_qe,mean_qh,mean_nee
 
      enddo ! do t = 1,Nt 
 
      close(fid)

    endif ! If data_assim
   
    !ensemble output files
    if (write_ensemble_output) then

      do e = 1,min(9,Ne)
       
        if (if_state_perturb) then

          if (data_assim) then
            fname = 'enks_'
          else
            fname = 'open_'
          endif
 
          if (e.lt.10) then ! SY: NOTE: only this option is active when we 
               !have done the enveloping "do e = 1,min(9,Ne)" as is the case now
            write(es1,'(i1)') e
            fname = trim(fname)//es1
          elseif (e.lt.100) then
            write(es2,'(i2)') e
            fname = trim(fname)//es2
          elseif (e.lt.1000) then
            write(es3,'(i3)') e
            fname = trim(fname)//es3
          endif ! if (e.lt.10) then
          fname = trim(fname)//'.out'

          ! SY: Begin section that I moved to inside of the containing "if (if_state_perturb) then"
          open(fid,file=trim(fname),status='replace')
        
          do t = 1,Nt
            write(fid,'(i7,i7,i7,i7,i7,         & 
                       f17.6,f17.6,f17.6, &
                       f17.6,f17.6,f17.6, &
                       f17.6,f17.6,f17.6,f17.6,f17.6, &
                       f17.6,f17.6, &
                       f17.6,f17.6, &
                       f17.6,f17.6,f17.6)')               & 
              date(t,:),time(t,:), &
              forcing(t,e)%sfctmp, forcing(t,e)%sfcspd, forcing(t,e)%sfcprs,   &
              forcing(t,e)%q2, forcing(t,e)%lwrad, forcing(t,e)%swrad,         &
              forcing(t,e)%prcprate, state(t,e)%smc(1:4),                      &
              state(t,e)%rtmass, state(t,e)%wood,                              &
              state(t,e)%lfmass, state(t,e)%stmass,                            &
              output(t,e)%qe, output(t,e)%qh, output(t,e)%nee
          enddo ! do t = 1,Nt
  
          close(fid)
          ! SY: End section that I moved to inside of the containing "if (if_state_perturb) then"

        endif ! if (if_state_perturb) then

        if (data_assim) then

          fname = 'back_'
          if (e.lt.10) then ! SY: NOTE: only this option is active when we have 
                  ! done the enveloping "do e = 1,min(9,Ne)" as is the case now
            write(es1,'(i1)') e
            fname = trim(fname)//es1
          elseif (e.lt.100) then
            write(es2,'(i2)') e
            fname = trim(fname)//es2
          elseif (e.lt.1000) then
            write(es3,'(i3)') e
            fname = trim(fname)//es3
          endif ! if (e.lt.10) then
          fname = trim(fname)//'.out'

          open(fid,file=trim(fname),status='replace')
       
          do t = 1,Nt
            write(fid,'(i7,i7,i7,i7,i7,         & 
                       f17.6,f17.6,f17.6, &
                       f17.6,f17.6,f17.6, &
                       f17.6,f17.6,f17.6,f17.6,f17.6, &
                       f17.6,f17.6, &
                       f17.6,f17.6, &
                       f17.6,f17.6,f17.6)')               & 
              date(t,:),time(t,:), &
              forcing(t,e)%sfctmp, forcing(t,e)%sfcspd, forcing(t,e)%sfcprs,   &
              forcing(t,e)%q2, forcing(t,e)%lwrad, forcing(t,e)%swrad,         &
              forcing(t,e)%prcprate, state(t,e)%smc(1:4),                      &
              state(t,e)%rtmass, state(t,e)%wood,                              &
              state(t,e)%lfmass, state(t,e)%stmass,                            &
              output(t,e)%qe, output(t,e)%qh, output(t,e)%nee
          enddo

          close(fid)

        endif ! if (data_assim) then

      enddo ! do e = 1,min(9,Ne)

    endif ! if (write_ensemble_output) then

  endif ! if (write_model_output) then
 
! SY: Commenting out this rmse_obj_fun call. I don't think that we need this. 
!     We're already doing this for every iteration of equilibriation loop 
!     further above. 
!  if (calc_obj_fun) then
!    call rmse_obj_fun(date,time,state(:,1),output(:,1),Nt)
! !  call transfer_entropy_obj_fun(date,time,forcing(:,1),state(:,1),output(:,1),Nt)
!  endif

  ! SY: begin deallocations

  deallocate(setup%sldpth)

  do ee = 1, Ne 
    do tt = 1, Nt

      deallocate( state( tt, ee )%stc )
      deallocate( state( tt, ee )%zsnso )
      deallocate( state( tt, ee )%tsno )
      deallocate( state( tt, ee )%snice )
      deallocate( state( tt, ee )%snliq )
      deallocate( state( tt, ee )%sh2o )
      deallocate( state( tt, ee )%smc )

      if (tt .gt. 1) then
        deallocate( background( tt, ee )%stc )
        deallocate( background( tt, ee )%zsnso )
        deallocate( background( tt, ee )%tsno )
        deallocate( background( tt, ee )%snice )
        deallocate( background( tt, ee )%snliq )
        deallocate( background( tt, ee )%sh2o )
        deallocate( background( tt, ee )%smc )
      endif

    enddo
  enddo

  deallocate(state)
  deallocate(background)
  deallocate(output)

  deallocate(forcing)

  deallocate(date)
  deallocate(time)

  if ( (data_assim) .or. (isOneStep) .or. (doEQ) .or. (calc_obj_fun) ) then
    deallocate(obs)
    deallocate(NEE)
    !deallocate(GPP)
    deallocate(Qle)
    deallocate(Qh)
!    deallocate(firstObs) ! SY: REMEMBER to uncomment this if states (smc etc)
                            !  are available as observations
  endif
  if( (data_assim) .or. (if_state_perturb) ) then
    deallocate(zcov)
  endif

  ! SY: end deallocations

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

subroutine sim_init(setup, state, data_assim, sig_sm, sig_veg, Nlag, Threshold_on_CC, &
     doEQ, EQs, Dz, write_model_output, write_ensemble_output, calc_obj_fun, fixensemble)

  use type_decs
 
  integer, parameter :: fid = 14
  type(state_data)   :: state
  type(setup_data)   :: setup
  logical, intent(in) :: data_assim
  !jmframe added sig_sm & Nlag here, instead of hard coded. July 2019
  real, intent(out) :: sig_sm, sig_veg
  integer, intent(out) :: Nlag, EQs, Dz
  logical, intent(out) :: doEQ, write_model_output, write_ensemble_output, &
          calc_obj_fun, fixensemble
  real, intent(out) :: Threshold_on_CC
  logical :: fexists

 ! simulation setup
  open(fid, file='init.txt', action='read')

  read(fid,*) setup%nsoil     
  read(fid,*) setup%nsnow     
 
! allocate state and setup component dimensions

  allocate( state%stc( -setup%nsnow+1: setup%nsoil ) )
  allocate( state%zsnso( -setup%nsnow+1: setup%nsoil ) )
  allocate( state%tsno( setup%nsnow ) )
  allocate( state%snice( setup%nsnow ) )
  allocate( state%snliq( setup%nsnow ) )
  allocate( state%sh2o( setup%nsoil ) )
  allocate( state%smc( setup%nsoil ) )

  allocate( setup%sldpth( setup%nsoil ) )

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

! initial states
  read(fid,*) state%stc( 1: setup%nsoil ) 
  read(fid,*) state%snowh   
  read(fid,*) state%sneqv   
  read(fid,*) state%canliq  
  read(fid,*) state%rtmass  
  read(fid,*) state%albold  
  read(fid,*) state%lai     
  read(fid,*) state%tv      
  read(fid,*) state%tg      

  close(fid)
 
  open(fid, file='startdate.txt')
  read(fid,*) setup%startdate
  close(fid)

  inquire(file='init_flag.txt',exist=fexists)
  if (fexists) then
    open(fid,file='init_flag.txt')
     read(fid,*) EQs
    close(fid)
    if (EQs.gt.0) then
      doEQ = .true.
    else
      doEQ = .false.
      EQs = 0 ! SY
    endif
  else
    doEQ = .false.
    EQs = 0 ! SY
  endif

  Dz = 3 ! commenting out next 3 lines for now
!  open(fid,file='Num_ObsVariables.txt')
!  read(fid,*) Dz ! should be 3 for NEE, Qle, Qh
!  close(fid)

  ! values needed for data assimilation
  if (data_assim) then
    open(fid,file='sig_sm.txt')
    read(fid,*) sig_sm
    close(fid)
    open(fid,file='sig_veg.txt')
    read(fid,*) sig_veg
    close(fid)
    open(fid,file='Nlag.txt')
    read(fid,*) Nlag
    close(fid)
    inquire(file='Threshold_on_CC.txt',exist=fexists)
    if (fexists) then
      open(fid,file='Threshold_on_CC.txt')
       read(fid,*) Threshold_on_CC
      close(fid)
    else
      Threshold_on_CC = -9999.0
    endif
  endif !data_assim
 
 ! specifications for output
  open(fid, file='output_specifications.txt', action='read')
  read(fid,*) write_model_output
  read(fid,*) write_ensemble_output
  read(fid,*) calc_obj_fun
  read(fid,*) fixensemble
  close(fid)
 
  open(fid, file='lat_lon.txt')
  read(fid,*) setup%latitude
  read(fid,*) setup%longitude
  close(fid)
 
! initial plant-related carbon states
  open(fid, file='plant_init.txt')
  read(fid,*) state%rtmass
  read(fid,*) state%wood
  read(fid,*) state%lfmass
  read(fid,*) state%stmass
  close(fid)
 
! initial soil-related moisture states
  open(fid, file='soil_init.txt')
  read(fid,*) state%smc(1)
  read(fid,*) state%smc(2)
  read(fid,*) state%smc(3)
  read(fid,*) state%smc(4)
  state%sh2o = state%smc
  close(fid)
 
  inquire(file='shdfac.txt', exist=fexists)
  if (fexists) then
    open(fid,file='shdfac.txt')
      read(fid,*) setup%shdfac_monthly
    close(fid)
  else
    setup%shdfac_monthly = (/0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98,0.98/)
  endif

end subroutine sim_init

subroutine perturb_state(setup, state, sig_sm, sig_veg)
 
  use type_decs
  use noahmp_veg_parameters
  use noahmp_globals

  implicit none
 
  type(state_data), intent(inout) :: state
  type(setup_data), intent(in)    :: setup
  real, intent(in)    :: sig_sm, sig_veg
! real, parameter :: sig_sm = 0.005
! real, parameter :: sig_veg = 0.01

  real :: eta
  real :: random_normal
  integer :: d
 
  ! soil moisture
  do d = 1,setup%nsoil
    eta = random_normal()
    state%smc(d) = state%smc(d) + eta * sig_sm
    state%sh2o(d) = state%sh2o(d) + eta * sig_sm
    eta = random_normal()
    state%stc(d) = state%stc(d) + eta * sig_sm ! SY: Assuming it's the 
            ! same sigma for soil moisture and temperature. So sig_sm 
            ! should be more appropriately named sig_soil. If it's 
            ! different for soil moisture and soil temperature, then 
            ! we perhaps need 2 separate sigmas of sig_sm and sig_st respectively
  enddo

  ! plant stores

  eta = random_normal()
  state%lfmass = state%lfmass + &
               eta * sig_veg * state%lfmass

  eta = random_normal()
  state%stmass = state%stmass + &
               eta * sig_veg * state%stmass
 
  eta = random_normal()
  state%rtmass = state%rtmass + &
               eta * sig_veg * state%rtmass

! eta = random_normal()
! state%lfmass = state%lfmass+eta*(sig_veg*(state%lfmass-50/SLA))
! if (state%lfmass.le.50/SLA) state%lfmass = 50/SLA+0.01
! if (state%lfmass.ge.5000/SLA) state%lfmass = 5000/SLA
!
! eta = random_normal()
! state%stmass = state%stmass+eta*(sig_veg*(state%stmass-0.05/0.003))
! if (state%stmass.le.0.05/0.003) state%stmass = 0.05/0.003+0.01 
!
! eta = random_normal()
! state%rtmass = state%rtmass+eta*(sig_veg*state%rtmass)
! if (state%rtmass.lt.5) state%rtmass = 5.01
!
! state%lai = max(state%lfmass*SLA/1000,0.05)
! state%sai = max(state%stmass*3*0.001,0.05)

  eta = random_normal()
  state%wood = state%wood + &
               eta * sig_veg * state%wood

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

subroutine printstatevalues(t, e, state, Nt, Ne)
  use type_decs
  implicit none
  integer, intent(in) :: t, e
  integer, intent(in) :: Nt, Ne
  type(state_data),   intent(in), dimension(Nt, Ne) :: state
  print*, 'States at time=', t, ', ensemble member=', e 
  print*, 'stc'
  print*, state(t, e)%stc
  print*, 'zsnso'
  print*, state(t, e)%zsnso
  print*, 'tsno'
  print*, state(t, e)%tsno
  print*, 'snice'
  print*, state(t, e)%snice
  print*, 'snliq'
  print*, state(t, e)%snliq
  print*, 'sh2o'
  print*, state(t, e)%sh2o
  print*, 'smc'
  print*, state(t, e)%smc
  print*, 'snowh'   
  print*, state(t, e)%snowh   
  print*, 'sneqv'   
  print*, state(t, e)%sneqv   
  print*, 'canliq'  
  print*, state(t, e)%canliq  
  print*, 'rtmass'  
  print*, state(t, e)%rtmass  
  print*, 'albold'  
  print*, state(t, e)%albold  
  print*, 'lai'     
  print*, state(t, e)%lai     
  print*, 'tv'      
  print*, state(t, e)%tv      
  print*, 'tg'      
  print*, state(t, e)%tg      
  print*, 'rtmass'
  print*, state(t, e)%rtmass
  print*, 'wood'
  print*, state(t, e)%wood
  print*, 'lfmass'
  print*, state(t, e)%lfmass
  print*, 'stmass'
  print*, state(t, e)%stmass
end subroutine printstatevalues


