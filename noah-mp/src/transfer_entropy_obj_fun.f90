subroutine transfer_entropy_obj_fun(mdate,mtime,forcing,state,output,Nt)!ofval,errcode)

 use type_decs
 implicit none

 ! ------------------------------------------------------------------
 ! variable declarations
 ! dimensions
 integer, intent(in) :: Nt
 integer, parameter :: Dx = 4, Du = 6
 integer, parameter :: Nof = Dx + Dx**2 + Dx*Du 
 integer :: ti, xi, ui, yi, ii
 integer, dimension(Dx) :: Gx

 ! dates
 integer, dimension(Nt,2) :: mdate, odate
 real, dimension(Nt) :: mtime, otime
 real :: dummy

 ! model vars
 type(forcing_data), intent(in), dimension(Nt) :: forcing
 type(state_data),   intent(in), dimension(Nt) :: state
 type(output_data),  intent(in), dimension(Nt) :: output

 ! obs/mod vars
 real, dimension(Nt,Dx) :: modl, obsv
 real, dimension(Dx) :: mean_modl, mean_obsv, std_modl, std_obsv
 real, dimension(Du) :: mean_forc, std_forc
 real, dimension(Nt,Du) :: forc

 ! objective function variables
 real, dimension(Nof)   :: ofval
 real, dimension(Dx)    :: rmse
 real, dimension(Du,Dx) :: TEuxObs, TEuxMdl, MIuxObs, MIuxMdl
 real, dimension(Dx,Dx) :: TExxObs, TExxMdl, MIxxObs, MIxxMdl
 real, dimension(Dx)    :: HHuxObs, HHuxMdl, HHxxObs, HHxxMdl 
 integer :: ofcounter

 ! transfer entropy calculations
 real, dimension(Nt-1)  :: Xt, Ys, Yt
 real, allocatable, dimension(:)  :: gXt, gYs, gYt
 integer, parameter :: Nbins = 35
 real, dimension(Du) :: Bu
 real, dimension(Dx) :: Bo, Bm

 ! I/O parms
 integer, parameter :: fid = 136

 ! init objective function counter
 ofcounter = 0
 ofval = 0.

 ! ------------------------------------------------------------------
 ! read observation file
 open(fid,file='obs.txt')
  do ti = 1,Nt
   read(fid,*) odate(ti,:),otime(ti),obsv(ti,1:3),dummy,obsv(ti,4)
  enddo ! times
 close(fid)

 ! ------------------------------------------------------------------
 ! check dates
 do ti = 1,Nt
  if ((odate(ti,1).ne.mdate(ti,1)).or.  &
      (odate(ti,2).ne.mdate(ti,2)).or.  &
      (otime(ti)  .ne.mtime(ti))) then
    print*, 'mismatch_obs_dates'
    stop 111
   return
  endif
 enddo ! times 

 ! ------------------------------------------------------------------
 ! extract pertinent modeled dimensions
 do ti = 1,Nt
  modl(ti,1) = output(ti)%qe 
  modl(ti,2) = output(ti)%qh
  modl(ti,3) = state(ti)%smc(1)
!  modl(ti,4) = state(ti)%smc(2)
  modl(ti,4) = output(ti)%nee
 enddo ! times 

 do ti = 1,Nt
  forc(ti,1) = forcing(ti)%sfcspd
  forc(ti,2) = forcing(ti)%lwrad + forcing(ti)%swrad
  forc(ti,3) = forcing(ti)%q2
  forc(ti,4) = forcing(ti)%prcprate
  forc(ti,5) = forcing(ti)%sfctmp
  forc(ti,6) = forcing(ti)%sfcprs
 enddo

 ! ------------------------------------------------------------------
 ! standardize all output dimensions

 ! calculate expected values
 mean_modl = 0.
 mean_obsv = 0.
 mean_forc = 0.
 Gx = 0

 do ti = 1,Nt

  do ui = 1,Du
   if (forc(ti,ui).gt.-9990) then
    mean_forc(ui) = mean_forc(ui) + forc(ti,ui)
   endif
  enddo

  do xi = 1,Dx
   if ((obsv(ti,xi).gt.-9990).and.(modl(ti,xi).gt.-9990)) then
    Gx(xi) = Gx(xi) + 1
    mean_modl(xi)  = mean_modl(xi)  + modl(ti,xi) 
    mean_obsv(xi)  = mean_obsv(xi)  + obsv(ti,xi) 
   endif
  enddo

 enddo

 do ui = 1,Du
  mean_forc(ui)  = mean_forc(ui)/Nt 
 enddo

 do xi = 1,Dx
  mean_modl(xi)  = mean_modl(xi)/Gx(xi) 
  mean_obsv(xi)  = mean_obsv(xi)/Gx(xi) 
 enddo

 ! calculate standard deviations
 std_modl = 0.
 std_obsv = 0.
 std_forc = 0.

 do ti = 1,Nt

  do ui = 1,Du
   if (forc(ti,ui).gt.-9990) then
    std_forc(ui)  = std_forc(ui)  + (forc(ti,ui) - mean_forc(ui))**2 
   endif
  enddo

  do xi = 1,Dx
   if ((obsv(ti,xi).gt.-9990).and.(modl(ti,xi).gt.-9990)) then
    std_modl(xi)  = std_modl(xi)  + (modl(ti,xi) - mean_modl(xi))**2 
    std_obsv(xi)  = std_obsv(xi)  + (obsv(ti,xi) - mean_obsv(xi))**2 
   endif
  enddo

 enddo ! times 

 do ui = 1,Du
  std_forc(ui)  = sqrt(std_forc(ui)/(Nt-1)) 
 enddo

 do xi = 1,Dx
  std_modl(xi)  = sqrt(std_modl(xi)/(Gx(xi)-1)) 
  std_obsv(xi)  = sqrt(std_obsv(xi)/(Gx(xi)-1)) 
 enddo

 ! ------------------------------------------------------------------
 ! calcualte squared errors for each output dimension
 rmse = 0.
 do ti = 1,Nt
  do xi = 1,Dx
   if ((obsv(ti,xi).gt.-9990).and.(modl(ti,xi).gt.-9990).and.Gx(xi).gt.1000) then
    rmse(xi) = rmse(xi) + ((obsv(ti,xi) - modl(ti,xi))**2)/(Gx(xi)-1) 
   endif 
  enddo 
 enddo

 ! standardize to R^2 (coefficeint of determination)
 do xi = 1,Dx
  if (Gx(xi).gt.1000) then
   rmse(xi) = rmse(xi)/std_obsv(xi)**2
  endif
  if (abs(std_obsv(xi)).lt.1e-6) then
   rmse(xi) = 2.
  endif
  ofcounter = ofcounter +1
  ofval(ofcounter) = rmse(xi)
 enddo

 ! ------------------------------------------------------------------
 ! calcualte bin widths
 do ui = 1,Du
  Bu(ui) = (maxval(forc(:,ui),forc(:,ui).gt.-9990) - &
            minval(forc(:,ui),forc(:,ui).gt.-9990) + 2e-6) / Nbins
 enddo

 do xi = 1,Dx
  Bo(xi) = (maxval(obsv(:,xi),obsv(:,xi).gt.-9990) - &
            minval(obsv(:,xi),obsv(:,xi).gt.-9990) + 2e-6) / Nbins
 enddo

 do xi = 1,Dx
  Bm(xi) = (maxval(modl(:,xi),modl(:,xi).gt.-9990) - &
            minval(modl(:,xi),modl(:,xi).gt.-9990) + 2e-6) / Nbins
 enddo

 ! ------------------------------------------------------------------
 ! transfer entropy pathways: forcings -> states/fluxes
 TEuxMdl = -1.
 MIuxMdl = -1.
 HHuxMdl = -1.
 TEuxObs = -1.
 MIuxObs = -1.
 HHuxObs = -1.
 do ui = 1,Du                  ! loop through senders
  do yi = 1,Dx                 ! loop through targets
   if (Gx(yi).gt.1000) then      ! if enough good data
    allocate(gXt(Gx(yi)-1))    ! allocate sender vector
    allocate(gYt(Gx(yi)-1))    ! allocate conditional vector
    allocate(gYs(Gx(yi)-1))    ! allocate target vector

    ! ------------
    Xt = forc(2:Nt,ui)         ! pull sender data
    Ys = obsv(2:Nt,yi)         ! pull target data
    Yt = obsv(1:Nt-1,yi)       ! pull conditional data

    ii = 0
    gXt = mean_forc(ui)
    gYt = mean_obsv(yi)
    gYs = mean_obsv(yi)
    do ti = 1,Nt-1             ! pull good data
     if ((Xt(ti).gt.-9990).and.(Ys(ti).gt.-9990).and.(Yt(ti).gt.-9990)) then
      ii = ii + 1
      gXt(ii) = Xt(ti)
      gYs(ii) = Ys(ti)
      gYt(ii) = Yt(ti)
      if (ii.eq.Gx(yi)-1) then; exit; endif
     endif
    enddo

    ! calcualte entropy metrics - observation
!print*,'---------------------------------------'
!print*,'te 1',ui,yi,ii,Gx(yi)-1
!print*,'te 1',Bu(ui),maxval(forc(:,ui),forc(:,ui).gt.-9990),minval(forc(:,ui),forc(:,ui).gt.-9990)
!print*,'te 1',Bo(yi),maxval(obsv(:,yi),obsv(:,yi).gt.-9990),minval(obsv(:,yi),obsv(:,yi).gt.-9990)
!print*,'te 1',maxval(Ys),minval(Ys),maxval(gYs),minval(gYs)
!print*,'te 1',maxval(Yt),minval(Yt),maxval(gYt),minval(gYt)
    call transfer_entropy(gYs,gXt,gYt,Gx(yi)-1,Bu(ui),Bo(yi),TEuxObs(ui,yi),HHuxObs(yi),MIuxObs(ui,yi))

    ! ------------
    Xt = forc(2:Nt,ui)         ! pull sender data
    Ys = modl(2:Nt,yi)         ! pull target data
    Yt = modl(1:Nt-1,yi)       ! pull conditional data

    ii = 0
    gXt = mean_forc(ui)
    gYt = mean_modl(yi)
    gYs = mean_modl(yi)
    do ti = 1,Nt-1             ! pull good data
     if ((Xt(ti).gt.-9990).and.(Ys(ti).gt.-9990).and.(Yt(ti).gt.-9990)) then
      ii = ii + 1
      gXt(ii) = Xt(ti)
      gYs(ii) = Ys(ti)
      gYt(ii) = Yt(ti)
      if (ii.eq.Gx(yi)-1) then; exit; endif
     endif
    enddo

    ! calcualte entropy metrics - model
!print*,'---------------------------------------'
!print*,'te 2',ui,yi,ii,Gx(yi)-1
!print*,'te 2',Bu(ui),maxval(forc(:,ui),forc(:,ui).gt.-9990),minval(forc(:,ui),forc(:,ui).gt.-9990)
!print*,'te 2',Bm(yi),maxval(modl(:,yi),modl(:,yi).gt.-9990),minval(modl(:,yi),modl(:,yi).gt.-9990)
    call transfer_entropy(gYs,gXt,gYt,Gx(yi)-1,Bu(ui),Bm(yi),TEuxMdl(ui,yi),HHuxMdl(yi),MIuxMdl(ui,yi))

    ! ------------
    ! store in of output variable
    ofcounter = ofcounter +1
    ofval(ofcounter) = abs((TEuxMdl(ui,yi)-TEuxObs(ui,yi))/HHuxObs(yi))

    deallocate(gYt)               ! deallocate conditional storage
    deallocate(gYs)               ! deallocate target storage
    deallocate(gXt)               ! deallocate sender storage
   endif                          ! enough good data
  enddo ! u-loop
 enddo ! x-loop

 ! transfer entropy pathways: states/fluxes -> states/fluxes
 TExxMdl = -1.
 MIxxMdl = -1.
 HHxxMdl = -1.
 TExxObs = -1.
 MIxxObs = -1.
 HHxxObs = -1.
 do xi = 1,Dx                             ! loop through senders
  do yi = 1,Dx                            ! loop through targets
   if ((Gx(xi).gt.1000).and.(Gx(yi).gt.1000)) then
    allocate(gXt(min(Gx(xi),Gx(yi))-1))   ! allocate sender vector
    allocate(gYt(min(Gx(xi),Gx(yi))-1))   ! allocate conditional vector
    allocate(gYs(min(Gx(xi),Gx(yi))-1))   ! allocate target vector

    ! ------------- 
    Xt = obsv(1:Nt-1,xi)                  ! pull sender data
    Ys = obsv(2:Nt,yi)                    ! pull target data
    Yt = obsv(1:Nt-1,yi)                  ! pull conditional data

    if (xi.eq.yi) then
     Yt = mean_obsv(yi)
    endif

    ii = 0
    gXt = mean_obsv(xi)
    gYt = mean_obsv(yi)
    gYs = mean_obsv(yi)
    do ti = 1,Nt-1                        ! pull good data
     if ((Xt(ti).gt.-9990).and.(Ys(ti).gt.-9990).and.(Yt(ti).gt.-9990)) then
      ii = ii + 1
      gXt(ii) = Xt(ti)
      gYs(ii) = Ys(ti)
      gYt(ii) = Yt(ti)
      if (ii.eq.min(Gx(xi),Gx(yi))-1) then; exit; endif
     endif
    enddo

    ! calculate transfer entropy metrics - observation data
    call transfer_entropy(gYs,gXt,gYt,min(Gx(xi),Gx(yi))-1,    &
          Bo(xi),Bo(yi),TExxObs(xi,yi),HHxxObs(yi),MIxxObs(xi,yi))

    ! ------------- 
    Xt = modl(1:Nt-1,xi)                  ! pull sender data
    Ys = modl(2:Nt,yi)                    ! pull target data
    Yt = modl(1:Nt-1,yi)                  ! pull conditional data

    if (xi.eq.yi) then
     Yt = mean_modl(yi)
    endif

    ii = 0
    gXt = mean_modl(xi)
    gYt = mean_modl(yi)
    gYs = mean_modl(yi)
    do ti = 1,Nt-1                        ! pull good data
     if ((Xt(ti).gt.-9990).and.(Ys(ti).gt.-9990).and.(Yt(ti).gt.-9990)) then
      ii = ii + 1
      gXt(ii) = Xt(ti)
      gYs(ii) = Ys(ti)
      gYt(ii) = Yt(ti)
      if (ii.eq.min(Gx(xi),Gx(yi))-1) then; exit; endif
     endif
    enddo

    ! calculate transfer entropy metrics - modeled data
    call transfer_entropy(gYs,gXt,gYt,min(Gx(xi),Gx(yi))-1,    &
          Bm(xi),Bm(yi),TExxMdl(xi,yi),HHxxMdl(yi),MIxxMdl(xi,yi))

    ! ------------
    ! store in of output variable
    ofcounter = ofcounter +1
    ofval(ofcounter) = abs((TExxMdl(xi,yi)-TExxObs(xi,yi))/HHxxObs(yi))

    deallocate(gYt)               ! deallocate conditional storage
    deallocate(gYs)               ! deallocate target storage
    deallocate(gXt)               ! deallocate sender storage
   endif                          ! enough good data
  enddo ! x-loop
 enddo ! x-loop

 ! ------------------------------------------------------------------
 ! write to output file
 open(fid,file='noahmp_objfun.out')
 do xi = 1,Nof
  write(fid,*) ofval(xi),','
 enddo
 close(fid)

! write(fid,'(A15,F20.10)') 'QeMSE = ',  rmse(1)
! write(fid,'(A15,F20.10)') 'QhMSE = ',  rmse(2)
! write(fid,'(A15,F20.10)') 'SM1MSE = ', rmse(3)
! write(fid,'(A15,F20.10)') 'SM2MSE = ', rmse(4)
! write(fid,'(A15,F20.10)') 'NEEMSE = ', rmse(5)
!
! write(fid,'(A15,F20.10)') 'TEuxDelta11 = ', TEuxDelta(1,1)
! write(fid,'(A15,F20.10)') 'TEuxDelta12 = ', TEuxDelta(1,2)
! write(fid,'(A15,F20.10)') 'TEuxDelta13 = ', TEuxDelta(1,3)
! write(fid,'(A15,F20.10)') 'TEuxDelta14 = ', TEuxDelta(1,4)
! write(fid,'(A15,F20.10)') 'TEuxDelta15 = ', TEuxDelta(1,5)
! write(fid,'(A15,F20.10)') 'TEuxDelta21 = ', TEuxDelta(2,1)
! write(fid,'(A15,F20.10)') 'TEuxDelta22 = ', TEuxDelta(2,2)
! write(fid,'(A15,F20.10)') 'TEuxDelta23 = ', TEuxDelta(2,3)
! write(fid,'(A15,F20.10)') 'TEuxDelta24 = ', TEuxDelta(2,4)
! write(fid,'(A15,F20.10)') 'TEuxDelta25 = ', TEuxDelta(2,5)
! write(fid,'(A15,F20.10)') 'TEuxDelta31 = ', TEuxDelta(3,1)
! write(fid,'(A15,F20.10)') 'TEuxDelta32 = ', TEuxDelta(3,2)
! write(fid,'(A15,F20.10)') 'TEuxDelta33 = ', TEuxDelta(3,3)
! write(fid,'(A15,F20.10)') 'TEuxDelta34 = ', TEuxDelta(3,4)
! write(fid,'(A15,F20.10)') 'TEuxDelta35 = ', TEuxDelta(3,5)
! write(fid,'(A15,F20.10)') 'TEuxDelta41 = ', TEuxDelta(4,1)
! write(fid,'(A15,F20.10)') 'TEuxDelta42 = ', TEuxDelta(4,2)
! write(fid,'(A15,F20.10)') 'TEuxDelta43 = ', TEuxDelta(4,3)
! write(fid,'(A15,F20.10)') 'TEuxDelta44 = ', TEuxDelta(4,4)
! write(fid,'(A15,F20.10)') 'TEuxDelta45 = ', TEuxDelta(4,5)
! write(fid,'(A15,F20.10)') 'TEuxDelta51 = ', TEuxDelta(5,1)
! write(fid,'(A15,F20.10)') 'TEuxDelta52 = ', TEuxDelta(5,2)
! write(fid,'(A15,F20.10)') 'TEuxDelta53 = ', TEuxDelta(5,3)
! write(fid,'(A15,F20.10)') 'TEuxDelta54 = ', TEuxDelta(5,4)
! write(fid,'(A15,F20.10)') 'TEuxDelta55 = ', TEuxDelta(5,5)
! write(fid,'(A15,F20.10)') 'TEuxDelta61 = ', TEuxDelta(6,1)
! write(fid,'(A15,F20.10)') 'TEuxDelta62 = ', TEuxDelta(6,2)
! write(fid,'(A15,F20.10)') 'TEuxDelta63 = ', TEuxDelta(6,3)
! write(fid,'(A15,F20.10)') 'TEuxDelta64 = ', TEuxDelta(6,4)
! write(fid,'(A15,F20.10)') 'TEuxDelta65 = ', TEuxDelta(6,5)
!
! write(fid,'(A15,F20.10)') 'TExxDelta11 = ', TExxDelta(1,1)
! write(fid,'(A15,F20.10)') 'TExxDelta12 = ', TExxDelta(1,2)
! write(fid,'(A15,F20.10)') 'TExxDelta13 = ', TExxDelta(1,3)
! write(fid,'(A15,F20.10)') 'TExxDelta14 = ', TExxDelta(1,4)
! write(fid,'(A15,F20.10)') 'TExxDelta15 = ', TExxDelta(1,5)
! write(fid,'(A15,F20.10)') 'TExxDelta21 = ', TExxDelta(2,1)
! write(fid,'(A15,F20.10)') 'TExxDelta22 = ', TExxDelta(2,2)
! write(fid,'(A15,F20.10)') 'TExxDelta23 = ', TExxDelta(2,3)
! write(fid,'(A15,F20.10)') 'TExxDelta24 = ', TExxDelta(2,4)
! write(fid,'(A15,F20.10)') 'TExxDelta25 = ', TExxDelta(2,5)
! write(fid,'(A15,F20.10)') 'TExxDelta31 = ', TExxDelta(3,1)
! write(fid,'(A15,F20.10)') 'TExxDelta32 = ', TExxDelta(3,2)
! write(fid,'(A15,F20.10)') 'TExxDelta33 = ', TExxDelta(3,3)
! write(fid,'(A15,F20.10)') 'TExxDelta34 = ', TExxDelta(3,4)
! write(fid,'(A15,F20.10)') 'TExxDelta35 = ', TExxDelta(3,5)
! write(fid,'(A15,F20.10)') 'TExxDelta41 = ', TExxDelta(4,1)
! write(fid,'(A15,F20.10)') 'TExxDelta42 = ', TExxDelta(4,2)
! write(fid,'(A15,F20.10)') 'TExxDelta43 = ', TExxDelta(4,3)
! write(fid,'(A15,F20.10)') 'TExxDelta44 = ', TExxDelta(4,4)
! write(fid,'(A15,F20.10)') 'TExxDelta45 = ', TExxDelta(4,5)
! write(fid,'(A15,F20.10)') 'TExxDelta51 = ', TExxDelta(5,1)
! write(fid,'(A15,F20.10)') 'TExxDelta52 = ', TExxDelta(5,2)
! write(fid,'(A15,F20.10)') 'TExxDelta53 = ', TExxDelta(5,3)
! write(fid,'(A15,F20.10)') 'TExxDelta54 = ', TExxDelta(5,4)
! write(fid,'(A15,F20.10)') 'TExxDelta55 = ', TExxDelta(5,5)
!
! close(fid)

 ! ------------------------------------------------------------------
 return
end subroutine
