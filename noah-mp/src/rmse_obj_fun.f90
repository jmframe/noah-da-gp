subroutine rmse_obj_fun(mdate,mtime,state,output,Nt)!ofval,errcode)

 use type_decs
 implicit none

 ! ------------------------------------------------------------------
 ! variable declarations
 ! dimensions
 integer, intent(in) :: Nt
 integer, parameter :: Dx = 3
 integer :: ti, xi
 integer, dimension(Dx) :: Gx

 ! dates
 integer, dimension(Nt,3) :: mdate, odate
 integer, dimension(Nt,2) :: mtime, otime
 real :: dummy

 ! obs/mod vars
 real, dimension(Nt,Dx) :: modl, obsv

 ! model vars
 type(state_data),   intent(in), dimension(Nt) :: state
 type(output_data),  intent(in), dimension(Nt) :: output

 ! objective function variables
 real, dimension(Dx) :: sse, nsmean, nsmeandiff

 ! I/O parms
 integer, parameter :: fid = 136

 ! ------------------------------------------------------------------
 ! read observation file
 ! jmframe: PLUMBER-2 flux data Ordered from text file make from NetCDF data.
 ! year,month,day,hour,minute,NEE,GPP,Qle,Qh
 open(fid,file='obs.txt')
  do ti = 1,Nt
   read(fid,*) odate(ti,:),otime(ti,:),obsv(ti,1),dummy,obsv(ti,2:3)
  enddo ! times
 close(fid)

 ! ------------------------------------------------------------------
 ! check dates
 do ti = 1,Nt
  if ((odate(ti,1).ne.mdate(ti,1)).or.  &
      (odate(ti,2).ne.mdate(ti,2)).or.  &
      (odate(ti,3).ne.mdate(ti,3)).or.  &
      (otime(ti,1).ne.mtime(ti,1)).or.  &
      (otime(ti,2).ne.mtime(ti,2))) then
    print*, 'mismatch_obs_dates'
    stop 111
   return
  endif
 enddo ! times 

 ! ------------------------------------------------------------------
 ! extract pertinent modeled dimensions
 do ti = 1,Nt
  modl(ti,1) = output(ti)%nee
  modl(ti,2) = output(ti)%qe 
  modl(ti,3) = output(ti)%qh
 enddo ! times 

 ! ------------------------------------------------------------------
 ! calcualte squared errors for each output dimension
 ! And calc the mean for NSE
 sse = 0.
 Gx = 0
 nsmean = 0. 
 do ti = 1,Nt
   do xi = 1,Dx
     if ((obsv(ti,xi).gt.-9990).and.(modl(ti,xi).gt.-9990)) then
       sse(xi) = sse(xi) + ((obsv(ti,xi) - modl(ti,xi))**2)
       Gx(xi) = Gx(xi) + 1
       ! Calculate the record mean for a Nash Sutcliffe efficiency
       ! will devide by number of good records in next loop
       nsmean(xi) = nsmean(xi) + obsv(ti,xi) 
     endif
   enddo 
 enddo
 ! Divide by good records to get mean
 do xi = 1,Dx
   nsmean(xi) = nsmean(xi)/Gx(xi)
 enddo
 ! calcualte anomaly (diff from mean) for each output dimension
 ! for the Nash-Sutcliffe efficiency
 nsmeandiff = 0.
 do ti = 1,Nt
   do xi = 1,Dx
     if ((obsv(ti,xi).gt.-9990).and.(modl(ti,xi).gt.-9990)) then
       nsmeandiff(xi) = nsmeandiff(xi) + ((obsv(ti,xi) - nsmean(xi))**2)
     endif
   enddo 
 enddo
 ! ------------------------------------------------------------------
 ! write to output file - Root Mean Sum of Squared Error
 open(fid,file='noahmp_objfun.out')
   write(fid,*) "Root Mean Sum of Squared Error"
   write(fid,*) "NEE:", sqrt(sse(1)/Gx(1))
   write(fid,*) "Qle:", sqrt(sse(2)/Gx(2))
   write(fid,*) "Qh:", sqrt(sse(2)/Gx(2))
   write(fid,*) "Nash-Sutcliffe efficiency"
   write(fid,*) "NEE", 1-(sse(1)/nsmeandiff(1))
   write(fid,*) "Qle", 1-(sse(2)/nsmeandiff(2))
   write(fid,*) "Qh", 1-(sse(3)/nsmeandiff(3))
 close(fid)
 print*, "Nash-Sutcliffe efficiency"
 print*, "NEE", 1-(sse(1)/nsmeandiff(1))
 print*, "Qle", 1-(sse(2)/nsmeandiff(2))
 print*, "Qh", 1-(sse(3)/nsmeandiff(3))

 ! ------------------------------------------------------------------
 return
end subroutine
