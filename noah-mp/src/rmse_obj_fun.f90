subroutine rmse_obj_fun(mdate,mtime,state,output,Nt)!ofval,errcode)

 use type_decs
 implicit none

 ! ------------------------------------------------------------------
 ! variable declarations
 ! dimensions
 integer, intent(in) :: Nt
 integer, parameter :: Dx = 4
 integer :: ti, xi
 integer :: Gx

 ! dates
 integer, dimension(Nt,3) :: mdate, odate
 real, dimension(Nt) :: mtime, otime
 real :: dummy

 ! obs/mod vars
 real, dimension(Nt) :: modl, obsv

 ! model vars
 type(state_data),   intent(in), dimension(Nt) :: state
 type(output_data),  intent(in), dimension(Nt) :: output

 ! objective function variables
 real    :: sse, nsmean, nsmeandiff

 ! I/O parms
 integer, parameter :: fid = 136

 ! ------------------------------------------------------------------
 ! read observation file
 open(fid,file='obs.txt')
  do ti = 1,Nt
   read(fid,*) odate(ti,:),otime(ti),obsv(ti),dummy
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
  modl(ti) = state(ti)%smc(1)
 enddo ! times 

 ! ------------------------------------------------------------------
 ! calcualte squared errors for each output dimension
 sse = 0.
 Gx = 0
 nsmean = 0. 
do ti = 1,Nt
  if ((obsv(ti).gt.-9990).and.(modl(ti).gt.-9990)) then
   sse = sse + ((obsv(ti) - modl(ti))**2)
   Gx = Gx + 1
   ! Calculate the record mean for a Nash Sutcliffe efficiency
   nsmean = nsmean + obsv(ti)
   ! calcualte anomaly (diff from mean) for each output dimension
   ! for the Nash-Sutcliffe efficiency
  endif 
 enddo

 ! ------------------------------------------------------------------
 ! calcualte anomaly (diff from mean) for each output dimension
 ! for the Nash-Sutcliffe efficiency
 nsmean = nsmean/Gx
 nsmeandiff = 0.
 do ti = 1,Nt
  if ((obsv(ti).gt.-9990).and.(modl(ti).gt.-9990)) then
   nsmeandiff = nsmeandiff + ((obsv(ti) - nsmean)**2)
  endif 
 enddo

 ! ------------------------------------------------------------------
 ! write to output file - Mean Sum of Squared Error
 !open(fid,file='noahmp_objfun.out')
 !do xi = 1,Dx
 ! write(fid,*) sse(xi)/Gx(xi),','
 !enddo
 !close(fid)

 ! ------------------------------------------------------------------
 ! write to output file - Root Mean Sum of Squared Error
 open(fid,file='noahmp_objfun.out')
 write(fid,*) sqrt(sse/Gx),','
  print*, 'RMSE = ', sqrt(sse/Gx)
 close(fid)

 ! ------------------------------------------------------------------
 ! write to output file - AS NASH-SUTCLIFFE EFFICIENCY
 ! NOTE: can't do 1-(sse/nsmeandiff) because ostrich only minimizes the obj.
 ! So I just removed the 1-
 !open(fid,file='noahmp_objfun.out')
 !do xi = 1,Dx
 ! write(fid,*) (sse(xi)/nsmeandiff(xi)),','
 !enddo
 !close(fid)

 ! ------------------------------------------------------------------
 return
end subroutine
