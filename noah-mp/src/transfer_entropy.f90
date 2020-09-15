subroutine transfer_entropy(Ys,Xt,Yt,N,dx,dy,TE,HY,MI)
 implicit none

 ! input data vectors
 integer, intent(in) :: N
 real, dimension(N), intent(in) :: Xt, Ys, Yt

 ! histogram bins
 real, intent(in) :: dx, dy
 integer :: ni, ii, jj, kk
 integer :: Mx, Ms, Mt
 real, allocatable, dimension(:) :: Bx, Bt, Bs

 ! histograms
 real, allocatable, dimension(:,:,:) :: Pxst 
 real, allocatable, dimension(:,:)   :: Pxs, Pxt, Pst 
 real, allocatable, dimension(:)     :: Px, Ps, Pt 

 ! entropies
 real :: Hx, Hs, Ht, Hxs, Hxt, Hst, Hxst, Ixs, Ixt, Ist, Ixst
 real, intent(out) :: TE, HY, MI

 ! check for missing data
 do ni = 1,N
  if (Xt(ni) < -9990) then
   print*, 'missing values in Xt'
   stop 1000
  endif
  if (Ys(ni) < -9990) then
   print*, 'missing values in Ys'
   stop 2000
  endif
  if (Xt(ni) < -9990) then
   print*, 'missing values in Yt'
   stop 3000
  endif
 enddo

 ! set histogram bins
 Mx = ceiling((maxval(Xt) - minval(Xt) + 2e-6)/dx) 
!print*,'Xt',dx,maxval(Xt),minval(Xt),Mx
 allocate(Bx(Mx))
 do ni = 1,Mx
  Bx(ni) = minval(Xt) - 1e-6 + (ni-1)*dx
 enddo

 Ms = ceiling((maxval(Ys) - minval(Ys) + 2e-6)/dy)  
!print*,'Ys',dy,maxval(Ys),minval(Ys),Ms
 allocate(Bs(Ms))
 do ni = 1,Ms
  Bs(ni) = minval(Ys) - 1e-6 + (ni-1)*dy
 enddo

 Mt = ceiling((maxval(Yt) - minval(Yt) + 2e-6)/dy)  
!print*,'Yt',dy,maxval(Yt),minval(Yt),Mt
 allocate(Bt(Mt))
 do ni = 1,Mt
  Bt(ni) = minval(Yt) - 1e-6 + (ni-1)*dy
 enddo

 ! create joint (and marginal) histogram(s)
 allocate(Pxst(Mx,Ms,Mt))
 allocate(Pxs(Mx,Ms))
 allocate(Pxt(Mx,Mt))
 allocate(Pst(Ms,Mt))
 allocate(Px(Mx))
 allocate(Ps(Ms))
 allocate(Pt(Mt))

 call hist3(Xt,Ys,Yt,Bx,Bs,Bt,N,Mx,Ms,Mt,Pxst,Pxs,Pxt,Pst,Px,Ps,Pt)

 ! error checking
 if (abs(sum(Px)-1).ge.1./N) then 
  print*, 'Px does not sum to 1'
  print*,sum(Px),abs(sum(Px)-1)
  stop 1111
 endif

 if (abs(sum(Ps)-1).ge.1./N) then 
  print*, 'Ps does not sum to 1'
  print*,sum(Ps),abs(sum(Ps)-1)
  stop 2111
 endif

 if (abs(sum(Pt)-1).ge.1./N) then 
  print*, 'Pt does not sum to 1'
  print*,sum(Pt),abs(sum(Pt)-1)
  stop 3111
 endif

 ! calculate entropies
 Hx = 0.
 do ii = 1,Mx
  if (Px(ii).gt.0.) then
   Hx = Hx -Px(ii)*log(Px(ii))
  endif
 enddo

 Hs = 0.
 do ii = 1,Ms
  if (Ps(ii).gt.0.) then
   Hs = Hs -Ps(ii)*log(Ps(ii))
  endif
 enddo

 Ht = 0.
 do ii = 1,Mt
  if (Pt(ii).gt.0.) then
   Ht = Ht -Pt(ii)*log(Pt(ii))
  endif
 enddo

 Hxs = 0.
 do ii = 1,Mx
  do jj = 1,Ms
   if (Pxs(ii,jj).gt.0.) then
    Hxs = Hxs -Pxs(ii,jj)*log(Pxs(ii,jj))
   endif
  enddo
 enddo

 Hxt = 0.
 do ii = 1,Mx
  do jj = 1,Mt
   if (Pxt(ii,jj).gt.0.) then
    Hxt = Hxt -Pxt(ii,jj)*log(Pxt(ii,jj))
   endif
  enddo
 enddo

 Hst = 0.
 do ii = 1,Ms
  do jj = 1,Mt
   if (Pst(ii,jj).gt.0.) then
    Hst = Hst -Pst(ii,jj)*log(Pst(ii,jj))
   endif
  enddo
 enddo

 Hxst = 0.
 do ii = 1,Mx
  do jj = 1,Ms
   do kk = 1,Mt
    if (Pxst(ii,jj,kk).gt.0.) then
     Hxst = Hxst -Pxst(ii,jj,kk)*log(Pxst(ii,jj,kk))
    endif
   enddo
  enddo
 enddo

 ! mutual information values
 Ixs = Hs+Hx-Hxs
 Ist = Hs+Ht-Hst
 Ixt = Hx+Ht-Hxt
 Ixst = -(Hs+Hx+Ht-Hxst-Ixs-Ist-Ixt);

 ! output values
 TE = Ixs-Ixst
 HY = Hs
 MI = Ixs

 ! deallocate everything
 deallocate(Bx)
 deallocate(Bs)
 deallocate(Bt)
 deallocate(Pxst)
 deallocate(Pxs)
 deallocate(Pxt)
 deallocate(Pst)
 deallocate(Px)
 deallocate(Ps)
 deallocate(Pt)

end subroutine
