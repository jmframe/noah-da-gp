subroutine hist3(X,Y,Z,Bx,By,Bz,N,Mx,My,Mz,Pxyz,Pxy,Pxz,Pyz,Px,Py,Pz)
 implicit none

 integer, intent(in) :: N, Mx, My, Mz
 real, dimension(Mx), intent(in) :: Bx
 real, dimension(My), intent(in) :: By
 real, dimension(Mz), intent(in) :: Bz

 real, dimension(N) :: X, Y, Z

 real, dimension(Mx,My,Mz) :: Pxyz 
 real, dimension(Mx,My)    :: Pxy 
 real, dimension(Mx,Mz)    :: Pxz 
 real, dimension(My,Mz)    :: Pyz 
 real, dimension(Mx)       :: Px 
 real, dimension(My)       :: Py 
 real, dimension(Mz)       :: Pz 

 integer :: ni, xi, yi, zi

 real, dimension(My) :: B

 Pxyz = 0.
 do ni = 1,N
  xi = maxloc(Bx-X(ni),1,mask=(X(ni).gt.Bx))
  yi = maxloc(By-Y(ni),1,mask=(Y(ni).gt.By))
  zi = maxloc(Bz-Z(ni),1,mask=(Z(ni).gt.Bz))
  !print*, Bx, '---', X(ni), xi
  !print*, By, '---', Y(ni), yi
  !print*, Bz, '---', Z(ni), zi
  Pxyz(xi,yi,zi) = Pxyz(xi,yi,zi) + 1.
 enddo

 ! check summation
 if (sum(sum(sum(Pxyz,3),2),1).ne.N) then
  print*,sum(sum(sum(Pxyz,3),2),1)
  stop 10031
 endif 

 ! marginalise for error checking
 Pxy = sum(Pxyz,3)
 Pxz = sum(Pxyz,2)
 Pyz = sum(Pxyz,1)
 Px  = sum(Pxy,2)
 Py  = sum(Pxy,1)
 Pz  = sum(Pxz,1)

 ! check summation
 if (sum(sum(Pxy,2),1).ne.N) then
  print*,sum(sum(Pxy,2),1)
  stop 10021
 endif 

 ! check summation
 if (sum(sum(Pxz,2),1).ne.N) then
  print*,sum(sum(Pxz,2),1)
  stop 10022
 endif 

 ! check summation
 if (sum(sum(Pyz,2),1).ne.N) then
  print*,sum(sum(Pyz,2),1)
  stop 10023
 endif 

 ! check summation
 if (sum(Px,1).ne.N) then
  print*,sum(Px,1)
  stop 10011
 endif 

 ! check summation
 if (sum(Py,1).ne.N) then
  print*,sum(Py,1)
  stop 10012
 endif 

 ! check summation
 if (sum(Pz,1).ne.N) then
  print*,sum(Pz,1)
  stop 10013
 endif 

 ! normalize
 Pxyz = Pxyz/N
 Pxy  = Pxy/N
 Pxz  = Pxz/N
 Pyz  = Pyz/N
 Px   = Px/N
 Py   = Py/N
 Pz   = Pz/N

! ! ensure numerical consistency
! Pxy = Pxy/sum(sum(Pxy,2),1)
! Pxz = Pxz/sum(sum(Pxz,2),1)
! Pyz = Pyz/sum(sum(Pyz,2),1)
! Px = Px/sum(Px)
! Py = Py/sum(Py)
! Pz = Pz/sum(Pz)

end subroutine
