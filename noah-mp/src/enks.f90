subroutine enks(X,Y,Zbar,Zsig,Nt,Ne,Dx,Dz,Threshold_on_CC)

  implicit none
 
 ! in/out 
  integer, intent(in)                      :: Nt ! This is actually the number of lagged time steps considered
  integer, intent(in)                      :: Ne ! Number of ensemble numbers
  integer, intent(in)                      :: Dz ! Number of observed variables at current time step
  integer, intent(in)                      :: Dx ! Number of modeled states at lagged time steps that need to be adjusted based on the information from the observed variables at current time step
  real, dimension(Dz), intent(in)          :: Zbar ! Vector of observed variables at current time step
  real, dimension(Dz), intent(in)          :: Zsig ! Vector of standard deviations corresponding to the observed variables at current time step. This is read in from obs_cov.txt, but the 'cov' apparently denoting a covariance matrix is a misnomer since this Zsig is only a standard deviation related to a variance vector (instead of covariance). Its values are used for the diagonal matrix R below where covariances are 0 (hence only variance-related terms remain and using the term 'variance' is probably ok).
                                                   ! Alternately Zsig can be vector of relative standard deviations. Comment out and uncomment the relevant 2 lines further below as per what is actually used. 
  real, dimension(Nt,Dx,Ne), intent(inout) :: X ! Model states that need to be adjusted
  real, dimension(Dz,Ne), intent(in)       :: Y ! Modeled variables (states of fluxes) at this time step that directly correspond to the observed variables
 
 ! constants
  real, dimension(Dz,Dz) :: eyeZ ! Identity matrix 
 
 ! manipulation vectors
  real, dimension(Dx,Ne) :: Xbar ! Deviations of X at the considered time lag from respective ensemble averages
  real, dimension(Dx,Ne) :: KK 
  real, dimension(Dz,Ne) :: Ybar ! Deviation of Y from respective ensemble averages
  real, dimension(Dz,Ne) :: Z ! Ensemble created from Zbar
  real, dimension(Dz,Dz) :: R ! Covariance matrix version of the Zsig standard deviations vector
  real, dimension(Dz,Dz) :: Cyy ! Matrix of covariances among observed variables
  real, dimension(Dz,Dz) :: Qinv
  real, dimension(Dx,Dz) :: Cxy, K
  real, dimension(Ne,Ne) :: X5
  real, dimension(Ne,Dz) :: X4
 
  real, dimension(Dx,Dx) :: Cxx
  real, dimension(Dx,Dz) :: CC
 
  real, intent(in)                      :: Threshold_on_CC ! If this is not -9999, then for CC terms less than this value, set corresponding K terms to 0 
 
 ! indexes
  integer e, t, d, i
  integer, dimension(Dz) :: ipiv
  integer :: info
  real :: eta, random_normal
 
 ! identity matrix
  eyeZ = 0.
  do d = 1,Dz
    eyeZ(d,d) = 1.
  enddo
 
 ! creating ensembles for the observation variables at this time step
  do e = 1,Ne 
    do d = 1,Dz
      eta = random_normal()
   !   Z(d,e) = Zbar(d) * (1 + eta*Zsig(d)) ! For Zsig being relative standard deviations. If we want to use Zsig being standard deviations, then use the line immediately below instead 
      Z(d,e) = Zbar(d) + eta*Zsig(d) ! For Zsig being standard deviations. If we want to use Zsig being relative standard deviations, then use the line immediately above instead
    enddo
  enddo
 
 ! compute bar matrices (i.e. deviations from ensemble means)
  do d = 1,Dz
    Ybar(d,:) = Y(d,:) - sum(Y(d,:))/Ne ! Ybar dimensions are (Dz,Ne) 
  enddo
 
 ! Begin computing cross-covariances
 
 !Between modeled variables that correspond to observed ones
  Cyy = 0.
  do d = 1,Dz
    do e = 1,Dz
      do i = 1,Ne 
        Cyy(d,e) = Cyy(d,e) + Ybar(d,i)*Ybar(e,i) ! Cyy dimensions are (Dz,Dz) 
      enddo
      Cyy(d,e) = Cyy(d,e)/(Ne-1) 
    enddo
  enddo
 
 !Between observed variables (but cross- terms are assumed zero so that this will have variances only) 
  R = 0.
  do d = 1,Dz
    R(d,d) = Zsig(d)**2
  enddo
 
 ! End computing cross-covariances
 
 ! invert covariance matrix
  Qinv = Cyy+R ! dimensions are (Dz,Dz)
  call sgesv(Dz,Dz,Qinv,Dz,ipiv,eyeZ,Dz,info)
  ! SY: Begin code lines to check if sgesv returned correct value, else halt program
  IF (info .NE. 0) THEN
    PRINT*,'******sgesv returns an info value of ', info, '******'
    STOP '*****sgesv did not exit successfully*****'
  END IF
  ! SY: End code lines to check if sgesv returned correct value, else halt program
  Qinv = eyeZ ! dimensions are (Dz,Dz)
 
 ! compute Kalman Gain
  Z = Z-Y
 
  IF (Threshold_on_CC .LT. -9998) THEN  
  
    X4 = matmul(transpose(Ybar),Qinv) ! dimensions of X4 and transpose(Ybar) are (Ne, Dz), and of Qinv are (Dz,Dz)
    X5 = matmul(X4,Z)/(Ne-1) ! dimensions of X5 are (Ne, Ne)

  ELSE IF (Threshold_on_CC .LT. -1) THEN  

    PRINT*,'*******Threshold_on_CC needs to be >= -1, currently it is ', Threshold_on_CC,'******'
    STOP '*******Threshold_on_CC needs to be >= -1******'

  END IF ! IF (Threshold_on_CC .LT. -9998) THEN  

  do t = 1,Nt ! Computations at each time lag
 
    IF (Threshold_on_CC .LT. -9998) THEN  
  
      ! perform update 
      X(t,:,:) = X(t,:,:) + matmul(Xbar,X5) 
    
    ELSE ! Threshold_on_CC >= -1 case 
  
      ! compute bar matrices (i.e. deviations from ensemble means) at the considered time lag
      do d = 1,Dx
        Xbar(d,:) = X(t,d,:) - sum(X(t,d,:))/Ne ! dimensions of Xbar are (Dx,Ne)
      enddo
    
      ! Begin computing cross-covariances
    
      !Between modeled states (at the considered time lag) and modeled variables (that correspond to observed ones)
      Cxy = 0.
      do d = 1,Dx
        do e = 1,Dz
          do i = 1,Ne
            Cxy(d,e) = Cxy(d,e) + Xbar(d,i)*Ybar(e,i) ! Cxy dimensions are (Dx,Dz)
          enddo
          Cxy(d,e) = Cxy(d,e)/(Ne-1)
        enddo
      enddo
    
      !Between modeled states (at the considered time lag) 
      Cxx = 0.0
      do d = 1,Dx
        do e = 1,Dx
          do i = 1,Ne
            Cxx(d,e) = Cxx(d,e) + Xbar(d,i)*Xbar(e,i) ! Cxx dimensions are (Dx,Dx)
          enddo
          Cxx(d,e) = Cxx(d,e)/(Ne-1)
        enddo
      enddo
    
      ! End computing cross-covariances
    
      ! compute cross-correlations
      do d=1,Dx
        do e = 1,Dz
          CC(d,e) = Cxy(d,e)/sqrt(Cxx(d,d))/sqrt(Cyy(e,e)) ! CC dimensions are (Dx,Dz)
        enddo
      enddo
    
      K = matmul(Cxy,Qinv) ! K dimensions are (Dx,Dz)    
      do d = 1,Dx
        do e = 1,Dz
          if (abs(CC(d,e)) .lt. Threshold_on_CC) K(d,e) = 0.0
        enddo
      enddo
      KK = matmul(K,Z) 
  
      ! perform update
      X(t,:,:) = X(t,:,:) + KK
 
    END IF ! IF (Threshold_on_CC .LT. -9998) THEN  
   
  enddo ! do t = 1,Nt ! Computations at each time lag
 
  return

end subroutine enks




