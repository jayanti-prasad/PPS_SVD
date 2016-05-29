program cl2pk_svd 
  implicit none 
  integer, parameter :: M = 2499, N = 3914 
  double precision, parameter :: pi =  4.0d0*atan(1.0d0)
  double precision, parameter :: outputscale=7.42835025d12
  double precision, dimension(N,M) :: A_IN,SS_IN,VSS_IN
  double precision, dimension(M,N) :: SS,TMP1,B
  double precision, dimension(M,M) :: U,UT
  double precision, dimension(N,N) :: VT,V
  double precision, dimension(M) :: S,ll, Cl_TT,Cl_EE,Cl_TE,cl1,cl   
  double precision, dimension(N) :: kk,dk,Pk,fk,fk1 
  character(len=70) :: string,svd_file,pk_file,cl_file,ckf,pkf   
  integer :: i,j,l,nsvd,opt,numarg  
  double precision :: coeff,ctnorm,fact     

  if ( iargc() < 1) then
     write(*,*)"./cl2pk_svd  # svd components"
     stop
  end if

  call getarg(1,string)
  read(string,'(I4)') nsvd
  
  ! user input goes here 

  opt=1  ! 1=TT, 2=EE, 3=TE 
  open(11, status='unknown',file='../../input_data/SVD1.DAT',form='unformatted')
  open(22,file='../../input_data/wmap9_Pk.dat') ! this is read for reference only 
  open(33,file='../../input_data/wmap9_Cls.dat')
  
  ! do not change anything below   
 
  write(pkf,'(A4,I4,A4)')"pps_",1000*opt+nsvd,".dat"
  write(ckf,'(A4,I4,A4)')"cls_",1000*opt+nsvd,".dat"
 
  write(*,*)"# svd=",nsvd  
  write(*,*)"P(k) file : ",pkf 
  write(*,*)"Cl   file : ",ckf 
 
  read(11)S
  read(11)U
  read(11)VT
  close(11)
  
  do  i=1, N
     read(22,*) kk(i),dk(i),pk(i)
  end do
  close(22)
  
  do l  =  1, M 
     read(33,*) ll(l),Cl_TT(l),Cl_EE(l),CL_TE(l)
     coeff = 2.0d0 * ll(l) *(ll(l)+1.0)* outputscale
     ctnorm=(ll(l)*ll(l)-1)*(ll(l)+2)*ll(l)
     CL_TT(l) = Cl_TT(l)/coeff 
     CL_EE(l) = Cl_EE(l)/(coeff*ctnorm) 
     Cl_TE(l) = Cl_TE(l)/(coeff*sqrt(ctnorm))
  end do
  close(33)
  
  SELECT CASE (opt)
  CASE (1)
     Cl = Cl_TT 
  CASE (2)
     Cl = Cl_EE 
  CASE (3)
     Cl = Cl_TE 
  END SELECT
  
  B(:,:) = 0.d0
  SS(:,:) = 0.D0
  
  do  i =  1, M
     if ( i .lt. nsvd) then
        SS(i,i) = S(i)
     end if
  end do
  
  TMP1 = MATMUL(U,SS)
  B    = MATMUL(TMP1,VT)
  
  fk = pk *dk 
  
  cl1 = MATMUL(B,fk)
  
  open(12,file=trim(ckf))
  do l  =  1, M 
     coeff = 2.0d0 * ll(l) *(ll(l)+1.0)* outputscale
     ctnorm=(ll(l)*ll(l)-1)*(ll(l)+2)*ll(l)
     SELECT CASE (opt)
     CASE (1)
        fact = coeff 
     CASE (2)
        fact =  coeff * ctnorm
     CASE (3)
        fact = coeff * sqrt(ctnorm)
     END SELECT
     write(12,*) ll(l), fact * cl1(l)  
  end do
  close(12)
  
  SS_IN(:,:) = 0.D0
  
  do i = 1, M
     IF (S(I) .gt. S(nsvd)) then
        SS_IN(i,i) = 1.0D0/S(I)
     END IF
  end do
  
  V      =transpose(VT)
  UT     =transpose(U)
  VSS_IN =MATMUL(V,SS_IN)
  
  A_IN = MATMUL(VSS_IN,UT)
  
  fk1 = MATMUL(A_IN,Cl)
  
  open(11,file=pkf)
  do i  = 1, N
     write(11,*)kk(I),pk(i),fk1(i)/dk(i)
  end do
  close(11)
  
end program cl2pk_svd

