
  !   real   rQmatrix(MaxNumbDegrees,MaxNumbDegrees),  rRmatrix(MaxNumbDegrees,MaxNumbDegrees)
  !   real   rInvMatrix(MaxNumbDegrees,MaxNumbDegrees), rMa(MaxNumbDegrees,MaxNumbDegrees)

  !   ! perform the QR decomposition of the matrix Ma
  !   CALL QRDecomposition(rMA, rQmatrix, rRmatrix, MaxNumbDegrees)
	! ! invert the QR matrix: rInvMatrix = (QR)^(-1)
  !   Call InvertQR(rQmatrix, rRmatrix, rInvMatrix,MaxNumbDegrees)  
  !   !    test = matmul(rQmatrix,rRmatrix) -rMA


 !> \brief L2 norm of a vector
 !! \param[in] x input vector
 !! \param[in] n length of the vector x
module qr_decomposition

  contains

  real  function l2norm(x,n)
   implicit none 
   integer, intent(in) :: n
   real, intent(in) ::  x(n)
   l2norm = sqrt(dot_product(x,x))
  end function
 
 !\brief
  real  function l2norm_vec3(x)
   real, intent(in) ::  x(1:3)    
   l2norm_vec3 = sqrt(dot_product(x,x))   
 end function
   
 
 !%%% QR decomposition stuff %%%
!
!!>\brief invert square matrix
!Subroutine InvertMaxtrix(N,A,InvA,deg_num)
! integer,intent(in) :: N,deg_num
! real,dimension(1:n,1:n) :: A,InvA,R,Q
! 
! ! perform the QR decomposition of the matrix a
!  CALL QRDecomposition(A, Q, R, deg_num)
! ! invert the QR matrix: rInvMatrix = (QR)^(-1)
!  Call InvertQR(Q, R, InvA,deg_num)  
!End subroutine
    
  !> \brief construct a householder vector v that  annihilates all but the first component of x
  !! \param[in] x input x
  !! \param[in] n dimensions
  !! \param[out] v resulting vector
subroutine house(x,v,n)
    integer, intent(in) :: n
    real,intent(in)  ::  x(n)
    real,intent(out) ::  v(n)
    v = x
    v(1) = x(1) + sign(1.0,x(1))*l2norm(x,n)    
End subroutine


  !> \brief construct a householder reflection matrix from a Householder vector v 
  !! \param[in] n dimension
  !! \param[in] v input vector
  !! \param[out] P resulting matrix
  subroutine ComputeHouseMatrix(P,v,n)
    integer, intent(in) :: n
    real, intent(in) ::  v(n)
    real,intent(out) ::  P(n,n)
    real vnorm
    integer i,j
    P = 0.0
    do i=1,n
      P(i,i) = 1.0
    enddo   
    vnorm = l2norm(v,N)
    vnorm = 1./vnorm**2
    do i=1,n
    do j=1,n
      p(i,j) = p(i,j) - 2*v(i)*v(j)*vnorm
    enddo
    enddo    
  End subroutine

 !> \brief Compute the inverse of QR decomposition
  !! \param[in]  Q Q  matrix
 !! \param[in]  R R  matrix
 !! \param[out] InvertedQR inverted QR matrix
 !! \param[in]  N dimension
  Subroutine  InvertQR(Q,R,InvertedQR,N)
    integer, intent(in) :: n
    real, intent(in) ::  Q(N,N),R(N,N)
    real, intent(out) :: InvertedQR(N,N)
    real  InvR(N,N),QT(N,N)
 
    ! invert R
    CALL InvertRmaxtrix(R,InvR,N)
    ! transpose Q
    Qt = transpose(Q)
    !  final inverted R^(-1)*Q^(-1)
    InvertedQR = matmul(InvR,QT)
 End subroutine 

 !> \brief QR decomposition
 !! \param[in]  A input matrix
 !! \param[in]  Q Q  matrix
 !! \param[in]  R R  matrix
 !! \param[out] InvertedQR inverted QR matrix
 !! \param[in]  N dimension
  subroutine  QRDecomposition_slow(A,Q,R,N)
    integer, intent(in) :: N
    real , intent(in) :: A(N,N)
    real , intent(out) ::  Q(N,N),R(N,N)
    real  Identity(N,N),V(N),P(N,N)
    integer i,j,l,pdim
    real, allocatable :: v1(:),x(:)
 
    Identity = 0.D0
    do i=1,n
      Identity(i,i) = 1.d0
    enddo

    Q = Identity
    R = A

    do l=1,n
     ! allocate vector and reflection matrix
     pdim = n-l+1
     allocate(v1(pdim),x(pdim)) 
     x = R(L:N,L)
     ! compute the partial vector    
     CALL House(x,v1,pdim)
     v = 0.
     v(L:N) = v1
     ! compute the reflection matrix
     CALL ComputeHouseMatrix(P,v,N)
     ! construct the Q(l) matrix
     R = MATMUL(P,R)
     Q = MATMUL(Q,P)
     deallocate(x,v1)
    enddo

    ! just in case - make 
    do i=1,N
      do j=1,i-1
	   R(i,j) =0.
	 enddo
   enddo 
 End subroutine

 !>\brief
  subroutine  QRDecomposition(A,Q,R,N) ! _fast
    integer, intent(in) :: N
    real, intent(in) :: A(N,N)
    real, intent(out) ::  Q(N,N),R(N,N)
    real   V(N),P(N,N),v1(1:N),x(1:N)
    integer i,j,l,pdim,m
 
    Q = 0.D0
    do i=1,n
      Q(i,i) = 1.d0
    enddo    
    v1 = 0.
    x  = 0.
    R = A
    do l=1,n
     ! allocate vector and reflection matrix
     pdim = n-l+1
     !allocate(v1(pdim),x(pdim)) 
     do m=1,pdim
       !x(1:pdim) = R(L:N,L)
        x(m) = R(L+m-1,L)
     enddo
     ! compute the partial vector    
     CALL House(x(1:pdim),v1(1:pdim),pdim)
     v = 0.
     do m=1,pdim
       !v(L:N) = v1(1:pdim)
       v(L+m-1) = v1(m)
     enddo
     ! compute the reflection matrix
     CALL ComputeHouseMatrix(P,v,N)
     ! construct the Q(l) matrix
     R = MATMUL(P,R)
     Q = MATMUL(Q,P)
     !deallocate(x,v1)
    enddo

    ! just in case - make 
    do i=1,N
      do j=1,i-1
	   R(i,j) =0.
	 enddo
   enddo
 
  End subroutine

 !> \brief Invert R matrix
 !! \param[in] R input R matrix  
 !! \param[out] InvR  inverse of R matrix
 !! \param[in] N  dimension of the matrix
  Subroutine InvertRmaxtrix(R,InvR,N)
  integer,intent(in) :: N
  real , intent(in) :: R(N,N)
  real , intent(out):: InvR(N,N)
  integer i,j,k
 !  implicit real  (t)
 
  InvR = 0.d0

  do i=N,1,-1
    invr(i,i) = 1./r(i,i)
	do j=i+1,N
	  invr(i,j) = 0.
	  do k= 1,j-1
	   invr(i,j) = invr(i,j) - r(k,j)*invr(i,k)
	  enddo
 	  invr(i,j) =invr(i,j) /r(j,j)
	enddo
  enddo
 End subroutine
 
end module