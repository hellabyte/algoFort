subroutine cheb(L,D,x)
! =====================================================================
! - cheb -------------------------------------------------------------- 
!  Computes Chebyshev-Gauss-Lobatto collocation points 
!    (extrema of Chebyshev polynomial of first kind, degree L), and 
!    corresponding differentation matrix D.
!  
!  ASSUMES DOUBLE PRECISION
!
! - References --------------------------------------------------------
!  Created from Trefethen "Spectral Methods in MATLAB", 2000, SIAM.
!    See Chapter 6, specifically pg 54 algo cheb.m
! =====================================================================
  implicit none
  integer,parameter :: SP=kind(0e0),DP=kind(0d0)
  integer,intent(in)                      :: L
  real(DP),dimension(0:L,0:L),intent(out) :: D
  real(DP),intent(out)                    :: x(0:L)
  real(DP) :: c(0:L),ic(0:L),dXI(0:L,0:L),sDT(0:L)
  integer  :: k
  ! = Zero initialization =============================================
  c=0d0; ic=0d0; D=0d0
  sDT = 0d0; dXI = 0d0
  ! = Compute domain, x ===============================================
  forall(k=0:L) x(k) = dcos(pi/L*k)
  ! = Compute derivative operator, D ==================================
  ! - Get barycentric weights, c --------------------------------------
  c=1d0; c(0)=2d0 ; c(L)=2d0;
  do k=1,L,2
    c(k) = -c(k)
  end do
  ! - Get inverse barycentric weights, ic -----------------------------
  ic = 1d0/c
  ! - compute dX + I --------------------------------------------------
  do k = 0,L
    ! - Get spatial metric (antisymmetric) ----------------------------
    dXI(:,k) = x - x(k)
    ! - Add identity to spatial metric --------------------------------
    dXI(k,k) = 1d0
  end do
  ! - compute first pass of D -----------------------------------------
  do k = 0,L 
    D(:,k) = (c(:)*ic(k)) / dXI(:,k)
  end do
  ! - sum along the rows of current D transpose (columns of D) --------
  sDT = sum(D,2)
  ! - Update diagonal of D with SDT result ----------------------------
  do k = 0,L
    D(k,k) = D(k,k) - sDT(k)
  end do
  return
end subroutine cheb
