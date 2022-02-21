module smatrix_core_coupled_class
  use math
  use bessel
  use smatrix_core_class
  implicit none
  private

  type,public,extends(smatrix_core) :: smatrix_core_coupled
     !> backgroundの流体の密度
     real(8) :: rho
     !> backgroundの流体の体積弾性率
     real(8) :: kappa
     !> backgroundの流体の波数
     real(8) :: k
     !> この層のLame定数
     real(8) :: dl
     !> この層のLame定数
     real(8) :: dm
     !> この層の密度
     real(8) :: rho_e
     !> この層の縦波の波数
     real(8) :: kL
     !> この層の横波の波数
     real(8) :: kT
     
     ! real(8) :: w     
     ! ! 半径
     ! real(8) :: R
     ! ! S行列のサイズ
     ! integer :: m_min
     ! integer :: m_max
     ! ! S行列
     ! complex(8),allocatable :: S(:,:)
   contains
     procedure :: del
     procedure :: new

     procedure :: T11
     procedure :: T12
     procedure :: T41
     procedure :: T42
     procedure :: Tt11
     procedure :: Tt12
     procedure :: Tt41
     procedure :: Tt42
     procedure :: W11
     procedure :: W12
     procedure :: W21
     procedure :: W22
     procedure :: Wt11
     procedure :: Wt12
     procedure :: Wt21
     procedure :: Wt22
  end type smatrix_core_coupled

  

contains
  elemental subroutine del(self)
    class(smatrix_core_coupled),intent(inout) :: self

    if(allocated(self%S)) deallocate(self%S)

    
    
  end subroutine del

  subroutine new(self, w, R, rho, kappa, rho_e, dl, dm, m_min, m_max)
    class(smatrix_core_coupled),intent(inout) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: R
    real(8),intent(in) :: rho
    real(8),intent(in) :: kappa
    real(8),intent(in) :: rho_e
    real(8),intent(in) :: dl
    real(8),intent(in) :: dm
    integer,intent(in) :: m_min
    integer,intent(in) :: m_max
    
    integer :: n
    complex(8) :: X(3,3), X_inv(3,3)

    call self%del

    self%w = w
    self%R = R
    self%rho = rho
    self%kappa = kappa
    self%rho_e = rho_e
    self%dl = dl
    self%dm = dm
    
    self%m_min = m_min
    self%m_max = m_max

    self%k = w * sqrt(self%rho/self%kappa)
    self%kL = w * sqrt(self%rho_e/(self%dl+2*self%dm))
    self%kT = w * sqrt(self%rho_e/(          self%dm))

    ! self%w = w
    
    ! self%rho_int = rho_int
    ! self%kappa_int = kappa_int
    ! self%rho_ext = one
    ! self%kappa_ext = one
    ! self%R = R
    ! self%m_min = m_min
    ! self%m_max = m_max
    
    ! self%k_int = w * sqrt(self%rho_int / self%kappa_int)
    ! self%k_ext = w * sqrt(self%rho_ext / self%kappa_ext)


    allocate(self%S(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S = zero

    do n=self%m_min,self%m_max

       X(1,:) = [self%T11(n,self%R), self%T12(n,self%R), besh(n,self%k*self%R)]
       X(2,:) = [self%W11(n,self%R), self%W12(n,self%R), -self%k/(self%rho*self%w**2)*besh_d(n,self%k*self%R)]
       X(3,:) = [self%T41(n,self%R), self%T42(n,self%R), zero]

       X_inv = mat_inv(3,X)

       self%S(n:n,n) = matmul(X_inv(3:3,1:2),[-besj(n,self%k*self%R),self%k/(self%rho*self%w**2)*besj_d(n,self%k*self%R)])
    end do
    
  end subroutine new

  complex(8) function T11(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    T11 = (n**2+n-0.5d0*(self%kT*r)**2) * besj(n,self%kL*r) - self%kL*r*besj(n-1,self%kL*r)

    T11 = T11 * 2*self%dm/r**2
  end function T11

  complex(8) function T12(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    T12 = -ione*n*( (n+1)*besj(n,self%kT*r) - self%kT*r*besj(n-1,self%kT*r)  )

    T12 = T12 * 2*self%dm/r**2
  end function T12

  complex(8) function T41(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    T41 = -ione*n*( (n+1)*besj(n,self%kL*r) - self%kL*r*besj(n-1,self%kL*r)  )
    
    T41 = T41 * 2*self%dm/r**2
  end function T41

  complex(8) function T42(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    T42 = -(n**2+n-0.5d0*(self%kT*r)**2) * besj(n,self%kT*r) + self%kT*r*besj(n-1,self%kT*r)

    T42 = T42 * 2*self%dm/r**2
  end function T42


  complex(8) function Tt11(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Tt11 = (n**2+n-0.5d0*(self%kT*r)**2) * besh(n,self%kL*r) - self%kL*r*besh(n-1,self%kL*r)

    Tt11 = Tt11 * 2*self%dm/r**2
  end function Tt11

  complex(8) function Tt12(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Tt12 = -ione*n*( (n+1)*besh(n,self%kT*r) - self%kT*r*besh(n-1,self%kT*r)  )

    Tt12 = Tt12 * 2*self%dm/r**2
  end function Tt12

  complex(8) function Tt41(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Tt41 = -ione*n*( (n+1)*besh(n,self%kL*r) - self%kL*r*besh(n-1,self%kL*r)  )

    Tt41 = Tt41 * 2*self%dm/r**2
  end function Tt41

  complex(8) function Tt42(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Tt42 = -(n**2+n-0.5d0*(self%kT*r)**2) * besh(n,self%kT*r) + self%kT*r*besh(n-1,self%kT*r)

    Tt42 = Tt42 * 2*self%dm/r**2
  end function Tt42

  complex(8) function W11(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W11 = self%kL * besj_d(n,self%kL*r)
  end function W11

  complex(8) function W12(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W12 = ione*n/r * besj(n,self%kT*r)
  end function W12

  complex(8) function W21(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W21 = ione*n/r * besj(n,self%kL*r)
  end function W21

  complex(8) function W22(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W22 = -self%kT * besj_d(n,self%kT*r)
  end function W22


  complex(8) function Wt11(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt11 = self%kL * besh_d(n,self%kL*r)
  end function Wt11

  complex(8) function Wt12(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt12 = ione*n/r * besh(n,self%kT*r)
  end function Wt12

  complex(8) function Wt21(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt21 = ione*n/r * besh(n,self%kL*r)
  end function Wt21

  complex(8) function Wt22(self, n, r)
    class(smatrix_core_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt22 = -self%kT * besh_d(n,self%kT*r)
  end function Wt22

  ! A^{-1}
  function mat_inv(n, A) result(out)
    integer,intent(in) :: n
    complex(8),intent(in) :: A(n,n)
    complex(8) :: out(n,n)

    complex(8) :: A_copy(n,n)

    integer :: info
    integer,allocatable :: ipiv(:)
    integer :: i

    out(:,:) = zero
    do i=1,n
       out(i,i) = one
    end do

    A_copy = A

    allocate(ipiv(n))
    ! LU分解
    call zgetrf(n, n, A_copy, n, ipiv, info)
    if(info /= 0) then
       stop "info /= 0"       
    end if
    call zgetrs('N', n, n, A_copy, n, ipiv, out, n, info)
    if(info /= 0) stop "info /= 0"
    
  end function mat_inv
  
end module smatrix_core_coupled_class
