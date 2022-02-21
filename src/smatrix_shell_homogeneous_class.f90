module smatrix_shell_homogeneous_class
  use smatrix_shell_class
  use math
  use misc
  use bessel
  implicit none
  private

  type,public,extends(smatrix_shell) :: smatrix_shell_homogeneous
     !> この層の比密度
     real(8) :: rho
     !> この層の比体積弾性率
     real(8) :: kappa
     !> この層の波数
     real(8) :: k
     !> この層の屈折率
     real(8) :: n
     !> この層のインピーダンス
     real(8) :: z

   contains
     procedure :: del
     procedure :: new
     procedure :: calc_derivative
  end type smatrix_shell_homogeneous

contains

  elemental subroutine del(self)
    class(smatrix_shell_homogeneous),intent(inout) :: self

    if(allocated(self%S11)) deallocate(self%S11)
    if(allocated(self%S21)) deallocate(self%S21)
    if(allocated(self%S12)) deallocate(self%S12)
    if(allocated(self%S22)) deallocate(self%S22)
    
  end subroutine del

  subroutine new(self, w, Rout, Rin, rho, kappa, m_min, m_max)
    class(smatrix_shell_homogeneous),intent(inout) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: Rout
    real(8),intent(in) :: Rin
    real(8),intent(in) :: rho
    real(8),intent(in) :: kappa
    integer,intent(in) :: m_min
    integer,intent(in) :: m_max

    integer :: n
    complex(8) :: L1(2,2), L2(2,2), T(2,2)
    
    call self%del

    self%w = w
    self%Rin = Rin
    self%Rout = Rout
    self%rho = rho
    self%kappa = kappa
    self%m_min = m_min
    self%m_max = m_max

    self%k = w * sqrt(self%rho/self%kappa)
    self%n = sqrt(self%rho/self%kappa)
    self%z = sqrt(self%rho*self%kappa)

    allocate(self%S11(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S11 = zero
    allocate(self%S21(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S21 = zero
    allocate(self%S12(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S12 = zero
    allocate(self%S22(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S22 = zero
    do n=self%m_min,self%m_max
       
       ! Lを計算
       L1(1,1) = +self%z*besh_d(n,self%w*self%Rin)*besj(n,self%k*self%Rin) - besh(n,self%w*self%Rin)*besj_d(n,self%k*self%Rin)
       L1(2,1) = -self%z*besj_d(n,self%w*self%Rin)*besj(n,self%k*self%Rin) + besj(n,self%w*self%Rin)*besj_d(n,self%k*self%Rin)
       L1(1,2) = +self%z*besh_d(n,self%w*self%Rin)*besh(n,self%k*self%Rin) - besh(n,self%w*self%Rin)*besh_d(n,self%k*self%Rin)
       L1(2,2) = -self%z*besj_d(n,self%w*self%Rin)*besh(n,self%k*self%Rin) + besj(n,self%w*self%Rin)*besh_d(n,self%k*self%Rin)

       L2(1,1) = +besh_d(n,self%k*self%Rout)*besj(n,self%w*self%Rout)/self%z - besh(n,self%k*self%Rout)*besj_d(n,self%w*self%Rout)
       L2(2,1) = -besj_d(n,self%k*self%Rout)*besj(n,self%w*self%Rout)/self%z + besj(n,self%k*self%Rout)*besj_d(n,self%w*self%Rout)
       L2(1,2) = +besh_d(n,self%k*self%Rout)*besh(n,self%w*self%Rout)/self%z - besh(n,self%k*self%Rout)*besh_d(n,self%w*self%Rout)
       L2(2,2) = -besj_d(n,self%k*self%Rout)*besh(n,self%w*self%Rout)/self%z + besj(n,self%k*self%Rout)*besh_d(n,self%w*self%Rout)

       ! T行列
       T = -pi**2/4*self%w*self%k*self%Rin*self%Rout * matmul(L1,L2)

       ! S行列の(n,n)成分
       self%S11(n,n) = -T(2,1)/T(2,2)
       self%S21(n,n) = T(1,1) - T(1,2)*T(2,1)/T(2,2)
       self%S12(n,n) = 1.d0 / T(2,2)
       self%S22(n,n) = T(1,2) / T(2,2)
    end do
       
  end subroutine new

  subroutine calc_derivative(self, idv, S11_d, S21_d, S12_d, S22_d)
    class(smatrix_shell_homogeneous),intent(in) :: self
    integer,intent(in) :: idv
    complex(8),intent(out) :: S11_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S21_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S12_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S22_d(self%m_min:self%m_max,self%m_min:self%m_max)

    stop "todo"
  end subroutine calc_derivative
end module smatrix_shell_homogeneous_class
