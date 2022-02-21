module smatrix_unit_donut_class
  use math
  use misc
  use bessel
  use smatrix_unit_class
  implicit none
  private

  type,public,extends(smatrix_unit) :: smatrix_unit_donut
     ! 層の数 >= 1
     integer :: nshell
     ! 半径 (外側から）
     real(8),allocatable :: Rs(:) ! (nshell)
     ! 密度
     real(8),allocatable :: rhos(:) ! (0:nshell)
     ! 体積弾性率
     real(8),allocatable :: kappas(:) ! (0:nshell)
     ! インピーダンスの逆数
     real(8),allocatable :: betas(:) ! (0:nshell)
   contains
     procedure :: del
     procedure :: new
     procedure :: calc_derivative
     procedure :: calc_derivative_radius
  end type smatrix_unit_donut

contains
  elemental subroutine del(self)
    class(smatrix_unit_donut),intent(inout) :: self

    if(allocated(self%S)) deallocate(self%S)

    if(allocated(self%Rs)) deallocate(self%Rs)
    if(allocated(self%rhos)) deallocate(self%rhos)
    if(allocated(self%kappas)) deallocate(self%kappas)
    if(allocated(self%betas)) deallocate(self%betas)
  end subroutine del

  subroutine new(self, w, nshell, Rs, rhos, kappas)
    class(smatrix_unit_donut),intent(inout) :: self
    real(8),intent(in) :: w
    integer,intent(in) :: nshell
    real(8),intent(in) :: Rs(1:nshell)
    real(8),intent(in) :: rhos(1:nshell)
    real(8),intent(in) :: kappas(1:nshell)

    integer :: n
    complex(8) :: T(2,2)
    complex(8) :: M1(2,2), M2(2,2)
    integer :: j
    real(8) :: k, R

    self%w = w
    self%nshell = nshell

    allocate(self%Rs(1:self%nshell))
    self%Rs(:) = Rs(:)
    
    allocate(self%rhos(0:self%nshell))
    self%rhos(0) = 1.d0
    self%rhos(1:self%nshell) = rhos(:)

    allocate(self%kappas(0:self%nshell))
    self%kappas(0) = 1.d0
    self%kappas(1:self%nshell) = kappas(:)

    allocate(self%betas(0:self%nshell))
    self%betas(:) = 1.d0 / sqrt(self%rhos(:)*self%kappas(:))
    
    self%R = self%Rs(1)

    ! 設計変数 = 半径
    self%ndv = self%nshell

    ! Rokhlin
    self%m_max = ceiling(2*w*self%R + 5*log(w*2*self%R+pi))
    self%m_min = -self%m_max

    self%size = self%m_max - self%m_min + 1

    allocate(self%S(self%m_min:self%m_max,self%m_min:self%m_max))
    self%S(:,:) = zero

    !
    ! 各nについてTを計算
    !
    do n=self%m_min,self%m_max

       T(1,1) = 1.d0
       T(2,1) = 0.d0
       T(1,2) = 0.d0
       T(2,2) = 1.d0

       do j=1,self%nshell
          k = self%w * sqrt(self%rhos(j)/self%kappas(j))
          R = self%Rs(j)
          M1(1,1) = +self%betas(j)*besh_d(n,k*R)
          M1(2,1) = -self%betas(j)*besj_d(n,k*R)
          M1(1,2) =               -besh(n,k*R)
          M1(2,2) =               +besj(n,k*R)

          k = self%w * sqrt(self%rhos(j-1)/self%kappas(j-1))
          R = self%Rs(j)
          M2(1,1) =                 +besj(n,k*R)
          M2(2,1) = +self%betas(j-1)*besj_d(n,k*R)
          M2(1,2) =                 +besh(n,k*R)
          M2(2,2) = +self%betas(j-1)*besh_d(n,k*R)

          T = pi*self%w * sqrt(self%rhos(j)/self%kappas(j))*self%Rs(j)/(2*ione*self%betas(j)) * matmul(M1,matmul(M2,T))
       end do

       self%S(n,n) = -T(2,1)/T(2,2)
       
    end do
    
    ! self%S = zero
    ! do m=self%m_min,self%m_max
    !    self%S(m,m) = -besj_d(m,self%w*self%R) / besh_d(m,self%w*self%R)
    ! end do
  end subroutine new

  ! idv番設計変数に関するSの微分
  subroutine calc_derivative(self, idv, out)
    class(smatrix_unit_donut),intent(in) :: self
    integer,intent(in) :: idv
    complex(8),intent(out) :: out(self%m_min:self%m_max,self%m_min:self%m_max)

    call self%calc_derivative_radius(idv, out)
  end subroutine calc_derivative

  ! i番半径に関するS行列の導関数
  subroutine calc_derivative_radius(self, i, out)
    class(smatrix_unit_donut),intent(in) :: self
    integer,intent(in) :: i
    complex(8),intent(out) :: out(self%m_min:self%m_max,self%m_min:self%m_max)

    integer :: n
    complex(8) :: T(2,2), Td(2,2)
    complex(8) :: M1(2,2), M2(2,2), X1(2,2), X2(2,2), X3(2,2)
    integer :: j
    real(8) :: k, R

    out = zero
    
    do n=self%m_min,self%m_max

       !
       ! T
       !
       T(1,1) = 1.d0
       T(2,1) = 0.d0
       T(1,2) = 0.d0
       T(2,2) = 1.d0

       do j=1,self%nshell
          k = self%w * sqrt(self%rhos(j)/self%kappas(j))
          R = self%Rs(j)
          M1(1,1) = +self%betas(j)*besh_d(n,k*R)
          M1(2,1) = -self%betas(j)*besj_d(n,k*R)
          M1(1,2) =               -besh(n,k*R)
          M1(2,2) =               +besj(n,k*R)

          k = self%w * sqrt(self%rhos(j-1)/self%kappas(j-1))
          R = self%Rs(j)
          M2(1,1) =                 +besj(n,k*R)
          M2(2,1) = +self%betas(j-1)*besj_d(n,k*R)
          M2(1,2) =                 +besh(n,k*R)
          M2(2,2) = +self%betas(j-1)*besh_d(n,k*R)

          T = pi*self%w * sqrt(self%rhos(j)/self%kappas(j))*self%Rs(j)/(2*ione*self%betas(j)) * matmul(M1,matmul(M2,T))
       end do

       !
       ! Tの導関数
       !
       Td(1,1) = 1.d0
       Td(2,1) = 0.d0
       Td(1,2) = 0.d0
       Td(2,2) = 1.d0

       do j=1,self%nshell

          if(j == i) then
             ! 第1項
             k = self%w * sqrt(self%rhos(j)/self%kappas(j))
             R = self%Rs(j)
             M1(1,1) = +self%betas(j)*besh_d(n,k*R)
             M1(2,1) = -self%betas(j)*besj_d(n,k*R)
             M1(1,2) =               -besh(n,k*R)
             M1(2,2) =               +besj(n,k*R)

             k = self%w * sqrt(self%rhos(j-1)/self%kappas(j-1))
             R = self%Rs(j)
             M2(1,1) =                 +besj(n,k*R)
             M2(2,1) = +self%betas(j-1)*besj_d(n,k*R)
             M2(1,2) =                 +besh(n,k*R)
             M2(2,2) = +self%betas(j-1)*besh_d(n,k*R)

             X1 = pi*self%w * sqrt(self%rhos(j)/self%kappas(j))/(2*ione*self%betas(j)) * matmul(M1,M2)

             ! 第2項
             k = self%w * sqrt(self%rhos(j)/self%kappas(j))
             R = self%Rs(j)
             M1(1,1) = +self%betas(j)*besh_dd(n,k*R)
             M1(2,1) = -self%betas(j)*besj_dd(n,k*R)
             M1(1,2) =               -besh_d(n,k*R)
             M1(2,2) =               +besj_d(n,k*R)

             k = self%w * sqrt(self%rhos(j-1)/self%kappas(j-1))
             R = self%Rs(j)
             M2(1,1) =                 +besj(n,k*R)
             M2(2,1) = +self%betas(j-1)*besj_d(n,k*R)
             M2(1,2) =                 +besh(n,k*R)
             M2(2,2) = +self%betas(j-1)*besh_d(n,k*R)

             X2 = pi*self%w * sqrt(self%rhos(j)/self%kappas(j))*self%Rs(j)/(2*ione*self%betas(j)) * &
                  self%w*sqrt(self%rhos(j)/self%kappas(j)) *matmul(M1,M2)

             ! 第3項
             k = self%w * sqrt(self%rhos(j)/self%kappas(j))
             R = self%Rs(j)
             M1(1,1) = +self%betas(j)*besh_d(n,k*R)
             M1(2,1) = -self%betas(j)*besj_d(n,k*R)
             M1(1,2) =               -besh(n,k*R)
             M1(2,2) =               +besj(n,k*R)

             k = self%w * sqrt(self%rhos(j-1)/self%kappas(j-1))
             R = self%Rs(j)
             M2(1,1) =                 +besj_d(n,k*R)
             M2(2,1) = +self%betas(j-1)*besj_dd(n,k*R)
             M2(1,2) =                 +besh_d(n,k*R)
             M2(2,2) = +self%betas(j-1)*besh_dd(n,k*R)

             X3 = pi*self%w * sqrt(self%rhos(j)/self%kappas(j))*self%Rs(j)/(2*ione*self%betas(j)) * &
                  self%w*sqrt(self%rhos(j-1)/self%kappas(j-1)) *matmul(M1,M2)

             Td = matmul(X1+X2+X3,Td)
             
          else          
             k = self%w * sqrt(self%rhos(j)/self%kappas(j))
             R = self%Rs(j)
             M1(1,1) = +self%betas(j)*besh_d(n,k*R)
             M1(2,1) = -self%betas(j)*besj_d(n,k*R)
             M1(1,2) =               -besh(n,k*R)
             M1(2,2) =               +besj(n,k*R)

             k = self%w * sqrt(self%rhos(j-1)/self%kappas(j-1))
             R = self%Rs(j)
             M2(1,1) =                 +besj(n,k*R)
             M2(2,1) = +self%betas(j-1)*besj_d(n,k*R)
             M2(1,2) =                 +besh(n,k*R)
             M2(2,2) = +self%betas(j-1)*besh_d(n,k*R)

             Td = pi*self%w * sqrt(self%rhos(j)/self%kappas(j))*self%Rs(j)/(2*ione*self%betas(j)) * matmul(M1,matmul(M2,Td))

          end if
          
       end do

       ! self%S(n,n) = -T(2,1)/T(2,2)
       out(n,n) = - (Td(2,1)*T(2,2) - T(2,1)*Td(2,2))/T(2,2)**2
       
    end do

    
  end subroutine calc_derivative_radius

  ! Hの2階微分
  complex(8) function besh_dd(n, z)
    integer :: n
    real(8) :: z

    besh_dd = -besh_d(n,z)/z + (n**2-z**2)*besh(n,z)/z**2
    
  end function besh_dd

  ! Jの2階微分
  real(8) function besj_dd(n, z)
    integer :: n
    real(8) :: z

    besj_dd = -besj_d(n,z)/z + (n**2-z**2)*besj(n,z)/z**2
    
  end function besj_dd
  
end module smatrix_unit_donut_class
