module smatrix_core_class
  use math
  use bessel
  implicit none
  private

  type,public,abstract :: smatrix_core
     real(8) :: w
     ! ! 外部のパラメータ
     ! complex(8) :: rho_ext
     ! complex(8) :: kappa_ext
     ! complex(8) :: k_ext
     ! ! 内部のパラメータ
     ! complex(8) :: rho_int
     ! complex(8) :: kappa_int
     ! complex(8) :: k_int
     ! 半径
     real(8) :: R
     ! S行列のサイズ
     integer :: m_min
     integer :: m_max
     ! S行列
     complex(8),allocatable :: S(:,:)
   contains
     ! procedure :: del
     ! procedure :: new
  end type smatrix_core

  

contains
  ! elemental subroutine del(self)
  !   class(smatrix_core),intent(inout) :: self

  !   if(allocated(self%S)) deallocate(self%S)

    
    
  ! end subroutine del

  ! subroutine new(self, w, rho_int, kappa_int, R, m_min, m_max)
  !   class(smatrix_core),intent(inout) :: self
  !   real(8),intent(in) :: w
  !   complex(8),intent(in) :: rho_int
  !   complex(8),intent(in) :: kappa_int
  !   real(8),intent(in) :: R
  !   integer,intent(in) :: m_min
  !   integer,intent(in) :: m_max

  !   complex(8),allocatable :: H_int(:), H_ext(:)
  !   complex(8),allocatable :: J_int(:), J_ext(:)
  !   integer :: n

  !   call self%del

  !   self%w = w
  !   self%rho_int = rho_int
  !   self%kappa_int = kappa_int
  !   self%rho_ext = one
  !   self%kappa_ext = one
  !   self%R = R
  !   self%m_min = m_min
  !   self%m_max = m_max
    
  !   self%k_int = w * sqrt(self%rho_int / self%kappa_int)
  !   self%k_ext = w * sqrt(self%rho_ext / self%kappa_ext)

  !   ! write(*,*) self%k_int, self%k_ext

  !   ! 波数の虚部は必ず正
  !   if(aimag(self%k_int) < 0.d0) self%k_int = conjg(self%k_int)
  !   if(aimag(self%k_ext) < 0.d0) self%k_ext = conjg(self%k_ext)

  !   !
  !   ! Sの計算
  !   !
  !   allocate(self%S(self%m_min:self%m_max,self%m_min:self%m_max))
  !   self%S(:,:) = zero

  !   ! bessel関数の計算
  !   ! n-1 ~ n次
  !   allocate(H_int(self%m_min-1:self%m_max),H_ext(self%m_min-1:self%m_max))
  !   call cbeshn(self%k_int*self%R, self%m_min-1, self%m_max, H_int)
  !   call cbeshn(self%k_ext*self%R, self%m_min-1, self%m_max, H_ext)
  !   allocate(J_int(self%m_min-1:self%m_max),J_ext(self%m_min-1:self%m_max))
  !   call cbesjn(self%k_int*self%R, self%m_min-1, self%m_max, J_int)
  !   call cbesjn(self%k_ext*self%R, self%m_min-1, self%m_max, J_ext)
    
  !   do n=self%m_min,self%m_max
  !      self%S(n,n) = -(self%rho_int*J_int(n)*(self%k_ext*J_ext(n-1) - n/self%R*J_ext(n)) &
  !           - self%rho_ext*J_ext(n)*(self%k_int*J_int(n-1) - n/self%R*J_int(n))) &
  !           / &
  !           (self%rho_int*J_int(n)*(self%k_ext*H_ext(n-1) - n/self%R*H_ext(n)) &
  !           - self%rho_ext*H_ext(n)*(self%k_int*J_int(n-1) - n/self%R*J_int(n)))

  !   end do
    
  ! end subroutine new
  
end module smatrix_core_class
