!> smatrix_coreと"1つ"のsmatrix_shellで構成される
module smatrix_shellcore_class
  use misc
  use math
  use smatrix_shell_class
  use smatrix_core_class
  implicit none
  private
  
  type,public :: smatrix_shellcore
     !> 角周波数
     real(8) :: w
     !> shell
     class(smatrix_shell),pointer :: shell => null()
     !> core
     class(smatrix_core),pointer :: core => null()
     !> 各shellのSのサイズ
     integer :: m_min
     integer :: m_max

     !> core + shellによる全体のS行列
     complex(8),allocatable :: S(:,:)
   contains
     procedure :: del
     procedure :: new
     procedure :: calc_derivative
  end type smatrix_shellcore

contains
  elemental subroutine del(self)
    class(smatrix_shellcore),intent(inout) :: self

    if(associated(self%shell)) then
       ! call self%shells%del
       ! deallocate(self%shells)
       nullify(self%shell)
    end if    

    if(associated(self%core)) nullify(self%core)

    if(allocated(self%S)) deallocate(self%S)
  end subroutine del

  subroutine new(self, shell, core)
    class(smatrix_shellcore),intent(inout) :: self
    class(smatrix_shell),intent(in),target :: shell
    class(smatrix_core),intent(in),target :: core

    integer :: ishell

    complex(8),allocatable :: S11(:,:)
    complex(8),allocatable :: S21(:,:)
    complex(8),allocatable :: S12(:,:)
    complex(8),allocatable :: S22(:,:)

    call self%del
    
    ! allocate(self%shells(self%nshell), source=shells(1))
    self%shell => shell
    self%core => core

    self%m_min = self%core%m_min
    self%m_max = self%core%m_max

    ! shellとcoreのSのサイズは全部同じであると仮定する
    call assert(self%m_min == self%shell%m_min)
    call assert(self%m_max == self%shell%m_max)

    ! 第1shellから順に積み重ねていき，そのshell S行列を計算する
    allocate(S11(self%m_min:self%m_max,self%m_min:self%m_max)) ; S11 = self%shell%S11
    allocate(S21(self%m_min:self%m_max,self%m_min:self%m_max)) ; S21 = self%shell%S21
    allocate(S12(self%m_min:self%m_max,self%m_min:self%m_max)) ; S12 = self%shell%S12    
    allocate(S22(self%m_min:self%m_max,self%m_min:self%m_max)) ; S22 = self%shell%S22
    ! do ishell=2,self%nshell
    !    call stack_shell(self%m_min,self%m_max, &
    !         S11, S21, S12, S22, &
    !         self%shells(ishell)%S11, self%shells(ishell)%S21, self%shells(ishell)%S12, self%shells(ishell)%S22)
    ! end do

    ! write(*,*) S12

    ! coreをつなげて，self%Sにする
    allocate(self%S(self%m_min:self%m_max,self%m_min:self%m_max))
    self%S = self%core%S
    call stack_core(self%m_min,self%m_max,S11, S21, S12, S22, self%S)

    ! ! tmp
    ! self%S = S11
    
  end subroutine new

  subroutine stack_core(m_min,m_max,S11, S21, S12, S22, Score)
    integer,intent(in) :: m_min, m_max    
    complex(8),intent(in) :: S11(m_min:m_max,m_min:m_max),S21(m_min:m_max,m_min:m_max),&
         S12(m_min:m_max,m_min:m_max),S22(m_min:m_max,m_min:m_max)
    complex(8),intent(inout) :: Score(m_min:m_max,m_min:m_max)

    complex(8),allocatable :: ident(:,:), inv(:,:), S(:,:)
    integer :: i
    
    ! identity matrix
    allocate(ident(m_min:m_max,m_min:m_max))
    ident(:,:) = zero
    do i=m_min,m_max
       ident(i,i) = one
    end do

    ! inv = I-S22*S_core
    allocate(inv(m_min:m_max,m_min:m_max))
    inv = mat_inv(m_max-m_min+1, ident - matmul(S22,Score))

    ! write(*,*) matmul(S12,matmul(inv,S21))
    allocate(S(m_min:m_max,m_min:m_max))
    S = S11 + matmul(S12,matmul(Score,matmul(inv,S21)))

    Score = S

    deallocate(ident,inv, S)
    
  end subroutine stack_core

  ! self%shellのidv番設計変数に関するself%Sの微分を計算
  subroutine calc_derivative(self, idv, Sd)
    class(smatrix_shellcore),intent(in) :: self
    integer,intent(in) :: idv
    complex(8),intent(out) :: Sd(self%m_min:self%m_max,self%m_min:self%m_max)
    
    complex(8),allocatable :: ident(:,:), inv(:,:)
    integer :: i, info
    integer,allocatable :: ipiv(:)

    complex(8),allocatable :: S11_d(:,:), S21_d(:,:), S12_d(:,:), S22_d(:,:)
    complex(8),allocatable :: tmp(:,:)

    ! self%shellの微分
    allocate(S11_d(self%shell%m_min:self%shell%m_max,self%shell%m_min:self%shell%m_max))
    allocate(S21_d(self%shell%m_min:self%shell%m_max,self%shell%m_min:self%shell%m_max))
    allocate(S12_d(self%shell%m_min:self%shell%m_max,self%shell%m_min:self%shell%m_max))
    allocate(S22_d(self%shell%m_min:self%shell%m_max,self%shell%m_min:self%shell%m_max))
    call self%shell%calc_derivative(idv, S11_d, S21_d, S12_d, S22_d)

    ! identity matrix
    allocate(ident(self%m_min:self%m_max,self%m_min:self%m_max))
    ident(:,:) = zero
    do i=self%m_min,self%m_max
       ident(i,i) = one
    end do

    ! (I-S22*Score)のLU分解
    allocate(inv(self%m_min:self%m_max,self%m_min:self%m_max))
    allocate(ipiv(self%m_min:self%m_max))
    inv = ident - matmul(self%shell%S22,self%core%S)
    call zgetrf(self%m_max-self%m_min+1,self%m_max-self%m_min+1, inv, self%m_max-self%m_min+1, ipiv, info)
    call assert(info == 0)

    allocate(tmp(self%m_min:self%m_max,self%m_min:self%m_max))

    !
    ! S11'
    !
    Sd = S11_d

    !
    ! S12' * Score * (I-S22*Score)^-1 * S21
    !
    tmp = self%shell%S21
    call zgetrs("N", self%m_max-self%m_min+1, self%m_max-self%m_min+1, inv, self%m_max-self%m_min+1, ipiv, tmp, &
         self%m_max-self%m_min+1, info)
    call assert(info == 0)
    tmp = matmul(self%core%S,tmp)
    tmp = matmul(S12_d,tmp)
    
    Sd = Sd + tmp

    !
    ! S22' * Score * (I-S22*Score)^-1 * S21 + S21'
    !
    tmp = self%shell%S21
    call zgetrs("N", self%m_max-self%m_min+1, self%m_max-self%m_min+1, inv, self%m_max-self%m_min+1, ipiv, tmp, &
         self%m_max-self%m_min+1, info)
    call assert(info == 0)
    tmp = matmul(self%core%S,tmp)
    tmp = matmul(S22_d,tmp)
    tmp = tmp + S21_d

    !
    ! S12 * Score * (I-S22*Score)^-1 * tmp
    !
    call zgetrs("N", self%m_max-self%m_min+1, self%m_max-self%m_min+1, inv, self%m_max-self%m_min+1, ipiv, tmp, &
         self%m_max-self%m_min+1, info)
    call assert(info == 0)
    tmp = matmul(self%core%S,tmp)
    tmp = matmul(self%shell%S12,tmp)

    Sd = Sd + tmp
    
  end subroutine calc_derivative

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
    if(info /= 0) stop "info /= 0"
    call zgetrs('N', n, n, A_copy, n, ipiv, out, n, info)
    if(info /= 0) stop "info /= 0"
    
  end function mat_inv
  
end module smatrix_shellcore_class
