!> smatrix_coreと複数のsmatrix_shellで構成される
module smatrix_class
  use misc
  use math
  use smatrix_shell_class
  use smatrix_core_class
  implicit none
  private
  
  type,public :: smatrix
     !> 角周波数
     real(8) :: w
     !> shellの数
     integer :: nshell
     !> 各shell
     class(smatrix_shell),pointer :: shells(:) => null()
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
  end type smatrix

contains
  elemental subroutine del(self)
    class(smatrix),intent(inout) :: self

    if(associated(self%shells)) then
       ! call self%shells%del
       ! deallocate(self%shells)
       nullify(self%shells)
    end if    

    if(associated(self%core)) nullify(self%core)

    if(allocated(self%S)) deallocate(self%S)
  end subroutine del

  subroutine new(self, nshell, shells, core)
    class(smatrix),intent(inout) :: self
    integer,intent(in) :: nshell
    class(smatrix_shell),intent(in),target :: shells(nshell)
    class(smatrix_core),intent(in),target :: core

    integer :: ishell

    complex(8),allocatable :: S11(:,:)
    complex(8),allocatable :: S21(:,:)
    complex(8),allocatable :: S12(:,:)
    complex(8),allocatable :: S22(:,:)

    call self%del

    self%nshell = nshell
    ! allocate(self%shells(self%nshell), source=shells(1))
    self%shells => shells
    self%core => core

    self%m_min = self%core%m_min
    self%m_max = self%core%m_max

    ! shellとcoreのSのサイズは全部同じであると仮定する
    do ishell=1,self%nshell
       call assert(self%m_min == self%shells(ishell)%m_min)
       call assert(self%m_max == self%shells(ishell)%m_max)
    end do

    ! 第1shellから順に積み重ねていき，そのshell S行列を計算する
    allocate(S11(self%m_min:self%m_max,self%m_min:self%m_max)) ; S11 = self%shells(1)%S11
    allocate(S21(self%m_min:self%m_max,self%m_min:self%m_max)) ; S21 = self%shells(1)%S21
    allocate(S12(self%m_min:self%m_max,self%m_min:self%m_max)) ; S12 = self%shells(1)%S12    
    allocate(S22(self%m_min:self%m_max,self%m_min:self%m_max)) ; S22 = self%shells(1)%S22
    do ishell=2,self%nshell
       call stack_shell(self%m_min,self%m_max, &
            S11, S21, S12, S22, &
            self%shells(ishell)%S11, self%shells(ishell)%S21, self%shells(ishell)%S12, self%shells(ishell)%S22)
    end do

    ! write(*,*) S12

    ! coreをつなげて，self%Sにする
    allocate(self%S(self%m_min:self%m_max,self%m_min:self%m_max))
    self%S = self%core%S
    call stack_core(self%m_min,self%m_max,S11, S21, S12, S22, self%S)

    ! ! tmp
    ! self%S = S11
    
  end subroutine new

  subroutine stack_shell(m_min,m_max,&
       S_btm_11, S_btm_21, S_btm_12, S_btm_22, S_top_11, S_top_21, S_top_12, S_top_22)
    integer,intent(in) :: m_min, m_max
    complex(8),intent(inout)  :: S_btm_11(m_min:m_max,m_min:m_max),S_btm_21(m_min:m_max,m_min:m_max),&
         S_btm_12(m_min:m_max,m_min:m_max),S_btm_22(m_min:m_max,m_min:m_max)
    complex(8),intent(in)  :: S_top_11(m_min:m_max,m_min:m_max),S_top_21(m_min:m_max,m_min:m_max),&
         S_top_12(m_min:m_max,m_min:m_max),S_top_22(m_min:m_max,m_min:m_max)
    ! complex(8),intent(out) :: S_all_11(m_min:m_max,m_min:m_max),S_all_21(m_min:m_max,m_min:m_max),&
    !      S_all_12(m_min:m_max,m_min:m_max),S_all_22(m_min:m_max,m_min:m_max)

    complex(8),allocatable :: S11(:,:), S21(:,:), S12(:,:), S22(:,:)

    complex(8),allocatable :: ident(:,:), inv1(:,:), inv2(:,:)
    integer :: i
    
    ! identity matrix
    allocate(ident(m_min:m_max,m_min:m_max))
    ident(:,:) = zero
    do i=m_min,m_max
       ident(i,i) = one
    end do

    ! inv1 = (I-S^(1)_22*S^(2)_11)^-1
    allocate(inv1(m_min:m_max,m_min:m_max))
    inv1 = mat_inv(m_max-m_min+1, ident - matmul(S_btm_22,S_top_11))
    ! inv2 = (I-S^(2)_11*S^(1)_22)^-1
    allocate(inv2(m_min:m_max,m_min:m_max))
    inv2 = mat_inv(m_max-m_min+1, ident - matmul(S_top_11,S_btm_22))

    allocate(S11(m_min:m_max,m_min:m_max))
    allocate(S21(m_min:m_max,m_min:m_max))
    allocate(S12(m_min:m_max,m_min:m_max))
    allocate(S22(m_min:m_max,m_min:m_max))
    
    ! S11 = S^(1)_12 * inv2 * S^(2)_11 * S^(1)_21 + S^(1)_11
    S11 = matmul(S_btm_12,matmul(inv2,matmul(S_top_11,S_btm_21))) + S_btm_11

    ! S21 = S^(2)_21 * inv1 * S^(1)_21
    S21 = matmul(S_top_21,matmul(inv1,S_btm_21))

    ! S12 = S^(1)_12 * inv2 * S^(2)_12
    S12 = matmul(S_btm_12,matmul(inv2,S_top_12))

    ! S22 = S^(2)_21 * inv1 * S^(1)_22 * S^(2)_12 + S^(2)_22
    S22 = matmul(S_top_21,matmul(inv1,matmul(S_btm_22,S_top_12))) + S_top_22

    ! 最後に代入
    S_btm_11 = S11
    S_btm_21 = S21
    S_btm_12 = S12
    S_btm_22 = S22

    deallocate(ident,inv1,inv2,S11,S21,S12,S22)
    
  end subroutine stack_shell

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
  
end module smatrix_class
