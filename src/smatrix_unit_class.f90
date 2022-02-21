!> 単位構造のS行列
module smatrix_unit_class
  use math
  use misc
  use bessel
  implicit none
  private

  type,public,abstract :: smatrix_unit
     !> 周波数
     real(8) :: w
     !> 最外半径 (|x|>Rは外部領域)
     real(8) :: R
     !> Sのサイズ
     integer :: m_min
     !> Sのサイズ
     integer :: m_max
     !> m_max-m_min+1
     integer :: size

     !> S行列
     complex(8),allocatable :: S(:,:)

     !> 設計変数の数
     integer :: ndv
   contains
     procedure(def_del),deferred :: del
     procedure(def_calc_derivative),deferred :: calc_derivative
  end type smatrix_unit

  interface
     elemental subroutine def_del(self)
       import smatrix_unit
       class(smatrix_unit),intent(inout) :: self
     end subroutine def_del

     subroutine def_calc_derivative(self, idv, out)
       import smatrix_unit
       class(smatrix_unit),intent(in) :: self
       integer,intent(in) :: idv
       complex(8),intent(out) :: out(self%m_min:self%m_max,self%m_min:self%m_max)
     end subroutine def_calc_derivative
  end interface
  
  
  
end module smatrix_unit_class
