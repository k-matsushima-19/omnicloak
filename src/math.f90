!> 定数とよく使う静的な関数
module math
  use misc
  implicit none
  
  complex(8),parameter :: zero = (0.d0,0.d0)
  complex(8),parameter :: ione = (0.d0,1.d0)
  complex(8),parameter ::  one = (1.d0,0.d0)
  real(8),parameter :: pi = 3.14159265358979323846d0

  interface rdot_product
     module procedure rdot_product1
     module procedure rdot_product2
     module procedure rdot_product3
     module procedure rdot_product4
  end interface rdot_product

  interface rotate_vec
     module procedure rotate_vec_real
     module procedure rotate_vec_complex
  end interface rotate_vec

contains
  ! 共役を取らないドット積
  complex(8) function rdot_product1(a, b)
    complex(8),intent(in) :: a(:)
    real(8),intent(in) :: b(:)

    rdot_product1 = dot_product(b, a)
  end function rdot_product1

  complex(8) function rdot_product2(a, b)
    real(8),intent(in) :: a(:)
    complex(8),intent(in) :: b(:)

    rdot_product2 = dot_product(a, b)
  end function rdot_product2

  complex(8) function rdot_product3(a, b)
    complex(8),intent(in) :: a(:)
    complex(8),intent(in) :: b(:)

    rdot_product3 = dot_product(conjg(a), b)
  end function rdot_product3

  complex(8) function rdot_product4(a, b)
    real(8),intent(in) :: a(:)
    real(8),intent(in) :: b(:)

    rdot_product4 = dot_product(a, b)
  end function rdot_product4

  !> \f$ ||x|| = \sqrt{x\cdot x} \f$
  real(8) function length_vector(x)
    real(8),intent(in) :: x(:)

    length_vector = sqrt(dot_product(x,x))
  end function length_vector

  !> 2次元のベクトルを回転
  function rotate_vec_real(x, theta)
    real(8),intent(in) :: x(2)
    real(8),intent(in) :: theta
    real(8) :: rotate_vec_real(2)

    real(8) :: mat(2,2)

    mat(:,1) = [ cos(theta), sin(theta)]
    mat(:,2) = [-sin(theta), cos(theta)]

    rotate_vec_real = matmul(mat,x)
    
  end function rotate_vec_real

  !> 2次元のベクトルを回転
  function rotate_vec_complex(x, theta)
    complex(8),intent(in) :: x(2)
    real(8),intent(in) :: theta
    complex(8) :: rotate_vec_complex(2)

    real(8) :: mat(2,2)

    mat(:,1) = [ cos(theta), sin(theta)]
    mat(:,2) = [-sin(theta), cos(theta)]

    rotate_vec_complex = matmul(mat,x)
    
  end function rotate_vec_complex

  !> Compute SVD
  ! A = U*S*V^H, where U: m x m, S: m x n, V^H: n x n
  ! S is diagonal and S(i,i)=0 if i>min(m,n)
  subroutine compute_SVD(m, n, mat, U, S, VH)
    integer,intent(in) :: m, n ! size of mat
    complex(8),intent(in) :: mat(m,n)
    complex(8),intent(out) :: U(m,m)
    real(8),intent(out) :: S(min(m,n))
    complex(8),intent(out) :: VH(n,n)
    

    complex(8),allocatable :: A(:,:)
    integer :: i, info
    complex(8),allocatable :: work(:)
    integer :: lwork
    real(8),allocatable :: rwork(:)
    
    ! copy
    allocate(A(n,n))
    A(:,:) = mat(:,:)

    ! SVD
    lwork = 2*MIN(M,N)+MAX(M,N)
    allocate(work(lwork))
    allocate(rwork(5*min(m,n)))
    call zgesvd(&
         "A", & ! JOBU
         "A", & ! JOBVT
         m, & ! M
         n, & ! N
         A, & ! A
         m, & ! LDA
         S, & ! S
         U, & ! U
         m, & ! LDU
         VH, & ! VT
         n, & ! LDVT
         work, & ! WORK
         lwork, & ! lwork
         rwork, & ! rwork
         info & ! info
         )
    call assert(info == 0)

    
  end subroutine compute_SVD
end module math
