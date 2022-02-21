module smatrix_shell_homogeneous_coupled_class
  use smatrix_shell_class
  use math
  use misc
  use bessel
  implicit none
  private

  type,public,extends(smatrix_shell) :: smatrix_shell_homogeneous_coupled
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

     !> Transfer行列の対角成分
     complex(8),allocatable :: T11(:)
     complex(8),allocatable :: T21(:)
     complex(8),allocatable :: T12(:)
     complex(8),allocatable :: T22(:)
   contains
     procedure :: del
     procedure :: new
     procedure :: calc_derivative

     procedure :: func_T11
     procedure :: func_T12
     procedure :: func_T41
     procedure :: func_T42
     procedure :: func_Tt11
     procedure :: func_Tt12
     procedure :: func_Tt41
     procedure :: func_Tt42
     procedure :: W11
     procedure :: W12
     procedure :: W21
     procedure :: W22
     procedure :: Wt11
     procedure :: Wt12
     procedure :: Wt21
     procedure :: Wt22
  end type smatrix_shell_homogeneous_coupled

contains

  elemental subroutine del(self)
    class(smatrix_shell_homogeneous_coupled),intent(inout) :: self

    if(allocated(self%S11)) deallocate(self%S11)
    if(allocated(self%S21)) deallocate(self%S21)
    if(allocated(self%S12)) deallocate(self%S12)
    if(allocated(self%S22)) deallocate(self%S22)
    
    if(allocated(self%T11)) deallocate(self%T11)
    if(allocated(self%T21)) deallocate(self%T21)
    if(allocated(self%T12)) deallocate(self%T12)
    if(allocated(self%T22)) deallocate(self%T22)
  end subroutine del

  subroutine new(self, w, Rout, Rin, rho, kappa, rho_e, dl, dm, m_min, m_max)
    class(smatrix_shell_homogeneous_coupled),intent(inout) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: Rout
    real(8),intent(in) :: Rin
    real(8),intent(in) :: rho
    real(8),intent(in) :: kappa
    real(8),intent(in) :: rho_e
    real(8),intent(in) :: dl
    real(8),intent(in) :: dm
    integer,intent(in) :: m_min
    integer,intent(in) :: m_max

    integer :: n

    complex(8) :: X11(4,4), X12(4,2)
    complex(8) :: X21(2,4), X22(2,2)
    complex(8) :: Y(4,2)

    complex(8) :: X(6,6), X_inv(6,6)

    complex(8) :: X11_inv(4,4), tmp(2,2), inv(2,2), T(2,2)

    call self%del

    self%w = w
    self%Rin = Rin
    self%Rout = Rout
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

    allocate(self%S11(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S11 = zero
    allocate(self%S21(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S21 = zero
    allocate(self%S12(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S12 = zero
    allocate(self%S22(self%m_min:self%m_max,self%m_min:self%m_max)) ; self%S22 = zero

    allocate(self%T11(self%m_min:self%m_max)) ; self%T11 = zero
    allocate(self%T12(self%m_min:self%m_max)) ; self%T12 = zero
    allocate(self%T21(self%m_min:self%m_max)) ; self%T21 = zero
    allocate(self%T22(self%m_min:self%m_max)) ; self%T22 = zero

    do n=self%m_min,self%m_max

       ! X11
       X11(1,:) = [self%func_T11(n,self%Rout), self%func_T12(n,self%Rout),&
        self%func_Tt11(n,self%Rout), self%func_Tt12(n,self%Rout)]
       X11(2,:) = [self%func_T41(n,self%Rout), self%func_T42(n,self%Rout), self%func_Tt41(n,self%Rout), self%func_Tt42(n,self%Rout)]
       X11(3,:) = [self%W11(n,self%Rout), self%W12(n,self%Rout), &
       self%Wt11(n,self%Rout), self%Wt12(n,self%Rout)]
       X11(4,:) = [self%func_T11(n,self%Rin), self%func_T12(n,self%Rin),&
        self%func_Tt11(n,self%Rin), self%func_Tt12(n,self%Rin)]

       ! X12
       X12(:,:) = zero
       X12(4,:) = [one*besj(n,self%k*self%Rin), besh(n,self%k*self%Rin)]

       ! X21
       X21(1,:) = [self%func_T41(n,self%Rin), self%func_T42(n,self%Rin), self%func_Tt41(n,self%Rin), self%func_Tt42(n,self%Rin)]
       X21(2,:) = [self%W11(n,self%Rin), self%W12(n,self%Rin), self%Wt11(n,self%Rin), self%Wt12(n,self%Rin)]

       ! X22
       X22(1,:) = [zero, zero]
       X22(2,:) = -[one*self%k/(self%rho*self%w**2)*besj_d(n,self%k*self%Rin),self%k/(self%rho*self%w**2)*besh_d(n,self%k*self%Rin)]

       ! Y
       Y(1,:) = [-one*besj(n,self%k*self%Rout), -besh(n,self%k*self%Rout)]
       Y(2,:) = [zero, zero]
       Y(3,:) = [one*self%k/(self%rho*self%w**2)*besj_d(n,self%k*self%Rout), self%k/(self%rho*self%w**2)*besh_d(n,self%k*self%Rout)]
       Y(4,:) = [zero, zero]


       ! X
       X(1:4,1:4) = X11
       X(5:6,1:4) = X21
       X(1:4,5:6) = X12
       X(5:6,5:6) = X22
       X_inv = mat_inv(6,X)

       T = matmul(X_inv(5:6,1:4),Y)

       self%T11(n) = T(1,1)
       self%T21(n) = T(2,1)
       self%T12(n) = T(1,2)
       self%T22(n) = T(2,2)

       
       
       ! ! write(*,*) n, X_inv
       
       ! ! write(*,*) n, X11(1,:)
       ! ! write(*,*) n, X11(2,:)
       ! ! write(*,*) n, X11(3,:)
       ! ! write(*,*) n, X11(4,:)

       ! ! X11のinverse
       ! X11_inv = mat_inv(4,X11)

       ! ! (X22-X21*X11^-1*X12)^1
       ! tmp = X22 - matmul(X21,matmul(X11_inv,X12))
       ! inv = mat_inv(2,tmp)

       ! ! T
       ! T = -matmul(inv,matmul(X21,matmul(X11_inv,Y)))

       ! S
       self%S11(n,n) = -T(2,1)/T(2,2)
       self%S12(n,n) = 1.d0 / T(2,2)
       self%S21(n,n) = T(1,1) - T(1,2)/T(2,2)*T(2,1)
       self%S22(n,n) = T(1,2)/T(2,2)
       
    end do

  end subroutine new

  subroutine calc_derivative(self, idv, S11_d, S21_d, S12_d, S22_d)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: idv
    complex(8),intent(out) :: S11_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S21_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S12_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S22_d(self%m_min:self%m_max,self%m_min:self%m_max)

    stop "todo"
  end subroutine calc_derivative

  complex(8) function func_T11(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_T11 = (n**2+n-0.5d0*(self%kT*r)**2) * besj(n,self%kL*r) - self%kL*r*besj(n-1,self%kL*r)

    func_T11 = func_T11 * 2*self%dm/r**2
  end function func_T11

  complex(8) function func_T12(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_T12 = -ione*n*( (n+1)*besj(n,self%kT*r) - self%kT*r*besj(n-1,self%kT*r)  )

    func_T12 = func_T12 * 2*self%dm/r**2
  end function func_T12

  complex(8) function func_T41(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_T41 = -ione*n*( (n+1)*besj(n,self%kL*r) - self%kL*r*besj(n-1,self%kL*r)  )
    
    func_T41 = func_T41 * 2*self%dm/r**2
  end function func_T41

  complex(8) function func_T42(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_T42 = -(n**2+n-0.5d0*(self%kT*r)**2) * besj(n,self%kT*r) + self%kT*r*besj(n-1,self%kT*r)

    func_T42 = func_T42 * 2*self%dm/r**2
  end function func_T42


  complex(8) function func_Tt11(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_Tt11 = (n**2+n-0.5d0*(self%kT*r)**2) * besh(n,self%kL*r) - self%kL*r*besh(n-1,self%kL*r)

    func_Tt11 = func_Tt11 * 2*self%dm/r**2
  end function func_Tt11

  complex(8) function func_Tt12(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_Tt12 = -ione*n*( (n+1)*besh(n,self%kT*r) - self%kT*r*besh(n-1,self%kT*r)  )

    func_Tt12 = func_Tt12 * 2*self%dm/r**2
  end function func_Tt12

  complex(8) function func_Tt41(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_Tt41 = -ione*n*( (n+1)*besh(n,self%kL*r) - self%kL*r*besh(n-1,self%kL*r)  )

    func_Tt41 = func_Tt41 * 2*self%dm/r**2
  end function func_Tt41

  complex(8) function func_Tt42(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    func_Tt42 = -(n**2+n-0.5d0*(self%kT*r)**2) * besh(n,self%kT*r) + self%kT*r*besh(n-1,self%kT*r)

    func_Tt42 = func_Tt42 * 2*self%dm/r**2
  end function func_Tt42

  complex(8) function W11(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W11 = self%kL * besj_d(n,self%kL*r)
  end function W11

  complex(8) function W12(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W12 = ione*n/r * besj(n,self%kT*r)
  end function W12

  complex(8) function W21(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W21 = ione*n/r * besj(n,self%kL*r)
  end function W21

  complex(8) function W22(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    W22 = -self%kT * besj_d(n,self%kT*r)
  end function W22


  complex(8) function Wt11(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt11 = self%kL * besh_d(n,self%kL*r)
  end function Wt11

  complex(8) function Wt12(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt12 = ione*n/r * besh(n,self%kT*r)
  end function Wt12

  complex(8) function Wt21(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r

    Wt21 = ione*n/r * besh(n,self%kL*r)
  end function Wt21

  complex(8) function Wt22(self, n, r)
    class(smatrix_shell_homogeneous_coupled),intent(in) :: self
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

end module smatrix_shell_homogeneous_coupled_class
