module nlopt_driver_homogenize_class
  use nlopt_driver_class
  use math
  use bessel
  use smatrix_shell_multiple_class
  use smatrix_unit_donut_class
  use smatrix_shell_homogeneous_class
  implicit none
  private

  include 'nlopt.f'

  type,public,extends(nlopt_driver) :: nlopt_driver_homogenize
     ! 散乱体の数
     integer :: nobj
     ! 各散乱体の座標
     real(8),allocatable :: centers(:,:)
     ! 比密度
     real(8) :: rhos(2)
     ! 比体積弾性率
     real(8) :: kappas(2)
     ! 角周波数
     real(8) :: w
     ! Sのサイズ
     integer :: m_min
     ! Sのサイズ
     integer :: m_max
     ! 目標とする均質なshell S行列
     type(smatrix_shell_homogeneous),pointer :: S_obj => null()
   contains
     procedure :: del
     procedure :: new
     procedure :: objective
     procedure :: objective_grad
     procedure :: constraint
     procedure :: constraint_grad
  end type nlopt_driver_homogenize


contains
  elemental subroutine del(self)
    class(nlopt_driver_homogenize),intent(inout) :: self

    if(allocated(self%lb)) deallocate(self%lb)
    if(allocated(self%ub)) deallocate(self%ub)
    if(allocated(self%x0)) deallocate(self%x0)
  end subroutine del

  subroutine new(self, w, rhos, kappas, m_min, m_max, Rmax, distfile, S_obj)
    class(nlopt_driver_homogenize),intent(inout) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: rhos(2)
    real(8),intent(in) :: kappas(2)
    integer,intent(in) :: m_min
    integer,intent(in) :: m_max
    real(8),intent(in) :: Rmax
    character(*),intent(in) :: distfile
    type(smatrix_shell_homogeneous),intent(in),target :: S_obj

    integer :: nobj, i, iobj
    real(8) :: angle

    call self%del

    self%w = w
    self%rhos = rhos
    self%kappas = kappas
    self%m_min = m_min
    self%m_max = m_max

    self%S_obj => S_obj

    !
    ! 配置を決定する
    !
    ! dist.txtから配置を読み込む
    open(10,file=distfile)
    read(10,*) self%nobj ; allocate(self%centers(2,self%nobj))
    do iobj=1,self%nobj
       read(10,*) i, self%centers(:,iobj), angle
    end do
    close(10)

    ! 設計変数は散乱体(ドーナツ)の外径と内径の2つ
    ! x(1): 外径，x(2): 内径
    self%n = 2

    allocate(self%lb(self%n))

    allocate(self%ub(self%n))
    ! 0 <= x(1) <= Rmax
    self%lb(1) = Rmax*1d-3
    self%ub(1) = Rmax
    ! 0 <= x(2) <= Rmax
    self%lb(2) = Rmax*1d-3
    self%ub(2) = Rmax
    
    allocate(self%x0(self%n))
    self%x0 = [Rmax, Rmax/2]
    ! self%x0 = [4.768145480444083E-002,  1.433872940251256E-002]

    self%n_constraint = 1

    ! self%ALGORITHM = NLOPT_GN_AGS
    self%ALGORITHM = NLOPT_LD_MMA
    
    
  end subroutine new

  ! 目的関数
  real(8) function objective(self, x)
    class(nlopt_driver_homogenize),intent(in) :: self
    real(8),intent(in) :: x(self%n)

    type(smatrix_unit_donut) :: smat_unit(1)
    type(smatrix_shell_multiple) :: smat_shell

    integer,allocatable :: dict(:)
    integer :: m

    ! 単位構造のS行列 (silicon - lead)
    call smat_unit(1)%new(self%w, 2, [x(1),x(2)], self%rhos, self%kappas)
    ! shell
    allocate(dict(self%nobj)) ; dict(:) = 1 ! 散乱体は1種類しかない
    call smat_shell%new(self%w, self%m_min, self%m_max, self%nobj, self%centers, 1, dict, smat_unit)

    ! 目的関数 = Sの対角成分の差の絶対値の2乗の総和
    objective = 0.d0
    do m=self%m_min,self%m_max
       objective = objective + abs(smat_shell%S11(m,m)-self%S_obj%S11(m,m))**2
       ! objective = objective + abs(smat_shell%S21(m,m)-self%S_obj%S21(m,m))**2
       ! objective = objective + abs(smat_shell%S12(m,m)-self%S_obj%S12(m,m))**2
       ! objective = objective + abs(smat_shell%S22(m,m)-self%S_obj%S22(m,m))**2

       ! write(*,*) x
       ! write(*,*) "11", abs(smat_shell%S11(m,m)-self%S_obj%S11(m,m))**2
       ! write(*,*) "21", abs(smat_shell%S21(m,m)-self%S_obj%S21(m,m))**2
       ! write(*,*) "12", abs(smat_shell%S12(m,m)-self%S_obj%S12(m,m))**2
       ! write(*,*) "22", abs(smat_shell%S22(m,m)-self%S_obj%S22(m,m))**2
    end do

    write(*,*) x, objective

    call smat_shell%del
    call smat_unit(1)%del

  end function objective

  ! 目的関数の勾配
  function objective_grad(self, x) result(out)
    class(nlopt_driver_homogenize),intent(in) :: self
    real(8),intent(in) :: x(self%n)
    real(8) :: out(self%n)

    type(smatrix_unit_donut) :: smat_unit(1)
    type(smatrix_shell_multiple) :: smat_shell

    integer,allocatable :: dict(:)
    integer :: m

    complex(8),allocatable :: S11_d(:,:)
    complex(8),allocatable :: S21_d(:,:)
    complex(8),allocatable :: S12_d(:,:)
    complex(8),allocatable :: S22_d(:,:)

    ! 単位構造のS行列 (silicon - lead)
    call smat_unit(1)%new(self%w, 2, [x(1),x(2)], self%rhos, self%kappas)
    ! shell
    allocate(dict(self%nobj)) ; dict(:) = 1 ! 散乱体は1種類しかない
    call smat_shell%new(self%w, self%m_min, self%m_max, self%nobj, self%centers, 1, dict, smat_unit)

    allocate(S11_d(smat_shell%m_min:smat_shell%m_max,smat_shell%m_min:smat_shell%m_max))
    allocate(S21_d(smat_shell%m_min:smat_shell%m_max,smat_shell%m_min:smat_shell%m_max))
    allocate(S12_d(smat_shell%m_min:smat_shell%m_max,smat_shell%m_min:smat_shell%m_max))
    allocate(S22_d(smat_shell%m_min:smat_shell%m_max,smat_shell%m_min:smat_shell%m_max))
     
    
    !
    ! x(1)に関する微分
    !
    call smat_shell%calc_derivative(1, S11_d, S21_d, S12_d, S22_d)
    out(1) = 0.d0
    do m=self%m_min,self%m_max
       out(1) = out(1) + 2*real(conjg(smat_shell%S11(m,m)-self%S_obj%S11(m,m))*S11_d(m,m))
       ! out(1) = out(1) + 2*real(conjg(smat_shell%S21(m,m)-self%S_obj%S21(m,m))*S21_d(m,m))
       ! out(1) = out(1) + 2*real(conjg(smat_shell%S12(m,m)-self%S_obj%S12(m,m))*S12_d(m,m))
       ! out(1) = out(1) + 2*real(conjg(smat_shell%S22(m,m)-self%S_obj%S22(m,m))*S22_d(m,m))
    end do

    !
    ! x(2)に関する微分
    !
    call smat_shell%calc_derivative(2, S11_d, S21_d, S12_d, S22_d)
    out(2) = 0.d0
    do m=self%m_min,self%m_max
       out(2) = out(2) + 2*real(conjg(smat_shell%S11(m,m)-self%S_obj%S11(m,m))*S11_d(m,m))
       ! out(2) = out(2) + 2*real(conjg(smat_shell%S21(m,m)-self%S_obj%S21(m,m))*S21_d(m,m))
       ! out(2) = out(2) + 2*real(conjg(smat_shell%S12(m,m)-self%S_obj%S12(m,m))*S12_d(m,m))
       ! out(2) = out(2) + 2*real(conjg(smat_shell%S22(m,m)-self%S_obj%S22(m,m))*S22_d(m,m))
    end do

    call smat_shell%del
    call smat_unit(1)%del
    
  end function objective_grad

  ! ic番目の制約関数の定義
  real(8) function constraint(self, ic, x)
    class(nlopt_driver_homogenize),intent(in) :: self
    integer,intent(in) :: ic
    real(8),intent(in) :: x(self%n)

    ! ic番目の制約 g(x) <= 0
    select case(ic)
    case(1)
       ! x(1) > x(2) => x(2)-x(1) < 0.0
       constraint = x(2) - x(1)
    end select
  end function constraint

  function constraint_grad(self, ic, x) result(out)
    class(nlopt_driver_homogenize),intent(in) :: self
    integer,intent(in) :: ic
    real(8),intent(in) :: x(self%n)
    real(8) :: out(self%n)

    ! ic番目の制約 g(x) <= 0
    select case(ic)
    case(1)
       ! x(1) > x(2) => x(2)-x(1) < 0.0
       out = [-1.d0, +1.d0]
    end select
  end function constraint_grad
  
end module nlopt_driver_homogenize_class
