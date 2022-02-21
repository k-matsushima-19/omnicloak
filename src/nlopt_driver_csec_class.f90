module nlopt_driver_csec_class
  use nlopt_driver_class
  use math
  use bessel
  use smatrix_shell_multiple_class
  use smatrix_unit_donut_coupled_class
  use smatrix_core_class
  use smatrix_shellcore_class
  implicit none
  private

  include 'nlopt.f'

  type,public,extends(nlopt_driver) :: nlopt_driver_csec
     ! 散乱体の数
     integer :: nobj
     ! 散乱体の種類
     integer :: nkind
     !> 散乱体の番号と種類を対応させる配列．iobj番散乱体はdict(iobj)種
     integer,allocatable :: dict(:)
     ! 各散乱体の座標
     real(8),allocatable :: centers(:,:)
     ! 比密度
     real(8) :: rhos(2)
     ! 比Lame定数
     real(8) :: dls(2)
     ! 比Lame定数
     real(8) :: dms(2)
     ! 角周波数
     real(8) :: w
     ! Sのサイズ
     integer :: m_min
     ! Sのサイズ
     integer :: m_max
     ! core
     class(smatrix_core),pointer :: S_core => null()
   contains
     procedure :: del
     procedure :: new
     procedure :: objective
     procedure :: objective_grad
     procedure :: constraint
     procedure :: constraint_grad
  end type nlopt_driver_csec


contains
  elemental subroutine del(self)
    class(nlopt_driver_csec),intent(inout) :: self

    if(allocated(self%lb)) deallocate(self%lb)
    if(allocated(self%ub)) deallocate(self%ub)
    if(allocated(self%x0)) deallocate(self%x0)
  end subroutine del

  subroutine new(self, w, rhos, dls, dms, m_min, m_max, Rmax, distfile, S_core)
    class(nlopt_driver_csec),intent(inout) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: rhos(2)
    real(8),intent(in) :: dls(2)
    real(8),intent(in) :: dms(2)
    integer,intent(in) :: m_min
    integer,intent(in) :: m_max
    real(8),intent(in) :: Rmax
    character(*),intent(in) :: distfile
    class(smatrix_core),intent(in),target :: S_core

    integer :: nobj, i, iobj
    real(8) :: angle
    integer :: ikind

    call self%del

    self%w = w
    self%rhos = rhos
    self%dls = dls
    self%dms = dms
    self%m_min = m_min
    self%m_max = m_max

    self%S_core => S_core

    !
    ! 配置を決定する
    !
    ! dist.txtから配置を読み込む
    open(10,file=distfile, status="old")
    read(10,*) self%nobj, self%nkind
    allocate(self%centers(2,self%nobj))
    allocate(self%dict(self%nobj))
    do iobj=1,self%nobj
       read(10,*) self%dict(iobj), self%centers(:,iobj), angle
    end do
    close(10)

    ! 設計変数は各散乱体(ドーナツ)の外径と内径の2つ
    ! x(1): 外径，x(2): 内径, ...
    self%n = self%nkind*2

    allocate(self%lb(self%n))

    allocate(self%ub(self%n))
    do ikind=1,self%nkind
       ! 外径の範囲
       self%lb(2*ikind-1) = Rmax*1d-5
       self%ub(2*ikind-1) = Rmax

       ! 内径の範囲
       self%lb(2*ikind  ) = Rmax*1d-5
       self%ub(2*ikind  ) = Rmax
    end do
        
    !    ! 0 <= x(2) <= Rmax
    ! self%lb(2) = Rmax*1d-3
    ! self%ub(2) = Rmax
    
    allocate(self%x0(self%n))

    !
    ! initial value
    !
    do ikind=1,self%nkind
       self%x0(2*ikind-1) = 0.18804571611927015
       self%x0(2*ikind  ) = 0.027306730722850392
    end do

    !   5.322431182861328E-002  7.127876281738277E-003   2.33118445462747     
    ! self%x0 = [Rmax, Rmax/2]
    ! self%x0 = [5.322431182861328E-002,  7.127876281738277E-003]
    ! self%x0 = [4.768145480444083E-002,  1.433872940251256E-002]

    self%n_constraint = self%nkind

    ! self%ALGORITHM = NLOPT_GN_AGS
    
    self%ALGORITHM = NLOPT_LD_SLSQP
    ! self%ALGORITHM = NLOPT_GN_AGS
    ! self%ALGORITHM = NLOPT_LD_MMA


    ! write(*,*) self%x0
    ! write(*,*) self%n_constraint
    
    
  end subroutine new

  ! 目的関数
  real(8) function objective(self, x)
    class(nlopt_driver_csec),intent(in) :: self
    real(8),intent(in) :: x(self%n)

    type(smatrix_unit_donut_coupled),allocatable :: smat_unit(:)
    type(smatrix_shell_multiple) :: smat_shell
    type(smatrix_shellcore) :: smat_shellcore

    ! integer,allocatable :: dict(:)
    integer :: m, ikind

    ! 単位構造のS行列をnkind個作る
    allocate(smat_unit(self%nkind))
    do ikind=1,self%nkind
       call smat_unit(ikind)%new(self%w, 1.d0, 1.d0, 2, [x(ikind*2-1),x(ikind*2)], self%rhos, self%dls, self%dms)
    end do
    ! shell
    ! allocate(dict(self%nobj)) ; dict(:) = 1 ! 散乱体は1種類しかない
    call smat_shell%new(self%w, self%m_min, self%m_max, self%nobj, self%centers, self%nkind, self%dict, smat_unit)
    call smat_shellcore%new(smat_shell, self%S_core)

    ! 目的関数 = Sの対角成分の差の絶対値の2乗の総和 ≒ 平面派入射の散乱断面積
    objective = 0.d0
    do m=self%m_min,self%m_max
       objective = objective + abs(smat_shellcore%S(m,m))**2
    end do

    write(*,*) objective
    
    call smat_shell%del
    call smat_unit%del
    call smat_shellcore%del

  end function objective

  ! 目的関数の勾配
  function objective_grad(self, x) result(out)
    class(nlopt_driver_csec),intent(in) :: self
    real(8),intent(in) :: x(self%n)
    real(8) :: out(self%n)

    type(smatrix_unit_donut_coupled),allocatable :: smat_unit(:)
    type(smatrix_shell_multiple) :: smat_shell
    type(smatrix_shellcore) :: smat_shellcore

    ! integer,allocatable :: dict(:)
    integer :: m, ikind, i

    complex(8),allocatable :: Sd(:,:)

    ! ! 単位構造のS行列 (silicon - lead)
    ! call smat_unit(1)%new(self%w, 2, [x(1),x(2)], self%rhos, self%kappas)
    ! 単位構造のS行列をnkind個作る
    allocate(smat_unit(self%nkind))
    do ikind=1,self%nkind
       call smat_unit(ikind)%new(self%w, 1.d0, 1.d0, 2, [x(ikind*2-1),x(ikind*2)], self%rhos, self%dls, self%dms)
    end do
    
    ! shell
    ! allocate(dict(self%nobj)) ; dict(:) = 1 ! 散乱体は1種類しかない
    call smat_shell%new(self%w, self%m_min, self%m_max, self%nobj, self%centers, self%nkind, self%dict, smat_unit)
    call smat_shellcore%new(smat_shell, self%S_core)

    allocate(Sd(smat_shellcore%m_min:smat_shellcore%m_max,smat_shellcore%m_min:smat_shellcore%m_max))
     
    
    !
    ! x(i)に関する微分
    !
    do i=1,self%n
       call smat_shellcore%calc_derivative(i, Sd)
       out(i) = 0.d0
       do m=self%m_min,self%m_max
          out(i) = out(i) + 2*real(conjg(smat_shellcore%S(m,m))*Sd(m,m))
       end do
    end do

    ! !
    ! ! x(2)に関する微分
    ! !
    ! call smat_shellcore%calc_derivative(2, Sd)
    ! out(2) = 0.d0
    ! do m=self%m_min,self%m_max
    !    out(2) = out(2) + 2*real(conjg(smat_shellcore%S(m,m))*Sd(m,m))
    ! end do

    call smat_shell%del
    call smat_unit%del
    call smat_shellcore%del
    
  end function objective_grad

  ! ic番目の制約関数の定義
  real(8) function constraint(self, ic, x)
    class(nlopt_driver_csec),intent(in) :: self
    integer,intent(in) :: ic
    real(8),intent(in) :: x(self%n)

    ! ic番目の制約 g(x) <= 0
    ! x(2*ic-1) > x(2*ic) => x(2*ic)-x(2*ic-1) < 0.0
    constraint = x(2*ic) - x(2*ic-1)
    
  end function constraint

  function constraint_grad(self, ic, x) result(out)
    class(nlopt_driver_csec),intent(in) :: self
    integer,intent(in) :: ic
    real(8),intent(in) :: x(self%n)
    real(8) :: out(self%n)

    ! ic番目の制約 g(x) <= 0
    ! x(2*ic-1) > x(2*ic) => x(2*ic)-x(2*ic-1) < 0.0
    out = 0.d0
    out(2*ic) = +1.d0
    out(2*ic-1) = -1.d0
  end function constraint_grad
  
end module nlopt_driver_csec_class
