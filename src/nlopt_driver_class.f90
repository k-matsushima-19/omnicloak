! libnlopt.aへのLIBRARY_PATHとnlopt.fへのINCLUDE PATHが必要
module nlopt_driver_class
  use iso_c_binding
  use misc
  implicit none  
  private

  include 'nlopt.f'

  integer :: step

  type,public,abstract :: nlopt_driver
     !> 変数の数
     integer :: n
     !> 各変数のlower bounds
     real(8),allocatable :: lb(:) ! (n)
     !> 各変数のupper bounds
     real(8),allocatable :: ub(:) ! (n)
     !> 変数の初期値
     real(8),allocatable :: x0(:) ! (n)
     !> 制約関数g_i(x) <= 0 の数 (0もok)
     integer :: n_constraint
     
     !> 解法 (局所最適: NLOPT_LD_SLSQP or NLOPT_LD_MMA, 大域最適: NLOPT_GN_ISRES)
     !> https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/ を参照
     integer :: ALGORITHM = NLOPT_GN_ISRES ! NLOPT_LD_MMA
     !> 最小化 or 最大化
     character(3) :: min_or_max = "min"
     !> torelence
     real(8) :: tol = 1d-4

     ! !> 今取り扱っている制約関数のindex (触らない)
     ! integer :: i_constraint
   contains
     procedure :: run
     procedure(def_del),deferred :: del
     procedure(def_objective),deferred :: objective
     procedure(def_objective_grad),deferred :: objective_grad
     procedure(def_constraint),deferred :: constraint
     procedure(def_constraint_grad),deferred :: constraint_grad
  end type nlopt_driver

  interface
     elemental subroutine def_del(self)
       import nlopt_driver
       class(nlopt_driver),intent(inout) :: self
     end subroutine def_del

     real(8) function def_objective(self, x)
       import nlopt_driver
       class(nlopt_driver),intent(in) :: self
       real(8),intent(in) :: x(self%n)
     end function def_objective

     function def_objective_grad(self, x) result(out)
       import nlopt_driver
       class(nlopt_driver),intent(in) :: self
       real(8),intent(in) :: x(self%n)
       real(8) :: out(self%n)
     end function def_objective_grad

     ! ic番目の制約関数の定義
     real(8) function def_constraint(self, ic, x)
       import nlopt_driver       
       class(nlopt_driver),intent(in) :: self
       integer,intent(in) :: ic
       real(8),intent(in) :: x(self%n)
     end function def_constraint

     function def_constraint_grad(self, ic, x) result(out)
       import nlopt_driver
       class(nlopt_driver),intent(in) :: self
       integer,intent(in) :: ic
       real(8),intent(in) :: x(self%n)
       real(8) :: out(self%n)
     end function def_constraint_grad
  end interface

  type :: data
     class(nlopt_driver),pointer :: self
     integer :: ic ! 制約関数の番号
  end type data

contains

  subroutine run(self, x_opt, f_opt)
    class(nlopt_driver),intent(in),target :: self
    real(8),intent(out) :: x_opt(self%n)
    real(8),intent(out) :: f_opt
    
    integer(8) :: opt ! work variable
    integer :: ires, ic
    type(data),allocatable :: ds(:)
    type(data) :: d
    
    !
    ! initialize
    !
    opt = 0
    call nlo_create(opt, self%ALGORITHM, self%n) 

    !
    ! upper and lower bounds
    !
    call nlo_set_lower_bounds(ires, opt, self%lb) ; call assert(ires >= 0)
    call nlo_set_upper_bounds(ires, opt, self%ub) ; call assert(ires >= 0)   

    !
    ! objective function
    !
    d%self => self
    if(self%min_or_max == "min") then
       call nlo_set_min_objective(ires, opt, myfunc, d) ; call assert(ires >= 0)
    else if(self%min_or_max == "max") then
       call nlo_set_max_objective(ires, opt, myfunc, d) ; call assert(ires >= 0)
    else
       stop "# invalid min_or_max at nlopt_driver"
    end if

    ! n_constraint個の不等式制約
    allocate(ds(self%n_constraint))
    do ic=1,self%n_constraint
       ds(ic)%self => self
       ds(ic)%ic = ic
       
       call nlo_add_inequality_constraint(ires, opt, myconstraint, ds(ic), 1d-8) ; call assert(ires >= 0)
    end do

    !
    ! tolerance
    !
    call nlo_set_xtol_rel(ires, opt, self%tol) ; call assert(ires >= 0)

    open(111,file="history.dat")

    !
    ! start
    !
    step = 0
    x_opt = self%x0
    call nlo_optimize(ires, opt, x_opt, f_opt)

    close(111)

    if (ires.lt.0) then
       write(*,*) 'nlopt failed, ires = ', ires
    end if


    ! 掃除
    do ic=1,self%n_constraint
       ds(ic)%self => null()
    end do
    deallocate(ds)
  end subroutine run

  subroutine myfunc(val, n, x, grad, need_gradient, d)
    real(8),intent(out) :: val
    integer,intent(in) :: n
    real(8),intent(in) :: x(n)
    real(8),intent(out) :: grad(n)
    integer,intent(in) :: need_gradient
    ! class(nlopt_driver),pointer :: self
    type(data) :: d

    integer :: i
    
    val = d%self%objective(x)

    if(need_gradient /= 0) then
       grad = d%self%objective_grad(x)
    end if

    ! 設計変数をスペース区切りで
    write(111,'(i6,a,1e24.16)',advance='no') step, " ", val
    do i=1,n
       write(111,'(1e24.16)',advance='no') x(i)
    end do
    ! write(111,'(1e24.16)',advance='no') val
    write(111,'(/)',advance='no') 
    ! write(111,*) ""
    flush(111)

    step = step + 1
  end subroutine myfunc

  subroutine myconstraint(val, n, x, grad, need_gradient, d)
    real(8),intent(out) :: val
    integer,intent(in) :: n
    real(8),intent(in) :: x(n)
    real(8),intent(out) :: grad(n)
    integer,intent(in) :: need_gradient
    ! class(nlopt_driver),pointer :: self
    type(data) :: d

    val = d%self%constraint(d%ic, x)

    if(need_gradient /= 0) then
       grad = d%self%constraint_grad(d%ic, x)
    end if
  end subroutine myconstraint
  
end module nlopt_driver_class
