module smatrix_unit_donut_coupled_class
  use math
  use misc
  use bessel
  use smatrix_unit_class
  ! use smatrix_shell_homogeneous_coupled_class
  ! use smatrix_core_coupled_class
  ! use smatrix_shellcore_class
  ! use smatrix_class
  implicit none
  private

  type,public,extends(smatrix_unit) :: smatrix_unit_donut_coupled
     ! 外部音響の質量密度
     real(8) :: rho
     ! 外部音響の体積弾性率
     real(8) :: kappa
     ! 外部音響の波数
     real(8) :: k
     ! 層の数 >= 1
     integer :: nshell
     ! 半径 (外側から）
     real(8),allocatable :: Rs(:) ! (1:nshell)
     ! 弾性体の密度
     real(8),allocatable :: rhos(:) ! (1:nshell)
     ! 弾性体のLame定数
     real(8),allocatable :: dls(:) ! (1:nshell)
     ! 弾性体のLame定数
     real(8),allocatable :: dms(:) ! (1:nshell)
     ! 弾性体の波数
     real(8),allocatable :: kLs(:)
     ! 弾性体の波数
     real(8),allocatable :: kTs(:)

     ! self%Sの対角成分のRsに関する勾配
     complex(8),allocatable :: S_grad(:,:) ! (nshell,m_min:m_max)
    
     ! 以下、後で消す
     ! S^elastic
     complex(8),allocatable :: S_elastic(:,:,:) ! (2,2,m_min:m_max)
     ! S^elasticのR微分
     complex(8),allocatable :: S_elastic_grad(:,:,:,:) ! (2,2,m_min:m_max,nshell)

     
   contains
     procedure :: del
     procedure :: new
     procedure :: calc_derivative
     procedure :: calc_derivative_radius

     procedure :: calc_T_layer
     procedure :: calc_T_layer_deri

     procedure :: func_T11
     procedure :: func_T12
     procedure :: func_T41
     procedure :: func_T42
     procedure :: func_Tt11
     procedure :: func_Tt12
     procedure :: func_Tt41
     procedure :: func_Tt42
     procedure :: func_W11
     procedure :: func_W12
     procedure :: func_W41
     procedure :: func_W42
     procedure :: func_Wt11
     procedure :: func_Wt12
     procedure :: func_Wt41
     procedure :: func_Wt42

     procedure :: func_T11_d
     procedure :: func_T12_d
     procedure :: func_T41_d
     procedure :: func_T42_d
     procedure :: func_Tt11_d
     procedure :: func_Tt12_d
     procedure :: func_Tt41_d
     procedure :: func_Tt42_d
     procedure :: func_W11_d
     procedure :: func_W12_d
     procedure :: func_W41_d
     procedure :: func_W42_d
     procedure :: func_Wt11_d
     procedure :: func_Wt12_d
     procedure :: func_Wt41_d
     procedure :: func_Wt42_d
  end type smatrix_unit_donut_coupled

contains
  elemental subroutine del(self)
    class(smatrix_unit_donut_coupled),intent(inout) :: self

    if(allocated(self%S)) deallocate(self%S)

    if(allocated(self%Rs)) deallocate(self%Rs)
    if(allocated(self%rhos)) deallocate(self%rhos)
    if(allocated(self%dls)) deallocate(self%dls)
    if(allocated(self%dms)) deallocate(self%dms)
    if(allocated(self%kLs)) deallocate(self%kLs)
    if(allocated(self%kTs)) deallocate(self%kTs)

    if(allocated(self%S_grad)) deallocate(self%S_grad)

    if(allocated(self%S_elastic)) deallocate(self%S_elastic)
    if(allocated(self%S_elastic_grad)) deallocate(self%S_elastic_grad)
  end subroutine del

  subroutine new(self, w, rho, kappa, nshell, Rs, rhos, dls, dms)
    class(smatrix_unit_donut_coupled),intent(inout) :: self
    real(8),intent(in) :: w
    real(8),intent(in) :: rho
    real(8),intent(in) :: kappa
    integer,intent(in) :: nshell
    real(8),intent(in) :: Rs(1:nshell)
    real(8),intent(in) :: rhos(1:nshell)
    real(8),intent(in) :: dls(1:nshell)
    real(8),intent(in) :: dms(1:nshell)

    integer :: n, ishell, jshell
    complex(8) :: T(4,4), T_elastic(4,4), S_elastic(2,2), X(5,5), X_inv(5,5), tmp(4,4), T_elastic_22_inv(2,2)

    complex(8),allocatable :: T_elastic_grad(:,:,:) ! (4,4,nshell)

    complex(8) :: adj(5) ! adjoint variable
    complex(8) :: bvec(5), bvec_d(5), Xd(5,5), fwd(5)


    call self%del

    self%w = w
    self%rho = rho
    self%kappa = kappa
    self%nshell = nshell

    allocate(self%Rs(1:self%nshell))
    self%Rs(:) = Rs(:)
    
    allocate(self%rhos(1:self%nshell))
    self%rhos(:) = rhos(:)

    allocate(self%dls(1:self%nshell))
    self%dls(:) = dls(:)

    allocate(self%dms(1:self%nshell))
    self%dms(:) = dms(:)
    
    self%R = self%Rs(1)

    self%k = w * sqrt(self%rho/self%kappa)
    allocate(self%kLs(1:self%nshell))
    allocate(self%kTs(1:self%nshell))
    self%kLs(:) = w * sqrt(self%rhos(:)/(self%dls(:)+2*self%dms(:)))
    self%kTs(:) = w * sqrt(self%rhos(:)/(              self%dms(:)))

    ! 設計変数 = 半径
    self%ndv = self%nshell

    ! Rokhlin
    self%m_max = ceiling(2*self%k*self%R + 5*log(self%k*2*self%R+pi))
    self%m_min = -self%m_max

    self%size = self%m_max - self%m_min + 1

    allocate(self%S(self%m_min:self%m_max,self%m_min:self%m_max))
    self%S(:,:) = zero

    allocate(self%S_elastic(2,2,self%m_min:self%m_max))
    allocate(T_elastic_grad(4,4,self%nshell))
    allocate(self%S_elastic_grad(2,2,self%m_min:self%m_max,self%nshell))

    allocate(self%S_grad(self%nshell,self%m_min:self%m_max))

    do n=self%m_min,self%m_max

       !
       ! S_elasticの計算
       !
       ! 単位行列
       T_elastic = 0.d0
       T_elastic(1,1) = 1.d0 ; T_elastic(2,2) = 1.d0 ; T_elastic(3,3) = 1.d0 ; T_elastic(4,4) = 1.d0
       do ishell=2,self%nshell
          ! T(ishell,ishell-1)
          call self%calc_T_layer(n, ishell, T)

          T_elastic = matmul(T,T_elastic)
          
       end do

       S_elastic = -matmul(mat_inv(2,T_elastic(3:4,3:4)),T_elastic(3:4,1:2))

       self%S_elastic(1:2,1:2,n) = S_elastic(:,:)

       X(:,1) = [&
       self%func_T11(n,self%Rs(1),1),self%func_T41(n,self%Rs(1),1),self%func_W11(n,self%Rs(1),1),&
       S_elastic(1,1), S_elastic(2,1) &
       ]

       X(:,2) = [&
       self%func_T12(n,self%Rs(1),1),self%func_T42(n,self%Rs(1),1),self%func_W12(n,self%Rs(1),1),&
       S_elastic(1,2), S_elastic(2,2) &
       ]

       X(:,3) = [&
       self%func_Tt11(n,self%Rs(1),1),self%func_Tt41(n,self%Rs(1),1),self%func_Wt11(n,self%Rs(1),1),&
       -one, zero &
       ]

       X(:,4) = [&
       self%func_Tt12(n,self%Rs(1),1),self%func_Tt42(n,self%Rs(1),1),self%func_Wt12(n,self%Rs(1),1),&
       zero, -one &
       ]

       X(:,5) = [besh(n,self%k*self%Rs(1)),zero, -self%k/(self%rho*self%w**2)*besh_d(n,self%k*self%Rs(1)),zero,zero]
      
      X_inv = mat_inv(5,X)

      self%S(n:n,n:n) = rdot_product(X_inv(5,1:4), &
      [-one*besj(n,self%k*self%Rs(1)),zero,one*self%k/(self%rho*self%w**2)*besj_d(n,self%k*self%Rs(1)),zero])
      
      !
      ! 半径Rに関する微分を計算
      !
      ! T^elasticのRs(jshell)に関する微分
      do jshell=1,self%nshell
        

        if(jshell == 1) then
          ! T^elasticはRs(1)に依存しない
          T_elastic_grad(:,:,1) = 0.d0

        else 
          ! T(jshell-1,jshell-2)*...*T(3,2)*T(2,1)を計算
          ! 単位行列
          tmp = 0.d0
          tmp(1,1) = 1.d0 ; tmp(2,2) = 1.d0 ; tmp(3,3) = 1.d0 ; tmp(4,4) = 1.d0
          do ishell=2,jshell-1
            ! T(ishell,ishell-1)
            call self%calc_T_layer(n, ishell, T)

            tmp = matmul(T,tmp)
          end do

          ! T(jshell,jshell-1)の微分をかける
          call self%calc_T_layer_deri(n, jshell, T)
          tmp = matmul(T,tmp)

          ! T(nshell,nshell-1)*...*T(jshell+1,jshell)をかける
          do ishell=jshell+1,self%nshell
            ! T(ishell,ishell-1)
            call self%calc_T_layer(n, ishell, T)

            tmp = matmul(T,tmp)
          end do

          T_elastic_grad(:,:,jshell) = tmp(:,:)

        end if

      end do

      ! S^elasticのRs(jshell)に関する微分
      do jshell=1,self%nshell
        T_elastic_22_inv = mat_inv(2,T_elastic(3:4,3:4))
        self%S_elastic_grad(:,:,n,jshell) = &
        matmul(T_elastic_22_inv, &
        matmul(T_elastic_grad(3:4,3:4,jshell),matmul(T_elastic_22_inv,T_elastic(3:4,1:2)))-T_elastic_grad(3:4,1:2,jshell) )
      end do

      !
      ! Solve the adjoint problem
      !
      adj = matmul(mat_inv(5,transpose(X)), [zero,zero,zero,zero,one])
      
      ! RHS vector
      bvec = [-one*besj(n,self%k*self%Rs(1)),zero,one*self%k/(self%rho*self%w**2)*besj_d(n,self%k*self%Rs(1)),zero,zero]

      ! X^-1*b
      fwd = matmul(X_inv,bvec)

      do jshell=1,self%nshell
        !
        ! XのRs(jshell)微分
        !
        if(jshell == 1) then
          Xd(:,1) = [&
          self%func_T11_d(n,self%Rs(1),1),self%func_T41_d(n,self%Rs(1),1),self%func_W11_d(n,self%Rs(1),1),&
          self%S_elastic_grad(1,1,n,jshell), self%S_elastic_grad(2,1,n,jshell) &
          ]

          Xd(:,2) = [&
          self%func_T12_d(n,self%Rs(1),1),self%func_T42_d(n,self%Rs(1),1),self%func_W12_d(n,self%Rs(1),1),&
          self%S_elastic_grad(1,2,n,jshell), self%S_elastic_grad(2,2,n,jshell) &
          ]

          Xd(:,3) = [&
          self%func_Tt11_d(n,self%Rs(1),1),self%func_Tt41_d(n,self%Rs(1),1),self%func_Wt11_d(n,self%Rs(1),1),&
          zero, zero &
          ]

          Xd(:,4) = [&
          self%func_Tt12_d(n,self%Rs(1),1),self%func_Tt42_d(n,self%Rs(1),1),self%func_Wt12_d(n,self%Rs(1),1),&
          zero, zero &
          ]

          Xd(:,5) = [self%k*besh_d(n,self%k*self%Rs(1)),zero,-self%k**2/(self%rho*self%w**2)*besh_dd(n,self%k*self%Rs(1)),zero,zero]

        else
          Xd(:,1) = [&
          zero, zero, zero,&
          self%S_elastic_grad(1,1,n,jshell), self%S_elastic_grad(2,1,n,jshell) &
          ]

          Xd(:,2) = [&
          zero, zero, zero,&
          self%S_elastic_grad(1,2,n,jshell), self%S_elastic_grad(2,2,n,jshell) &
          ]

          Xd(:,3) = [&
          zero, zero, zero,&
          zero, zero &
          ]

          Xd(:,4) = [&
          zero, zero, zero,&
          zero, zero &
          ]

          Xd(:,5) = [zero,zero, zero,zero,zero]
        end if

        !
        ! bvecのRs(jshell)微分
        !
        if(jshell == 1) then
          bvec_d = [-one*self%k*besj_d(n,self%k*self%Rs(1)),zero,&
          one*self%k**2/(self%rho*self%w**2)*besj_dd(n,self%k*self%Rs(1)),zero,zero]
        else
          bvec_d = [zero,zero,zero,zero,zero]
        end if

        
        self%S_grad(jshell,n) = -rdot_product(adj, matmul(Xd,fwd)-bvec_d)

        

      end do

    end do
    
  end subroutine new

  ! idv番設計変数に関するSの微分
  subroutine calc_derivative(self, idv, out)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: idv
    complex(8),intent(out) :: out(self%m_min:self%m_max,self%m_min:self%m_max)

    call self%calc_derivative_radius(idv, out)
  end subroutine calc_derivative

  ! i番半径に関するS行列の導関数
  subroutine calc_derivative_radius(self, i, out)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: i
    complex(8),intent(out) :: out(self%m_min:self%m_max,self%m_min:self%m_max)

    integer :: n

    ! stop "todo"
       
    out = zero
    do n=self%m_min,self%m_max
       out(n,n) = self%S_grad(i,n)
    end do

    
  end subroutine calc_derivative_radius

  ! T(ishell,ishell-1)を計算
  subroutine calc_T_layer(self, n, ishell, T)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    integer,intent(in) :: ishell
    complex(8),intent(out) :: T(4,4)
    
    complex(8) :: X(4,4), Y(4,4), X_inv(4,4)

    call assert(2 <= ishell .and. ishell <= self%nshell)

    ! T = X^-1 * Y

    X(:,1) = [ &
         self%func_W11(n,self%Rs(ishell),ishell),self%func_W41(n,self%Rs(ishell),ishell), &
         self%func_T11(n,self%Rs(ishell),ishell),self%func_T41(n,self%Rs(ishell),ishell)]
    X(:,2) = [ &
         self%func_W12(n,self%Rs(ishell),ishell),self%func_W42(n,self%Rs(ishell),ishell), &
         self%func_T12(n,self%Rs(ishell),ishell),self%func_T42(n,self%Rs(ishell),ishell)]
    X(:,3) = [ &
         self%func_Wt11(n,self%Rs(ishell),ishell),self%func_Wt41(n,self%Rs(ishell),ishell), &
         self%func_Tt11(n,self%Rs(ishell),ishell),self%func_Tt41(n,self%Rs(ishell),ishell)]
    X(:,4) = [ &
         self%func_Wt12(n,self%Rs(ishell),ishell),self%func_Wt42(n,self%Rs(ishell),ishell), &
         self%func_Tt12(n,self%Rs(ishell),ishell),self%func_Tt42(n,self%Rs(ishell),ishell)]

    
    Y(:,1) = [ &
         self%func_W11(n,self%Rs(ishell),ishell-1),self%func_W41(n,self%Rs(ishell),ishell-1), &
         self%func_T11(n,self%Rs(ishell),ishell-1),self%func_T41(n,self%Rs(ishell),ishell-1)]
    Y(:,2) = [ &
         self%func_W12(n,self%Rs(ishell),ishell-1),self%func_W42(n,self%Rs(ishell),ishell-1), &
         self%func_T12(n,self%Rs(ishell),ishell-1),self%func_T42(n,self%Rs(ishell),ishell-1)]
    Y(:,3) = [ &
         self%func_Wt11(n,self%Rs(ishell),ishell-1),self%func_Wt41(n,self%Rs(ishell),ishell-1), &
         self%func_Tt11(n,self%Rs(ishell),ishell-1),self%func_Tt41(n,self%Rs(ishell),ishell-1)]
    Y(:,4) = [ &
         self%func_Wt12(n,self%Rs(ishell),ishell-1),self%func_Wt42(n,self%Rs(ishell),ishell-1), &
         self%func_Tt12(n,self%Rs(ishell),ishell-1),self%func_Tt42(n,self%Rs(ishell),ishell-1)]

    X_inv = mat_inv(4,X)

    T = matmul(X_inv,Y)
    
    
  end subroutine calc_T_layer

  ! T(ishell,ishell-1)のRs(ishell)に関する微分を計算
  subroutine calc_T_layer_deri(self, n, ishell, T)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    integer,intent(in) :: ishell
    complex(8),intent(out) :: T(4,4)
    
    complex(8) :: X(4,4), Y(4,4), X_inv(4,4)
    complex(8) :: Xd(4,4), Yd(4,4)

    call assert(2 <= ishell .and. ishell <= self%nshell)

    ! T = X^-1 * Y

    X(:,1) = [ &
         self%func_W11(n,self%Rs(ishell),ishell),self%func_W41(n,self%Rs(ishell),ishell), &
         self%func_T11(n,self%Rs(ishell),ishell),self%func_T41(n,self%Rs(ishell),ishell)]
    X(:,2) = [ &
         self%func_W12(n,self%Rs(ishell),ishell),self%func_W42(n,self%Rs(ishell),ishell), &
         self%func_T12(n,self%Rs(ishell),ishell),self%func_T42(n,self%Rs(ishell),ishell)]
    X(:,3) = [ &
         self%func_Wt11(n,self%Rs(ishell),ishell),self%func_Wt41(n,self%Rs(ishell),ishell), &
         self%func_Tt11(n,self%Rs(ishell),ishell),self%func_Tt41(n,self%Rs(ishell),ishell)]
    X(:,4) = [ &
         self%func_Wt12(n,self%Rs(ishell),ishell),self%func_Wt42(n,self%Rs(ishell),ishell), &
         self%func_Tt12(n,self%Rs(ishell),ishell),self%func_Tt42(n,self%Rs(ishell),ishell)]

    
    Y(:,1) = [ &
         self%func_W11(n,self%Rs(ishell),ishell-1),self%func_W41(n,self%Rs(ishell),ishell-1), &
         self%func_T11(n,self%Rs(ishell),ishell-1),self%func_T41(n,self%Rs(ishell),ishell-1)]
    Y(:,2) = [ &
         self%func_W12(n,self%Rs(ishell),ishell-1),self%func_W42(n,self%Rs(ishell),ishell-1), &
         self%func_T12(n,self%Rs(ishell),ishell-1),self%func_T42(n,self%Rs(ishell),ishell-1)]
    Y(:,3) = [ &
         self%func_Wt11(n,self%Rs(ishell),ishell-1),self%func_Wt41(n,self%Rs(ishell),ishell-1), &
         self%func_Tt11(n,self%Rs(ishell),ishell-1),self%func_Tt41(n,self%Rs(ishell),ishell-1)]
    Y(:,4) = [ &
         self%func_Wt12(n,self%Rs(ishell),ishell-1),self%func_Wt42(n,self%Rs(ishell),ishell-1), &
         self%func_Tt12(n,self%Rs(ishell),ishell-1),self%func_Tt42(n,self%Rs(ishell),ishell-1)]

    X_inv = mat_inv(4,X)

    ! T = matmul(X_inv,Y)
    
    ! 微分
    Xd(:,1) = [ &
         self%func_W11_d(n,self%Rs(ishell),ishell),self%func_W41_d(n,self%Rs(ishell),ishell), &
         self%func_T11_d(n,self%Rs(ishell),ishell),self%func_T41_d(n,self%Rs(ishell),ishell)]
    Xd(:,2) = [ &
         self%func_W12_d(n,self%Rs(ishell),ishell),self%func_W42_d(n,self%Rs(ishell),ishell), &
         self%func_T12_d(n,self%Rs(ishell),ishell),self%func_T42_d(n,self%Rs(ishell),ishell)]
    Xd(:,3) = [ &
         self%func_Wt11_d(n,self%Rs(ishell),ishell),self%func_Wt41_d(n,self%Rs(ishell),ishell), &
         self%func_Tt11_d(n,self%Rs(ishell),ishell),self%func_Tt41_d(n,self%Rs(ishell),ishell)]
    Xd(:,4) = [ &
         self%func_Wt12_d(n,self%Rs(ishell),ishell),self%func_Wt42_d(n,self%Rs(ishell),ishell), &
         self%func_Tt12_d(n,self%Rs(ishell),ishell),self%func_Tt42_d(n,self%Rs(ishell),ishell)]

    
    Yd(:,1) = [ &
         self%func_W11_d(n,self%Rs(ishell),ishell-1),self%func_W41_d(n,self%Rs(ishell),ishell-1), &
         self%func_T11_d(n,self%Rs(ishell),ishell-1),self%func_T41_d(n,self%Rs(ishell),ishell-1)]
    Yd(:,2) = [ &
         self%func_W12_d(n,self%Rs(ishell),ishell-1),self%func_W42_d(n,self%Rs(ishell),ishell-1), &
         self%func_T12_d(n,self%Rs(ishell),ishell-1),self%func_T42_d(n,self%Rs(ishell),ishell-1)]
    Yd(:,3) = [ &
         self%func_Wt11_d(n,self%Rs(ishell),ishell-1),self%func_Wt41_d(n,self%Rs(ishell),ishell-1), &
         self%func_Tt11_d(n,self%Rs(ishell),ishell-1),self%func_Tt41_d(n,self%Rs(ishell),ishell-1)]
    Yd(:,4) = [ &
         self%func_Wt12_d(n,self%Rs(ishell),ishell-1),self%func_Wt42_d(n,self%Rs(ishell),ishell-1), &
         self%func_Tt12_d(n,self%Rs(ishell),ishell-1),self%func_Tt42_d(n,self%Rs(ishell),ishell-1)]

    T = matmul(X_inv, -matmul(Xd,matmul(X_inv,Y)) + Yd )
    
  end subroutine calc_T_layer_deri

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
  
  complex(8) function func_T11(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T11 = (n**2+n-0.5d0*(kT*r)**2) * besj(n,kL*r) - kL*r*besj(n-1,kL*r)

    func_T11 = func_T11 * 2*dm/r**2
  end function func_T11

  complex(8) function func_T11_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T11_d = -2*self%func_T11(n,r,i)/r + 2*dm/r**2* &
    ( -kT**2*r*besj(n,kL*r) + (n**2+n-0.5d0*(kT*r)**2)*kL*besj_d(n,kL*r) - kL*besj(n-1,kL*r) - kL**2*r*besj_d(n-1,kL*r)  )
  end function func_T11_d

  complex(8) function func_T12(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T12 = -ione*n*( (n+1)*besj(n,kT*r) - kT*r*besj(n-1,kT*r)  )

    func_T12 = func_T12 * 2*dm/r**2
  end function func_T12

  complex(8) function func_T12_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T12_d = -2*self%func_T12(n,r,i)/r - ione*n* 2*dm/r**2* &
    ( (n+1)*kT*besj_d(n,kT*r) - kT*besj(n-1,kT*r) - kT**2*r*besj_d(n-1,kT*r) )
  end function func_T12_d

  complex(8) function func_T41(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T41 = -ione*n*( (n+1)*besj(n,kL*r) - kL*r*besj(n-1,kL*r)  )
    
    func_T41 = func_T41 * 2*dm/r**2
  end function func_T41

  complex(8) function func_T41_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T41_d = -2*self%func_T41(n,r,i)/r - ione*n* 2*dm/r**2* &
    ( (n+1)*kL*besj_d(n,kL*r) - kL*besj(n-1,kL*r) - kL**2*r*besj_d(n-1,kL*r) )
  end function func_T41_d

  complex(8) function func_T42(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T42 = -(n**2+n-0.5d0*(kT*r)**2) * besj(n,kT*r) + kT*r*besj(n-1,kT*r)

    func_T42 = func_T42 * 2*dm/r**2
  end function func_T42

  complex(8) function func_T42_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_T42_d = -2*self%func_T42(n,r,i)/r - 2*dm/r**2* &
    ( -kT**2*r*besj(n,kT*r) + (n**2+n-0.5d0*(kT*r)**2)*kT*besj_d(n,kT*r) - kT*besj(n-1,kT*r) - kT**2*r*besj_d(n-1,kT*r)  )
  end function func_T42_d

  complex(8) function func_Tt11(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt11 = (n**2+n-0.5d0*(kT*r)**2) * besh(n,kL*r) - kL*r*besh(n-1,kL*r)

    func_Tt11 = func_Tt11 * 2*dm/r**2
  end function func_Tt11

  complex(8) function func_Tt12(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt12 = -ione*n*( (n+1)*besh(n,kT*r) - kT*r*besh(n-1,kT*r)  )

    func_Tt12 = func_Tt12 * 2*dm/r**2
  end function func_Tt12

  complex(8) function func_Tt41(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt41 = -ione*n*( (n+1)*besh(n,kL*r) - kL*r*besh(n-1,kL*r)  )

    func_Tt41 = func_Tt41 * 2*dm/r**2
  end function func_Tt41

  complex(8) function func_Tt42(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt42 = -(n**2+n-0.5d0*(kT*r)**2) * besh(n,kT*r) + kT*r*besh(n-1,kT*r)

    func_Tt42 = func_Tt42 * 2*dm/r**2
  end function func_Tt42

  complex(8) function func_Tt11_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt11_d = -2*self%func_Tt11(n,r,i)/r + 2*dm/r**2* &
    ( -kT**2*r*besh(n,kL*r) + (n**2+n-0.5d0*(kT*r)**2)*kL*besh_d(n,kL*r) - kL*besh(n-1,kL*r) - kL**2*r*besh_d(n-1,kL*r)  )
  end function func_Tt11_d

  complex(8) function func_Tt12_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt12_d = -2*self%func_Tt12(n,r,i)/r - ione*n* 2*dm/r**2* &
    ( (n+1)*kT*besh_d(n,kT*r) - kT*besh(n-1,kT*r) - kT**2*r*besh_d(n-1,kT*r) )
  end function func_Tt12_d

  complex(8) function func_Tt41_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt41_d = -2*self%func_Tt41(n,r,i)/r - ione*n* 2*dm/r**2* &
    ( (n+1)*kL*besh_d(n,kL*r) - kL*besh(n-1,kL*r) - kL**2*r*besh_d(n-1,kL*r) )
  end function func_Tt41_d

  complex(8) function func_Tt42_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Tt42_d = -2*self%func_Tt42(n,r,i)/r - 2*dm/r**2* &
    ( -kT**2*r*besh(n,kT*r) + (n**2+n-0.5d0*(kT*r)**2)*kT*besh_d(n,kT*r) - kT*besh(n-1,kT*r) - kT**2*r*besh_d(n-1,kT*r)  )
  end function func_Tt42_d

  complex(8) function func_W11(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W11 = kL * besj_d(n,kL*r)
  end function func_W11

  complex(8) function func_W11_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W11_d = kL**2 * besj_dd(n,kL*r)
  end function func_W11_d

  complex(8) function func_W12(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W12 = ione*n/r * besj(n,kT*r)
  end function func_W12

  complex(8) function func_W12_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W12_d = -ione*n/r**2*besj(n,kT*r) + ione*n/r * kT * besj_d(n,kT*r)
  end function func_W12_d

  complex(8) function func_W41(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W41 = ione*n/r * besj(n,kL*r)
  end function func_W41

  complex(8) function func_W41_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W41_d = -ione*n/r**2*besj(n,kL*r) + ione*n/r * kL * besj_d(n,kL*r)
  end function func_W41_d

  complex(8) function func_W42(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W42 = -kT * besj_d(n,kT*r)
  end function func_W42

  complex(8) function func_W42_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_W42_d = -kT**2 * besj_dd(n,kT*r)
  end function func_W42_d

  complex(8) function func_Wt11(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt11 = kL * besh_d(n,kL*r)
  end function func_Wt11

  complex(8) function func_Wt12(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt12 = ione*n/r * besh(n,kT*r)
  end function func_Wt12

  complex(8) function func_Wt41(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt41 = ione*n/r * besh(n,kL*r)
  end function func_Wt41

  complex(8) function func_Wt42(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt42 = -kT * besh_d(n,kT*r)
  end function func_Wt42

  complex(8) function func_Wt11_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt11_d = kL**2 * besh_dd(n,kL*r)
  end function func_Wt11_d

  complex(8) function func_Wt12_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt12_d = -ione*n/r**2*besh(n,kT*r) + ione*n/r * kT * besh_d(n,kT*r)
  end function func_Wt12_d

  complex(8) function func_Wt41_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt41_d = -ione*n/r**2*besh(n,kL*r) + ione*n/r * kL * besh_d(n,kL*r)
  end function func_Wt41_d

  complex(8) function func_Wt42_d(self, n, r, i)
    class(smatrix_unit_donut_coupled),intent(in) :: self
    integer,intent(in) :: n
    real(8),intent(in) :: r
    integer,intent(in) :: i

    real(8) :: dm, kL, kT

    dm = self%dms(i)
    kL = self%kLs(i)
    kT = self%kTs(i)

    func_Wt42_d = -kT**2 * besh_dd(n,kT*r)
  end function func_Wt42_d

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
  
end module smatrix_unit_donut_coupled_class
