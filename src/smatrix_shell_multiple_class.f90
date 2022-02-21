!> 多重散乱理論からshell S行列を構成する
module smatrix_shell_multiple_class
  use math
  use misc
  use bessel
  use smatrix_shell_class
  use smatrix_unit_class
  implicit none
  private

  type,public,extends(smatrix_shell) :: smatrix_shell_multiple
     !> このshellを構成する散乱体の数
     integer :: nobj
     !> 各散乱体が置かれる(中心)座標 (2,nobj)
     real(8),allocatable :: centers(:,:)
     !> 何種類の散乱体があるか
     integer :: nkind
     !> 散乱体の番号と種類を対応させる配列．iobj番散乱体はdict(iobj)種
     integer,allocatable :: dict(:)
     !> nkind個の散乱行列
     class(smatrix_unit),pointer :: smats(:) => null()

     !> I-STのサイズ
     integer :: size_all
     !> translation行列 (size_all,size_all)
     complex(8),allocatable :: T(:,:)
     !> 行列G (size,size_all)
     complex(8),allocatable :: G(:,:)
     !> 行列H (size,size_all)
     complex(8),allocatable :: H(:,:)
     !> 行列J (size_all,size)
     complex(8),allocatable :: J(:,:)
     !> 行列J (size_all,size)
     complex(8),allocatable :: K(:,:)
     !> (I-ST)のLU
     complex(8),allocatable :: amat_inv(:,:)
     !> (I-ST)のLUのipiv
     integer,allocatable :: ipiv(:)

     !> idv番設計変数は第dict_dv(1,idv)種散乱体のdict_dv(2,idv)番設計変数
     integer,allocatable :: dict_dv(:,:) ! (2,ndv)
   contains
     procedure :: del
     procedure :: new
     procedure :: calc_tmatrix
     procedure :: calc_G
     procedure :: calc_H
     procedure :: calc_J
     procedure :: calc_K
     procedure :: product_S
     procedure :: product_S_derivative
     procedure :: calc_derivative
  end type smatrix_shell_multiple

contains
  elemental subroutine del(self)
    class(smatrix_shell_multiple),intent(inout) :: self

    if(allocated(self%S11)) deallocate(self%S11)
    if(allocated(self%S21)) deallocate(self%S21)
    if(allocated(self%S12)) deallocate(self%S12)
    if(allocated(self%S22)) deallocate(self%S22)

    if(allocated(self%centers)) deallocate(self%centers)
    if(allocated(self%dict)) deallocate(self%dict)
    if(associated(self%smats)) nullify(self%smats)
    if(allocated(self%T)) deallocate(self%T)
    if(allocated(self%G)) deallocate(self%G)
    if(allocated(self%H)) deallocate(self%H)
    if(allocated(self%J)) deallocate(self%J)
    if(allocated(self%K)) deallocate(self%K)

    if(allocated(self%amat_inv)) deallocate(self%amat_inv)
    if(allocated(self%ipiv)) deallocate(self%ipiv)
    if(allocated(self%dict_dv)) deallocate(self%dict_dv)
    
  end subroutine del

  subroutine new(self, w, m_min, m_max, nobj, centers, nkind, dict, smats)
    class(smatrix_shell_multiple),intent(inout) :: self
    real(8),intent(in) :: w
    integer,intent(in) :: m_min
    integer,intent(in) :: m_max
    integer,intent(in) :: nobj
    real(8),intent(in) :: centers(2,nobj)
    integer,intent(in) :: nkind
    integer,intent(in) :: dict(nobj)
    class(smatrix_unit),intent(in),target :: smats(:)

    integer :: iobj, jobj, ioff, joff, i, j, m, ikind, idv
    complex(8),allocatable :: amat(:,:)
    integer :: info
    integer,allocatable :: ipiv(:)

    complex(8),allocatable :: tmp(:,:)

    call self%del

    self%w = w
    self%m_min = m_min ! self%smats(:)%m_minとは異なることに注意
    self%m_max = m_max ! self%smats(:)%m_maxとは異なることに注意
    self%nobj = nobj
    self%centers = centers
    self%nkind = nkind
    self%dict = dict
    self%smats => smats

    self%size = self%m_max - self%m_min + 1

    self%size_all = 0
    do iobj=1,self%nobj
       self%size_all = self%size_all + smats(self%dict(iobj))%size
    end do

    ! 内径
    self%Rin = sqrt(dot_product(self%centers(:,1),self%centers(:,1)))
    self%Rout = 0.0d0
    do iobj=1,self%nobj
       self%Rin = min(self%Rin,sqrt(dot_product(self%centers(:,iobj),self%centers(:,iobj)))-self%smats(self%dict(iobj))%R)
       self%Rout = max(self%Rout,sqrt(dot_product(self%centers(:,iobj),self%centers(:,iobj)))+self%smats(self%dict(iobj))%R)
    end do

    ! 設計変数の総数
    self%ndv = 0
    do ikind=1,self%nkind
       self%ndv = self%ndv + self%smats(ikind)%ndv

       ! write(*,*) self%smats(ikind)%ndv
    end do

    ! write(*,*) "ndv", self%ndv

    ! idv番設計変数は第dict_dv(1,idv)種散乱体のdict_dv(2,idv)番設計変数
    allocate(self%dict_dv(2,self%ndv))
    i = 0
    do ikind=1,self%nkind
       do idv=1,self%smats(ikind)%ndv
          i = i + 1
          self%dict_dv(1,i) = ikind
          self%dict_dv(2,i) = idv
       end do       
    end do

    ! write(*,*) self%dict_dv(:,1)
    ! write(*,*) self%dict_dv(:,2)
    
    !
    ! Tの計算
    !
    allocate(self%T(self%size_all,self%size_all))    
    joff = 0
    do jobj=1,self%nobj
       ioff = 0
       do iobj=1,self%nobj
          ! 部分行列T(ioff+1:ioff+self%smats(self%dict(iobj))%size,joff+1:joff+self%smats(self%dict(jobj))%size) = Tij
          if(iobj == jobj) then
             self%T(ioff+1:ioff+self%smats(self%dict(iobj))%size,joff+1:joff+self%smats(self%dict(jobj))%size) = zero
          else
             call self%calc_tmatrix(self%smats(self%dict(iobj))%m_min,self%smats(self%dict(iobj))%m_max,&
                  self%smats(self%dict(jobj))%m_min,self%smats(self%dict(jobj))%m_max,&
                  self%centers(:,iobj) - self%centers(:,jobj), &
                  self%T(ioff+1:ioff+self%smats(self%dict(iobj))%size,joff+1:joff+self%smats(self%dict(jobj))%size))
          end if
          
          
          ioff = ioff + self%smats(self%dict(iobj))%size   
       end do
       joff = joff + self%smats(self%dict(jobj))%size
    end do

    !
    ! 行列Gの計算
    !
    allocate(self%G(self%m_min:self%m_max,self%size_all))
    ioff = 0
    do iobj=1,self%nobj
       ! 部分行列G(:,ioff+1:ioff+self%smats(self%dict(iobj))%size) = G^(i)
       call self%calc_G(self%m_min, self%m_max, &
            self%smats(self%dict(iobj))%m_min, self%smats(self%dict(iobj))%m_max, &
            self%centers(:,iobj), &
            self%G(:,ioff+1:ioff+self%smats(self%dict(iobj))%size))

       ioff = ioff + self%smats(self%dict(iobj))%size   
    end do

    !
    ! 行列Hの計算
    !
    allocate(self%H(self%m_min:self%m_max,self%size_all))
    ioff = 0
    do iobj=1,self%nobj
       ! 部分行列H(:,ioff+1:ioff+self%smats(self%dict(iobj))%size) = H^(i)
       call self%calc_H(self%m_min, self%m_max, &
            self%smats(self%dict(iobj))%m_min, self%smats(self%dict(iobj))%m_max, &
            self%centers(:,iobj), &
            self%H(:,ioff+1:ioff+self%smats(self%dict(iobj))%size))

       ioff = ioff + self%smats(self%dict(iobj))%size   
    end do

    !
    ! 行列Jの計算
    !
    allocate(self%J(self%size_all,self%m_min:self%m_max))
    ioff = 0
    do iobj=1,self%nobj
       ! 部分行列J(ioff+1:ioff+self%smats(self%dict(iobj))%size,:) = J^(i)
       call self%calc_J(self%smats(self%dict(iobj))%m_min, self%smats(self%dict(iobj))%m_max, &
            self%m_min, self%m_max, &
            self%centers(:,iobj), &
            self%J(ioff+1:ioff+self%smats(self%dict(iobj))%size,:))

       ioff = ioff + self%smats(self%dict(iobj))%size   
    end do

    !
    ! 行列Kの計算
    !
    allocate(self%K(self%size_all,self%m_min:self%m_max))
    ioff = 0
    do iobj=1,self%nobj
       ! 部分行列K(ioff+1:ioff+self%smats(self%dict(iobj))%size,:) = K^(i)
       call self%calc_K(self%smats(self%dict(iobj))%m_min, self%smats(self%dict(iobj))%m_max, &
            self%m_min, self%m_max, &
            self%centers(:,iobj), &
            self%K(ioff+1:ioff+self%smats(self%dict(iobj))%size,:))

       ioff = ioff + self%smats(self%dict(iobj))%size   
    end do

    !
    ! self%amat_inv = (I-ST)^-1
    !
    allocate(self%amat_inv(self%size_all,self%size_all), self%ipiv(self%size_all))
    ! -ST
    do i=1,self%size_all
       call self%product_S(self%T(:,i), self%amat_inv(:,i))
    end do
    self%amat_inv = -self%amat_inv
    ! I - ST
    do i=1,self%size_all
       self%amat_inv(i,i) = self%amat_inv(i,i) + one
    end do
    ! LU分解
    call zgetrf(self%size_all,self%size_all,self%amat_inv,self%size_all,self%ipiv,info) ; call assert(info == 0)

    allocate(tmp(self%size_all,self%m_min:self%m_max))

    !
    ! S11
    !
    allocate(self%S11(self%m_min:self%m_max,self%m_min:self%m_max))
    ! tmp = S*J
    do m=self%m_min,self%m_max
       call self%product_S(self%J(:,m), tmp(:,m))       
    end do
    ! tmp = self%amat_inv*S*J
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp, self%size_all, info)
    call assert(info == 0)
    ! S11 = G*tmp
    self%S11 = matmul(self%G,tmp)

    !
    ! S21
    !
    allocate(self%S21(self%m_min:self%m_max,self%m_min:self%m_max))
    ! tmp = S*J
    do m=self%m_min,self%m_max
       call self%product_S(self%J(:,m), tmp(:,m))       
    end do
    ! tmp = self%amat_inv*S*J
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp, self%size_all, info)
    call assert(info == 0)
    ! S21 = H*tmp
    self%S21 = matmul(self%H,tmp)
    ! +I
    do m=self%m_min,self%m_max
       self%S21(m,m) = self%S21(m,m) + 1.d0
    end do

    !
    ! S12
    !
    allocate(self%S12(self%m_min:self%m_max,self%m_min:self%m_max))
    ! tmp = S*K
    do m=self%m_min,self%m_max
       call self%product_S(self%K(:,m), tmp(:,m))       
    end do
    ! tmp = self%amat_inv*S*K
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp, self%size_all, info)
    call assert(info == 0)
    ! S12 = G*tmp
    self%S12 = matmul(self%G,tmp)
    ! +I
    do m=self%m_min,self%m_max
       self%S12(m,m) = self%S12(m,m) + 1.d0
    end do

    !
    ! S22
    !
    allocate(self%S22(self%m_min:self%m_max,self%m_min:self%m_max))
    ! tmp = S*K
    do m=self%m_min,self%m_max
       call self%product_S(self%K(:,m), tmp(:,m))       
    end do
    ! tmp = self%amat_inv*S*K
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp, self%size_all, info)
    call assert(info == 0)
    ! S22 = H*tmp
    self%S22 = matmul(self%H,tmp)
    
  end subroutine new

  !> Translation matrix \f$ T_{mm^\prime} =  H^{(1)}_{m-m^\prime}(-ka)\exp(\mathrm{i}(m-m^\prime)\theta(-a)) \f$の計算
  !!\param[in] self
  !!\param[in] a ベクトル\f$ a \f$
  !!\param[out] out 行列\f$ T_{mm^\prime} \f$  
  subroutine calc_tmatrix(self, m_min_i, m_max_i, m_min_j, m_max_j, a, out)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: m_min_i
    integer,intent(in) :: m_max_i
    integer,intent(in) :: m_min_j
    integer,intent(in) :: m_max_j
    real(8),intent(in) :: a(2)
    complex(8),intent(out) :: out(m_min_i:m_max_i,m_min_j:m_max_j)
    
    complex(8) :: ctmp(min(m_min_i,m_min_j)-max(m_max_i,m_max_j):-min(m_min_i,m_min_j)+max(m_max_i,m_max_j))
    real(8) :: theta
    integer :: m, n

    ! Hankel関数の計算
    call cbeshn(one*self%w*length_vector(a), min(m_min_i,m_min_j)-max(m_max_i,m_max_j),&
         -min(m_min_i,m_min_j)+max(m_max_i,m_max_j), ctmp)
    theta = atan2(a(2),a(1))
    
    do m=min(m_min_i,m_min_j)-max(m_max_i,m_max_j),-min(m_min_i,m_min_j)+max(m_max_i,m_max_j)
       ctmp(m) = ctmp(m) * exp(ione*m*theta)
    end do

    do n=m_min_j,m_max_j
       do m=m_min_i,m_max_i
          out(m,n) = ctmp(n-m)
       end do
    end do
    
    
  end subroutine calc_tmatrix
  
  
  subroutine calc_G(self, m_min_i, m_max_i, m_min_j, m_max_j, a, out)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: m_min_i
    integer,intent(in) :: m_max_i
    integer,intent(in) :: m_min_j
    integer,intent(in) :: m_max_j
    real(8),intent(in) :: a(2)
    complex(8),intent(out) :: out(m_min_i:m_max_i,m_min_j:m_max_j)
    
    complex(8) :: ctmp(min(m_min_i,m_min_j)-max(m_max_i,m_max_j):-min(m_min_i,m_min_j)+max(m_max_i,m_max_j))
    real(8) :: theta
    integer :: m, n

    ! Hankel関数の計算
    call cbesjn(one*self%w*length_vector(a), min(m_min_i,m_min_j)-max(m_max_i,m_max_j),&
         -min(m_min_i,m_min_j)+max(m_max_i,m_max_j), ctmp)
    theta = atan2(-a(2),-a(1))
    
    do m=min(m_min_i,m_min_j)-max(m_max_i,m_max_j),-min(m_min_i,m_min_j)+max(m_max_i,m_max_j)
       ctmp(m) = ctmp(m) * exp(ione*m*theta)
    end do

    do n=m_min_j,m_max_j
       do m=m_min_i,m_max_i
          out(m,n) = ctmp(n-m)
       end do
    end do
    
    
  end subroutine calc_G

  subroutine calc_H(self, m_min_i, m_max_i, m_min_j, m_max_j, a, out)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: m_min_i
    integer,intent(in) :: m_max_i
    integer,intent(in) :: m_min_j
    integer,intent(in) :: m_max_j
    real(8),intent(in) :: a(2)
    complex(8),intent(out) :: out(m_min_i:m_max_i,m_min_j:m_max_j)
    
    complex(8) :: ctmp(min(m_min_i,m_min_j)-max(m_max_i,m_max_j):-min(m_min_i,m_min_j)+max(m_max_i,m_max_j))
    real(8) :: theta
    integer :: m, n

    ! Hankel関数の計算
    call cbeshn(one*self%w*length_vector(a), min(m_min_i,m_min_j)-max(m_max_i,m_max_j),&
         -min(m_min_i,m_min_j)+max(m_max_i,m_max_j), ctmp)
    theta = atan2(-a(2),-a(1))
    
    do m=min(m_min_i,m_min_j)-max(m_max_i,m_max_j),-min(m_min_i,m_min_j)+max(m_max_i,m_max_j)
       ctmp(m) = ctmp(m) * exp(ione*m*theta)
    end do

    do n=m_min_j,m_max_j
       do m=m_min_i,m_max_i
          out(m,n) = ctmp(n-m)
       end do
    end do
    
    
  end subroutine calc_H

  subroutine calc_J(self, m_min_i, m_max_i, m_min_j, m_max_j, a, out)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: m_min_i
    integer,intent(in) :: m_max_i
    integer,intent(in) :: m_min_j
    integer,intent(in) :: m_max_j
    real(8),intent(in) :: a(2)
    complex(8),intent(out) :: out(m_min_i:m_max_i,m_min_j:m_max_j)
    
    call self%calc_G(m_min_i, m_max_i, m_min_j, m_max_j, -a, out)
    
    
  end subroutine calc_J

  subroutine calc_K(self, m_min_i, m_max_i, m_min_j, m_max_j, a, out)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: m_min_i
    integer,intent(in) :: m_max_i
    integer,intent(in) :: m_min_j
    integer,intent(in) :: m_max_j
    real(8),intent(in) :: a(2)
    complex(8),intent(out) :: out(m_min_i:m_max_i,m_min_j:m_max_j)
    
    call self%calc_H(m_min_i, m_max_i, m_min_j, m_max_j, -a, out)
    
    
  end subroutine calc_K

  !> ブロック対角行列Sとベクトルの積
  subroutine product_S(self, vec, out)
    class(smatrix_shell_multiple),intent(in) :: self
    complex(8),intent(in) :: vec(self%size_all)
    complex(8),intent(out) :: out(self%size_all)

    integer :: iobj, off

    off = 0
    ! 並列化？
    do iobj=1,self%nobj

       out(off+1:off+self%smats(self%dict(iobj))%size) = &
            matmul(self%smats(self%dict(iobj))%S, vec(off+1:off+self%smats(self%dict(iobj))%size))

       off = off + self%smats(self%dict(iobj))%size
    end do
    
  end subroutine product_S

  !> ブロック対角行列Sのidv設計導関数とベクトルの積
  subroutine product_S_derivative(self, ncol, idv, vec, out)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: ncol
    integer,intent(in) :: idv
    complex(8),intent(in) :: vec(self%size_all,ncol)
    complex(8),intent(out) :: out(self%size_all,ncol)

    integer :: iobj, off
    complex(8),allocatable :: S_deri(:,:)

    ! 第dict_dv(1,idv)種散乱体S行列のdict_dv(2,idv)番設計変数に関する導関数
    allocate(S_deri(self%smats(self%dict_dv(1,idv))%m_min:self%smats(self%dict_dv(1,idv))%m_max,&
         self%smats(self%dict_dv(1,idv))%m_min:self%smats(self%dict_dv(1,idv))%m_max))
    call self%smats(self%dict_dv(1,idv))%calc_derivative(self%dict_dv(2,idv), S_deri)

    off = 0
    ! 並列化？
    do iobj=1,self%nobj

       ! 導関数に対応するブロックでなければ，零行列を掛ける
       if(self%dict(iobj) == self%dict_dv(1,idv)) then
          out(off+1:off+self%smats(self%dict(iobj))%size,:) = &
               matmul(S_deri, vec(off+1:off+self%smats(self%dict(iobj))%size,:))
       else
          out(off+1:off+self%smats(self%dict(iobj))%size,:) = zero
       end if

       off = off + self%smats(self%dict(iobj))%size
    end do
    
  end subroutine product_S_derivative

  !> idv番設計変数に関する自身の微分を計算
  subroutine calc_derivative(self, idv, S11_d, S21_d, S12_d, S22_d)
    class(smatrix_shell_multiple),intent(in) :: self
    integer,intent(in) :: idv
    complex(8),intent(out) :: S11_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S21_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S12_d(self%m_min:self%m_max,self%m_min:self%m_max)
    complex(8),intent(out) :: S22_d(self%m_min:self%m_max,self%m_min:self%m_max)

    complex(8),allocatable :: tmp(:,:), tmp2(:,:)
    integer :: m, info

    allocate(tmp(self%size_all,self%m_min:self%m_max))
    allocate(tmp2(self%size_all,self%m_min:self%m_max))

    !
    ! (I-ST)^1*S'*(T*(I-ST)^-1*S*J + J)
    !
    ! tmp = S*J
    do m=self%m_min,self%m_max
       call self%product_S(self%J(:,m), tmp(:,m))       
    end do
    ! tmp = (I-ST)^-1*S*J
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp, self%size_all, info)
    call assert(info == 0)
    ! tmp = T * (I-ST)^-1*S*J
    tmp = matmul(self%T,tmp)
    ! tmp = T * (I-ST)^-1*S*J + J
    tmp = tmp + self%J
    ! tmp2 = S'*(T*(I-ST)^-1*S*J + J)
    call self%product_S_derivative(-self%m_min+self%m_max+1, idv, tmp, tmp2)
    ! tmp2 = (I-ST)^1 * S'*(T*(I-ST)^-1*S*J + J)
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp2, self%size_all, info)
    call assert(info == 0)

    ! S11' = G * (I-ST)^1 * S'*(T*(I-ST)^-1*S*J + J)
    S11_d = matmul(self%G, tmp2)
    ! S21' = H * (I-ST)^1 * S'*(T*(I-ST)^-1*S*J + J)
    S21_d = matmul(self%H, tmp2)

    !
    ! (I-ST)^1*S'*(T*(I-ST)^-1*S*K + K)
    !
    ! tmp = S*K
    do m=self%m_min,self%m_max
       call self%product_S(self%K(:,m), tmp(:,m))       
    end do
    ! tmp = (I-ST)^-1*S*K
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp, self%size_all, info)
    call assert(info == 0)
    ! tmp = T * (I-ST)^-1*S*K
    tmp = matmul(self%T,tmp)
    ! tmp = T * (I-ST)^-1*S*K + K
    tmp = tmp + self%K
    ! tmp2 = S'*(T*(I-ST)^-1*S*K + K)
    call self%product_S_derivative(-self%m_min+self%m_max+1, idv, tmp, tmp2)
    ! tmp2 = (I-ST)^1 * S'*(T*(I-ST)^-1*S*K + K)
    call zgetrs("N", self%size_all, self%size, self%amat_inv, self%size_all, self%ipiv, tmp2, self%size_all, info)
    call assert(info == 0)

    ! S12' = G * (I-ST)^1 * S'*(T*(I-ST)^-1*S*K + K)
    S12_d = matmul(self%G, tmp2)
    ! S22' = H * (I-ST)^1 * S'*(T*(I-ST)^-1*S*K + K)
    S22_d = matmul(self%H, tmp2)
    
  end subroutine calc_derivative
end module smatrix_shell_multiple_class
