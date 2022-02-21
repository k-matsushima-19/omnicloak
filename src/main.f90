program main
  use math
  use misc
  use nlopt_driver_csec_class
  ! use smatrix_shell_homogeneous_class
  use smatrix_core_homogeneous_class
  implicit none

  

  ! w*Rmax/(2*pi) = 0.036くらいで共鳴  
  real(8),parameter :: Rmax = 0.25d0*0.99d0
  real(8),parameter :: w = 1.d0
  real(8),parameter :: Rout = 1.5d0
  real(8),parameter :: Rin  = 1.0d0
  real(8),parameter :: rho_obj = 0.27594d0
  real(8),parameter :: kappa_obj = 0.59873d0
  integer,parameter :: m_min = -20
  integer,parameter :: m_max = +20

  
  real(8),parameter :: rhos(2) = [1.3d3/1.2d0, 1.16d4/1.2d0]! [1.083d3,9.667d3]
  real(8),parameter :: dls(2) = [(6.d5)/1.4d5, (4.23d10)/1.4d5]
  real(8),parameter :: dms(2) = [(4.d4)/1.4d5, (1.49d10)/1.4d5]

  ! core
  type(smatrix_core_homogeneous) :: core

  type(nlopt_driver_csec) :: opt

  real(8) :: x_opt(2), f_opt

  ! core
  call core%new(w, 1d9*one, 1d9*one, Rin, m_min, m_max)
  
  call opt%new(w, rhos, dls, dms, m_min, m_max, Rmax, "../data/dist.txt", core)
  ! call opt%new(w, [1.2d-3,1.0d0], [0.63636363636*1d-4,1.0d0], m_min, m_max, Rmax, "dist.txt", S_obj)

  ! x= 4.768145480444083E-002  1.433872940251256E-002 付近に最適解？
  call opt%run(x_opt, f_opt)
 
  write(*,*) "# optimizer: ", x_opt
  write(*,*) "# value    : ", f_opt
  
  
end program main

