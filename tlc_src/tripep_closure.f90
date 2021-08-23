!----------------------------------------------------------------------
! Copyright (C) 2003, 2004, 2005
!      Chaok Seok, Evangelos Coutsias, Matthew Jacobson, and Ken Dill.
!      UCSF, Univeristy of New Mexico, Seoul National University
                                                                                                                         
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!----------------------------------------------------------------------------
! written by Chaok Seok, 2003
!----------------------------------------------------------------------------
!*************** Tripeptide Loop Closure Algorithm *****************
! files to be compiled with:
!    tripeptide_closure.f90
!    sturm.c
! 
!*******************************************************************
! subroutine  initialize_loop_closure(b_len, b_ang, t_ang)
!*******************************************************************
! connectivity of atoms:
!   N1-A1-C1-N2-A2-C2-N3-A3-C3
!
! input:
!
!  * b_len(1:6): bond lengths (A1-C1, C1-N2, ..., N3-A3)
!  * b_ang(1:7): bond angles (N1-A1-C1, A1-C1-N2, ..., N3-A3-C3)
!  * t_ang(1:2): torsion angles (A1-C1-N2-A2, A2-C2-N3-A3)
!*******************************************************************
!
!*******************************************************************
! subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
!     r_soln_n, r_soln_a, r_soln_c, n_soln)
!*******************************************************************
! input: 
!  * r_n1(3), r_a1(3), r_a3(3), r_c3(3): 
!       Cartesian coordinates of N and CA atoms of the first residue and
!        CA and C atoms of the last (third) residue.
! output:
!  * n_soln: number of alternative loop closure solutions.
!  * r_soln_n(3,3,8), r_soln_a(3,3,8), r_soln_c(3,3,8): 
!       Cartesian coordinates of loop closure solutions. 
!       first dimension: x, y, z component
!       second dim: residue number
!       third dim: solution number
!*******************************************************************
!----------------------------------------------------------------------------
MODULE tripep_closure
!----------------------------------------------------------------------------
  real(KIND=8), parameter :: pi = 3.141592653589793238462643383279502884197d0
  real(KIND=8), parameter :: two_pi=2.0d0*pi, deg2rad = pi/180.0d0, rad2deg = 180.0d0/pi
  integer, parameter :: max_soln = 16
  integer, parameter :: deg_pol = 16
  integer, parameter :: print_level = 0
  ! parameters for tripeptide loop (including bond lengths & angles)
  real(KIND=8) :: len0(6), b_ang0(7), t_ang0(2)
  real(KIND=8) :: aa13_min_sqr, aa13_max_sqr
  real(KIND=8) :: delta(0:3), xi(3), eta(3), alpha(3), theta(3)
  real(KIND=8) :: cos_alpha(3), sin_alpha(3), cos_theta(3), sin_theta(3)
  real(KIND=8) :: cos_delta(0:3), sin_delta(0:3)
  real(KIND=8) :: cos_xi(3), cos_eta(3), sin_xi(3), sin_eta(3)
  real(KIND=8) :: r_a1a3(3), r_a1n1(3), r_a3c3(3)
  real(KIND=8) :: b_a1a3(3), b_a1n1(3), b_a3c3(3)
  real(KIND=8) :: len_na(3), len_ac(3), len_aa(3)
  ! used for polynomial coefficients
  real(KIND=8) :: C0(0:2,3), C1(0:2,3), C2(0:2,3)
  real(KIND=8) :: Q(0:16,0:4), R(0:16,0:2)
CONTAINS
!-----------------------------------------------------------------------
subroutine solv_3pep_poly(r_n1, r_a1, r_a3, r_c3, &
     r_soln_n, r_soln_a, r_soln_c, n_soln)
  implicit none
  real(KIND=8), intent(in) :: r_n1(3), r_a1(3), r_a3(3), r_c3(3)
  integer, intent(out) :: n_soln
  real(KIND=8), intent(out) :: r_soln_n(:,:,:), r_soln_a(:,:,:), r_soln_c(:,:,:)
  real(KIND=8) :: poly_coeff(0:deg_pol), roots(max_soln)

  call get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3)

  if (n_soln == 0) then
     return
  end if

  call get_poly_coeff(poly_coeff)

  call solve_sturm(deg_pol, n_soln, poly_coeff, roots)

  if (n_soln == 0) then
!     print*, 'return 2'
     return
  end if

  call coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c)

end subroutine solv_3pep_poly
!-----------------------------------------------------------------------
subroutine initialize_loop_closure(b_len, b_ang, t_ang)
!-----------------------------------------------------------------------
! Input angles for the given bond lengths and angles
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), intent(in) :: b_len(6), b_ang(7), t_ang(2)
  real(KIND=8) :: len1, len2, a_min, a_max
  real(KIND=8), dimension(3) :: axis, rr_a1, rr_c1, rr_n2, rr_a2, rr_n2a2_ref, rr_c1a1
  real(KIND=8), dimension(3) :: rr_a1a2, dr, bb_c1a1, bb_a1a2, bb_a2n2
  real(KIND=8), dimension(3) :: r_a1, r_n1
  real(KIND=8) :: p(4), Us(3,3)
  real(KIND=8), parameter :: tol_secant = 1.0d-15
  integer, parameter :: max_iter_sturm = 100, max_iter_secant = 20
  integer :: i

  call initialize_sturm(tol_secant, max_iter_sturm, max_iter_secant)

  len0(1:6) = b_len(1:6)
  b_ang0(1:7) = b_ang(1:7)
  t_ang0(1:2) = t_ang(1:2)

  rr_c1(1:3) = 0.0d0
  axis(1:3) = (/ 1.0d0, 0.0d0, 0.0d0 /)

  do i = 0, 1
     rr_a1(1:3) = (/ cos(b_ang0(3*i+2))*len0(3*i+1), sin(b_ang0(3*i+2))*len0(3*i+1), 0.0d0 /)
     rr_n2(1:3) = (/ len0(3*i+2), 0.0d0, 0.0d0 /)
     rr_c1a1(:) = rr_a1(:) - rr_c1(:)
     rr_n2a2_ref(1:3) = (/ -cos(b_ang0(3*i+3))*len0(3*i+3), sin(b_ang0(3*i+3))*len0(3*i+3), 0.0d0 /)
     call quaternion(axis, t_ang0(i+1)*0.25d0, p)
     call rotation_matrix(p, Us)
     rr_a2(:) =  matmul(Us, rr_n2a2_ref) + rr_n2(:)
     rr_a1a2(:) = rr_a2(:) - rr_a1(:)
     dr(:) = rr_a1a2(:)
     len2 = dot_product(dr, dr)
     len1 = sqrt(len2)
     ! len_aa
     len_aa(i+2) = len1
     bb_c1a1(:) = rr_c1a1(:)/len0(3*i+1)
     bb_a1a2(:) = rr_a1a2(:)/len1
     bb_a2n2(:) = (rr_n2(:) - rr_a2(:))/len0(3*i+3)
     ! xi
     call calc_bnd_ang(-bb_a1a2, bb_a2n2, xi(i+2))
     ! eta
     call calc_bnd_ang(bb_a1a2, -bb_c1a1, eta(i+1))
     ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)
     call calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2, delta(i+1))
     delta(i+1) = pi - delta(i+1)
  end do

  a_min = b_ang(4) - (xi(2) + eta(2))
  a_max = min(b_ang(4) + (xi(2) + eta(2)), pi)

  ! min/max of base length
!  print*, 'len1, len3=', len_aa(2:3)
!  print*, 'a_min, a_max=', a_min*rad2deg, a_max*rad2deg
  aa13_min_sqr = len_aa(2)**2 + len_aa(3)**2 - 2.0d0*len_aa(2)*len_aa(3)*cos(a_min)
  aa13_max_sqr = len_aa(2)**2 + len_aa(3)**2 - 2.0d0*len_aa(2)*len_aa(3)*cos(a_max)
!  print*, 'aa13_min_sqr,aa13_max_sqr', aa13_min_sqr,aa13_max_sqr

end subroutine initialize_loop_closure
!-----------------------------------------------------------------------
subroutine get_input_angles(n_soln, r_n1, r_a1, r_a3, r_c3)
!-----------------------------------------------------------------------
! Input angles and vectors (later used in coordinates)
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), intent(in) :: r_n1(:), r_a1(:), r_a3(:), r_c3(:)
  integer, intent(out) :: n_soln
  real(KIND=8) :: dr_sqr
  integer :: i
  character(len=2) :: cone_type
!-----------------------------------------------------------------------

  n_soln = max_soln

  ! vertual bond
  r_a1a3(:) = r_a3(:) - r_a1(:)
  dr_sqr = dot_product(r_a1a3,r_a1a3)
  len_aa(1) = sqrt(dr_sqr)

  if (dr_sqr < aa13_min_sqr .or. dr_sqr > aa13_max_sqr) then
     n_soln = 0
!     print*, 'return 0'
!     print*, sqrt(dr_sqr), sqrt(aa13_min_sqr), sqrt(aa13_max_sqr)
     return
  end if

  ! bond lengths
  r_a1n1(:) = r_n1(:) - r_a1(:)
  len_na(1) = sqrt(dot_product(r_a1n1,r_a1n1))
  len_na(2) = len0(3)
  len_na(3) = len0(6)
  r_a3c3(:) = r_c3(:) - r_a3(:)
  len_ac(1) = len0(1)
  len_ac(2) = len0(4)
  len_ac(3) = sqrt(dot_product(r_a3c3,r_a3c3))

  ! unit vectors
  b_a1n1(:) = r_a1n1(:)/len_na(1)
  b_a3c3(:) = r_a3c3(:)/len_ac(3)
  b_a1a3(:) = r_a1a3(:)/len_aa(1)

  ! delta(3): dih of N(1)CA(1)CA(3)C(3)
  call calc_dih_ang(-b_a1n1, b_a1a3, b_a3c3, delta(3))
  delta(0) = delta(3)

  ! xi(1)
  call calc_bnd_ang(-b_a1a3, b_a1n1, xi(1))

  ! eta(3)
  call calc_bnd_ang(b_a1a3, b_a3c3, eta(3))

  do i = 1, 3
     cos_delta(i) = cos(delta(i))
     sin_delta(i) = sin(delta(i))
     cos_xi(i) = cos(xi(i))
     sin_xi(i) = sin(xi(i))
     cos_eta(i) = cos(eta(i))
     sin_eta(i) = sin(eta(i))
  end do
  cos_delta(0) = cos_delta(3)
  sin_delta(0) = sin_delta(3)

  ! theta (N, CA, C) bond angle
  theta(1) = b_ang0(1)
  theta(2) = b_ang0(4)
  theta(3) = b_ang0(7)
  do i = 1, 3
     cos_theta(i) = cos(theta(i))
  end do

  ! alpha
  cos_alpha(1) = -(len_aa(1)**2 + len_aa(2)**2 - len_aa(3)**2)/(2.0d0*len_aa(1)*len_aa(2))
  alpha(1) = acos(cos_alpha(1))
  sin_alpha(1) = sin(alpha(1))
  cos_alpha(2) = (len_aa(2)**2 + len_aa(3)**2 - len_aa(1)**2)/(2.0d0*len_aa(2)*len_aa(3))
  alpha(2) = acos(cos_alpha(2))
  sin_alpha(2) = sin(alpha(2))
  alpha(3) = pi - alpha(1) + alpha(2)
  cos_alpha(3) = cos(alpha(3))
  sin_alpha(3) = sin(alpha(3))

  if (print_level > 0) then
     write(*,'(a,3f9.4)') 'xi = ', xi(1:3)*rad2deg
     write(*,'(a,3f9.4)') 'eta = ', eta(1:3)*rad2deg
     write(*,'(a,3f9.4)') 'delta = ', delta(1:3)*rad2deg
     write(*,'(a,3f9.4)') 'theta = ', theta(1:3)*rad2deg
     write(*,'(a,3f9.4)') 'alpha = ', alpha(1:3)*rad2deg
  end if

  ! check for existence of soln
  do i = 1, 3
     call test_two_cone_existence_soln(theta(i), xi(i), eta(i), alpha(i), &
          n_soln, cone_type)
     if (n_soln == 0) then
        ! print*, 'return 1', i
        return
     end if
  end do

end subroutine get_input_angles
!-----------------------------------------------------------------------
subroutine test_two_cone_existence_soln(tt, kx, et, ap, n_soln, cone_type)
  implicit none
  real(KIND=8), intent(in) :: tt, kx, et, ap
  integer, intent(out) :: n_soln
  character(len=2), intent(out) :: cone_type
  character(len=2) :: case_type
  real(KIND=8) :: at, ex, abs_at, ap1, kx1, et1
  real(KIND=8) :: cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1, cos_ea2, cos_xa1, cos_xa2
  logical :: s1, s2, t1, t2, complicated = .false.
  real(KIND=8), parameter :: half_pi = 0.5d0*pi

  n_soln = max_soln

  ap1 = ap
  kx1 = kx
  et1 = et

  at = ap1 - tt
  ex = kx1 + et1
  abs_at = abs(at)

  ! case of no soln
  if (abs_at > ex) then
     n_soln = 0
     return
  end if

  if (complicated) then
     ! find type of intersection
     cos_tx1 = cos(tt+kx1)
     cos_tx2 = cos(tt-kx1)
     cos_te1 = cos(tt+et1)
     cos_te2 = cos(tt-et1)
     cos_ea1 = cos(et1+ap1)
     cos_ea2 = cos(et1-ap1)
     cos_xa1 = cos(kx1+ap1)
     cos_xa2 = cos(kx1-ap1)
     s1 = .false.; s2 = .false.; t1 = .false.; t2 = .false.
     if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) <= 0.0d0) s1 = .true.
     if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) <= 0.0d0) s2 = .true.
     if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) <= 0.0d0) t1 = .true.
     if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) <= 0.0d0) t2 = .true.

  end if

end subroutine test_two_cone_existence_soln
!-----------------------------------------------------------------------
subroutine get_poly_coeff(poly_coeff)
  implicit none
  real(KIND=8), intent(out) :: poly_coeff(0:deg_pol)
  integer :: i, j
  real(KIND=8) :: A0, A1, A2, A3, A4, A21, A22, A31, A32, A41, A42
  real(KIND=8) :: B0(3), B1(3), B2(3), B3(3), B4(3), B5(3), B6(3), B7(3), B8(3)
  real(KIND=8), dimension(0:4,0:4) :: u11, u12, u13, u31, u32, u33
  real(KIND=8), dimension(0:4,0:4) :: um1, um2, um3, um4, um5, um6, q_tmp
  integer, dimension(2) :: p1, p3, p_um1, p_um2, p_um3, p_um4, p_um5, p_um6, p_Q
  integer :: p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, &
       p_f8, p_f9, p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, &
       p_f18, p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26, p_f27
  integer :: p_final
  real(KIND=8), dimension(0:16) :: f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, &
       f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24, f25, f26, f27

  ! A0, B0
  do i = 1, 3
     A0 = cos_alpha(i)*cos_xi(i)*cos_eta(i) - cos_theta(i)
     A1 = -sin_alpha(i)*cos_xi(i)*sin_eta(i)
     A2 = sin_alpha(i)*sin_xi(i)*cos_eta(i)
     A3 = sin_xi(i)*sin_eta(i)
     A4 = A3*cos_alpha(i)
     j = i - 1
     A21 = A2*cos_delta(j)
     A22 = A2*sin_delta(j)
     A31 = A3*cos_delta(j)
     A32 = A3*sin_delta(j)
     A41 = A4*cos_delta(j)
     A42 = A4*sin_delta(j)
     B0(i) = A0 + A22 + A31
     B1(i) = 2.0d0*(A1 + A42)
     B2(i) = 2.0d0*(A32 - A21)
     B3(i) = -4.0d0*A41
     B4(i) = A0 + A22 - A31
     B5(i) = A0 - A22 - A31
     B6(i) = -2.0d0*(A21 + A32)
     B7(i) = 2.0d0*(A1 - A42)
     B8(i) = A0 - A22 + A31
  end do

  ! C0i
  i = 1
  C0(0:2,i) = (/ B0(i), B2(i), B5(i) /)
  C1(0:2,i) = (/ B1(i), B3(i), B7(i) /)
  C2(0:2,i) = (/ B4(i), B6(i), B8(i) /)
  do i = 2, 3
     C0(0:2,i) = (/ B0(i), B1(i), B4(i) /)
     C1(0:2,i) = (/ B2(i), B3(i), B6(i) /)
     C2(0:2,i) = (/ B5(i), B7(i), B8(i) /)
  end do

  ! first determinant
  do i = 0, 2
     u11(i,0) = C0(i,1)
     u12(i,0) = C1(i,1)
     u13(i,0) = C2(i,1)
     u31(0,i) = C0(i,2)
     u32(0,i) = C1(i,2)
     u33(0,i) = C2(i,2)
  end do

  p1(1:2) = (/ 2, 0 /)
  p3(1:2) = (/ 0, 2 /)

  call poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1)
  call poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2)
  call poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3)
  call poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4)
  call poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5)
  call poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6)
  call poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q)
  Q(0:4,0:4) = q_tmp(0:4,0:4)

  ! second determinant
  R(:,:) = 0.0d0
  R(0:2,0) = C0(0:2,3)
  R(0:2,1) = C1(0:2,3)
  R(0:2,2) = C2(0:2,3)
  p2 = 2
  p4 = 4

  call poly_mul_sub1(R(:,1), R(:,1), R(:,0), R(:,2), p2, p2, p2, p2, f1, p_f1)
  call poly_mul1(R(:,1), R(:,2), p2, p2, f2, p_f2)
  call poly_mul_sub1(R(:,1), f1, R(:,0), f2, p2, p_f1, p2, p_f2, f3, p_f3)
  call poly_mul1(R(:,2), f1, p2, p_f1, f4, p_f4)
  call poly_mul_sub1(R(:,1), f3, R(:,0), f4, p2, p_f3, p2, p_f4, f5, p_f5)

  call poly_mul_sub1(Q(:,1), R(:,1), Q(:,0), R(:,2), p4, p2, p4, p2, f6, p_f6)
  call poly_mul_sub1(Q(:,2), f1, R(:,2), f6, p4, p_f1, p2, p_f6, f7, p_f7)
  call poly_mul_sub1(Q(:,3), f3, R(:,2), f7, p4, p_f3, p2, p_f7, f8, p_f8)
  call poly_mul_sub1(Q(:,4), f5, R(:,2), f8, p4, p_f5, p2, p_f8, f9, p_f9)

  call poly_mul_sub1(Q(:,3), R(:,1), Q(:,4), R(:,0), p4, p2, p4, p2, f10, p_f10)
  call poly_mul_sub1(Q(:,2), f1, R(:,0), f10, p4, p_f1, p2, p_f10, f11, p_f11)
  call poly_mul_sub1(Q(:,1), f3, R(:,0), f11, p4, p_f3, p2, p_f11, f12, p_f12)

  call poly_mul_sub1(Q(:,2), R(:,1), Q(:,1), R(:,2), p4, p2, p4, p2, f13, p_f13)
  call poly_mul_sub1(Q(:,3), f1, R(:,2), f13, p4, p_f1, p2, p_f13, f14, p_f14)
  call poly_mul_sub1(Q(:,3), R(:,1), Q(:,2), R(:,2), p4, p2, p4, p2, f15, p_f15)
  call poly_mul_sub1(Q(:,4), f1, R(:,2), f15, p4, p_f1, p2, p_f15, f16, p_f16)
  call poly_mul_sub1(Q(:,1), f14, Q(:,0), f16, p4, p_f14, p4, p_f16, f17, p_f17)

  call poly_mul_sub1(Q(:,2), R(:,2), Q(:,3), R(:,1), p4, p2, p4, p2, f18, p_f18)
  call poly_mul_sub1(Q(:,1), R(:,2), Q(:,3), R(:,0), p4, p2, p4, p2, f19, p_f19)
  call poly_mul_sub1(Q(:,3), f19, Q(:,2), f18, p4, p_f19, p4, p_f18, f20, p_f20)
  call poly_mul_sub1(Q(:,1), R(:,1), Q(:,2), R(:,0), p4, p2, p4, p2, f21, p_f21)
  call poly_mul1(Q(:,4), f21, p4, p_f21, f22, p_f22)
  call poly_sub1(f20, f22, p_f20, p_f22, f23, p_f23)
  call poly_mul1(R(:,0), f23, p2, p_f23, f24, p_f24)
  call poly_sub1(f17, f24, p_f17, p_f24, f25, p_f25)
  call poly_mul_sub1(Q(:,4), f12, R(:,2), f25, p4, p_f12, p2, p_f25, f26, p_f26)
  call poly_mul_sub1(Q(:,0), f9, R(:,0), f26, p4, p_f9, p2, p_f26, poly_coeff, p_final)

  if (p_final /= 16) then
     print*, 'Error. Degree of polynomial is not 16!'
     stop
  end if

  if (poly_coeff(16) < 0.0d0) then
     poly_coeff(0:16) = -poly_coeff(0:16)
  end if

  if (print_level > 0) then
     ! print*, 'poly_coeff'
     do i = 0, 16
        write(*,"(i5,e15.6)") i, poly_coeff(i)
     end do
  end if

end subroutine get_poly_coeff
!----------------------------------------------------------------------------
subroutine poly_mul_sub2(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
  implicit none
  real(KIND=8), dimension(0:4,0:4), intent(in) :: u1, u2, u3, u4
  integer, dimension(2), intent(in) :: p1, p2, p3, p4
  real(KIND=8), dimension(0:4,0:4), intent(out) :: u5
  integer, dimension(2), intent(out) :: p5
  real(KIND=8), dimension(0:4,0:4) :: d1, d2
  integer, dimension(2) :: pd1, pd2

  call poly_mul2(u1, u2, p1, p2, d1, pd1)
  call poly_mul2(u3, u4, p3, p4, d2, pd2)
  call poly_sub2(d1, d2, pd1, pd2, u5, p5)

end subroutine poly_mul_sub2
!----------------------------------------------------------------------------
subroutine poly_mul2(u1, u2, p1, p2, u3, p3)
  implicit none
  real(KIND=8), dimension(0:4,0:4), intent(in) :: u1, u2
  integer, dimension(2), intent(in) :: p1, p2
  real(KIND=8), dimension(0:4,0:4), intent(out) :: u3
  integer, intent(out) :: p3(2)
  integer :: i1, j1, i2, j2, i3, j3, p11, p12, p21, p22
  real(KIND=8) :: u1ij

  p3(:) = p1(:) + p2(:)
  u3(:,:) = 0.0d0

  p11 = p1(1)
  p12 = p1(2)
  p21 = p2(1)
  p22 = p2(2)

  do i1 = 0, p12
     do j1 = 0, p11
        u1ij = u1(j1,i1)
        do i2 = 0, p22
           i3 = i1 + i2
           do j2 = 0, p21
              j3 = j1 + j2
              u3(j3,i3) = u3(j3,i3) + u1ij*u2(j2,i2)
           end do
        end do
     end do
  end do

end subroutine poly_mul2
!----------------------------------------------------------------------------
subroutine poly_sub2(u1, u2, p1, p2, u3, p3)
  implicit none
  real(KIND=8), dimension(0:4,0:4), intent(in) :: u1, u2
  integer, intent(in) :: p1(2), p2(2)
  real(KIND=8), dimension(0:4,0:4), intent(out) :: u3
  integer, intent(out) :: p3(2)
  integer :: i, j, p11, p12, p21, p22
  logical :: i1_ok, i2_ok

  p11 = p1(1)
  p12 = p1(2)
  p21 = p2(1)
  p22 = p2(2)
  p3(1) = max(p11,p21)
  p3(2) = max(p12,p22)

  do i = 0, p3(2)
     i1_ok = (i > p12)
     i2_ok = (i > p22)
     do j = 0, p3(1)
        if (i2_ok .or. (j > p21)) then
           u3(j,i) = u1(j,i)
        else if (i1_ok .or. (j > p11)) then
           u3(j,i) = -u2(j,i)
        else
           u3(j,i) = u1(j,i) - u2(j,i)
        end if
     end do
  end do

end subroutine poly_sub2
!----------------------------------------------------------------------------
subroutine poly_mul_sub1(u1, u2, u3, u4, p1, p2, p3, p4, u5, p5)
  implicit none
  real(KIND=8), dimension(0:16), intent(in) :: u1, u2, u3, u4
  integer, intent(in) :: p1, p2, p3, p4
  real(KIND=8), dimension(0:16), intent(out) :: u5
  integer, intent(out) :: p5
  real(KIND=8), dimension(0:16) :: d1, d2
  integer :: pd1, pd2

  call poly_mul1(u1, u2, p1, p2, d1, pd1)
  call poly_mul1(u3, u4, p3, p4, d2, pd2)
  call poly_sub1(d1, d2, pd1, pd2, u5, p5)

end subroutine poly_mul_sub1
!----------------------------------------------------------------------------
subroutine poly_mul1(u1, u2, p1, p2, u3, p3)
  implicit none
  real(KIND=8), dimension(0:16), intent(in) :: u1, u2
  integer, intent(in) :: p1, p2
  real(KIND=8), dimension(0:16), intent(out) :: u3
  integer, intent(out) :: p3
  integer :: i1, i2, i3
  real(KIND=8) :: u1i

  p3 = p1 + p2
  u3(:) = 0.0d0

  do i1 = 0, p1
     u1i = u1(i1)
     do i2 = 0, p2
        i3 = i1 + i2
        u3(i3) = u3(i3) + u1i*u2(i2)
     end do
  end do

end subroutine poly_mul1
!----------------------------------------------------------------------------
subroutine poly_sub1(u1, u2, p1, p2, u3, p3)
  implicit none
  real(KIND=8), dimension(0:16), intent(in) :: u1, u2
  integer, intent(in) :: p1, p2
  real(KIND=8), dimension(0:16), intent(out) :: u3
  integer, intent(out) :: p3
  integer :: i

  p3 = max(p1, p2)

  do i = 0, p3
     if (i > p2) then
        u3(i) = u1(i)
     else if (i > p1) then
        u3(i) = -u2(i)
     else
        u3(i) = u1(i) - u2(i)
     end if
  end do

end subroutine poly_sub1
!----------------------------------------------------------------------------
subroutine coord_from_poly_roots(n_soln, roots, r_n1, r_a1, r_a3, r_c3, r_soln_n, r_soln_a, r_soln_c)
  implicit none
  integer, intent(in) :: n_soln
  real(KIND=8), intent(in) :: r_n1(:), r_a1(:), r_a3(:), r_c3(:), roots(n_soln)
  real(KIND=8), intent(out) :: r_soln_n(:,:,:), r_soln_a(:,:,:), r_soln_c(:,:,:)
  real(KIND=8) :: ex(3), ey(3), ez(3), b_a1a2(3), b_a3a2(3), r_tmp(3)
  real(KIND=8) :: p_s(3,3), s1(3,3), s2(3,3), p_t(3,3), t1(3,3), t2(3,3)
  real(KIND=8) :: p_s_c(3,3), s1_s(3,3), s2_s(3,3), p_t_c(3,3), t1_s(3,3), t2_s(3,3)
  real(KIND=8) :: angle, sig1_init, half_tan(3)
  real(KIND=8) :: cos_tau(0:3), sin_tau(0:3), cos_sig(3), sin_sig(3), ht, tmp, sig1
  real(KIND=8) :: r_s(3), r_t(3), r0(3), r_n(3,3), r_a(3,3), r_c(3,3), p(4), Us(3,3)
  integer :: i_soln, i, j
  real(KIND=8) :: a1c1, c1n2, n2a2, a2c2, c2n3, n3a3, a1a2, a2a3
  real(KIND=8) :: rr_a1c1(3), rr_c1n2(3), rr_n2a2(3), rr_a2c2(3), rr_c2n3(3), rr_n3a3(3), rr_a1a2(3), rr_a2a3(3)
  real(KIND=8) :: a3a1a2, a2a3a1, n1a1c1, n2a2c2, n3a3c3, a1c1n2a2, a2c2n3a3

  if (n_soln == 0) return

  ! Define body frame (ex, ey, ez)
  ex(:) = b_a1a3(:)
  call cross(r_a1n1, ex, ez)
  ez(:) = ez(:)/sqrt(dot_product(ez,ez))
  call cross(ez, ex, ey)
  ! vertual bond vectors in the reference plane
  b_a1a2(:) = -cos_alpha(1)*ex(:) + sin_alpha(1)*ey(:)
  b_a3a2(:) = cos_alpha(3)*ex(:) + sin_alpha(3)*ey(:)
  !! Define cone coordinates for each angle joint.
  ! (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems
  ! residue 1
  p_s(:,1) = -ex(:)
  s1(:,1)  = ez(:)  ! (p_s)X(p_t)/||(p_s)X(p_t)||
  s2(:,1)  = ey(:)  ! p_s X s1
  p_t(:,1) = b_a1a2(:)
  t1(:,1)  = ez(:)  ! s1
  t2(:,1)  = sin_alpha(1)*ex(:) + cos_alpha(1)*ey(:) ! p_t X t1
  ! residue 2
  p_s(:,2) = -b_a1a2(:)
  s1(:,2)  = -ez(:)
  s2(:,2)  = t2(:,1)  ! sina1*ex(:) + cosa1*ey(:)
  p_t(:,2) = -b_a3a2(:)
  t1(:,2)  = -ez(:)
  t2(:,2)  = sin_alpha(3)*ex(:) - cos_alpha(3)*ey(:)
  ! residue 3
  p_s(:,3) = b_a3a2(:)
  s2(:,3)  = t2(:,2)   ! sina3*ex(:) + cosa3*ey(:)
  s1(:,3)  = ez(:)
  p_t(:,3) = ex(:)
  t1(:,3) =  ez(:)
  t2(:,3) = -ey(:)
  ! scale vectors
  do i = 1, 3
     p_s_c(:,i) = p_s(:,i)*cos_xi(i)
     s1_s(:,i)  =  s1(:,i)*sin_xi(i)
     s2_s(:,i)  =  s2(:,i)*sin_xi(i)
     p_t_c(:,i) = p_t(:,i)*cos_eta(i)
     t1_s(:,i)  =  t1(:,i)*sin_eta(i)
     t2_s(:,i)  =  t2(:,i)*sin_eta(i)
  end do

  ! initial sig(1)
  r_tmp(:) = (r_a1n1(:)/len_na(1) - p_s_c(:,1))/sin_xi(1)
  call calc_bnd_ang(s1(:,1), r_tmp, angle)
  sig1_init = sign(angle, dot_product(r_tmp(:),s2(:,1)))

  ! CA
  r_a(:,1) = r_a1(:)
  r_a(:,2) = r_a1(:) + len_aa(2)*b_a1a2(:)
  r_a(:,3) = r_a3(:)
  r0(:) = r_a1(:)

  do i_soln = 1, n_soln
     half_tan(3) = roots(i_soln)
     half_tan(2) = calc_t2(half_tan(3))
     half_tan(1) = calc_t1(half_tan(3), half_tan(2))
     do i = 1, 3
        ht = half_tan(i)
        tmp = 1.0d0 + ht*ht
        cos_tau(i) = (1.0d0 - ht*ht)/tmp
        sin_tau(i) = 2.0d0*ht/tmp
     end do
     cos_tau(0) = cos_tau(3)
     sin_tau(0) = sin_tau(3)
     do i = 1, 3
        j = i - 1
        cos_sig(i) = cos_delta(j)*cos_tau(j) + sin_delta(j)*sin_tau(j)
        sin_sig(i) = sin_delta(j)*cos_tau(j) - cos_delta(j)*sin_tau(j)
     end do
     do i = 1, 3
        r_s(:) = p_s_c(:,i) + cos_sig(i)*s1_s(:,i) + sin_sig(i)*s2_s(:,i)
        r_t(:) = p_t_c(:,i) + cos_tau(i)*t1_s(:,i) + sin_tau(i)*t2_s(:,i)
        r_n(:,i) = r_s(:)*len_na(i) + r_a(:,i)
        r_c(:,i) = r_t(:)*len_ac(i) + r_a(:,i)
     end do

     ! rotate back atoms by -(sig(1) - sig1_init) around -ex
     sig1 = atan2(sin_sig(1), cos_sig(1))
     call quaternion(-ex, -(sig1 - sig1_init)*0.25d0, p)
     call rotation_matrix(p, Us)

     r_soln_n(:,1,i_soln) = r_n1(:)
     r_soln_a(:,1,i_soln) = r_a1(:)
     r_soln_c(:,1,i_soln) = matmul(Us, r_c(:,1) - r0(:)) + r0(:)
     r_soln_n(:,2,i_soln) = matmul(Us, r_n(:,2) - r0(:)) + r0(:)
     r_soln_a(:,2,i_soln) = matmul(Us, r_a(:,2) - r0(:)) + r0(:)
     r_soln_c(:,2,i_soln) = matmul(Us, r_c(:,2) - r0(:)) + r0(:)
     r_soln_n(:,3,i_soln) = matmul(Us, r_n(:,3) - r0(:)) + r0(:)
     r_soln_a(:,3,i_soln) = r_a3(:)
     r_soln_c(:,3,i_soln) = r_c3(:)


     if (print_level > 0) then
        ! print*, 'roots: t0, t2, t1', i_soln
        write(*,"(3f15.6)") half_tan(3), half_tan(2), half_tan(1)

        rr_a1c1(:) = r_soln_c(:,1,i_soln) - r_soln_a(:,1,i_soln)
        rr_c1n2(:) = r_soln_n(:,2,i_soln) - r_soln_c(:,1,i_soln)
        rr_n2a2(:) = r_soln_a(:,2,i_soln) - r_soln_n(:,2,i_soln)
        rr_a2c2(:) = r_soln_c(:,2,i_soln) - r_soln_a(:,2,i_soln)
        rr_c2n3(:) = r_soln_n(:,3,i_soln) - r_soln_c(:,2,i_soln)
        rr_n3a3(:) = r_soln_a(:,3,i_soln) - r_soln_n(:,3,i_soln)
        rr_a1a2(:) = r_soln_a(:,2,i_soln) - r_soln_a(:,1,i_soln)
        rr_a2a3(:) = r_soln_a(:,3,i_soln) - r_soln_a(:,2,i_soln)

        a1c1 = sqrt(dot_product(rr_a1c1, rr_a1c1))
        c1n2 = sqrt(dot_product(rr_c1n2, rr_c1n2))
        n2a2 = sqrt(dot_product(rr_n2a2, rr_n2a2))
        a2c2 = sqrt(dot_product(rr_a2c2, rr_a2c2))
        c2n3 = sqrt(dot_product(rr_c2n3, rr_c2n3))
        n3a3 = sqrt(dot_product(rr_n3a3, rr_n3a3))
        a1a2 = sqrt(dot_product(rr_a1a2, rr_a1a2))
        a2a3 = sqrt(dot_product(rr_a2a3, rr_a2a3))

        write(*,"('na: n2a2, n3a3 = ', 4f9.3)") len0(3), n2a2, len0(6), n3a3
        write(*,"('ac: a1c1, a2c2 = ', 4f9.3)") len0(1), a1c1, len0(4), a2c2
        write(*,"('cn: c1n2, c2n3 = ', 4f9.3)") len0(2), c1n2, len0(5), c2n3
        write(*,"('aa: a1a2, a2a3 = ', 4f9.3)") len_aa(2), a1a2, len_aa(3), a2a3

        call calc_bnd_ang(b_a1a3, rr_a1a2/a1a2, a3a1a2)
        call calc_bnd_ang(rr_a2a3/a2a3, b_a1a3, a2a3a1)
        write(*,"('alpha1, alpha3 = ', 2f9.3)") (pi-a3a1a2)*RAD2DEG, (pi-a2a3a1)*RAD2DEG
        call calc_bnd_ang(b_a1n1, rr_a1c1/a1c1, n1a1c1)
        call calc_bnd_ang(b_a3c3, -rr_n3a3/n3a3, n3a3c3)
        call calc_bnd_ang(rr_a2c2/a2c2, -rr_n2a2/n2a2, n2a2c2)
        write(*,"('ang_nac = ', 2f9.3)") b_ang0(1)*RAD2DEG, n1a1c1*RAD2DEG
        write(*,"('ang_nac = ', 2f9.3)") b_ang0(4)*RAD2DEG, n2a2c2*RAD2DEG
        write(*,"('ang_nac = ', 2f9.3)") b_ang0(7)*RAD2DEG, n3a3c3*RAD2DEG

        call calc_dih_ang(rr_a1c1/a1c1, rr_c1n2/c1n2, rr_n2a2/n2a2, a1c1n2a2)
        call calc_dih_ang(rr_a2c2/a2c2, rr_c2n3/c2n3, rr_n3a3/n3a3, a2c2n3a3)
        write(*,"('t_ang1 = ', 2f9.3)") t_ang0(1)*RAD2DEG, a1c1n2a2*RAD2DEG
        write(*,"('t_ang2 = ', 2f9.3)") t_ang0(2)*RAD2DEG, a2c2n3a3*RAD2DEG
     end if

  end do

end subroutine coord_from_poly_roots
!-----------------------------------------------------------------------
function calc_t2(t0)
  implicit none
  real(KIND=8), intent(in) :: t0
  real(KIND=8) :: calc_t2
  real(KIND=8) :: B0, B1, B2, A0, A1, A2, A3, A4, B2_2, B2_3
  real(KIND=8) :: K0, K1, K2, K3, t0_2, t0_3, t0_4

  t0_2 = t0*t0
  t0_3 = t0_2*t0
  t0_4 = t0_3*t0

  A0 = Q(0,0) + Q(1,0)*t0 + Q(2,0)*t0_2 + Q(3,0)*t0_3 + Q(4,0)*t0_4
  A1 = Q(0,1) + Q(1,1)*t0 + Q(2,1)*t0_2 + Q(3,1)*t0_3 + Q(4,1)*t0_4
  A2 = Q(0,2) + Q(1,2)*t0 + Q(2,2)*t0_2 + Q(3,2)*t0_3 + Q(4,2)*t0_4
  A3 = Q(0,3) + Q(1,3)*t0 + Q(2,3)*t0_2 + Q(3,3)*t0_3 + Q(4,3)*t0_4
  A4 = Q(0,4) + Q(1,4)*t0 + Q(2,4)*t0_2 + Q(3,4)*t0_3 + Q(4,4)*t0_4

  B0 = R(0,0) + R(1,0)*t0 + R(2,0)*t0_2
  B1 = R(0,1) + R(1,1)*t0 + R(2,1)*t0_2
  B2 = R(0,2) + R(1,2)*t0 + R(2,2)*t0_2

  B2_2 = B2*B2
  B2_3 = B2_2*B2

  K0 = A2*B2 - A4*B0
  K1 = A3*B2 - A4*B1
  K2 = A1*B2_2 - K1*B0
  K3 = K0*B2 - K1*B1

  calc_t2 = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1)

end function calc_t2
!-----------------------------------------------------------------------
function calc_t1(t0, t2)
  implicit none
  real(KIND=8), intent(in) :: t0, t2
  real(KIND=8) :: calc_t1
  real(KIND=8) :: U11, U12, U13, U31, U32, U33
  real(KIND=8) :: t0_2, t2_2

  t0_2 = t0*t0
  t2_2 = t2*t2

  U11 = C0(0,1) + C0(1,1)*t0 + C0(2,1)*t0_2
  U12 = C1(0,1) + C1(1,1)*t0 + C1(2,1)*t0_2
  U13 = C2(0,1) + C2(1,1)*t0 + C2(2,1)*t0_2
  U31 = C0(0,2) + C0(1,2)*t2 + C0(2,2)*t2_2
  U32 = C1(0,2) + C1(1,2)*t2 + C1(2,2)*t2_2
  U33 = C2(0,2) + C2(1,2)*t2 + C2(2,2)*t2_2

  calc_t1 = (U31*U13-U11*U33)/(U12*U33-U13*U32)

end function calc_t1
!-----------------------------------------------------------------------
subroutine calc_dih_ang(r1, r2, r3, angle)
!-----------------------------------------------------------------------
! r1=Rab, r2=Rbc, r3=Rcd : angle between planes abc and bcd
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), intent(in) :: r1(3), r2(3), r3(3)
  real(KIND=8), intent(out) :: angle
  real(KIND=8), dimension(3) :: p, q, s
  real(KIND=8) :: arg

  call cross(r1, r2, p)
  call cross(r2, r3, q)
  call cross(r3, r1, s)
  arg = dot_product(p,q)/sqrt(dot_product(p,p)*dot_product(q,q))
  arg = sign(min(abs(arg),1.0d0),arg) ! to be sure abs(arg)<=1
  angle = sign(acos(arg), dot_product(s,r2))

end subroutine calc_dih_ang
!-----------------------------------------------------------------------
subroutine calc_bnd_ang(r1, r2, angle)
!-----------------------------------------------------------------------
! assume that each vector is normalized
! r1=Rba, r2=Rbc: angle between r1 and r2
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), intent(in) :: r1(3), r2(3)
  real(KIND=8), intent(out) :: angle
  real(KIND=8) :: arg

  arg = dot_product(r1, r2)
  arg = sign(min(abs(arg),1.0d0),arg) ! to be sure abs(arg)<=1
  angle = acos(arg)

end subroutine calc_bnd_ang
!-----------------------------------------------------------------------
subroutine cross(p, q, s)
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), dimension(:), intent(in) :: p, q
  real(KIND=8), dimension(:), intent(out) :: s

  s(1) = p(2)*q(3) - p(3)*q(2)
  s(2) = p(3)*q(1) - p(1)*q(3)
  s(3) = p(1)*q(2) - p(2)*q(1)

end subroutine cross
!-----------------------------------------------------------------------
subroutine quaternion(axis, quarter_ang, p)
!-----------------------------------------------------------------------
! calculate quaternion, given rotation axis and angle.
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), dimension(:), intent(in) :: axis
  real(KIND=8), intent(in) :: quarter_ang
  real(KIND=8), dimension(:), intent(out) :: p
  real(KIND=8) :: tan_w, tan_sqr, tan1, cosine, sine

  tan_w = tan(quarter_ang)
  tan_sqr = tan_w * tan_w
  tan1 = 1.0d0 + tan_sqr
  cosine = (1.0d0 - tan_sqr)/tan1
  sine = 2.0d0*tan_w/tan1
  p(1) = cosine
  p(2:4) = axis(1:3) * sine

end subroutine quaternion
!-----------------------------------------------------------------------
subroutine rotation_matrix(q, U)
!-----------------------------------------------------------------------
! constructs rotation matrix U from quaternion q.
!-----------------------------------------------------------------------
  implicit none
  real(KIND=8), dimension(:), intent(in) :: q
  real(KIND=8), dimension(:,:), intent(out) :: U
  real(KIND=8) :: q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12,q13,q22,q23,q33

  q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4)
  b0 = 2.0d0*q0; b1 = 2.0d0*q1
  q00 = b0*q0-1.0d0; q02 = b0*q2; q03 = b0*q3
  q11 = b1*q1;       q12 = b1*q2; q13 = b1*q3  
  b2 = 2.0d0*q2; b3 = 2.0d0*q3
  q01 = b0*q1; q22 = b2*q2; q23 = b2*q3; q33 = b3*q3 
  U(1,1) = q00+q11; U(1,2) = q12-q03; U(1,3) = q13+q02
  U(2,1) = q12+q03; U(2,2) = q00+q22; U(2,3) = q23-q01
  U(3,1) = q13-q02; U(3,2) = q23+q01; U(3,3) = q00+q33

end subroutine rotation_matrix
!----------------------------------------------------------------------------
END MODULE tripep_closure
!----------------------------------------------------------------------------
