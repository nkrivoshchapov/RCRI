!----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
MODULE main_pep
!-----------------------------------------------------------------------
! This is a sample driver routine to reconstruct tripeptide loops
! from the coordinates in a pdb file using a canonical bond lengths and
! angles.
!-----------------------------------------------------------------------
  use tripep_closure

CONTAINS
!-----------------------------------------------------------------------
subroutine runtlc(b_len,b_ang,t_ang,r0_n,r0_a,r0_c,r_soln_n,r_soln_a,r_soln_c,nsol)
  implicit none
  character(len=15) :: res_name(5)
  character(len=10) :: tor1char, tor2char,confcountstr
  integer :: n0, n_soln, unit = 17, i, j, k, confint, n1, n2, order(16)
  integer, intent(out) :: nsol;
  real(KIND=8) :: r_n(3,5), r_a(3,5), r_c(3,5)
  real(KIND=8), intent(out) :: r_soln_n(3,3,16), r_soln_a(3,3,16), r_soln_c(3,3,16)
  real(KIND=8) :: rmsd, dr(3), tor1fl,tor2fl
  real(KIND=8), intent(in) :: r0_n(3,3), r0_a(3,3), r0_c(3,3)

  real(KIND=8), intent(in) :: b_len(6), b_ang(7), t_ang(2)

  ! input parameters (lengths)
  real(KIND=8), parameter :: bl = 1.56d0
  ! input parameters (angles in radians)
  !do i = 1,6
!	write(*,*) 'i = ',i,' bond = ',b_len(i)
 ! end do
  !do i = 1,7
!	write(*,*) 'i = ',i,' va = ',b_ang(i)*rad2deg
 ! end do
  !do i = 1,2
!	write(*,*) 'i = ',i,' ta = ',t_ang(i)*rad2deg
 ! end do

  call initialize_loop_closure(b_len, b_ang, t_ang)

  call solv_3pep_poly(r0_n(:,1), r0_a(:,1), r0_a(:,3), r0_c(:,3), &
           r_soln_n, r_soln_a, r_soln_c, n_soln)
   nsol = n_soln
900 format(' ---------------------------------------------------')

end subroutine runtlc
!-----------------------------------------------------------------------
END MODULE main_pep
