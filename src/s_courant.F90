!$Id$
!***********************************************************************
SUBROUTINE courant(n_r,dtrkc,dthkc,vr,vt,vp,br,bt,bp, &
     n_theta_min,n_theta_block)
  !***********************************************************************

  !    !------------ This is release 2 level 1  --------------!
  !    !------------ Created on 1/17/02  by JW. --------------!

  !  +-------------+----------------+------------------------------------+
  !  |                                                                   |
  !  |  courant condition check: calculates Courant                      |
  !  |  advection lengths in radial direction dtrkc                      |
  !  |  and in horizontal direction dthkc                                |
  !  |  on the local radial level n_r                                    |
  !  |                                                                   |
  !  |  for the effective velocity, the abs. sum of fluid                |
  !  |  velocity and Alfven velocity is taken                            |
  !  |                                                                   |
  !  |  instead of the full Alfven velocity                              |
  !  |  a modified Alfven velocity is employed that takes                |
  !  |  viscous and Joule damping into account. Different                |
  !  |  Courant factors are used for the fluid velocity and              |
  !  |  the such modified Alfven velocity                                |
  !  |                                                                   |
  !  +-------------------------------------------------------------------+
  !-----------------------------------------------------------------------

  USE truncation
  USE radial_functions
  USE physical_parameters
  USE num_param
  USE blocking
  USE horizontal_data
  USE logic

  IMPLICIT NONE

  !-- input:

  INTEGER,intent(IN) :: n_r           ! radial level
  integer,intent(IN) :: n_theta_min   ! first theta in block stored in fields
  integer,intent(IN) :: n_theta_block ! size of theta block

  real(kind=8),intent(IN) :: vr(nrp,nfs)      ! radial velocity
  real(kind=8),intent(IN) :: vt(nrp,nfs)      ! longitudinal velocity
  real(kind=8),intent(IN) :: vp(nrp,nfs)      ! azimuthal velocity
  real(kind=8),intent(IN) :: br(nrp,nfs)      ! radial magnetic field
  real(kind=8),intent(IN) :: bt(nrp,nfs)      ! longitudinal magnetic field
  real(kind=8),intent(IN) :: bp(nrp,nfs)      ! azimuthal magnetic field

  !-- output:
  real(kind=8),intent(INOUT) :: dtrkc          ! Courant step (based on radial advection)
  ! for the range of points covered
  real(kind=8),intent(INOUT) :: dthkc          ! Courant step based on horizontal advection

  !-- local:
  integer :: n_theta       ! absolut no of theta
  integer :: n_theta_rel   ! no of theta in block
  integer :: n_theta_nhs   ! no of theta in NHS
  integer :: n_phi         ! no of longitude

  real(kind=8) :: valri2,valhi2,valh2,valh2m
  real(kind=8) :: vr2max,vh2max
  real(kind=8) :: valr,valr2,vflr2,vflh2
  real(kind=8) :: O_r_E_2,O_r_E_4
  real(kind=8) :: cf2,af2

  !-- end of declaration
  !-----------------------------------------------------------------------


  valri2=(0.5d0*(1.d0+opm))**2/delxr2(n_r)
  valhi2=(0.5d0*(1.d0+opm))**2/delxh2(n_r)

  vr2max=0.d0
  vh2max=0.d0
  cf2=courfac*courfac
  O_r_E_4=or4(n_r)
  O_r_E_2=or2(n_r)

  n_theta=n_theta_min-1

  if ( l_mag .AND. l_mag_LF .AND. .NOT. l_mag_kin ) then

     af2=alffac*alffac

     do n_theta_rel=1,n_theta_block

        n_theta=n_theta+1
        n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

        do n_phi=1,n_phi_max

           vflr2=orho2(n_r)*vr(n_phi,n_theta_rel) &
                *vr(n_phi,n_theta_rel)
           valr =br(n_phi,n_theta_rel)*br(n_phi,n_theta_rel) * &
                LFfac*orho1(n_r)
           valr2=valr*valr/(valr+valri2)
           vr2max=max(vr2max,O_r_e_4*(cf2*vflr2+af2*valr2))

           vflh2= ( vt(n_phi,n_theta_rel)*vt(n_phi,n_theta_rel) + &
                vp(n_phi,n_theta_rel)*vp(n_phi,n_theta_rel) )* &
                osn2(n_theta_nhs)*orho2(n_r)
           valh2= ( bt(n_phi,n_theta_rel)*bt(n_phi,n_theta_rel) + &
                bp(n_phi,n_theta_rel)*bp(n_phi,n_theta_rel) )* &
                LFfac*osn2(n_theta_nhs)*orho1(n_r)
           valh2m=valh2*valh2/(valh2+valhi2)
           vh2max=MAX(vh2max,O_r_E_2*(cf2*vflh2+af2*valh2m))

        end do

     end do

  else   ! Magnetic field ?

     do n_theta_rel=1,n_theta_block

        n_theta=n_theta+1
        n_theta_nhs=(n_theta+1)/2 ! northern hemisphere=odd n_theta

        do n_phi=1,n_phi_max

           vflr2=orho2(n_r)*vr(n_phi,n_theta_rel)* &
                vr(n_phi,n_theta_rel)
           vr2max=max(vr2max,cf2*O_r_E_4*vflr2)

           vflh2= ( vt(n_phi,n_theta_rel)*vt(n_phi,n_theta_rel) + &
                vp(n_phi,n_theta_rel)*vp(n_phi,n_theta_rel) )* &
                osn2(n_theta_nhs)*orho2(n_r)
           vh2max=max(vh2max,cf2*O_r_E_2*vflh2)

        end do

     end do

  end if   ! Magnetic field ?


  !$OMP CRITICAL
  if ( vr2max /= 0.d0 ) dtrkc=min(dtrkc,dsqrt(delxr2(n_r)/vr2max))
  if ( vh2max /= 0.d0 ) dthkc=min(dthkc,dsqrt(delxh2(n_r)/vh2max))
  !$OMP END CRITICAL

  return
end SUBROUTINE courant

!-----------------------------------------------------------------------
