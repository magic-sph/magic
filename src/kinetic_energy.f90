!$Id$
module kinetic_energy

   use parallel_mod
   use precision_mod
   use truncation, only: n_r_max, l_max
   use radial_functions, only: r, or1, drx, i_costf_init, d_costf_init, &
                               or2, r_cmb, r_icb, orho1, orho2, sigma
   use physical_parameters, only: prmag, ek, nVarCond
   use num_param, only: tScale, eScale
   use blocking, only: lo_map, st_map
   use horizontal_data, only: dLh
   use logic, only: l_save_out, l_non_rot
   use output_data, only: n_e_kin_file, e_kin_file, tag, n_u_square_file, &
                          u_square_file
   use const, only: pi, vol_oc, one, two, half, four
   use LMLoop_data, only: llm,ulm
   use communications, only: get_global_sum
 
   use integration, only: rInt_R
   use useful, only: cc2real
 
   implicit none
 
   private
 
   real(cp), allocatable :: e_pA(:),e_p_asA(:)
   real(cp), allocatable :: e_tA(:),e_t_asA(:)
 
   public :: get_e_kin, get_u_square, initialize_kinetic_energy
 
contains

   subroutine initialize_kinetic_energy
 
      allocate( e_pA(n_r_max),e_p_asA(n_r_max) )
      allocate( e_tA(n_r_max),e_t_asA(n_r_max) )
     
   end subroutine initialize_kinetic_energy
!-----------------------------------------------------------------------------
   subroutine get_e_kin(time,l_write,l_stop_time,n_e_sets, &
       &               w,dw,z,e_p,e_t,e_p_as,e_t_as,       &
       &               ekinR,ekinRave)

      !--------------------------------------------------------------------
      !
      !  calculates kinetic energy  = 1/2 Integral (v^2 dV)
      !  integration in theta,phi by summation of spherical harmonics
      !  integration in r by using Chebycheff integrals
      !
      !  Output:
      !  e_p: Total poloidal     e_p_as: Total toroidal
      !  e_t: Axisym. poloidal   e_t_as: Axisym. toroidal
      !
      !--------------------------------------------------------------------

      !-- Input variables:
      integer,     intent(in) :: n_e_sets
      real(cp),    intent(in) :: time
      logical,     intent(in) :: l_write
      logical,     intent(in) :: l_stop_time
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)

      !-- Output variables:
      real(cp), intent(out), optional :: ekinR(n_r_max)   
      real(cp), intent(out), optional :: ekinRave(n_r_max) 
      real(cp), intent(out) :: e_p     ! poloidal energy
      real(cp), intent(out) :: e_t     ! toroidal energy
      real(cp), intent(out) :: e_p_as  ! axisymmetric poloidal energy
      real(cp), intent(out) :: e_t_as  ! axisymmetric toroidal energy

      !-- Local variables:
      real(cp) :: e_p_es  ! equatorially symmetric poloidal energy 
      real(cp) :: e_t_es  ! equatorially symmetric toroidal energy
      real(cp) :: e_p_eas ! equator. & axially symmetric poloidal energy
      real(cp) :: e_t_eas ! equator. & axially symmetric toroidal energy

      real(cp) :: e_p_temp,e_t_temp
      real(cp) :: e_p_r(n_r_max)
      real(cp) :: e_t_r(n_r_max)
      real(cp) :: e_p_as_r(n_r_max)
      real(cp) :: e_t_as_r(n_r_max)
      real(cp) :: e_p_es_r(n_r_max)
      real(cp) :: e_t_es_r(n_r_max)
      real(cp) :: e_p_eas_r(n_r_max)
      real(cp) :: e_t_eas_r(n_r_max)

      real(cp) :: e_p_r_global(n_r_max),e_t_r_global(n_r_max)
      real(cp) :: e_p_as_r_global(n_r_max),e_t_as_r_global(n_r_max)
      real(cp) :: e_p_es_r_global(n_r_max),e_t_es_r_global(n_r_max)
      real(cp) :: e_p_eas_r_global(n_r_max),e_t_eas_r_global(n_r_max)

      integer nR,lm,l,m
      real(cp) :: fac
      real(cp) :: O_rho ! 1/rho (anelastic)

      !-- time averaging of e(r):
      character(len=80) :: filename
      real(cp) :: dt,surf
      real(cp), save :: timeLast,timeTot

      !write(*,"(A,6ES22.14)") "ekin: w,dw,z = ",get_global_sum( w(llm:ulm,:) ),&
      !     & get_global_sum( dw(llm:ulm,:) ), get_global_sum( z(llm:ulm,:) )

      do nR=1,n_r_max
         e_p_r(nR)    =0.0_cp
         e_t_r(nR)    =0.0_cp
         e_p_as_r(nR) =0.0_cp
         e_t_as_r(nR) =0.0_cp
         e_p_es_r(nR) =0.0_cp
         e_t_es_r(nR) =0.0_cp
         e_p_eas_r(nR)=0.0_cp
         e_t_eas_r(nR)=0.0_cp
         O_rho        =orho1(nR)
         !do lm=2,lm_max
         do lm=max(2,llm),ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)

            e_p_temp= O_rho*dLh(st_map%lm2(l,m)) * ( &
                 &      dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                 &      + cc2real(dw(lm,nR),m) )
            e_t_temp= O_rho*dLh(st_map%lm2(l,m)) * cc2real(z(lm,nR),m)
            !write(*,"(A,3I4,ES22.14)") "e_p_temp = ",nR,l,m,e_p_temp
            if ( m == 0 ) then  ! axisymmetric part
               e_p_as_r(nR) = e_p_as_r(nR) + e_p_temp
               e_t_as_r(nR) = e_t_as_r(nR) + e_t_temp
               if ( mod(l,2) == 0 ) then
                  e_p_eas_r(nR)=e_p_eas_r(nR)+e_p_temp
               else
                  e_t_eas_r(nR)=e_t_eas_r(nR)+e_t_temp
               end if
            else
               e_p_r(nR)=e_p_r(nR) + e_p_temp
               e_t_r(nR)=e_t_r(nR) + e_t_temp
            end if
            if ( mod(l+m,2) == 0 ) then
               e_p_es_r(nR)=e_p_es_r(nR)+e_p_temp
            else
               e_t_es_r(nR)=e_t_es_r(nR)+e_t_temp
            end if

            !write(*,"(8X,A,4I4,ES22.14)") "e_p_r: ",lm,l,m,nR,e_p_r(nR)
            
         end do    ! do loop over lms in block
         e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
         e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)
         !write(*,"(4X,A,I4,2ES22.14)") "e_p_r: ",nR,e_p_r(nR),e_p_as_r(nR)
      end do    ! radial grid points

      ! reduce over the ranks
#ifdef WITH_MPI
      call MPI_Reduce(e_p_r,    e_p_r_global,     n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_r,    e_t_r_global,     n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_as_r, e_p_as_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_as_r, e_t_as_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_es_r, e_p_es_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_es_r, e_t_es_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_eas_r,e_p_eas_r_global, n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_eas_r,e_t_eas_r_global, n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      e_p_r_global    =e_p_r
      e_t_r_global    =e_t_r
      e_p_as_r_global =e_p_as_r
      e_t_as_r_global =e_t_as_r
      e_p_es_r_global =e_p_es_r
      e_t_es_r_global =e_t_es_r
      e_p_eas_r_global=e_p_eas_r
      e_t_eas_r_global=e_t_eas_r
#endif

      if ( rank == 0 ) then
         !do nR=1,n_r_max
         !   write(*,"(4X,A,I4,ES22.14)") "e_p_r_global: ",nR,e_p_r_global(nR)
         !end do
         !-- Radial Integrals:
         e_p    = rInt_R(e_p_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_t    = rInt_R(e_t_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_p_as = rInt_R(e_p_as_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_t_as = rInt_R(e_t_as_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_p_es = rInt_R(e_p_es_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_t_es = rInt_R(e_t_es_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_p_eas= rInt_R(e_p_eas_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
         e_t_eas= rInt_R(e_t_eas_r_global,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)

         fac    =half*eScale
         e_p    =fac*e_p
         e_t    =fac*e_t
         e_p_as =fac*e_p_as
         e_t_as =fac*e_t_as
         e_p_es =fac*e_p_es
         e_t_es =fac*e_t_es
         e_p_eas=fac*e_p_eas
         e_t_eas=fac*e_t_eas

         !-- OUTPUT:
         if ( present(ekinR) ) then
            do nR=1,n_r_max
               ekinR(nR)=fac*(e_p_r_global(nR)+e_t_r_global(nR))
            end do
         end if
         if ( l_write ) then
            if ( l_save_out ) then
               open(n_e_kin_file, file=e_kin_file, status='unknown', position='append')
            end if
            write(n_e_kin_file,'(1P,D20.12,8D16.8)')    &
                 & time*tScale, &  ! 1
                 & e_p,e_t,       &! 2,3
                 & e_p_as,e_t_as, &! 4,5
                 & e_p_es,e_t_es, &! 6,7
                 & e_p_eas,e_t_eas ! 8,9
            if ( l_save_out ) close(n_e_kin_file)
         end if

         ! NOTE: n_e_sets=0 prevents averaging
         if ( n_e_sets == 1 ) then
            timeTot=one
            e_pA    = e_p_r_global
            e_p_asA = e_p_r_global
            e_tA    = e_t_r_global
            e_t_asA = e_t_r_global
         else if ( n_e_sets == 2 ) then
            dt=time-timeLast
            timeTot=two*dt
            e_pA    = dt*(e_pA   +e_p_r_global   )
            e_p_asA = dt*(e_p_asA+e_p_as_r_global)
            e_tA    = dt*(e_tA   +e_t_r_global   )
            e_t_asA = dt*(e_t_asA+e_t_as_r_global)
         else
            dt=time-timeLast
            timeTot=timeTot+dt
            e_pA    = e_pA    + dt*e_p_r_global
            e_p_asA = e_p_asA + dt*e_p_as_r_global
            e_tA    = e_tA    + dt*e_t_r_global
            e_t_asA = e_t_asA + dt*e_t_as_r_global
         end if

         !write(*,"(A,2ES22.14)") "e_pA, e_tA = ",SUM( e_pA ),SUM( e_tA )
         if ( l_stop_time .and. (n_e_sets > 1) ) then
            fac=half*eScale
            filename='eKinR.'//tag
            open(99, file=filename, status='unknown')
            if ( present(ekinRave) ) then
               ekinRave(1)      =fac*(e_pA(1)+e_tA(1))/timetot
               ekinRave(n_r_max)=fac*(e_pA(n_r_max)+e_tA(n_r_max))/timetot
               do nR=1,n_r_max
                  ekinRave(nR)  =fac*e_pA(nR)/timetot+fac*e_tA(nR)/timetot
               end do
            end if
            do nR=1,n_r_max
               surf=four*pi*r(nR)**2
               write(99,'(2x,9D12.4)',advance='no') r(nR),        &
                    &               fac*e_pA(nR)/timetot,         &
                    &               fac*e_p_asA(nR)/timetot,      &
                    &               fac*e_tA(nR)/timetot,         &
                    &               fac*e_t_asA(nR)/timetot,      &
                    &               fac*e_pA(nR)/timetot/surf,    &
                    &               fac*e_p_asA(nR)/timetot/surf, &
                    &               fac*e_tA(nR)/timetot/surf,    &
                    &               fac*e_t_asA(nR)/timetot/surf

               if ( present(ekinR) ) then
                  write(99,'(D12.4)',advance='no') ekinR(nR)
               else
                  write(99,'(A)') ""
               end if
               if ( present(ekinRave) ) then
                  write(99,'(D12.4)') ekinRave(nR)
               else
                  write(99,'(A)') ""
               end if
            end do
            close(99)
         end if
         timeLast=time
      end if

      ! broadcast the output arguments of the function to have them on all ranks
      ! e_p,e_t,e_p_as,e_t_as
#ifdef WITH_MPI
      call MPI_Bcast(e_p,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(e_t,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(e_p_as,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(e_t_as,1,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)

      if ( present(ekinR) ) then
         call MPI_Bcast(ekinR,n_r_max,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      end if
      if ( present(ekinRave) ) then
         call MPI_Bcast(ekinRave,n_r_max,MPI_DEF_REAL,0,MPI_COMM_WORLD,ierr)
      end if
#endif

   end subroutine get_e_kin
!-----------------------------------------------------------------------------
   subroutine get_u_square(time,w,dw,z,RolR,dlR,dlRc)
      !--------------------------------------------------------------------

      !  calculates square velocity  = 1/2 Integral (v^2 dV)
      !  integration in theta,phi by summation of spherical harmonics
      !  integration in r by using Chebychef integrals

      !  Write the different contributions in u_square.TAG file

      !--------------------------------------------------------------------

      !-- Input of scalar fields:
      real(cp),    intent(in) :: time
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(in) :: dw(llm:ulm,n_r_max)
      complex(cp), intent(in) :: z(llm:ulm,n_r_max)

      !-- Output variables:
      real(cp), intent(out) :: dlR(n_r_max)
      real(cp), intent(out) :: dlRc(n_r_max)
      real(cp), intent(out) :: RolR(n_r_max)

      !-- Local variables:
      real(cp) :: e_p     ! poloidal u**2
      real(cp) :: e_t     ! toroidal u**2
      real(cp) :: e_p_as  ! axisymmetric poloidal u**2
      real(cp) :: e_t_as  ! axisymmetric toroidal u**2
      real(cp) :: e_kin   ! total u**2

      real(cp) :: e_p_temp,e_t_temp
      real(cp) :: e_p_r(n_r_max),e_p_r_global(n_r_max)
      real(cp) :: e_t_r(n_r_max),e_t_r_global(n_r_max)
      real(cp) :: e_p_as_r(n_r_max),e_p_as_r_global(n_r_max)
      real(cp) :: e_t_as_r(n_r_max),e_t_as_r_global(n_r_max)
      real(cp) :: e_lr(n_r_max,l_max), e_lr_global(n_r_max,l_max)
      real(cp) :: e_lr_c(n_r_max,l_max), e_lr_c_global(n_r_max,l_max)
      real(cp) :: ER(n_r_max),ELR(n_r_max),ReR(n_r_max),RoR(n_r_max)
      real(cp) :: ERc(n_r_max),ELRc(n_r_max)
      real(cp) :: ekinR(n_r_max)
      real(cp) :: RmR(n_r_max)
      real(cp) :: e_l,E,EL,Ec,ELc

      integer :: nR,lm,l,m
      real(cp) :: fac
      real(cp) :: O_rho ! 1/rho**2 (anelastic)

      !-- property parameters
      real(cp) :: Re,Rm,Ro,Rol,dl,dlc
      real(cp) :: ReConv,RoConv,RolC

      do nR=1,n_r_max
         e_p_r(nR)    =0.0_cp
         e_t_r(nR)    =0.0_cp
         e_p_as_r(nR) =0.0_cp
         e_t_as_r(nR) =0.0_cp
         O_rho        =orho2(nR) ! divided by rho**2
         do l=1,l_max
            e_lr(nR,l)=0.0_cp
            e_lr_c(nR,l)=0.0_cp
         end do

         do lm=max(2,llm),ulm
            l=lo_map%lm2l(lm)
            m=lo_map%lm2m(lm)
            !do lm=2,lm_max
            !  l=lm2l(lm)
            !  m=lm2m(lm)

            e_p_temp= O_rho*dLh(st_map%lm2(l,m)) * ( &
                 &      dLh(st_map%lm2(l,m))*or2(nR)*cc2real(w(lm,nR),m) &
                 &      + cc2real(dw(lm,nR),m) )
            e_t_temp= O_rho*dLh(st_map%lm2(l,m))*cc2real(z(lm,nR),m)
            if ( m == 0 ) then  ! axisymmetric part
               e_p_as_r(nR)=e_p_as_r(nR)+ e_p_temp
               e_t_as_r(nR)=e_t_as_r(nR)+ e_t_temp
            else
               e_p_r(nR)=e_p_r(nR) + e_p_temp
               e_t_r(nR)=e_t_r(nR) + e_t_temp
               e_lr_c(nR,l)=e_lr_c(nR,l) + e_p_temp + e_t_temp
            end if

            e_lr(nR,l)=e_lr(nR,l) + e_p_temp + e_t_temp
         end do    ! do loop over lms in block
         e_p_r(nR)=e_p_r(nR)+e_p_as_r(nR)
         e_t_r(nR)=e_t_r(nR)+e_t_as_r(nR)
      end do    ! radial grid points

      ! reduce over the ranks
#ifdef WITH_MPI
      call MPI_Reduce(e_p_r,    e_p_r_global,     n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_r,    e_t_r_global,     n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_p_as_r, e_p_as_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_t_as_r, e_t_as_r_global,  n_r_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_lr_c, e_lr_c_global,  n_r_max*l_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(e_lr, e_lr_global,  n_r_max*l_max, &
           & MPI_DEF_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#else
      e_p_r_global   =e_p_r
      e_t_r_global   =e_t_r
      e_p_as_r_global=e_p_as_r
      e_t_as_r_global=e_t_as_r
      e_lr_c_global  =e_lr_c
      e_lr_global    =e_lr
#endif

      if ( rank == 0 ) then
         !-- Radial Integrals:
         e_p    =rInt_R(e_p_r_global,    n_r_max,n_r_max,drx, &
              i_costf_init,d_costf_init)
         e_t    =rInt_R(e_t_r_global,    n_r_max,n_r_max,drx, &
              i_costf_init,d_costf_init)
         e_p_as =rInt_R(e_p_as_r_global, n_r_max,n_r_max,drx, &
              i_costf_init,d_costf_init)
         e_t_as =rInt_R(e_t_as_r_global, n_r_max,n_r_max,drx, &
              i_costf_init,d_costf_init)
         fac    =half*eScale
         e_p    =fac*e_p
         e_t    =fac*e_t
         e_p_as =fac*e_p_as
         e_t_as =fac*e_t_as

         e_kin  =e_t+e_p
         do nR=1,n_r_max
            ekinR(nR)=fac*(e_p_r_global(nR)+e_t_r_global(nR))
         end do

         !-- Rossby number
         Re=sqrt(two*e_kin/vol_oc)
         ReConv=sqrt(two*(e_kin-e_p_as-e_t_as)/vol_oc)
         if ( l_non_rot ) then
            Ro=0.0_cp
            RoConv=0.0_cp
         else
            Ro=Re*ek
            RoConv=ReConv*ek
         end if

         !-- Length Scale
         E  =0.0_cp
         EL =0.0_cp
         Ec =0.0_cp
         ELc=0.0_cp
         do l=1,l_max
            e_l=fac*rInt_R(e_lr_global(1,l),n_r_max,n_r_max,drx, &
                           i_costf_init,d_costf_init)
            E =E+e_l
            EL=EL+real(l,cp)*e_l
            e_l=fac*rInt_R(e_lr_c_global(1,l),n_r_max,n_r_max,drx, &
                           i_costf_init,d_costf_init)
            Ec =Ec+e_l
            ELc=ELc+real(l,cp)*e_l
         end do
         if ( EL /= 0.0_cp ) then
            dl=pi*E/EL
            dlc=pi*Ec/ELc
         else
            dl=0.0_cp
            dlc=0.0_cp
         end if
         do nR=1,n_r_max
            ER(nR)  =0.0_cp
            ELR(nR) =0.0_cp
            ERc(nR) =0.0_cp
            ELRc(nR)=0.0_cp
            do l=1,l_max
               e_l=fac*e_lr_global(nR,l)
               ER(nR) =ER(nR)+e_l
               ELR(nR)=ELR(nR)+real(l,cp)*e_l
               e_l=fac*e_lr_c_global(nR,l)
               ERc(nR) =ERc(nR)+e_l
               ELRc(nR)=ELRc(nR)+real(l,cp)*e_l
            end do
            if ( ELR(nR) /= 0.0_cp ) then
               dlR(nR)=pi*ER(nR)/ELR(nR)
               dlRc(nR)=pi*ERc(nR)/ELRc(nR)
            else
               dlR(nR)=0.0_cp
               dlRc(nR)=0.0_cp
            end if
         end do

         !-- Local Rossby number
         if ( dl/=0.0_cp ) then
            Rol = Ro/dl
            RolC = RoConv/dlc
         else
            Rol = Ro
            RolC = RoConv
         end if
         do nR=1,n_r_max
            ReR(nR)=sqrt(two*ekinR(nR)*or2(nR)/(4*pi))
            RoR(nR)=ReR(nR)*ek
            if ( dlR(nR) /= 0.0_cp ) then
               RolR(nR)=RoR(nR)/dlR(nR)
            else
               RolR(nR)=RoR(nR)
            end if
            RmR(nR)=ReR(nR)*prmag*sigma(nR)*r(nR)*r(nR)
         end do

         !-- Magnetic reynolds number
         if ( prmag /= 0 .and. nVarCond > 0 ) then
            Rm=0.0_cp
            Rm=rInt_R(RmR,n_r_max,n_r_max,drx,i_costf_init,d_costf_init)
            Rm=Rm*3/(r_cmb**3-r_icb**3)
         elseif ( prmag /= 0 ) then
            Rm=Re*prmag
         else
            Rm=Re
         end if

         !-- Output
         if ( l_save_out ) then
            open(n_u_square_file, file=u_square_file, status='unknown', &
                 position='append')
         end if
         write(n_u_square_file,'(1P,D20.12,10D16.8)') &
              &  time*tScale,     & ! 1
              &      e_p,e_t,     & ! 2,3
              &e_p_as,e_t_as,     & ! 4,5
              &        Ro,Rm,     & ! 6,7
              &       Rol,dl,     & ! 8,9
              &     RolC,dlc        ! 10,11
         if ( l_save_out ) close(n_u_square_file)
      end if

   end subroutine get_u_square
!-----------------------------------------------------------------------------
end module kinetic_energy
