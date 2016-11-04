module init_fields

   use truncation
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use blocking, only: nfs, nThetaBs, sizeThetaB, st_map, lmP2lmPS
   use horizontal_data, only: sinTheta, dLh, dTheta1S, dTheta1A, D_l, &
       &                      phi, cosTheta
   use logic, only: l_rot_ic, l_rot_ma, l_SRIC, l_SRMA, l_anelastic_liquid, &
       &            l_cond_ic, l_temperature_diff, l_single_matrix,    &
       &            l_chemical_conv, l_non_adia
   use radial_functions, only: r_icb, r, r_cmb, r_ic, or1, jVarCon,    &
       &                       cheb_norm, lambda, or2, d2cheb, dcheb,  &
       &                       cheb, dLlambda, or3, cheb_ic, dcheb_ic, &
       &                       d2cheb_ic, cheb_norm_ic, or1, r_ic,     &
       &                       orho1, chebt_oc, chebt_ic, temp0,       &
       &                       dLtemp0, kappa, dLkappa, beta, dbeta,   &
       &                       epscProf, ddLtemp0, ddLalpha0, rgrav,   &
       &                       rho0, dLalpha0, alpha0
   use radial_data, only: n_r_icb, n_r_cmb
   use constants, only: pi, y10_norm, c_z10_omega_ic, c_z10_omega_ma, osq4pi, &
       &            zero, one, two, three, four, third, half
   use useful, only: random
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
#ifdef WITH_SHTNS
   use shtns
#else
   use fft
#endif
   use physical_parameters, only: impS, n_impS_max, n_impS, phiS, thetaS, &
       &                          peakS, widthS, radratio, imagcon, opm,  &
       &                          sigma_ratio, O_sr, kbots, ktops, opr,   &
       &                          epsc, ViscHeatFac, ThExpNb, ogrun,      &
       &                          impXi, n_impXi_max, n_impXi, phiXi,     &
       &                          thetaXi, peakXi, widthXi, osc, epscxi,  &
       &                          kbotxi, ktopxi, BuoFac
   use algebra, only: sgesl, sgefa, cgesl
   use horizontal_data, only: D_lP1, hdif_B, dLh
   use matrices, only: jMat, jPivot, s0Mat, s0Pivot, ps0mat, ps0pivot, &
       &               ps0Mat_fac, xi0Mat, xi0Pivot
   use legendre_grid_to_spec, only: legTF1
   use cosine_transform_odd

   implicit none

   private

   !-- Initialisation of fields:
   integer, public :: init_s1,init_s2
   integer, public :: init_xi1,init_xi2
   integer, public :: init_b1,init_v1

   !----- Entropy amplitudes for initialisation:
   real(cp), public :: amp_s1,amp_s2,amp_v1,amp_b1,amp_xi1,amp_xi2

   !----- Entropy at CMB and ICB (input):
   integer, public, parameter :: n_s_bounds=20
   real(cp), public :: s_bot(4*n_s_bounds)  ! input variables for tops,bots
   real(cp), public :: s_top(4*n_s_bounds)
   complex(cp), public, allocatable :: tops(:,:)
   complex(cp), public, allocatable :: bots(:,:)

   !----- Chemical composition
   integer, public, parameter :: n_xi_bounds=20
   real(cp), public :: xi_bot(4*n_xi_bounds)  ! input variables for topxi,botxi
   real(cp), public :: xi_top(4*n_xi_bounds)
   complex(cp), public, allocatable :: topxi(:,:)
   complex(cp), public, allocatable :: botxi(:,:)

   !----- Peak values for magnetic field:
   real(cp), public :: bpeakbot,bpeaktop

   !----- Initialised IC and mantle rotation rates:
   integer, public :: nRotMa,nRotIc
   real(cp), public :: omega_ma1,omegaOsz_ma1,tShift_ma1,tOmega_ma1
   real(cp), public :: omega_ma2,omegaOsz_ma2,tShift_ma2,tOmega_ma2
   real(cp), public :: omega_ic1,omegaOsz_ic1,tShift_ic1,tOmega_ic1
   real(cp), public :: omega_ic2,omegaOsz_ic2,tShift_ic2,tOmega_ic2

   !----- About start-file:
   logical, public :: l_start_file     ! taking fields from startfile ?
   logical, public :: l_reset_t        ! reset time from startfile ?
   integer, public :: inform           ! format of start_file
   integer, public :: n_start_file     ! I/O unit of start_file
   character(len=72), public :: start_file  ! name of start_file           

   !-- Scales for input field:
   real(cp), public :: scale_s
   real(cp), public :: scale_xi
   real(cp), public :: scale_v
   real(cp), public :: scale_b
   real(cp), public :: tipdipole       ! adding to symetric field

   public :: initialize_init_fields, initV, initS, initB, s_cond, &
             ps_cond, initXi, xi_cond

contains

   subroutine initialize_init_fields
      !
      ! Memory allocation
      !

      n_start_file=8
      
      allocate( tops(0:l_max,0:m_max) )
      allocate( bots(0:l_max,0:m_max) )
      tops(:,:)=zero
      bots(:,:)=zero
      bytes_allocated = bytes_allocated+2*(l_max+1)*(m_max+1)*SIZEOF_DEF_COMPLEX

      if ( l_chemical_conv ) then
         allocate( topxi(0:l_max,0:m_max) )
         allocate( botxi(0:l_max,0:m_max) )
         topxi(:,:)=zero
         botxi(:,:)=zero
         bytes_allocated = bytes_allocated+2*(l_max+1)*(m_max+1)*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_init_fields
!-----------------------------------------------------------------------
   subroutine initV(w,z,omega_ic,omega_ma,lmStart,lmStop)
      !
      ! Purpose of this subroutine is to initialize the velocity field   
      ! So far it is only rudimentary and will be expanded later.        
      ! Because s is needed for dwdt init_s has to be called before.     
      !                                                                   

      !-- Input variables
      integer, intent(in) :: lmStart,lmStop
    
      !-- Output variables
      complex(cp), intent(inout) :: w(lm_max,n_r_max)
      complex(cp), intent(inout) :: z(lm_max,n_r_max)
      real(cp), intent(out) :: omega_ic,omega_ma
    
      !-- Local variables
      integer :: lm,l,m,n,lmStart00,st_lmP
      integer :: nR,nTheta,nThetaB,nThetaStart,nPhi
      real(cp) :: ra1,ra2,c_r,c_i
      real(cp) :: amp_r,rExp
      real(cp) :: rDep(n_r_max)
      integer :: l1m0
    
      real(cp) :: ss,ome(nrp,nfs)
      complex(cp) :: omeLM(lmP_max)
    
      ! This routine is called from a rank-0 region and is operating
      ! on the full fields in standard order.
    
      lmStart00=max(lmStart,1)
      l1m0=st_map%lm2(1,0)
    
      !-- Initialize rotation according to
      !   given inner core and mantel rotation rate:
      if ( init_v1 == 1 .and. ( omega_ic1 /= 0.0_cp .or. omega_ma1 /= 0.0_cp ) ) then
    
         !-- Approximating the Stewardson solution:
         do nR=1,n_r_max
    
            nTheta=0
            do n=1,nThetaBs ! loop over the theta blocks
    
               nThetaStart=(n-1)*sizeThetaB+1
               do nThetaB=1,sizeThetaB
                  nTheta=nTheta+1
                  ss=r(nR)*sinTheta(nTheta)
                  !------------ start with constructing rotation rate ome:
                  do nPhi=1,n_phi_max
                     if ( ss <= r_icb ) then
                        ome(nPhi,nThetaB)=omega_ma1+half*omega_ic1
                     else
                        ome(nPhi,nThetaB)=omega_ma1
                     end if
                  end do
               end do
               !------------ Transform to spherical hamonic space for each theta block
#ifndef WITH_SHTNS
               call fft_thetab(ome,-1)
               call legTF1(nThetaStart,omeLM,ome)
#endif
            end do ! End of loop over theta blocks
#ifdef WITH_SHTNS
               call spat_to_SH(ome, omeLM)
#endif
    
            !------------ ome now in spherical harmonic space,
            !             apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
            !             additional application of r**2/(l*(l+1)) then yields
            !             the axisymmetric toriodal flow contribution:
            do lm=lmStart00,lmStop
               l   =st_map%lm2l(lm)
               m   =st_map%lm2m(lm)
               st_lmP=st_map%lm2lmP(st_map%lm2(l,m))
               if ( l > m ) then
                  z(lm,nR)=z(lm,nR) + &
                       r(nR)**2/dLh(lm) * ( &
                       dTheta1S(lm)*omeLM(lmP2lmPS(st_lmP)) &
                       & - dTheta1A(st_map%lm2(l,m))*omeLM(st_map%lmP2lmPA(st_lmP)) )
               else if ( l == m ) then
                  if ( dLh(st_map%lm2(l,m)) /= 0.0_cp ) then 
                     z(lm,nR)=z(lm,nR) - &
                          r(nR)**2/dLh(st_map%lm2(l,m)) * &
                          dTheta1A(st_map%lm2(l,m))*omeLM(st_map%lmP2lmPA(st_lmP))
                  end if
               end if
            end do
    
         end do ! close loop over radial grid points
    
      else if ( init_v1 == 2 ) then
    
         !-- Approximating the Stewardson solution:
         do nR=1,n_r_max
    
            nTheta=0
            do n=1,nThetaBs ! loop over the theta blocks
    
               nThetaStart=(n-1)*sizeThetaB+1
               do nThetaB=1,sizeThetaB
                  nTheta=nTheta+1
                  ss=r(nR)*sinTheta(nTheta)
                  !------------ start with constructing rotation rate ome:
                  do nPhi=1,n_phi_max
                     !ome(nPhi,nThetaB)=amp_v1*(one-(r(nR)-r(n_r_max))**2/r(nR)**3)
                     !ome(nPhi,nThetaB)=amp_v1*r_icb/r(nR)
                     ome(nPhi,nThetaB)=amp_v1/sqrt(one+ss**4)
                  end do
               end do
               !------------ Transform to spherical hamonic space for each theta block
#ifndef WITH_SHTNS
               call fft_thetab(ome,-1)
               call legTF1(nThetaStart,omeLM,ome)
#endif
            end do ! End of loop over theta blocks
#ifdef WITH_SHTNS
            call spat_to_SH(ome, omeLM)
#endif
    
            !------------ ome now in spherical harmonic space,
            !             apply operator dTheta1=1/(r sinTheta) d/ d theta sinTheta**2,
            !             additional application of r**2/(l*(l+1)) then yields
            !             the axisymmetric toriodal flow contribution:
            do lm=lmStart00,lmStop
               l   =st_map%lm2l(lm)
               m   =st_map%lm2m(lm)
               st_lmP=st_map%lm2lmP(st_map%lm2(l,m))
               if ( l > m ) then
                  z(lm,nR)=z(lm,nR) + &
                       r(nR)**2/dLh(st_map%lm2(l,m)) * ( &
                       dTheta1S(st_map%lm2(l,m))*omeLM(st_map%lmP2lmPS(st_lmP)) &
                       & - dTheta1A(st_map%lm2(l,m))*omeLM(st_map%lmP2lmPA(st_lmP)) )
               else if ( l == m ) then
                  if ( dLh(st_map%lm2(l,m)) /= 0.0_cp ) then 
                      z(lm,nR)=z(lm,nR) - &
                           r(nR)**2/dLh(st_map%lm2(l,m)) * &
                           dTheta1A(st_map%lm2(l,m))*omeLM(st_map%lmP2lmPA(st_lmP))
                  end if
               end if
            end do
    
         end do ! close loop over radial grid points
    
    
      else if ( init_v1 > 2 ) then
    
         !--- Add random noise toroidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(init_v1-1)
         !    Amplitude is chosen so that the (1,0) term resembles amp_v1 *
         !    the 'solid body' rotation set by inner core and mantle rotation.
    
         rExp=4.
         if ( omega_ic1 /= 0 ) then
            amp_r=amp_v1*omega_ic1*r_ICB**(rExp+1.)/y10_norm
         else
            amp_r=amp_v1*r_ICB**(rExp+1.)/y10_norm
         end if
         do nR=1,n_r_max
            rDep(nR)=amp_r/r(nR)**(rExp-1.)
            !write(*,"(A,I3,A,ES20.12)") "rDep(",nR,") = ",rDep(nR)
            do lm=lmStart00,lmStop
               l=st_map%lm2l(lm)
               m=st_map%lm2m(lm)
               if ( D_l(st_map%lm2(l,m)) /= 0.0_cp ) then
                  ra1=(-one+two*random(0.0_cp))/D_l(st_map%lm2(l,m))**(init_v1-1)
                  ra2=(-one+two*random(0.0_cp))/D_l(st_map%lm2(l,m))**(init_v1-1)
                  c_r=ra1*rDep(nR)
                  c_i=ra2*rDep(nR)
                  if ( m == 0 ) then  ! non axisymmetric modes
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,0.0_cp,kind=cp)
                  else
                     z(lm,nR)=z(lm,nR)+cmplx(c_r,c_i,kind=cp)
                  end if
               end if
               write(*,"(A,4I4,2ES20.12)") "z = ",nR,lm,l,m,z(lm,nR)
            end do
         end do
    
      else if ( init_v1 < -1 ) then
    
         !--- Add random noise poloidal field of all (l,m) modes exept (l=0,m=0):
         !    It decays likes l**(init_v1-1)
         !    Amplitude is chosen to be comparable to amp * inner core roation speed
         !    at inner core boundary...
         if ( omega_ic1 /= 0.0_cp ) then
            amp_r=amp_v1*omega_ic1*r_icb*r_icb/(y10_norm*PI)
         else
            amp_r=amp_v1
         end if
         do nR=1,n_r_max
            rDep(nR)=-amp_r*sin( (r(nR)-r_ICB)*PI )
            do lm=lmStart00,lmStop
               l=st_map%lm2l(lm)
               m=st_map%lm2m(lm)
               ra1=(-one+two*random(0.0_cp))/D_l(st_map%lm2(l,m))**(-init_v1-1)
               ra2=(-one+two*random(0.0_cp))/D_l(st_map%lm2(l,m))**(-init_v1-1)
               c_r=ra1*rDep(nR)
               c_i=ra2*rDep(nR)
               if ( m > 0 ) then  ! no axisymmetric modes
                  w(lm,nR)=w(lm,nR)+cmplx(c_r,c_i,kind=cp)
                  z(lm,nR)=z(lm,nR)+cmplx(c_r,c_i,kind=cp)
               end if
            end do
         end do
    
      end if
    
      !----- Caring for IC and mantle rotation rates if this
      !      has not been done already in read_start_file.f:
      if ( lmStart == 1 ) then ! Selects one processor !
         if ( ( .not. l_start_file ) .and.( llm <= st_map%lm2(1,0) ) &
              & .and.( st_map%lm2(1,0) <= ulm ) ) then
            write(*,*) '! NO STARTFILE READ, SETTING Z10!'
            if ( l_SRIC .or. l_rot_ic .and. omega_ic1 /= 0.0_cp ) then
               omega_ic=omega_ic1*cos(omegaOsz_ic1*tShift_ic1) + &
                    omega_ic2*cos(omegaOsz_ic2*tShift_ic2)
               write(*,*)
               write(*,*) '! I use prescribed inner core rotation rate:'
               write(*,*) '! omega_ic=',omega_ic
               z(st_map%lm2(1,0),n_r_icb)=cmplx(omega_ic/c_z10_omega_ic,kind=cp)
            else if ( l_rot_ic .and. omega_ic1 == 0.0_cp ) then
               omega_ic=c_z10_omega_ic*real(z(st_map%lm2(1,0),n_r_icb))
            else
               omega_ic=0.0_cp
            end if
            if ( l_SRMA .or. l_rot_ma .and. omega_ma1 /= 0.0_cp ) then
               omega_ma=omega_ma1*cos(omegaOsz_ma1*tShift_ma1) + &
                        omega_ma2*cos(omegaOsz_ma2*tShift_ma2)
               write(*,*)
               write(*,*) '! I use prescribed mantle rotation rate:'
               write(*,*) '! omega_ma=',omega_ma
               z(st_map%lm2(1,0),n_r_cmb)=cmplx(omega_ma/c_z10_omega_ma,kind=cp)
            else if ( l_rot_ma .and. omega_ma1 == 0.0_cp ) then
               omega_ma=c_z10_omega_ma*real(z(st_map%lm2(1,0),n_r_cmb))
            else
               omega_ma=0.0_cp
            end if
         else
            if ( nRotIc == 2 ) omega_ic=omega_ic1
            if ( nRotMa == 2 ) omega_ma=omega_ma1
         end if
      end if ! lmStart=1 ?
    
   end subroutine initV
!--------------------------------------------------------------------
   subroutine initS(s,p,lmStart,lmStop)
      !
      ! Purpose of this subroutine is to initialize the entropy field    
      ! according to the input control parameters.                       
      !
      ! +-----------------+---------------------------------------------+
      ! | Input           | value                                       |
      ! +=================+=============================================+
      ! | init_s1 < 100:  | random noise initialized                    | 
      ! |                 | the noise spectrum decays as l ^ (init_s1-1)|    
      ! |                 | with peak amplitude amp_s1  for l=1         |    
      ! +-----------------+---------------------------------------------+
      ! | init_s1 >=100:  | a specific harmonic mode initialized        |  
      ! |                 | with amplitude amp_s1.                      |        
      ! |                 | init_s1 is interpreted as number llmm       |         
      ! |                 | where ll: harmonic degree,                  | 
      ! |                 | mm: harmonic order.                         | 
      ! +-----------------+---------------------------------------------+
      ! | init_s2 >100 :  | a second harmonic mode initialized          |  
      ! |                 | with amplitude amp_s2.                      |         
      ! |                 | init_s2 is again interpreted as number llmm |         
      ! |                 | where ll: harmonic degree,                  |
      ! |                 | mm: harmonic order.                         |
      ! +-----------------+---------------------------------------------+
      !                                                                   

      !-- Input variables:
      integer, intent(in) :: lmStart,lmStop

      !-- Output variables:
      complex(cp), intent(inout) :: s(lm_max,n_r_max)
      complex(cp), intent(inout) :: p(lm_max,n_r_max)

      !-- Local variables:
      integer :: n_r,lm,l,m,lm00,lmMin
      real(cp) :: x,rr,c_r,c_i,s_r,s_i
      real(cp) :: ra1,ra2
      real(cp) :: s0(n_r_max),p0(n_r_max),s1(n_r_max)

      integer :: nTheta,n,nThetaStart,nThetaB,nPhi,nS
      real(cp) :: xL,yL,zL,rH,angleL,s00,s00P
      real(cp) :: mata(n_impS_max,n_impS_max)
      real(cp) :: amp(n_impS_max)
      integer :: pivot(n_impS_max)
      real(cp) :: xS(n_impS_max),yS(n_impS_max)
      real(cp) :: zS(n_impS_max),sFac(n_impS_max)
      real(cp) :: sCMB(nrp,nfs)
      complex(cp) :: sLM(lmP_max)
      integer :: info,i,j,l1,m1


      lm00=st_map%lm2(0,0)
      lmMin=max(lmStart,2)

      if ( .not. l_start_file ) then

         if ( lmStart <= lm00 .and. lmStop >= lm00 ) then
            if ( l_single_matrix ) then ! ps cond is required
               if ( .not. l_non_adia ) then
                  call ps_cond(s0,p0)
                  open(unit=999, file='pscond.dat')
                  do n_r=1,n_r_max
                     s(lm00,n_r)=s0(n_r)
                     p(lm00,n_r)=p0(n_r)
                     write(999,*) r(n_r), s0(n_r)*osq4pi, p0(n_r)*osq4pi,  &
                     &            osq4pi*temp0(n_r)*(s0(n_r)+alpha0(n_r)*  &
                     &            orho1(n_r)*p0(n_r)*ThExpNb*ViscHeatFac), &
                     &            osq4pi*alpha0(n_r)*ThExpNb*(-rho0(n_r)*  &
                     &            temp0(n_r)*s0(n_r)+ViscHeatFac*ogrun*    &
                     &            p0(n_r))
                  end do
                  close(999)
               end if
            else ! s cond is enough
               if ( .not. l_non_adia ) then
                  call s_cond(s0)
                  open(unit=999, file='scond.dat')
                  do n_r=1,n_r_max
                     s(lm00,n_r)=s0(n_r)
                     write(999,*) r(n_r), s0(n_r)*osq4pi
                  end do
                  close(999)
               end if
            end if
         end if

      end if

      !-- Radial dependence of perturbation in s1:
      do n_r=1,n_r_max
         x=two*r(n_r)-r_cmb-r_icb
         s1(n_r)=one-three*x**2+three*x**4-x**6
      end do

      if ( init_s1 < 100 .and. init_s1 > 0 ) then

      !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
           
         rr=random(one)
         do lm=lmMin,lmStop
            m1 = st_map%lm2m(lm)
            l1 = st_map%lm2l(lm)
            ra1=(-one+two*random(0.0_cp))*amp_s1/D_l(st_map%lm2(l1,m1))**(init_s1-1)
            ra2=(-one+two*random(0.0_cp))*amp_s1/D_l(st_map%lm2(l1,m1))**(init_s1-1)
            do n_r=1,n_r_max
               c_r=ra1*s1(n_r)
               c_i=ra2*s1(n_r)
               if ( m1 > 0 ) then  ! non axisymmetric modes
                  s(lm,n_r)=s(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               else
                  s(lm,n_r)=s(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
               end if
            end do
         end do
           
      else  if ( init_s1 >= 100 ) then

      !-- Initialize one or two modes specifically

      !----- Initialize first mode:
         l=init_s1/100
         if ( l > 99 ) l=init_s1/1000
         m=mod(init_s1,100)
         if ( l > 99 ) m=mod(init_s1,1000)
         if ( mod(m,minc) /= 0 ) then
            write(*,*) '! Wave number of mode for entropy initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m
            stop
         end if
         if ( l > l_max .or. l < m ) then
            write(*,*) '! Degree of mode for entropy initialisation'
            write(*,*) '! > l_max or < m !',l
            stop
         end if
         lm=st_map%lm2(l,m)

         if ( lmMin <= lm .and. lmStop >= lm ) then
            do n_r=1,n_r_max
               c_r=s1(n_r)*amp_s1
               s(lm,n_r)=s(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do

            write(*,'(/'' ! Entropy initialized at mode:'', &
                &  '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_s1
         end if

      !----- Initialize second mode:
         if ( init_s2 > 99 ) then
            m=mod(init_s2,100)
            if ( mod(m,minc) /= 0 ) then
               write(*,*) '! Wave number of mode for entropy initialisation'
               write(*,*) '! not compatible with phi-symmetry:',m
               stop
            end if
            l=init_s2/100
            if ( l > l_max .or. l < m ) then
               write(*,*) '! Degree of mode for entropy initialisation'
               write(*,*) '! > l_max or < m !',l
               stop
            end if

            lm=st_map%lm2(l,m)
            if ( lmMin <= lm .and. lmStop >= lm ) then
               s_r=amp_s2
               s_i=0.0_cp
               if ( amp_s2 < 0.0_cp .and. m /= 0 ) then
               !-------- Sin(phi)-mode initialized for amp_s2<0
                  s_r = 0.0_cp
                  s_i = amp_s2
               end if
               do n_r=1,n_r_max
                  c_r=s1(n_r)*s_r
                  c_i=s1(n_r)*s_i
                  s(lm,n_r)=s(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               end do
               write(6,'('' ! Second mode:'', &
                    &  '' l='',i3,'' m='',i3,'' Ampl='',f8.5/)') l,m,amp_s2
            end if

         end if

      end if

      if ( lmStart > lm00 .or. impS == 0 ) then
         !sumVal = SUM( s(lmStart:lmStop,:) )
         !write(*,"(A,2I5,2(I4,F20.16))") "END   initS: ",lmStart,lmStop, EXPONENT(REAL(sumVal)),FRACTION(REAL(sumVal)),&
         !     &EXPONENT(aimag(sumVal)),FRACTION(aimag(sumVal))
         return
      end if

      !-- Now care for the prescribed boundary condition:

      if ( minc /= 1 ) then
         write(*,*) '! impS doesnt work for minc /= 1'
         stop
      end if

      if ( abs(impS) == 1 ) then
         n_impS=2
         peakS(2)=-peakS(1)
         thetaS(2)=pi-thetaS(1)
         phiS(2)  =pi+phiS(1)
         if ( phiS(2) > 2*pi ) phiS(2)=phiS(2)-2*pi
         widthS(2)=widthS(1)
      end if

      !-- Determine the peak value vector in (xS,yS,zS) space.
      !       Then get the proportionality factors for the linear dependence
      !       of the mean (l=0,m=0) contribution on the total peak amplitude
      !       amp:
      do nS=1,n_impS

         xS(nS)=sin(thetaS(nS))*cos(phiS(nS))
         yS(nS)=sin(thetaS(nS))*sin(phiS(nS))
         zS(nS)=cos(thetaS(nS))

         nTheta=0
         do n=1,nThetaBs ! loop over the theta blocks

            nThetaStart=(n-1)*sizeThetaB+1
            do nThetaB=1,sizeThetaB
               nTheta=nTheta+1
               do nPhi=1,n_phi_max
                  xL=sinTheta(nTheta)*cos(phi(nPhi))
                  yL=sinTheta(nTheta)*sin(phi(nPhi))
                  zL=cosTheta(nTheta)
                  rH=sqrt((xS(nS)-xL)**2 + (yS(nS)-yL)**2+(zS(nS)-zL)**2)
                  !------ Opening angleL with peak value vector:
                  angleL=two*abs(asin(rH/2))
                  if ( angleL <= widthS(nS) ) then
                     sCMB(nPhi,nThetaB) = (cos(angleL/widthS(nS)*pi)+1)/2
                  else
                     sCMB(nPhi,nThetaB)=0.0_cp
                  end if
               end do
            end do
         !------ Transform to spherical hamonic space for each theta block
#ifndef WITH_SHTNS
            call fft_thetab(sCMB,-1)
            call legTF1(nThetaStart,sLM,sCMB)
#endif

         end do ! Loop over theta blocks
#ifdef WITH_SHTNS
            call spat_to_SH(sCMB, sLM)
#endif

      !--- sFac describes the linear dependence of the (l=0,m=0) mode
      !    on the amplitude peakS, SQRT(4*pi) is a normalisation factor
      !    according to the spherical harmonic function form chosen here.
         sFac(nS)=real(sLM(st_map%lm2(0,0)))*osq4pi

      end do ! Loop over peak

      !-- Value due to prescribed (l=0,m=0) contribution
      s00P=real(tops(0,0))*osq4pi
      if ( s00P == 0.0_cp .and. impS < 0 ) then
         write(*,*) '! No relative amplitudes possible!'
         write(*,*) '! for impS<0 because the mean value!'
         write(*,*) '! is zero! Refince s_top?'
         stop
      end if
      if ( impS > 0 ) s00P=one

      !-- Determine the true amplitudes amp for the peaks by solving linear system:
      !    These amplitudes guarantee that the peak as an ampliture peakS
      !    above or below the mean (l=0,m=0)
      if ( n_impS == 1 ) then
         amp(1)=peakS(1)/(s00P*(one-sFac(1)))
      else
         do j=1,n_impS
            amp(j)=-peakS(j)/s00P
            do i=1,n_impS
               if ( i == j ) then
                  mata(i,j)=sFac(i)-1
               else
                  mata(i,j)=sFac(i)
               end if
            end do
         end do
        call sgefa(mata,n_impS_max,n_impS,pivot,info)
        call sgesl(mata,n_impS_max,n_impS,pivot,amp)
      end if
      s00=0.0_cp
      do nS=1,n_impS
         s00=s00+sFac(nS)*amp(nS)
      end do

      !--- Now get the total thing so that the mean (l=0,m=0) due
      !    to the peaks is zero. The (l=0,m=0) contribution is
      !    determined (prescribed) by other means.
      nTheta=0
      do n=1,nThetaBs ! loop over the theta blocks

         nThetaStart=(n-1)*sizeThetaB+1
         do nThetaB=1,sizeThetaB
            nTheta=nTheta+1
            do nPhi=1,n_phi_max
               xL=sinTheta(nTheta)*cos(phi(nPhi))
               yL=sinTheta(nTheta)*sin(phi(nPhi))
               zL=cosTheta(nTheta)
               sCMB(nPhi,nThetaB)=-s00
               do nS=1,n_impS
                  rH=sqrt((xS(nS)-xL)**2 + (yS(nS)-yL)**2+(zS(nS)-zL)**2)
                  !------ Opening angle with peak value vector:
                  angleL=two*abs(asin(rH/2))
                  if ( angleL <= widthS(nS) )                &
                     sCMB(nPhi,nThetaB)=sCMB(nPhi,nThetaB) + &
                                        amp(nS)*(cos(angleL/widthS(nS)*pi)+1)/2
               end do
            end do
         end do
      !------ Transform to spherical hamonic space for each theta block
#ifndef WITH_SHTNS
         call fft_thetab(sCMB,-1)
         call legTF1(nThetaStart,sLM,sCMB)
#endif

      end do ! Loop over theta blocks
#ifdef WITH_SHTNS
         call spat_to_SH(sCMB, sLM)
#endif


      !--- Finally store the boundary condition and care for
      !    the fact that peakS provides the relative amplitudes
      !    in comparison to the (l=0,m=0) contribution when impS<0:
      !    Note that the (l=0,m=0) has to be determined by other means
      !    for example by setting: s_top= 0 0 -1 0
      do m=0,l_max,minc
         do l=m,l_max
            lm=st_map%lmP2(l,m)
            if ( l <= l_max .and. l > 0 ) tops(l,m)=tops(l,m)+sLM(lm)
         end do
      end do

   end subroutine initS
!---------------------------------------------------------------------------
   subroutine initXi(xi,lmStart,lmStop)
      !
      ! Purpose of this subroutine is to initialize the chemical composition
      ! according to the input control parameters.                       
      !
      ! +-----------------+---------------------------------------------+
      ! | Input           | value                                       |
      ! +=================+=============================================+
      ! | init_xi1 < 100: | random noise initialized                    | 
      ! |                 | the noise spectrum decays as l^ (init_xi1-1)|    
      ! |                 | with peak amplitude amp_xi1  for l=1        |    
      ! +-----------------+---------------------------------------------+
      ! | init_xi1 >=100: | a specific harmonic mode initialized        |  
      ! |                 | with amplitude amp_xi1.                     |        
      ! |                 | init_xi1 is interpreted as number llmm      |         
      ! |                 | where ll: harmonic degree,                  | 
      ! |                 | mm: harmonic order.                         | 
      ! +-----------------+---------------------------------------------+
      ! | init_xi2 >100 : | a second harmonic mode initialized          |  
      ! |                 | with amplitude amp_xi2.                     |         
      ! |                 | init_xi2 is again interpreted as number llmm|         
      ! |                 | where ll: harmonic degree,                  |
      ! |                 | mm: harmonic order.                         |
      ! +-----------------+---------------------------------------------+
      !                                                                   

      !-- Input variables:
      integer, intent(in) :: lmStart,lmStop

      !-- Output variables:
      complex(cp), intent(inout) :: xi(lm_max,n_r_max)

      !-- Local variables:
      integer :: n_r,lm,l,m,lm00,lmMin
      real(cp) :: x,rr,c_r,c_i,xi_r,xi_i
      real(cp) :: ra1,ra2
      real(cp) :: xi0(n_r_max),xi1(n_r_max)

      integer :: nTheta,n,nThetaStart,nThetaB,nPhi,nXi
      real(cp) :: xL,yL,zL,rH,angleL,xi00,xi00P
      real(cp) :: mata(n_impXi_max,n_impXi_max)
      real(cp) :: amp(n_impXi_max)
      integer :: pivot(n_impXi_max)
      real(cp) :: xXi(n_impXi_max),yXi(n_impXi_max)
      real(cp) :: zXi(n_impXi_max),xiFac(n_impXi_max)
      real(cp) :: xiCMB(nrp,nfs)
      complex(cp) :: xiLM(lmP_max)
      integer :: info,i,j,l1,m1


      lm00=st_map%lm2(0,0)
      lmMin=max(lmStart,2)

      if ( .not. l_start_file ) then

         if ( lmStart <= lm00 .and. lmStop >= lm00 ) then
            call xi_cond(xi0)
            open(unit=999, file='xicond.dat')
            do n_r=1,n_r_max
               xi(lm00,n_r)=xi0(n_r)
               write(999,*) r(n_r), xi0(n_r)*osq4pi
            end do
            close(999)
         end if

      end if

      !-- Radial dependence of perturbation in xi1:
      do n_r=1,n_r_max
         x=two*r(n_r)-r_cmb-r_icb
         xi1(n_r)=one-three*x**2+three*x**4-x**6
      end do

      if ( init_xi1 < 100 .and. init_xi1 > 0 ) then

      !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
           
         rr=random(one)
         do lm=lmMin,lmStop
            m1 = st_map%lm2m(lm)
            l1 = st_map%lm2l(lm)
            ra1=(-one+two*random(0.0_cp))*amp_xi1/D_l(st_map%lm2(l1,m1))**(init_xi1-1)
            ra2=(-one+two*random(0.0_cp))*amp_xi1/D_l(st_map%lm2(l1,m1))**(init_xi1-1)
            do n_r=1,n_r_max
               c_r=ra1*xi1(n_r)
               c_i=ra2*xi1(n_r)
               if ( m1 > 0 ) then  ! non axisymmetric modes
                  xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               else
                  xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
               end if
            end do
         end do
           
      else if ( init_xi1 >= 100 ) then

      !-- Initialize one or two modes specifically

      !----- Initialize first mode:
         l=init_xi1/100
         if ( l > 99 ) l=init_xi1/1000
         m=mod(init_xi1,100)
         if ( l > 99 ) m=mod(init_xi1,1000)
         if ( mod(m,minc) /= 0 ) then
            write(*,*) '! Wave number of mode for chemical composition initialisation'
            write(*,*) '! not compatible with phi-symmetry:',m
            stop
         end if
         if ( l > l_max .or. l < m ) then
            write(*,*) '! Degree of mode for chemical composition initialisation'
            write(*,*) '! > l_max or < m !',l
            stop
         end if
         lm=st_map%lm2(l,m)

         if ( lmMin <= lm .and. lmStop >= lm ) then
            do n_r=1,n_r_max
               c_r=xi1(n_r)*amp_xi1
               xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,0.0_cp,kind=cp)
            end do

            write(*,'(/'' ! Chemical composition initialized at mode:'', &
                &  '' l='',i4,'' m='',i4,'' Ampl='',f8.5)') l,m,amp_s1
         end if

      !----- Initialize second mode:
         if ( init_xi2 > 99 ) then
            m=mod(init_xi2,100)
            if ( mod(m,minc) /= 0 ) then
               write(*,*) '! Wave number of mode for chemical composition initialisation'
               write(*,*) '! not compatible with phi-symmetry:',m
               stop
            end if
            l=init_xi2/100
            if ( l > l_max .or. l < m ) then
               write(*,*) '! Degree of mode for chemical composition initialisation'
               write(*,*) '! > l_max or < m !',l
               stop
            end if

            lm=st_map%lm2(l,m)
            if ( lmMin <= lm .and. lmStop >= lm ) then
               xi_r=amp_s2
               xi_i=0.0_cp
               if ( amp_s2 < 0.0_cp .and. m /= 0 ) then
               !-------- Sin(phi)-mode initialized for amp_xi2<0
                  xi_r = 0.0_cp
                  xi_i = amp_s2
               end if
               do n_r=1,n_r_max
                  c_r=xi1(n_r)*xi_r
                  c_i=xi1(n_r)*xi_i
                  xi(lm,n_r)=xi(lm,n_r)+cmplx(c_r,c_i,kind=cp)
               end do
               write(6,'('' ! Second mode:'', &
                    &  '' l='',i3,'' m='',i3,'' Ampl='',f8.5/)') l,m,amp_xi2
            end if

         end if

      end if

      if ( lmStart > lm00 .or. impXi == 0 ) then
         return
      end if

      !-- Now care for the prescribed boundary condition:
      if ( minc /= 1 ) then
         write(*,*) '! impXi doesnt work for minc /= 1'
         stop
      end if

      if ( abs(impXi) == 1 ) then
         n_impXi=2
         peakXi(2)=-peakXi(1)
         thetaXi(2)=pi-thetaXi(1)
         phiXi(2)  =pi+phiXi(1)
         if ( phiXi(2) > 2*pi ) phiXi(2)=phiXi(2)-2*pi
         widthXi(2)=widthXi(1)
      end if

      !-- Determine the peak value vector in (xXi,yXi,zXi) space.
      !     Then get the proportionality factors for the linear dependence
      !     of the mean (l=0,m=0) contribution on the total peak amplitude
      !     amp:
      do nXi=1,n_impXi

         xXi(nXi)=sin(thetaXi(nXi))*cos(phiXi(nXi))
         yXi(nXi)=sin(thetaXi(nXi))*sin(phiXi(nXi))
         zXi(nXi)=cos(thetaXi(nXi))

         nTheta=0
         do n=1,nThetaBs ! loop over the theta blocks

            nThetaStart=(n-1)*sizeThetaB+1
            do nThetaB=1,sizeThetaB
               nTheta=nTheta+1
               do nPhi=1,n_phi_max
                  xL=sinTheta(nTheta)*cos(phi(nPhi))
                  yL=sinTheta(nTheta)*sin(phi(nPhi))
                  zL=cosTheta(nTheta)
                  rH=sqrt((xXi(nXi)-xL)**2 + (yXi(nXi)-yL)**2+(zXi(nXi)-zL)**2)
                  !------ Opening angleL with peak value vector:
                  angleL=two*abs(asin(rH/2))
                  if ( angleL <= widthXi(nXi) ) then
                     xiCMB(nPhi,nThetaB) = half*(cos(angleL/widthXi(nXi)*pi)+1)
                  else
                     xiCMB(nPhi,nThetaB)=0.0_cp
                  end if
               end do
            end do
         !------ Transform to spherical hamonic space for each theta block
#ifndef WITH_SHTNS
            call fft_thetab(xiCMB,-1)
            call legTF1(nThetaStart,xiLM,xiCMB)
#endif

         end do ! Loop over theta blocks
#ifdef WITH_SHTNS
         call spat_to_SH(xiCMB, xiLM)
#endif

      !--- xiFac describes the linear dependence of the (l=0,m=0) mode
      !    on the amplitude peakXi, sqrt(4*pi) is a normalisation factor
      !    according to the spherical harmonic function form chosen here.
         xiFac(nXi)=real(xiLM(st_map%lm2(0,0)))*osq4pi

      end do ! Loop over peak

      !-- Value due to prescribed (l=0,m=0) contribution
      xi00P=real(topxi(0,0))*osq4pi
      if ( xi00P == 0.0_cp .and. impXi < 0 ) then
         write(*,*) '! No relative amplitudes possible!'
         write(*,*) '! for impXi<0 because the mean value!'
         write(*,*) '! is zero! Refince xi_top?'
         stop
      end if
      if ( impXi > 0 ) xi00P=one

      !-- Determine the true amplitudes amp for the peaks by solving linear system:
      !    These amplitudes guarantee that the peak as an ampliture peakXi
      !    above or below the mean (l=0,m=0)
      if ( n_impXi == 1 ) then
         amp(1)=peakXi(1)/(xi00P*(one-xiFac(1)))
      else
         do j=1,n_impXi
            amp(j)=-peakXi(j)/xi00P
            do i=1,n_impXi
               if ( i == j ) then
                  mata(i,j)=xiFac(i)-1
               else
                  mata(i,j)=xiFac(i)
               end if
            end do
         end do
        call sgefa(mata,n_impXi_max,n_impXi,pivot,info)
        call sgesl(mata,n_impXi_max,n_impXi,pivot,amp)
      end if
      xi00=0.0_cp
      do nXi=1,n_impXi
         xi00=xi00+xiFac(nXi)*amp(nXi)
      end do

      !--- Now get the total thing so that the mean (l=0,m=0) due
      !    to the peaks is zero. The (l=0,m=0) contribution is
      !    determined (prescribed) by other means.
      nTheta=0
      do n=1,nThetaBs ! loop over the theta blocks

         nThetaStart=(n-1)*sizeThetaB+1
         do nThetaB=1,sizeThetaB
            nTheta=nTheta+1
            do nPhi=1,n_phi_max
               xL=sinTheta(nTheta)*cos(phi(nPhi))
               yL=sinTheta(nTheta)*sin(phi(nPhi))
               zL=cosTheta(nTheta)
               xiCMB(nPhi,nThetaB)=-xi00
               do nXi=1,n_impXi
                  rH=sqrt((xXi(nXi)-xL)**2 + (yXi(nXi)-yL)**2+(zXi(nXi)-zL)**2)
                  !------ Opening angle with peak value vector:
                  angleL=two*abs(asin(rH/2))
                  if ( angleL <= widthXi(nXi) )              &
                     xiCMB(nPhi,nThetaB)=xiCMB(nPhi,nThetaB) + &
                                        amp(nXi)*half*(cos(angleL/widthXi(nXi)*pi)+1)
               end do
            end do
         end do
      !------ Transform to spherical hamonic space for each theta block
#ifndef WITH_SHTNS
         call fft_thetab(xiCMB,-1)
         call legTF1(nThetaStart,xiLM,xiCMB)
#endif

      end do ! Loop over theta blocks
#ifdef WITH_SHTNS
      call spat_to_SH(xiCMB, xiLM)
#endif

      !--- Finally store the boundary condition and care for
      !    the fact that peakS provides the relative amplitudes
      !    in comparison to the (l=0,m=0) contribution when impS<0:
      !    Note that the (l=0,m=0) has to be determined by other means
      !    for example by setting: s_top= 0 0 -1 0
      do m=0,l_max,minc
         do l=m,l_max
            lm=st_map%lmP2(l,m)
            if ( l <= l_max .and. l > 0 ) topxi(l,m)=topxi(l,m)+xiLM(lm)
         end do
      end do

   end subroutine initXi
!---------------------------------------------------------------------------
   subroutine initB(b,aj,b_ic,aj_ic,lorentz_torque_ic,lorentz_torque_ma, &
                    lmStart,lmStop)
      !
      ! Purpose of this subroutine is to initialize the magnetic field  
      ! according to the control parameters imagcon and init_b1/2.     
      ! In addition CMB and ICB peak values are calculated for        
      ! magneto convection.                                          
      !

      !-- Input variables:
      integer, intent(in) :: lmStart,lmStop

      !-- Output variables:
      real(cp), intent(out) :: lorentz_torque_ic
      real(cp), intent(out) :: lorentz_torque_ma

      complex(cp), intent(inout) :: b(lm_maxMag,n_r_maxMag)
      complex(cp), intent(inout) :: aj(lm_maxMag,n_r_maxMag)
      complex(cp), intent(inout) :: b_ic(lm_maxMag,n_r_ic_max)
      complex(cp), intent(inout) :: aj_ic(lm_maxMag,n_r_ic_max)

      !-- Local variables:
      integer :: lm,lm0,lmStart00,l1,m1
      integer :: n_r
      real(cp) :: b_pol,b_tor
      complex(cp) :: aj0(n_r_max+1)
      complex(cp) :: aj0_ic(n_r_ic_max)
      real(cp) :: arg,aj_ic1,aj_ic2

      real(cp) :: b1(n_r_max)
      real(cp) :: b1_ic(n_r_ic_max)
      real(cp) :: bR,bI,rr
      real(cp) :: aVarCon,bVarCon
      integer :: bExp

      integer :: l1m0,l2m0,l3m0,l1m1

      lmStart00 = max(lmStart,1)

      l1m0 = st_map%lm2(1,0)
      l2m0 = st_map%lm2(2,0)
      l3m0 = st_map%lm2(3,0)
      l1m1 = st_map%lm2(1,1)

      lm0=l2m0 ! Default quadrupole field

      if ( imagcon == -1 ) then

      !----- impose l=1,m=0 poloidal field at ICB:
         lm0 = l1m0
         bpeakbot = -sqrt(third*pi)*r_icb**2*amp_b1
         bpeaktop = 0.0_cp

      else if ( imagcon == -2 ) then

      !----- impose l=1,m=0 poloidal field at CMB:
         lm0 = l1m0
         bpeakbot = 0.0_cp
         bpeaktop = -sqrt(third*pi)*r_cmb**2*amp_b1

      else if ( imagcon == 1 ) then

      !----- impose l=2,m=0 toroidal field at ICB:
         lm0 = l2m0
         bpeakbot = four*third*sqrt(pi/5.0_cp)*r_icb*amp_b1
         bpeaktop = 0.0_cp

      else if ( imagcon == 10 ) then

      !----- impose l=2,m=0 toroidal field at ICB and CMB:
         lm0 = l2m0
         bpeakbot = four*third*sqrt(pi/5.0_cp)*r_icb*amp_b1
         bpeaktop = four*third*sqrt(pi/5.0_cp)*r_cmb*amp_b1

      else if ( imagcon == 11 ) then

      !----- same as imagcon == 10 but opposite sign at CMB:
         lm0 = l2m0
         bpeakbot = four*third*sqrt(pi/5.0_cp)*r_icb*amp_b1
         bpeaktop = -four*third*sqrt(pi/5.0_cp)*r_cmb*amp_b1

      else if ( imagcon == 12 ) then

      !----- impose l=1,m=0 toroidal field at ICB and CMB:
         lm0 = l1m0
         bpeakbot = two*sqrt(third*pi)*r_icb*amp_b1
         bpeaktop = two*sqrt(third*pi)*r_cmb*amp_b1

      else if ( imagcon == 0 ) then

         lm0 = l2m0
         bpeakbot = 0.0_cp
         bpeaktop = 0.0_cp

      else if ( imagcon == -10 ) then

      !----- Test of variable conductivity case with analytical solution:
      !      Assume the magnetic diffusivity is lambda=r**5, that the aspect ratio
      !      is 0.5, and that there is no flow.
      !      The analytical stationary solution for the (l=3,m=0) toroidal field
      !      with bounday condition aj(r=r_ICB)=1, aj(r=r_CMB)=0 is then
      !      given by jVarCon(r)!
      !      A disturbed solution is used to initialize aj,
      !      the disturbance should decay with time.
      !      The solution is stored in file testVarCond.TAG at the end of the run,
      !      where the first column denotes radius, the second is aj(l=3,m=0,r) and
      !      the third is jVarCon(r). Second and third should be identical when
      !      the stationary solution has been reached.
         lm0=l3m0  ! This is l=3,m=0
         bpeakbot=one
         bpeaktop=0.0_cp
         aVarCon =-one/255.0_cp
         bVarCon =256.0_cp/255.0_cp
         if ( lmStart <= lm0 .and. lmStop >= lm0 ) then ! select processor
            do n_r=1,n_r_max             ! Diffusive toroidal field
               jVarCon(n_r)=aVarCon*r(n_r)**2 + bVarCon/(r(n_r)**6)
               aj(lm0,n_r) =jVarCon(n_r) + 0.1_cp*sin((r(n_r)-r_ICB)*pi)
            end do
         end if

      end if

      if ( init_b1 == 1 .or. imagcon > 0 ) then
      !----- Conductive solution for toroidal field,
      !      diffusion equation solved in j_cond, amplitude defined
      !      by bpeaktop and bpeakbot respectively.
      !      bpeakbot is only used for insulating inner core !
         if ( lmStart <= lm0 .and. lmStop >= lm0 ) then ! select processor
            call j_cond(lm0,aj0,aj0_ic)
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(lm0,n_r)=aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(lm0,n_r)=aj0_ic(n_r)
               end do
            end if
         end if
           
      else if ( init_b1 == 2 ) then  ! l=1,m=0 analytical toroidal field
      ! with a maximum of amp_b1 at mid-radius
      ! between r_icb and r_cmb for an insulating
      ! inner core and at r_cmb/2 for a conducting
      ! inner core

         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            b_tor=-two*amp_b1*sqrt(third*pi)  ! minus sign makes phi comp. > 0
            if ( l_cond_ic ) then
               do n_r=1,n_r_max
                  aj(l1m0,n_r)=aj(l1m0,n_r) + b_tor*r(n_r)*sin(pi*r(n_r)/r_cmb)
               end do
               do n_r=1,n_r_ic_max
                  aj_ic(l1m0,n_r)=aj_ic(l1m0,n_r) + &
                                  b_tor*r_ic(n_r)*sin(pi*r_ic(n_r)/r_cmb)
              end do
            else
               do n_r=1,n_r_max
                  aj(l1m0,n_r)=aj(l1m0,n_r) + b_tor*r(n_r)*sin(pi*(r(n_r)-r_icb))
               end do
            end if
         end if
           
      else if ( init_b1 == 3 ) then  
         ! l=2,m=0 toroidal field and l=1,m=0 poloidal field
         ! toroidal field has again its maximum of amp_b1
         ! at mid-radius between r_icb and r_cmb for an
         ! insulating inner core and at r_cmb/2 for a
         ! conducting inner core
         ! The outer core poloidal field is defined by
         ! a homogeneous  current density, its maximum at
         ! the ICB is set to amp_b1.
         ! The inner core poloidal field is chosen accordingly.
         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            b_tor=-four*third*amp_b1*sqrt(pi/5.0_cp)
            if ( l_cond_ic ) then
               b_pol=amp_b1*sqrt(three*pi)/(three+r_cmb)
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r) + &
                              b_pol*(r(n_r)**3 - four*third*r_cmb*r(n_r)**2)
               end do
               arg=pi*r_icb/r_cmb
               aj_ic1=(arg-two*sin(arg)*cos(arg)) / (arg+sin(arg)*cos(arg))
               aj_ic2=(one-aj_ic1)*r_icb*sin(arg)/cos(arg)
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2 * (    &
                                            half*r_ic(n_r)**2/r_icb + &
                                                         half*r_icb - &
                                                    four*third*r_cmb )
               end do
            else
               b_pol=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  b(l1m0,n_r)=b(l1m0,n_r)+ b_pol *     ( &
                                             r(n_r)**3 - &
                             four*third*r_cmb*r(n_r)**2 + &
                          third*r_icb**4/r(n_r)    )
               end do
            end if
         end if

         if ( lmStart <= l2m0 .and. lmStop >= l2m0 ) then ! select processor
            b_tor=-four*third*amp_b1*sqrt(pi/5.0_cp)
            if ( l_cond_ic ) then
               b_pol=amp_b1*sqrt(three*pi)/(three+r_cmb)
               do n_r=1,n_r_max
                  aj(l2m0,n_r)=aj(l2m0,n_r) + b_tor*r(n_r)*sin(pi*(r(n_r)/r_cmb))
               end do
               arg=pi*r_icb/r_cmb
               aj_ic1=(arg-two*sin(arg)*cos(arg)) / (arg+sin(arg)*cos(arg))
               aj_ic2=(one-aj_ic1)*r_icb*sin(arg)/cos(arg)
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_tor*             ( &
                           aj_ic1*r_ic(n_r)*sin(pi*r_ic(n_r)/r_cmb) + &
                                     aj_ic2*cos(pi*r_ic(n_r)/r_cmb) )
               end do
            else
               b_pol=amp_b1*sqrt(three*pi)/four
               do n_r=1,n_r_max
                  aj(l2m0,n_r)=aj(l2m0,n_r) + b_tor*r(n_r)*sin(pi*(r(n_r)-r_icb))
               end do
            end if
         end if

      else if ( init_b1 == 4 .or. imagcon == -1 ) then  ! l=1,m0 poloidal field
      ! with max field amplitude amp_b1 at r_icb
       if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
          b_pol=-amp_b1*r_icb**3*sqrt(third*pi)
          do n_r=1,n_r_max
             b(l1m0,n_r)=b(l1m0,n_r)+b_pol*or1(n_r)
          end do
          if ( l_cond_ic ) then
             do n_r=1,n_r_ic_max
                b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol/r_icb* &
                               ( -three*half + half*(r_ic(n_r)/r_icb)**2 )
             end do
          end if
       end if

      else if ( init_b1 == 5 ) then  ! l=1,m0 poloidal field
      ! constant j density, defined max field value amp_v1 at r_cmb
       if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
          if ( l_cond_ic ) then
             b_pol=amp_b1*sqrt(three*pi)/r_cmb
             do n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol* (   r(n_r)**3 - &
                            four*third*r_cmb * r(n_r)**2 )
             end do
             do n_r=1,n_r_ic_max
                b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2 * &
                   (-5.0_cp/6.0_cp*r_icb-four*third+half*r_ic(n_r)**2/r_icb)
             end do
          else
             b_pol=amp_b1*sqrt(three*pi)/(r_cmb*(one-radratio**4))
             do n_r=1,n_r_max
                b(l1m0,n_r)=b(l1m0,n_r)+b_pol* (   r(n_r)**3 - &
                                 four*third*r_cmb * r(n_r)**2 + &
                                 third*r_icb**4 / r(n_r)    )
             end do
          end if
       end if

      else if ( init_b1 == 6 ) then  ! l=1,m=0 poloidal field , constant in r !
      ! no potential at r_cmb but simple
         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            b_pol=amp_b1
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*r(n_r)**2
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol*r_icb**2
               end do
            end if
         end if

      else if ( init_b1 == 7 .or. imagcon == -2 ) then  ! l=1,m0 poloidal field
      ! which is potential field at r_cmb
         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_icb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 8 ) then  ! l=1,m0 pol. field, l=2,m=0 toroidal field
      ! which is potential field at r_cmb
         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_cmb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

         if ( lmStart <= l2m0 .and. lmStop >= l2m0 ) then ! select processor
            b_tor=amp_b1*three*half*sqrt(pi/5.0_cp)*r_icb**2*radratio
            do n_r=1,n_r_max
               aj(l2m0,n_r)=aj(l2m0,n_r)+b_tor*(r(n_r)/r_icb)**3 * &
                            ( one - (r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+b_tor * &
                                  ( one - (r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 9 ) then  ! l=2,m0 poloidal field
      ! which is potential field at r_cmb
         if ( lmStart <= l2m0 .and. lmStop >= l2m0 ) then ! select processor
            b_pol=amp_b1*7.0_cp/6.0_cp*sqrt(pi/5.0_cp)*r_icb**2*radratio
            do n_r=1,n_r_max
               b(l2m0,n_r)=b(l2m0,n_r)+b_pol*(r(n_r)/r_icb)**3 * &
                           ( one - 5.0_cp/7.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l2m0,n_r)=b_ic(l2m0,n_r)+b_pol * &
                                 ( one - 5.0_cp/7.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 10 ) then  ! only equatorial dipole

       if ( l1m1 <= 0 ) then
          write(*,*) '! Can not initialize l=1,m=1 !'
          stop
       end if

       if ( lmStart <= l1m1 .and. lmStop >= l1m1 ) then ! select processor
          b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
          do n_r=1,n_r_max
             b(l1m1,n_r)=b(l1m1,n_r)+b_pol*(r(n_r)/r_icb)**2 * &
                          ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
          end do
          if ( l_cond_ic ) then
             do n_r=1,n_r_ic_max
                b_ic(l1m1,n_r)=b_ic(l1m1,n_r)+b_pol * &
                               ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
             end do
          end if
       end if

      else if ( init_b1 < 0 ) then  ! l,m mixture, random init

         bExp=abs(init_b1)
         do n_r=1,n_r_max
            b1(n_r)=(r(n_r)/r_cmb)**2 * ( one-three/5.0_cp*(r(n_r)/r_cmb)**2 )
         end do
         if ( l_cond_ic ) then
            do n_r=1,n_r_ic_max
               b1_ic(n_r)= ( one-three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
            end do
         end if

     !-- Random noise initialization of all (l,m) modes exept (l=0,m=0):
         rr=random(one)
         do lm=lmStart00,lmStop
            l1=st_map%lm2l(lm)
            m1=st_map%lm2m(lm)
            bR=(-one+two*random(0.0_cp))*amp_b1/D_l(st_map%lm2(l1,m1))**(bExp-1)
            bI=(-one+two*random(0.0_cp))*amp_b1/D_l(st_map%lm2(l1,m1))**(bExp-1)
            if ( m1 == 0 ) bI=0.0_cp
            do n_r=1,n_r_max
               b(lm,n_r)=b(lm,n_r) + cmplx(bR*b1(n_r),bI*b1(n_r),kind=cp)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(lm,n_r)=b_ic(lm,n_r) + cmplx(bR*b1_ic(n_r), bI*b1_ic(n_r), cp)
               end do
            end if
         end do

      else if ( init_b1 == 11 ) then  ! axial and equatorial dipole

         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m0,n_r)=b(l1m0,n_r)+b_pol*(r(n_r)/r_cmb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m0,n_r)=b_ic(l1m0,n_r)+b_pol * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

         if ( l1m1 <= 0 ) then
            write(*,*) '! Cannot initialize l=1,m=1 !'
            stop
         end if

         if ( lmStart <= l1m1 .and. lmStop >= l1m1 ) then ! select processor
            b_pol=amp_b1*5.0_cp*half*sqrt(third*pi)*r_icb**2
            do n_r=1,n_r_max
               b(l1m1,n_r)=b(l1m1,n_r) +                   &
                           b_pol/10.0_cp*(r(n_r)/r_icb)**2 * &
                           ( one - three/5.0_cp*(r(n_r)/r_cmb)**2 )
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  b_ic(l1m1,n_r)=b_ic(l1m1,n_r)+b_pol/5.0_cp * &
                                 ( one - three/5.0_cp*(r_ic(n_r)/r_cmb)**2 )
               end do
            end if
         end if

      else if ( init_b1 == 21 ) then ! toroidal field created by inner core rotation
      ! equatorialy symmetric
         if ( lmStart <= l1m0 .and. lmStop >= l1m0 ) then ! select processor
            do n_r=1,n_r_max
               aj0(n_r)=amp_b1*(r_icb/r(n_r))**6
            end do
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(l1m0,n_r)=aj(l1m0,n_r)+aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l1m0,n_r)=aj_ic(l1m0,n_r)+aj0(n_r_icb)
               end do
            end if
         end if

      else if ( init_b1 == 22 ) then ! toroidal field created by inner core rotation
      ! equatorialy asymmetric
         if ( lmStart <= l2m0 .and. lmStop >= l2m0 ) then ! select processor
            do n_r=1,n_r_max
               aj0(n_r)=amp_b1*(r_icb/r(n_r))**6
            end do
            do n_r=1,n_r_max             ! Diffusive toroidal field
               aj(l2m0,n_r)=aj(l2m0,n_r)+aj0(n_r)
            end do
            if ( l_cond_ic ) then
               do n_r=1,n_r_ic_max
                  aj_ic(l2m0,n_r)=aj_ic(l2m0,n_r)+aj0(n_r_icb)
               end do
            end if
         end if

      end if

      !-- Too lazy to calculate these:
      lorentz_torque_ic=0.0_cp
      lorentz_torque_ma=0.0_cp

   end subroutine initB
!-----------------------------------------------------------------------
   subroutine j_cond(lm0, aj0, aj0_ic)
      !
      ! Purpose of this subroutine is to solve the diffusion equation    
      ! for an initial toroidal magnetic field.                          
      !                                                                   

      !-- Input variable:
      integer, intent(in) :: lm0

      !-- Output variables:
      complex(cp), intent(out) :: aj0(:)    ! aj(l=0,m=0) in the outer core
      complex(cp), intent(out) :: aj0_ic(:) ! aj(l=0,m=0) in the inner core

      !-- Local variables
      integer :: n_cheb,n_r,info,n_r_real
      complex(cp) :: rhs(n_r_tot)
      complex(cp) :: work_l(n_r_max)
      complex(cp) :: work_l_ic(n_r_ic_max)

      n_r_real = n_r_max
      if ( l_cond_ic ) n_r_real = n_r_real+n_r_ic_max

      !----- Outer core:
      do n_cheb=1,n_cheb_max
         do n_r=2,n_r_max-1
            jMat(n_r,n_cheb,1)= cheb_norm * &
              hdif_B(lm0)*dLh(lm0)*opm*lambda(n_r)*or2(n_r) * &
                      (                  d2cheb(n_cheb,n_r) + &
                            dLlambda(n_r)*dcheb(n_cheb,n_r) - &
                         dLh(lm0)*or2(n_r)*cheb(n_cheb,n_r) )
         end do
      end do
       
      !----- boundary conditions:
      !----- CMB:
      do n_cheb=1,n_cheb_max     ! should be bpeaktop at CMB
         jMat(1,n_cheb,1)=cheb_norm*cheb(n_cheb,1)
      end do
      do n_cheb=n_cheb_max+1,n_r_max
         jMat(1,n_cheb,1)=0.0_cp
      end do

      !----- ICB:
      if ( l_cond_ic ) then  ! matching condition at inner core:

         do n_cheb=1,n_cheb_max
            jMat(n_r_max,n_cheb,1)=cheb_norm*cheb(n_cheb,n_r_max)
            jMat(n_r_max+1,n_cheb,1)= cheb_norm*sigma_ratio*dcheb(n_cheb,n_r_max)
         end do
         do n_cheb=n_cheb_max+1,n_r_max
            jMat(n_r_max,n_cheb,1)=0.0_cp
            jMat(n_r_max+1,n_cheb,1)=0.0_cp
         end do

      else

         do n_cheb=1,n_cheb_max
            jMat(n_r_max,n_cheb,1)=cheb(n_cheb,n_r_max)*cheb_norm
         end do
         do n_cheb=n_cheb_max+1,n_r_max
            jMat(n_r_max,n_cheb,1)=0.0_cp
         end do

      end if
       
      do n_r=1,n_r_max
         jMat(n_r,1,1)      =half*jMat(n_r,1,1)
         jMat(n_r,n_r_max,1)=half*jMat(n_r,n_r_max,1)
      end do

      !----- Inner core:
      if ( l_cond_ic ) then

         do n_cheb=1,n_r_ic_max ! counts even IC cheb modes
            do n_r=2,n_r_ic_max ! counts IC radial grid point
               jMat(n_r_max+n_r,n_r_max+n_cheb,1) =               &
                  cheb_norm_ic*dLh(lm0)*or3(n_r_max)*opm*O_sr * ( &
                                r_ic(n_r)*d2cheb_ic(n_cheb,n_r) + &
                            two*D_lP1(lm0)*dcheb_ic(n_cheb,n_r) )
            end do
         end do

      !-------- boundary conditions:
         do n_cheb=1,n_cheb_ic_max
            jMat(n_r_max,n_r_max+n_cheb,1)=-cheb_norm_ic*cheb_ic(n_cheb,1)
            jMat(n_r_max+1,n_r_max+n_cheb,1)= -cheb_norm_ic * (   &
                                             dcheb_ic(n_cheb,1) + &
                      D_lP1(lm0)*or1(n_r_max)*cheb_ic(n_cheb,1) )
         end do
         do n_cheb=n_r_max+n_cheb_ic_max+1,n_r_tot
            jMat(n_r_max,n_cheb,1)=0.0_cp
            jMat(n_r_max+1,n_cheb,1)=0.0_cp
         end do

      !-------- normalization for lowest Cheb mode:
         do n_r=n_r_max+1,n_r_tot
            jMat(n_r,n_r_max+1,1)=half*jMat(n_r,n_r_max+1,1)
         end do

      !-------- fill matrix up with zeros:
         do n_cheb=n_r_max+1,n_r_tot
            do n_r=1,n_r_max-1
               jMat(n_r,n_cheb,1)=0.0_cp
            end do
         end do
         do n_cheb=1,n_r_max
            do n_r=n_r_max+2,n_r_tot
               jMat(n_r,n_cheb,1)=0.0_cp
            end do
         end do

      end if ! conducting inner core ?

      !----- invert matrix:
      call sgefa(jMat(:,:,1),n_r_tot,n_r_real,jPivot(:,1),info)
      if ( info /= 0 ) then
         write(*,*) 'Singular matrix jMat in j_cond.'
         stop
      end if
       
      !----- zero RHS, except BC's
      do n_r=2,n_r_real-1
         rhs(n_r)=zero
      end do
      rhs(1)= bpeaktop                             ! Outer boundary
      if ( .not. l_cond_ic ) rhs(n_r_max)=bpeakbot  ! Inner boundary
       
      !----- solve linear system:
      call cgesl(jMat(1,1,1),n_r_tot,n_r_real,jPivot(1,1),rhs)

      !----- copy result for OC:
      do n_cheb=1,n_cheb_max
         aj0(n_cheb)=rhs(n_cheb)
      end do
      do n_cheb=n_cheb_max+1,n_r_max
         aj0(n_cheb)=zero
      end do

      !----- transform to radial space:
      call chebt_oc%costf1(aj0,work_l)
       
      if ( l_cond_ic ) then
           
      !----- copy result for IC:
         do n_cheb=1,n_cheb_ic_max
            aj0_ic(n_cheb)=rhs(n_r_max+n_cheb)
         end do
         do n_cheb=n_cheb_ic_max+1,n_r_ic_max
            aj0_ic(n_cheb)=zero
         end do
           
      !----- transform to radial space:
      !  Note: this is assuming that aj0_ic is an even function !
         call chebt_ic%costf1(aj0_ic,work_l_ic)

      end if
    
   end subroutine j_cond
!--------------------------------------------------------------------------------
   subroutine s_cond(s0)
      !
      ! Purpose of this subroutine is to solve the entropy equation      
      ! for an the conductive (l=0,m=0)-mode.                            
      ! Output is the radial dependence of the solution in s0.           
      !

      real(cp), intent(out) :: s0(:) ! spherically-symmetric part

      !-- local variables:
      integer :: n_cheb,n_r,info
      real(cp) :: rhs(n_r_max)
      real(cp) :: work(n_r_max)

      !-- Set Matrix:
      if ( l_anelastic_liquid ) then
         do n_cheb=1,n_r_max
            do n_r=2,n_r_max-1
               s0Mat(n_r,n_cheb)=cheb_norm*opr*kappa(n_r)* ( d2cheb(n_cheb,n_r) + &
              ( two*or1(n_r)+beta(n_r)+dLkappa(n_r) )                             &
                                                            * dcheb(n_cheb,n_r)  )
            end do
         end do
      else
         do n_cheb=1,n_r_max
            do n_r=2,n_r_max-1
               s0Mat(n_r,n_cheb)=cheb_norm*opr*kappa(n_r)* ( d2cheb(n_cheb,n_r) + &
              ( two*or1(n_r)+beta(n_r)+dLtemp0(n_r)+                              &
                                              dLkappa(n_r) )* dcheb(n_cheb,n_r)  )
            end do
         end do
      end if
       

      !-- Set boundary conditions:
      do n_cheb=1,n_cheb_max
         if ( ktops == 1 .or. kbots == 2 ) then
            s0Mat(1,n_cheb)=cheb_norm
         else
            s0Mat(1,n_cheb)=dcheb(n_cheb,1)*cheb_norm
         end if
         if ( kbots == 1 ) then   ! Constant entropy at ICB
            s0Mat(n_r_max,n_cheb)=cheb(n_cheb,n_r_max)*cheb_norm
         else                     ! Constant heat flux at ICB
            s0Mat(n_r_max,n_cheb)=dcheb(n_cheb,n_r_max)*cheb_norm
         end if
      end do
       
      !-- Fill with zeros:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            s0Mat(1,n_cheb)      =0.0_cp
            s0Mat(n_r_max,n_cheb)=0.0_cp
         end do
      end if
       
      !-- Renormalize:
      do n_r=1,n_r_max
         s0Mat(n_r,1)=half*s0Mat(n_r,1)
         s0Mat(n_r,n_r_max)=half*s0Mat(n_r,n_r_max)
      end do
       
      !-- Invert matrix:
      call sgefa(s0Mat,n_r_max,n_r_max,s0Pivot,info)
      if ( info /= 0 ) then
         write(*,*) '! Singular Matrix s0Mat in init_s!'
         stop
      end if
       
      !-- Set source terms in RHS:
      do n_r=2,n_r_max-1
         rhs(n_r)=-epsc*epscProf(n_r)*orho1(n_r)
      end do
       
      !-- Set boundary values:
      if ( ktops == 2 .and. kbots == 2 ) then
         rhs(1)=0.0_cp
      else
         rhs(1)=real(tops(0,0))
      end if
      rhs(n_r_max)=real(bots(0,0))
       
      !-- Solve for s0:
      call sgesl(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)
       
      !-- Copy result to s0:
      do n_r=1,n_r_max
         s0(n_r)=rhs(n_r)
      end do

      !-- Set cheb-modes > n_cheb_max to zero:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            s0(n_cheb)=0.0_cp
         end do
      end if
       
      !-- Transform to radial space:
      call chebt_oc%costf1(s0,work)

   end subroutine s_cond
!--------------------------------------------------------------------------------
   subroutine xi_cond(xi0)
      !
      ! Purpose of this subroutine is to solve the chemical composition equation      
      ! for an the conductive (l=0,m=0)-mode.                            
      ! Output is the radial dependence of the solution in s0.           
      !

      real(cp), intent(out) :: xi0(:) ! spherically-symmetric part

      !-- local variables:
      integer :: n_cheb,n_r,info
      real(cp) :: rhs(n_r_max)
      real(cp) :: work(n_r_max)

      !-- Set Matrix:
      do n_cheb=1,n_r_max
         do n_r=2,n_r_max-1
            xi0Mat(n_r,n_cheb)=cheb_norm*osc*( d2cheb(n_cheb,n_r) + &
           ( two*or1(n_r)+beta(n_r) )*          dcheb(n_cheb,n_r)  )
         end do
      end do
       

      !-- Set boundary conditions:
      do n_cheb=1,n_cheb_max
         if ( ktopxi == 1 .or. kbotxi == 2 ) then
            xi0Mat(1,n_cheb)=cheb_norm
         else
            xi0Mat(1,n_cheb)=dcheb(n_cheb,1)*cheb_norm
         end if
         if ( kbotxi == 1 ) then
            xi0Mat(n_r_max,n_cheb)=cheb(n_cheb,n_r_max)*cheb_norm
         else
            xi0Mat(n_r_max,n_cheb)=dcheb(n_cheb,n_r_max)*cheb_norm
         end if
      end do
       
      !-- Fill with zeros:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            xi0Mat(1,n_cheb)      =0.0_cp
            xi0Mat(n_r_max,n_cheb)=0.0_cp
         end do
      end if
       
      !-- Renormalize:
      do n_r=1,n_r_max
         xi0Mat(n_r,1)=half*xi0Mat(n_r,1)
         xi0Mat(n_r,n_r_max)=half*xi0Mat(n_r,n_r_max)
      end do
       
      !-- Invert matrix:
      call sgefa(xi0Mat,n_r_max,n_r_max,xi0Pivot,info)
      if ( info /= 0 ) then
         write(*,*) '! Singular Matrix xi0Mat in init_xi!'
         stop
      end if
       
      !-- Set source terms in RHS:
      do n_r=2,n_r_max-1
         rhs(n_r)=-epscxi
      end do
       
      !-- Set boundary values:
      if ( ktopxi == 2 .and. kbotxi == 2 ) then
         rhs(1)=0.0_cp
      else
         rhs(1)=real(topxi(0,0))
      end if
      rhs(n_r_max)=real(botxi(0,0))
       
      !-- Solve for s0:
      call sgesl(xi0Mat,n_r_max,n_r_max,xi0Pivot,rhs)
       
      !-- Copy result to s0:
      do n_r=1,n_r_max
         xi0(n_r)=rhs(n_r)
      end do

      !-- Set cheb-modes > n_cheb_max to zero:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            xi0(n_cheb)=0.0_cp
         end do
      end if
       
      !-- Transform to radial space:
      call chebt_oc%costf1(xi0,work)

   end subroutine xi_cond
!--------------------------------------------------------------------------------
   subroutine ps_cond(s0,p0)
      !
      ! Purpose of this subroutine is to solve the entropy equation      
      ! for an the conductive (l=0,m=0)-mode.                            
      ! Output is the radial dependence of the solution in s0 and p0.
      !

      real(cp), intent(out) :: s0(:) ! spherically-symmetric part
      real(cp), intent(out) :: p0(:) ! spherically-symmetric part

      !-- local variables:
      integer :: n_cheb,nCheb_p,n_r,n_r_p,info,n_cheb_in
      real(cp) :: rhs(2*n_r_max)
      real(cp) :: work(n_r_max), work2(n_r_max), work3(n_r_max)

      if ( l_temperature_diff ) then

         do n_cheb=1,n_r_max
            nCheb_p=n_cheb+n_r_max
            do n_r=1,n_r_max
               n_r_p=n_r+n_r_max

               ! Delta T = epsc
               ps0Mat(n_r,n_cheb)=cheb_norm*opr*kappa(n_r)* (                &
                  &                                    d2cheb(n_cheb,n_r) +  &
                  &      ( beta(n_r)+two*dLtemp0(n_r)+                       &
                  &        two*or1(n_r)+dLkappa(n_r) )* dcheb(n_cheb,n_r) +  &
                  &      ( ddLtemp0(n_r)+dLtemp0(n_r)*(                      &
                  &  two*or1(n_r)+dLkappa(n_r)+dLtemp0(n_r)+beta(n_r) ) ) *  &
                  &                                      cheb(n_cheb,n_r) ) 

               ps0Mat(n_r,nCheb_p)=cheb_norm*opr*kappa(n_r)*                  &
                  &       alpha0(n_r)*orho1(n_r)*ViscHeatFac*ThExpNb*(        &
                  &                                    d2cheb(n_cheb,n_r) +   &
                  &      ( dLkappa(n_r)+two*(dLalpha0(n_r)+dLtemp0(n_r)) -    &
                  &        beta(n_r) +two*or1(n_r) ) *  dcheb(n_cheb,n_r) +   &
                  & ( (dLkappa(n_r)+dLalpha0(n_r)+dLtemp0(n_r)+two*or1(n_r)) *&
                  &        (dLalpha0(n_r)+dLtemp0(n_r)-beta(n_r)) +           &
                  &        ddLalpha0(n_r)+ddLtemp0(n_r)-dbeta(n_r) ) *        &
                  &                                      cheb(n_cheb,n_r) )

               ! Hydrostatic equilibrium
               ps0Mat(n_r_p,n_cheb) = -cheb_norm*rho0(n_r)*BuoFac*rgrav(n_r)*&
                  &                                 cheb(n_cheb,n_r)
               ps0Mat(n_r_p,nCheb_p)= cheb_norm *( dcheb(n_cheb,n_r)- &
                  &                       beta(n_r)*cheb(n_cheb,n_r) )

            end do
         end do

      else ! entropy diffusion

         do n_cheb=1,n_r_max
            nCheb_p=n_cheb+n_r_max
            do n_r=1,n_r_max
               n_r_p=n_r+n_r_max

               ! Delta T = epsc
               ps0Mat(n_r,n_cheb)=cheb_norm*opr*kappa(n_r)* (                &
                  &                                    d2cheb(n_cheb,n_r) +  &
                  &      ( beta(n_r)+dLtemp0(n_r)+                           &
                  &        two*or1(n_r)+dLkappa(n_r) )* dcheb(n_cheb,n_r) )
               ps0Mat(n_r,nCheb_p)  =0.0_cp

               ! Hydrostatic equilibrium
               ps0Mat(n_r_p,n_cheb) = -cheb_norm*rho0(n_r)*BuoFac*rgrav(n_r)* &
                                                     cheb(n_cheb,n_r)
               ps0Mat(n_r_p,nCheb_p)= cheb_norm *(  dcheb(n_cheb,n_r)- &
               &                           beta(n_r)*cheb(n_cheb,n_r) )

            end do
         end do

      end if


      !-- Set boundary conditions:
      do n_cheb=1,n_cheb_max
         nCheb_p=n_cheb+n_r_max
         if ( ktops == 1 .or. kbots == 2 .or. kbots == 4 ) then
            ps0Mat(1,n_cheb)  =cheb_norm
            ps0Mat(1,nCheb_p) =0.0_cp
         else if ( ktops == 2) then ! constant entropy flux at CMB
            ps0Mat(1,n_cheb) =dcheb(n_cheb,1)*cheb_norm
            ps0Mat(1,nCheb_p)=0.0_cp
         else if ( ktops == 3) then ! constant temperature at CMB
            ps0Mat(1,n_cheb) =cheb_norm*temp0(1)
            ps0Mat(1,nCheb_p)=cheb_norm*alpha0(1)*temp0(1)*orho1(1)* &
            &                 ViscHeatFac*ThExpNb
         else if ( ktops == 4) then ! constant temperature flux at CMB
            ps0Mat(1,n_cheb)  =cheb_norm*temp0(1)*( dcheb(n_cheb,1)+ &
              &                         dLtemp0(1)*cheb(n_cheb,1) )
            ps0Mat(1,nCheb_p)=cheb_norm*orho1(1)*alpha0(1)*   &
              &              temp0(1)*ViscHeatFac*ThExpNb*(   &
              &              dcheb(n_cheb,1)+(dLalpha0(1)+    &
              &              dLtemp0(1)-beta(1))*cheb(n_cheb,1) )
         end if

         if ( kbots == 1 ) then        ! Constant entropy at ICB
            ps0Mat(n_r_max,n_cheb) =cheb(n_cheb,n_r_max)*cheb_norm
            ps0Mat(n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 2 ) then   ! Constant entropy flux at ICB
            ps0Mat(n_r_max,n_cheb) =dcheb(n_cheb,n_r_max)*cheb_norm
            ps0Mat(n_r_max,nCheb_p)=0.0_cp
         else if ( kbots == 3 ) then   ! Constant temperature at ICB
            ps0Mat(n_r_max,n_cheb) =cheb_norm*cheb(n_cheb,n_r_max)*temp0(n_r_max)
            ps0Mat(n_r_max,nCheb_p)=cheb_norm*cheb(n_cheb,n_r_max)*  &
              &                     alpha0(n_r_max)*temp0(n_r_max)*  &
              &                     orho1(n_r_max)*ViscHeatFac*ThExpNb
         else if ( kbots == 4 ) then   ! Constant temperature flux at ICB
            ps0Mat(n_r_max,n_cheb)  =cheb_norm*temp0(n_r_max)*(           &
              &                                    dcheb(n_cheb,n_r_max)+ &
              &                     dLtemp0(n_r_max)*cheb(n_cheb,n_r_max) )
            ps0Mat(n_r_max,nCheb_p)=cheb_norm*orho1(n_r_max)*alpha0(n_r_max)* &
              &                      temp0(n_r_max)*ViscHeatFac*ThExpNb*(     &
              &                                    dcheb(n_cheb,n_r_max)+     &
              &                      (dLalpha0(n_r_max)+dLtemp0(n_r_max)-     &
              &                       beta(n_r_max))*cheb(n_cheb,n_r_max) )
         end if

         !-- Boundary condition on spherically-symmetric pressure
         if ( ViscHeatFac*ThExpNb == 0.0_cp ) then ! No feedback of density on pressure
            ps0Mat(n_r_max+1,nCheb_p)=cheb_norm
         else
            ps0Mat(n_r_max+1,nCheb_p)=0.0_cp
         end if
         ps0Mat(n_r_max+1,n_cheb) =0.0_cp
         ps0Mat(2*n_r_max,n_cheb) =0.0_cp
         ps0Mat(2*n_r_max,nCheb_p)=0.0_cp

      end do

      ! In case density perturbations feed back on pressure (non-Boussinesq)
      ! Impose that the integral of (rho' r^2) vanishes
      if ( ViscHeatFac*ThExpNb /= 0.0_cp ) then

         work(:)=ThExpNb*ViscHeatFac*ogrun*alpha0(:)*r(:)*r(:)
         call chebt_oc%costf1(work,work2)
         work         =work*cheb_norm
         work(1)      =half*work(1)
         work(n_r_max)=half*work(n_r_max)

         work2(:)=-ThExpNb*alpha0(:)*temp0(:)*rho0(:)*r(:)*r(:) 
         call chebt_oc%costf1(work2,work3)
         work2         =work2*cheb_norm
         work2(1)      =half*work2(1)
         work2(n_r_max)=half*work2(n_r_max)

         do n_cheb=1,n_cheb_max
            nCheb_p=n_cheb+n_r_max
            ps0Mat(n_r_max+1,nCheb_p)=0.0_cp
            ps0Mat(n_r_max+1,n_cheb)=0.0_cp
            do n_cheb_in=1,n_cheb_max
               if (mod(n_cheb+n_cheb_in-2,2)==0) then
               ps0Mat(n_r_max+1,nCheb_p)=ps0Mat(n_r_max+1,nCheb_p)+ &
               &                       (one/(one-real(n_cheb_in-n_cheb,cp)**2)+&
               &                       one/(one-real(n_cheb_in+n_cheb-2,cp)**2))*&
               &                       work(n_cheb_in)*half*cheb_norm
               ps0Mat(n_r_max+1,n_cheb)=ps0Mat(n_r_max+1,n_cheb)+ &
               &                       (one/(one-real(n_cheb_in-n_cheb,cp)**2)+&
               &                       one/(one-real(n_cheb_in+n_cheb-2,cp)**2))*&
               &                       work2(n_cheb_in)*half*cheb_norm
               end if
            end do
         end do

      end if
       
      !-- Fill with zeros:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            nCheb_p=n_cheb+n_r_max
            ps0Mat(1,n_cheb)         =0.0_cp
            ps0Mat(n_r_max,n_cheb)   =0.0_cp
            ps0Mat(n_r_max+1,n_cheb) =0.0_cp
            ps0Mat(2*n_r_max,n_cheb) =0.0_cp
            ps0Mat(1,nCheb_p)        =0.0_cp
            ps0Mat(n_r_max,nCheb_p)  =0.0_cp
            ps0Mat(n_r_max+1,nCheb_p)=0.0_cp
         end do
      end if
       
      !-- Renormalize:
      do n_r=1,n_r_max
         n_r_p=n_r+n_r_max
         ps0Mat(n_r,1)          =half*ps0Mat(n_r,1)
         ps0Mat(n_r,n_r_max)    =half*ps0Mat(n_r,n_r_max)
         ps0Mat(n_r,n_r_max+1)  =half*ps0Mat(n_r,n_r_max+1)
         ps0Mat(n_r,2*n_r_max)  =half*ps0Mat(n_r,2*n_r_max)
         ps0Mat(n_r_p,1)        =half*ps0Mat(n_r_p,1)
         ps0Mat(n_r_p,n_r_max)  =half*ps0Mat(n_r_p,n_r_max)
         ps0Mat(n_r_p,n_r_max+1)=half*ps0Mat(n_r_p,n_r_max+1)
         ps0Mat(n_r_p,2*n_r_max)=half*ps0Mat(n_r_p,2*n_r_max)
      end do

      !-- Invert matrix:
      call sgefa(ps0Mat,2*n_r_max,2*n_r_max,ps0Pivot,info)
      if ( info /= 0 ) then
         write(*,*) '! Singular Matrix ps0Mat in ps_cond!'
         stop
      end if
       
      !-- Set source terms in RHS:
      do n_r=1,n_r_max
         rhs(n_r)          =-epsc*epscProf(n_r)*orho1(n_r)
         rhs(n_r+n_r_max)  =0.0_cp
      end do
       
      !-- Set boundary values:
      if ( (ktops==2 .and. kbots==2) .or. (ktops==4 .and. kbots==4) ) then
         rhs(1)=0.0_cp
      else
         rhs(1)=real(tops(0,0))
      end if
      rhs(n_r_max)=real(bots(0,0))

      !-- Pressure at the top boundary
      rhs(n_r_max+1)=0.0_cp

      !-- Solve for s0:
      call sgesl(ps0Mat,2*n_r_max,2*n_r_max,ps0Pivot,rhs)
       
      !-- Copy result to s0:
      do n_r=1,n_r_max
         s0(n_r)=rhs(n_r)
         p0(n_r)=rhs(n_r+n_r_max)
      end do

      !-- Set cheb-modes > n_cheb_max to zero:
      if ( n_cheb_max < n_r_max ) then
         do n_cheb=n_cheb_max+1,n_r_max
            s0(n_cheb)=0.0_cp
            p0(n_cheb)=0.0_cp
         end do
      end if
       
      !-- Transform to radial space:
      call chebt_oc%costf1(s0,work)
      call chebt_oc%costf1(p0,work)

   end subroutine ps_cond
!--------------------------------------------------------------------------------
end module init_fields
