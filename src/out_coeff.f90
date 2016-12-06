module out_coeff
   !
   ! This module contains the subroutines that calculate the Bcmb files
   ! and the [B|V|T]_coeff_r files
   !
  
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use logic, only: l_r_field, l_cmb_field, l_save_out, l_average
   use parallel_mod, only: rank
   use blocking, only: lm2
   use truncation, only: lm_max, l_max, minc
   use communications, only: gather_from_lo_to_rank0
   use LMLoop_data, only: llm, ulm

   implicit none

   private

   complex(cp), allocatable :: work(:) ! work array for r_field

   public :: write_Bcmb, write_coeff_r, initialize_coeffs, finalize_coeffs

contains

   subroutine initialize_coeffs()

      if ( l_r_field .or. l_cmb_field .or. l_average ) then
         allocate ( work(lm_max) )
         bytes_allocated = bytes_allocated+lm_max*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_coeffs
!-------------------------------------------------------------------------------
   subroutine finalize_coeffs()

      if ( l_r_field .or. l_cmb_field .or. l_average ) then
         deallocate( work )
      end if

   end subroutine finalize_coeffs
!-------------------------------------------------------------------------------
   subroutine write_Bcmb(time,b_LMloc,l_max_cmb,n_cmb_sets,cmb_file,n_cmb_file)
      !
      ! Each call of this subroutine writes time and the poloidal magnetic
      ! potential coefficients b at the CMB up to degree and order        
      ! l_max_cmb at the end of output file $cmb_file.                    
      ! The parameters l_max_cmb, minc and the number of stored coeffs.  
      ! are written into the first line of $cmb_file.                     
      ! Each further set contains:                                        
      !    
      !     .. code-block:: fortran
      !
      !                    time,
      !                    real(b(l=0,m=0)),imag(b(l=0,m=0)),                  
      !                    real(b(l=1,m=0)),imag(b(l=1,m=0)),    
      !                                                                   
      ! Real and imaginary part of b(*) for all orders m<=l are written   
      ! for a specific degree l, then for the degrees l+1, l+2, l_max_cmb.
      !

      !-- Input variables:
      real(cp),         intent(in) :: time             ! Time
      complex(cp),      intent(in) :: b_LMloc(llm:ulm) ! Poloidal field potential
      character(len=*), intent(in) :: cmb_file         ! Name of output file

      !-- Output variables:
      integer, intent(inout) :: n_cmb_file    ! Output unit for $cmb_file
      integer, intent(inout) :: l_max_cmb     ! Max degree of output
      integer, intent(inout) :: n_cmb_sets    ! Total no. of cmb sets,
      ! should be set to zero before first call

      !-- Local variables:
      integer :: n,n_out        ! counters
      integer :: l,m            ! degree and order
      integer :: lm             ! position of (l,m) in b(*,n_r_cmb)
      integer :: m_max_cmb      ! Max order of output
      integer :: lm_max_cmb     ! Max no of combinations l,m for output
      integer :: n_data         ! No of output data
      integer :: n_r_cmb        ! Position of cmb-radius on grid

      real(cp), allocatable ::  out(:) ! Output array

      call gather_from_lo_to_rank0(b_LMloc, work)

      if ( rank == 0 ) then

         !--- Definition of max degree for output
         if ( l_max < l_max_cmb ) l_max_cmb=l_max

         !--- Define postition of CMB on radial grid:
         n_r_cmb=1

         !--- Calculate no of data for l_max_cmb:
         m_max_cmb=(l_max_cmb/minc)*minc
         lm_max_cmb= m_max_cmb*(l_max_cmb+1)/minc &
              &     -m_max_cmb*(m_max_cmb-minc)/(2*minc) &
              &     +l_max_cmb-m_max_cmb+1
         n_data=2*lm_max_cmb-l_max_cmb-2

         allocate(out(n_data))

         !--- Increase no. of cmb_sets:
         n_cmb_sets=n_cmb_sets+1

         !--- Open output file name:
         if ( l_save_out .or. n_cmb_sets == 0 ) then
            open(newunit=n_cmb_file, file=cmb_file, position='append', &
            &    form='unformatted')
         end if

         !--- If this is the first set write l_max_cmb and minc into file:
         if ( n_cmb_sets <= 1 ) write(n_cmb_file) l_max_cmb,minc,n_data

         !--- Write b(*) into output array out(*):
         n_out=0

         !--- Axisymmetric part: (m=0) only real part stored
         do l=1,l_max_cmb
            lm=lm2(l,0)
            n_out=n_out+1
            out(n_out)=real(work(lm))
         end do

         !--- Non-axisymmetric part: store real and imag part
         do m=minc,l_max_cmb,minc
            do l=m,l_max_cmb
               lm=lm2(l,m)
               n_out=n_out+1
               out(n_out)=real(work(lm))
               n_out=n_out+1
               out(n_out)=aimag(work(lm))
            end do
         end do

         !--- Finally write output array out(*) into cmb_file:
         write(n_cmb_file) time,(out(n),n=1,n_out)

         !--- Close cmb_file
         if ( l_save_out .or. n_cmb_sets == 0 ) then
            close(n_cmb_file)
         end if

         deallocate(out)

      end if

   end subroutine write_Bcmb
!----------------------------------------------------------------------
   subroutine write_coeff_r(time,w_LMloc,dw_LMloc,ddw_LMloc,z_LMloc,r,  &
      &                     l_max_r,n_sets,file,n_file,nVBS)
      !
      ! Each call of this subroutine writes time and the poloidal and     
      ! toroidal coeffitients w,dw,z at a specific radius up to degree    
      ! and order l_max_r at the end of output file $file.                
      ! If the input is magnetic field (nVBS=2) ddw is stored as well.    
      ! If the input is entropy field (nVBS=3) only ddw=s is stored.      
      ! The parameters l_max_r, minc, the number of stored coeffs and     
      ! radius in the outer core are written into the first line of $file.
      ! Each further set contains:                                        
      !
      !     .. code-block :: fortran
      !
      !               time,
      !               real(w(l=0,m=0)),imag(w(l=0,m=0)),                  
      !               real(w(l=1,m=0)),imag(w(l=1,m=0)),                  
      !               ...
      !               real(dw(l=0,m=0)),imag(dw(l=0,m=0)),                
      !               real(dw(l=1,m=0)),imag(dw(l=1,m=0)),                
      !               ...
      !               real(z(l=0,m=0)),imag(z(l=0,m=0)),                  
      !               real(z(l=1,m=0)),imag(z(l=1,m=0)),                  
      !               ...
      !               real(ddw(l=0,m=0)),imag(ddw(l=0,m=0)),              
      !               real(ddw(l=1,m=0)),imag(ddw(l=1,m=0)),              
      !                                                                   
      ! Real and imaginary part of w(*) for all orders m<=l are written   
      ! for a specific degree l, then for the degrees l+1, l+2, l_max_r.  
      !                                                                   

      !-- Input variables:
      real(cp),         intent(in) :: r                 ! radius of coeffs
      real(cp),         intent(in) :: time              ! Time
      complex(cp),      intent(in) :: w_LMloc(llm:ulm)  ! Poloidal field potential
      complex(cp),      intent(in) :: dw_LMloc(llm:ulm) ! dr of Poloidal field potential
      complex(cp),      intent(in) :: ddw_LMloc(llm:ulm)! dr^2 of Poloidal field potential
      complex(cp),      intent(in) :: z_LMloc(llm:ulm)   ! Toroidal field potential
      character(len=*), intent(in) :: file         ! Name of output file
      integer,          intent(in) :: nVBS         ! True if output is flow

      !-- Output:
      integer, intent(inout) :: n_file      ! Output unit for $file
      integer, intent(inout) :: l_max_r     ! Max degree of output
      integer, intent(inout) :: n_sets      ! Total no. of cmb sets,
      ! should be set to zero before first call

      !-- Local variables:
      integer :: n,n_out        ! counter
      integer :: l,m            ! degree and order
      integer :: lm             ! position of (l,m) in v(*),..
      integer :: m_max_r        ! Max order of output
      integer :: lm_max_r       ! Max no of combinations l,m for output
      integer :: n_data         ! No of output data

      real(cp), allocatable ::  out(:)! Output array

      !--- Definition of max degree for output
      if ( l_max < l_max_r ) l_max_r=l_max

      !--- Calculate no of data for l_max_r:
      m_max_r=(l_max_r/minc)*minc
      lm_max_r=m_max_r*(l_max_r+1)/minc- &
           m_max_r*(m_max_r-minc)/(2*minc) + &
           l_max_r-m_max_r+1
      n_data=2*lm_max_r-l_max_r-2

      if ( rank == 0 ) then
         if ( nVBS == 1 ) then
            allocate(out(3*n_data))
         else if ( nVBS == 2 ) then
            allocate(out(4*n_data))
         else if ( nVBS == 3 ) then
            allocate(out(n_data+1))
         end if
      end if

      !--- Increase no. of sets:
      n_sets=n_sets+1

      !--- Write b(*) into output array out(*):
      n_out=0

      call gather_from_lo_to_rank0(w_LMloc, work)

      if ( rank == 0 ) then
         if ( nVBS == 3 ) then
            !--- Axisymmetric part of s: (m=0) only real part stored
            do l=0,l_max ! start with l=0
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(work(lm))
               end if
            end do
         else
            !--- Axisymmetric part of w: (m=0) only real part stored
            do l=1,l_max ! start with l=1
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(work(lm))
               end if
            end do
         end if

         !--- Non-axisymmetric part of w: store real and imag part
         do m=minc,l_max_r,minc
            do l=m,l_max
               lm=lm2(l,m)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(work(lm))
                  n_out=n_out+1
                  out(n_out)=aimag(work(lm))
               end if
            end do
         end do
      end if

      if ( nVBS /= 3 ) then

         call gather_from_lo_to_rank0(dw_LMloc, work)

         if ( rank == 0 ) then
            !-- Now output for flow or magnetic field only:
            !--- Axisymmetric part of dw: (m=0) only real part stored
            do l=1,l_max
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(work(lm))
               end if
            end do
            !--- Non-axisymmetric part of dv: store real and imag part
            do m=minc,l_max_r,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  if ( l <= l_max_r ) then
                     n_out=n_out+1
                     out(n_out)=real(work(lm))
                     n_out=n_out+1
                     out(n_out)=aimag(work(lm))
                  end if
               end do
            end do
         end if

         call gather_from_lo_to_rank0(z_LMloc, work)

         if ( rank == 0 ) then
            !--- Axisymmetric part of z: (m=0) only real part stored
            do l=1,l_max
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(work(lm))
               end if
            end do
            !--- Non-axisymmetric part of z: store real and imag part
            do m=minc,l_max_r,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  if ( l <= l_max_r ) then
                     n_out=n_out+1
                     out(n_out)=real(work(lm))
                     n_out=n_out+1
                     out(n_out)=aimag(work(lm))
                  end if
               end do
            end do
         end if
      end if

      !--- If this is a magnetic field I also store the second radial derivative
      !    of the poloidal potential to caluclate diffusion:
      if ( nVBS == 2 ) then

         call gather_from_lo_to_rank0(ddw_LMloc, work)

         if ( rank == 0 ) then
            !--- Axisymmetric part of ddw: (m=0) only real part stored
            do l=1,l_max
               lm=lm2(l,0)
               if ( l <= l_max_r ) then
                  n_out=n_out+1
                  out(n_out)=real(work(lm))
               end if
            end do
            !--- Non-axisymmetric part of ddw: store real and imag part
            do m=minc,l_max_r,minc
               do l=m,l_max
                  lm=lm2(l,m)
                  if ( l <= l_max_r ) then
                     n_out=n_out+1
                     out(n_out)=real(work(lm))
                     n_out=n_out+1
                     out(n_out)=aimag(work(lm))
                  end if
               end do
            end do
         end if 
      end if

      if ( rank == 0 ) then

         !--- Open output file with name $file:
         if ( l_save_out ) then
            open(newunit=n_file, file=file, form='unformatted', status='unknown', &
                 position='append')
         end if

         !--- If this is the first set write, l_max_r and minc into first line:
         if ( n_sets == 1 ) then
            write(n_file) l_max_r,minc,n_data,r
         end if


         !--- Finally write output array out(*) into file:
         write(n_file) time,(out(n),n=1,n_out)

         !--- Close file
         if ( l_save_out ) then
            close(n_file)
         end if

         deallocate(out)

      end if


   end subroutine write_coeff_r
!-------------------------------------------------------------------------------
end module out_coeff
