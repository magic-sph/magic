!$Id$
module fft_fac_mod

   use precision_mod, only: cp
   use const, only: sin36, sin60, sin72, cos36, cos72

   implicit none

   private

   public :: fft_fac_complex, fft_fac_real

contains

   subroutine fft_fac_real(a,b,c,d,trigs,nv,l1,l2,n,ifac,la)
      !---------------------------------------------------------------------------
      !     main part of Fourier / Chebychev transform
      !     called in costf1, costf2
      !---------------------------------------------------------------------------
    
      !-- Input variables:
      integer,  intent(in) :: nv,l1,l2,n,ifac,la
      real(cp), intent(in) :: a(*),b(*)
      real(cp), intent(in) :: trigs(2*n)
    
      !-- Output variables
      real(cp), intent(out) :: c(*),d(*)
    
      !-- Local variables:
      integer :: i,ia,ib,ic,id,ie
      integer :: j,ja,jb,jc,jd,je
      integer :: istart,istop,iink
      integer :: jump,jadd,jink
      integer :: k,kb,kc,kd,ke
      integer :: nv2
      integer :: l,m,lm1,lm2,ll,la1
    
      real(cp) :: c1,c2,c3,c4
      real(cp) :: s1,s2,s3,s4
    
      m    =n/ifac
      iink =m*2*nv
      jink =la*2*nv
      jump =(ifac-1)*jink
      nv2  =2*nv
      lm1  =l1-1
      lm2  =l2-l1
    
      if ( ifac == 2 ) then
    
         ia  =1
         ja  =1
         ib  =ia+iink
         jb  =ja+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop=istart+lm2
            do i=istart,istop
               c(ja+i)=a(ia+i)+a(ib+i)
               d(ja+i)=b(ia+i)+b(ib+i)
               c(jb+i)=a(ia+i)-a(ib+i)
               d(jb+i)=b(ia+i)-b(ib+i)
            end do
         end do
    
         if ( la == m ) return
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb  =k+k-2
            c1  =trigs(kb+1)
            s1  =trigs(kb+2)
            do l=1,la
               ll    =k+l-1
               istart=(ll-1)*nv2+lm1
               istop =istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=a(ia+i)+a(ib+i)
                  d(ja+j)=b(ia+i)+b(ib+i)
                  c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
                  d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
               end do
            end do
         end do
    
      else if ( ifac == 3 ) then
    
         ia  =1
         ja  =1
         ib  =ia+iink
         jb  =ja+jink
         ic  =ib+iink
         jc  =jb+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            do i=istart,istop
               c(ja+i)=a(ia+i)+(a(ib+i)+a(ic+i))
               d(ja+i)=b(ia+i)+(b(ib+i)+b(ic+i))
               c(jb+i)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
               c(jc+i)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
               d(jb+i)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
               d(jc+i)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
            end do
         end do
    
         if ( la == m ) RETURN
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb  =k+k-2
            kc  =kb+kb
            c1  =trigs(kb+1)
            s1  =trigs(kb+2)
            c2  =trigs(kc+1)
            s2  =trigs(kc+2)
            do l=1,la
               ll    =k+l-1
               istart=(ll-1)*nv2+lm1
               istop =istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                  d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                  c(jb+j)= c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    -(sin60*(b(ib+i)-b(ic+i))) ) &
                          -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    +(sin60*(a(ib+i)-a(ic+i))) )
                  d(jb+j)= s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    -(sin60*(b(ib+i)-b(ic+i))) ) &
                          +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    +(sin60*(a(ib+i)-a(ic+i))) )
                  c(jc+j)= c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    +(sin60*(b(ib+i)-b(ic+i))) ) &
                          -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    -(sin60*(a(ib+i)-a(ic+i))) )
                  d(jc+j)= s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    +(sin60*(b(ib+i)-b(ic+i))) ) &
                          +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    -(sin60*(a(ib+i)-a(ic+i))) )
               end do
            end do
         end do
    
      else if ( ifac == 4 ) then
    
         ia  =1
         ja  =1
         ib  =ia+iink
         jb  =ja+jink
         ic  =ib+iink
         jc  =jb+jink
         id  =ic+iink
         jd  =jc+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            do i=istart,istop
               c(ja+i)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
               c(jc+i)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
               d(ja+i)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
               d(jc+i)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
               c(jb+i)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
               c(jd+i)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
               d(jb+i)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
               d(jd+i)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
            end do
         end do
    
         if ( la == m ) RETURN
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb=k+k-2
            kc=kb+kb
            kd=kc+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            do l=1,la
               ll    =k+l-1
               istart=(ll-1)*nv2+lm1
               istop =istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                  d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                  c(jc+j)= c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                          -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                  d(jc+j)= s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                          +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                  c(jb+j)= c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                          -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                  d(jb+j)= s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                          +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                  c(jd+j)= c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                          -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                  d(jd+j)= s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                          +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
               end do
            end do
         end do
    
      else if ( ifac == 5 ) then
    
         ia=1
         ja=1
         ib=ia+iink
         jb=ja+jink
         ic=ib+iink
         jc=jb+jink
         id=ic+iink
         jd=jc+jink
         ie=id+iink
         je=jd+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            do i=istart,istop
               c(ja+i)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
               d(ja+i)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
               c(jb+i)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))   &
                               -cos36*(a(ic+i)+a(id+i)) ) &
                             -( sin72*(b(ib+i)-b(ie+i)) + &
                                sin36*(b(ic+i)-b(id+i)) )
               c(je+i)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))   &
                               -cos36*(a(ic+i)+a(id+i)) ) &
                             +( sin72*(b(ib+i)-b(ie+i)) + &
                                sin36*(b(ic+i)-b(id+i)) )
               d(jb+i)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))   &
                               -cos36*(b(ic+i)+b(id+i)) ) &
                             +( sin72*(a(ib+i)-a(ie+i)) + &
                                sin36*(a(ic+i)-a(id+i)) )
               d(je+i)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))   & 
                               -cos36*(b(ic+i)+b(id+i)) ) &
                             -( sin72*(a(ib+i)-a(ie+i)) + &
                                sin36*(a(ic+i)-a(id+i)) )
               c(jc+i)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))   &
                               +cos72*(a(ic+i)+a(id+i)) ) &
                             -( sin36*(b(ib+i)-b(ie+i)) - &
                                sin72*(b(ic+i)-b(id+i)) )
               c(jd+i)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))   &
                               +cos72*(a(ic+i)+a(id+i)) ) &
                             +( sin36*(b(ib+i)-b(ie+i)) - &
                                sin72*(b(ic+i)-b(id+i)) )
               d(jc+i)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))   &
                               +cos72*(b(ic+i)+b(id+i)) ) &
                             +( sin36*(a(ib+i)-a(ie+i)) - &
                    sin72*(a(ic+i)-a(id+i)) )
               d(jd+i)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))   &
                               +cos72*(b(ic+i)+b(id+i)) ) &
                             -( sin36*(a(ib+i)-a(ie+i)) - &
                                sin72*(a(ic+i)-a(id+i)) )
            end do
         end do
    
         if ( la == m ) return
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb=k+k-2
            kc=kb+kb
            kd=kc+kb
            ke=kd+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            c4=trigs(ke+1)
            s4=trigs(ke+2)
            do l=1,la
               ll=k+l-1
               istart=(ll-1)*nv2+lm1
               istop=istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=a(ia+i) + (a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                  d(ja+j)=b(ia+i) + (b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                  c(jb+j)=c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    -( sin72*(b(ib+i)-b(ie+i))     &
                                     + sin36*(b(ic+i)-b(id+i)) ) ) &
                         -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    +( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  d(jb+j)=s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    -( sin72*(b(ib+i)-b(ie+i))     &
                                      +sin36*(b(ic+i)-b(id+i)) ) ) &
                         +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    +( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  c(je+j)=c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    +( sin72*(b(ib+i)-b(ie+i))     &
                                      +sin36*(b(ic+i)-b(id+i)) ) ) &
                         -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    -( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  d(je+j)=s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    +( sin72*(b(ib+i)-b(ie+i))     &
                                      +sin36*(b(ic+i)-b(id+i)) ) ) &
                         +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    -( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  c(jc+j)=c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    -( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    +( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
                  d(jc+j)=s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    -( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    +( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
                  c(jd+j)=c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    +( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    -( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
                  d(jd+j)=s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    +( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    -( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
               end do
            end do
         end do
    
      end if

   end subroutine fft_fac_real
!-----------------------------------------------------------------------------------
   subroutine fft_fac_complex(a,b,c,d,trigs,nv,l1,l2,n,ifac,la)
      !---------------------------------------------------------------------------
      !     main part of Fourier / Chebychev transform
      !     called in costf1, costf2
      !---------------------------------------------------------------------------
    
      !-- Input variables:
      integer,     intent(in) :: nv,l1,l2,n,ifac,la
      complex(cp), intent(in) :: a(*),b(*)
      real(cp),    intent(in) :: trigs(2*n)
    
      !-- Output variables
      complex(cp), intent(out) :: c(*),d(*)
    
      !-- Local variables:
      integer :: i,ia,ib,ic,id,ie
      integer :: j,ja,jb,jc,jd,je
      integer :: istart,istop,iink
      integer :: jump,jadd,jink
      integer :: k,kb,kc,kd,ke
      integer :: nv2
      integer :: l,m,lm1,lm2,ll,la1
    
      real(cp) :: c1,c2,c3,c4
      real(cp) :: s1,s2,s3,s4
    
      m    =n/ifac
      iink =m*2*nv
      jink =la*2*nv
      jump =(ifac-1)*jink
      nv2  =2*nv
      lm1  =l1-1
      lm2  =l2-l1
    
      if ( ifac == 2 ) then
    
         ia  =1
         ja  =1
         ib  =ia+iink
         jb  =ja+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop=istart+lm2
            do i=istart,istop
               c(ja+i)=a(ia+i)+a(ib+i)
               d(ja+i)=b(ia+i)+b(ib+i)
               c(jb+i)=a(ia+i)-a(ib+i)
               d(jb+i)=b(ia+i)-b(ib+i)
            end do
         end do
    
         if ( la == m ) return
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb  =k+k-2
            c1  =trigs(kb+1)
            s1  =trigs(kb+2)
            do l=1,la
               ll    =k+l-1
               istart=(ll-1)*nv2+lm1
               istop =istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=a(ia+i)+a(ib+i)
                  d(ja+j)=b(ia+i)+b(ib+i)
                  c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
                  d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
               end do
            end do
         end do
    
      else if ( ifac == 3 ) then
    
         ia  =1
         ja  =1
         ib  =ia+iink
         jb  =ja+jink
         ic  =ib+iink
         jc  =jb+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            do i=istart,istop
               c(ja+i)=a(ia+i)+(a(ib+i)+a(ic+i))
               d(ja+i)=b(ia+i)+(b(ib+i)+b(ic+i))
               c(jb+i)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
               c(jc+i)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
               d(jb+i)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
               d(jc+i)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
            end do
         end do
    
         if ( la == m ) RETURN
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb  =k+k-2
            kc  =kb+kb
            c1  =trigs(kb+1)
            s1  =trigs(kb+2)
            c2  =trigs(kc+1)
            s2  =trigs(kc+2)
            do l=1,la
               ll    =k+l-1
               istart=(ll-1)*nv2+lm1
               istop =istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                  d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                  c(jb+j)= c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    -(sin60*(b(ib+i)-b(ic+i))) ) &
                          -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    +(sin60*(a(ib+i)-a(ic+i))) )
                  d(jb+j)= s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    -(sin60*(b(ib+i)-b(ic+i))) ) &
                          +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    +(sin60*(a(ib+i)-a(ic+i))) )
                  c(jc+j)= c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    +(sin60*(b(ib+i)-b(ic+i))) ) &
                          -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    -(sin60*(a(ib+i)-a(ic+i))) )
                  d(jc+j)= s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))   &
                                    +(sin60*(b(ib+i)-b(ic+i))) ) &
                          +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))   &
                                    -(sin60*(a(ib+i)-a(ic+i))) )
               end do
            end do
         end do
    
      else if ( ifac == 4 ) then
    
         ia  =1
         ja  =1
         ib  =ia+iink
         jb  =ja+jink
         ic  =ib+iink
         jc  =jb+jink
         id  =ic+iink
         jd  =jc+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            do i=istart,istop
               c(ja+i)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
               c(jc+i)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
               d(ja+i)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
               d(jc+i)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
               c(jb+i)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
               c(jd+i)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
               d(jb+i)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
               d(jd+i)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
            end do
         end do
    
         if ( la == m ) RETURN
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb=k+k-2
            kc=kb+kb
            kd=kc+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            do l=1,la
               ll    =k+l-1
               istart=(ll-1)*nv2+lm1
               istop =istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                  d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                  c(jc+j)= c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                          -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                  d(jc+j)= s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                          +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                  c(jb+j)= c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                          -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                  d(jb+j)= s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                          +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                  c(jd+j)= c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                          -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                  d(jd+j)= s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                          +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
               end do
            end do
         end do
    
      else if ( ifac == 5 ) then
    
         ia=1
         ja=1
         ib=ia+iink
         jb=ja+jink
         ic=ib+iink
         jc=jb+jink
         id=ic+iink
         jd=jc+jink
         ie=id+iink
         je=jd+jink
         jadd=0
         do l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            do i=istart,istop
               c(ja+i)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
               d(ja+i)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
               c(jb+i)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))   &
                               -cos36*(a(ic+i)+a(id+i)) ) &
                             -( sin72*(b(ib+i)-b(ie+i)) + &
                                sin36*(b(ic+i)-b(id+i)) )
               c(je+i)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))   &
                               -cos36*(a(ic+i)+a(id+i)) ) &
                             +( sin72*(b(ib+i)-b(ie+i)) + &
                                sin36*(b(ic+i)-b(id+i)) )
               d(jb+i)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))   &
                               -cos36*(b(ic+i)+b(id+i)) ) &
                             +( sin72*(a(ib+i)-a(ie+i)) + &
                                sin36*(a(ic+i)-a(id+i)) )
               d(je+i)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))   & 
                               -cos36*(b(ic+i)+b(id+i)) ) &
                             -( sin72*(a(ib+i)-a(ie+i)) + &
                                sin36*(a(ic+i)-a(id+i)) )
               c(jc+i)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))   &
                               +cos72*(a(ic+i)+a(id+i)) ) &
                             -( sin36*(b(ib+i)-b(ie+i)) - &
                                sin72*(b(ic+i)-b(id+i)) )
               c(jd+i)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))   &
                               +cos72*(a(ic+i)+a(id+i)) ) &
                             +( sin36*(b(ib+i)-b(ie+i)) - &
                                sin72*(b(ic+i)-b(id+i)) )
               d(jc+i)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))   &
                               +cos72*(b(ic+i)+b(id+i)) ) &
                             +( sin36*(a(ib+i)-a(ie+i)) - &
                    sin72*(a(ic+i)-a(id+i)) )
               d(jd+i)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))   &
                               +cos72*(b(ic+i)+b(id+i)) ) &
                             -( sin36*(a(ib+i)-a(ie+i)) - &
                                sin72*(a(ic+i)-a(id+i)) )
            end do
         end do
    
         if ( la == m ) return
    
         la1=la+1
         do k=la1,m,la
            jadd=jadd+jump
            kb=k+k-2
            kc=kb+kb
            kd=kc+kb
            ke=kd+kb
            c1=trigs(kb+1)
            s1=trigs(kb+2)
            c2=trigs(kc+1)
            s2=trigs(kc+2)
            c3=trigs(kd+1)
            s3=trigs(kd+2)
            c4=trigs(ke+1)
            s4=trigs(ke+2)
            do l=1,la
               ll=k+l-1
               istart=(ll-1)*nv2+lm1
               istop=istart+lm2
               do i=istart,istop
                  j=i+jadd
                  c(ja+j)=a(ia+i) + (a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                  d(ja+j)=b(ia+i) + (b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                  c(jb+j)=c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    -( sin72*(b(ib+i)-b(ie+i))     &
                                     + sin36*(b(ic+i)-b(id+i)) ) ) &
                         -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    +( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  d(jb+j)=s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    -( sin72*(b(ib+i)-b(ie+i))     &
                                      +sin36*(b(ic+i)-b(id+i)) ) ) &
                         +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    +( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  c(je+j)=c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    +( sin72*(b(ib+i)-b(ie+i))     &
                                      +sin36*(b(ic+i)-b(id+i)) ) ) &
                         -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    -( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  d(je+j)=s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))     &
                                      -cos36*(a(ic+i)+a(id+i)) )   &
                                    +( sin72*(b(ib+i)-b(ie+i))     &
                                      +sin36*(b(ic+i)-b(id+i)) ) ) &
                         +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))     &
                                      -cos36*(b(ic+i)+b(id+i)) )   &
                                    -( sin72*(a(ib+i)-a(ie+i))     &
                                      +sin36*(a(ic+i)-a(id+i)) ) )
                  c(jc+j)=c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    -( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    +( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
                  d(jc+j)=s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    -( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    +( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
                  c(jd+j)=c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    +( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    -( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
                  d(jd+j)=s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))     &
                                      +cos72*(a(ic+i)+a(id+i)) )   &
                                    +( sin36*(b(ib+i)-b(ie+i))     &
                                      -sin72*(b(ic+i)-b(id+i)) ) ) &
                         +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))     &
                                      +cos72*(b(ic+i)+b(id+i)) )   &
                                    -( sin36*(a(ib+i)-a(ie+i))     &
                                      -sin72*(a(ic+i)-a(id+i)) ) )
               end do
            end do
         end do
    
      end if

   end subroutine fft_fac_complex
!------------------------------------------------------------------------------
end module fft_fac_mod
