    subroutine extrapolate(input,radratio,br,bt,bp,np,nt,azsym)

    implicit none

! input variables
    integer, intent(in) :: azsym
    integer :: np,nt
    real(kind=8), intent(in) :: radratio
    complex(kind=8), dimension(np,nt),intent(in):: input

! output variables
    complex(kind=8), dimension(np,nt),intent(out) :: br,bt,bp

! local variables
    integer :: nalias=20
    integer :: n_phi_max,n_theta_max,minc
    integer :: lmP_max,m,mca
    integer :: l_max,m_max,lm_max,n_m_max
    integer :: n_theta,l,lm,lmP
    real(kind=8) :: dpi
    real(kind=8) :: colat
    real(kind=8), dimension(nt) :: XX,W
    real(kind=8), dimension(:), allocatable :: D_l
    integer, dimension(:), allocatable :: lmP2m,lmP2l,lmP2lm
    integer, dimension(:), allocatable :: lm2m,lm2l,lStopP
    logical, dimension(:), allocatable :: lmOdd
    integer, dimension(:), allocatable :: lStart,lStop,lStartP
    complex(kind=8), dimension(:), allocatable :: inputLM,cs2LM
    real(kind=8), dimension(:,:), allocatable :: Plm,dPlm,wPlm
    complex(kind=8), dimension(:,:), allocatable :: dPhi
    real(kind=8), dimension(:), allocatable :: plma,dtheta_plma

    dpi=4.d0*datan(1.d0)
    n_phi_max=np
    n_theta_max=nt
    minc=azsym

    l_max=(nalias*n_theta_max)/30
    m_max=(l_max/minc)*minc
    lm_max=m_max*(l_max+1)/minc-m_max*(m_max-minc)/(2*minc)+ &
                 (l_max+1-m_max)
    n_m_max=m_max/minc+1
    lmP_max=lm_max+n_m_max

    allocate(lm2l(1:lm_max))
    allocate(lm2m(1:lm_max))
    allocate(D_l(1:lm_max))
    allocate(inputLM(1:lmP_max))
    allocate(cs2LM(1:lmP_max))
    allocate(plma(1:lmP_max))
    allocate(lmP2m(1:lmP_max))
    allocate(lmP2l(1:lmP_max))
    allocate(lmP2lm(1:lmP_max))
    allocate(dtheta_plma(1:lmP_max))
    allocate(dPhi(1:lm_max,1:n_theta_max/2))
    allocate(Plm(1:lm_max,1:n_theta_max/2))
    allocate(dPlm(1:lm_max,1:n_theta_max/2))
    allocate(wPlm(1:lmP_max,1:n_theta_max/2))
    allocate(lStart(1:n_m_max))
    allocate(lStop(1:n_m_max))
    allocate(lStartP(1:n_m_max))
    allocate(lStopP(1:n_m_max))
    allocate(lmOdd(1:n_m_max))


! gauleg.f
    call gauleg(XX,W,n_theta_max)

! lmP from s_blocking.f
    call get_blocks(l_max,lmP_max,n_m_max,minc,lm2l,lm2m, &
                    lmP2l,lmP2m,lmP2lm,lStartP,lStopP,lStart,lStop,lmOdd)

! Legendre polynomials from s_horizontal.f
    do n_theta=1,n_theta_max/2  ! Loop over colat in NHS

        colat=XX(n_theta)

        call plm_theta(colat,l_max+1,m_max,minc, &
                       plma,dtheta_plma,lmP_max)
        do lmP=1,lmP_max
            l=lmP2l(lmP)
            m=lmP2m(lmP)
            if ( l <= l_max ) then
                lm=lmP2lm(lmP)
                Plm(lm,n_theta) =plma(lmP)
            ! True theta derivative !!!
                dPlm(lm,n_theta)=dtheta_plma(lmP)/dsin(colat)
            ! Add the theta dependence in dPhi to simplify the output
                dPhi(lm,n_theta)=dcmplx(0.d0,dble(m))/dsin(colat)
            end if
            wPlm(lmP,n_theta)=2.d0*dpi*W(n_theta)*plma(lmP)
        end do
    end do

    do lm=1,lm_max
        l=lm2l(lm)
        m=lm2m(lm)
        D_l(lm)=DBLE(l+1)
    enddo

!       spatial to spectral
    call spatSpec(input,inputLM,wPlm,lStartP,lStopP,lStart, &
                  lStop,np,nt,lmP_max,n_m_max)

!       Brsurf
    lm=0
    do mca=0,n_m_max-1
        do l=azsym*mca,l_max
            lm=lm+1
            cs2LM(lm)  =-inputLM(lm)*radratio**(l+2)/D_l(lm)
            inputLM(lm)= inputLM(lm)*radratio**(l+2)
        enddo
    enddo

!       spectral to spatial
    call specSpat(inputLM,cs2LM,br,bt,bp,Plm,dPlm,dPhi, &
                                 lStart,lStop,lmOdd,np,nt, &
                                   lm_max,lmP_max,n_m_max)

    deallocate(lm2l)
    deallocate(lm2m)
    deallocate(D_l)
    deallocate(dPhi)
    deallocate(inputLM)
    deallocate(cs2LM)
    deallocate(plma)
    deallocate(lmP2m)
    deallocate(lmP2l)
    deallocate(lmP2lm)
    deallocate(dtheta_plma)
    deallocate(Plm)
    deallocate(dPlm)
    deallocate(wPlm)
    deallocate(lStart)
    deallocate(lStop)
    deallocate(lStartP)
    deallocate(lStopP)
    deallocate(lmOdd)

    end subroutine
!---------------------------------------------------------------------
    subroutine gauleg(XX,W,n_theta_max)

    implicit none

    integer :: n_theta_max
    integer ::N,M,I,J
    real(kind=8) :: dpi,XXM,XXL,EPS,ZZ,ZZ1
    real(kind=8) :: P1,P2,P3,PP
    real(kind=8), dimension(n_theta_max) :: XX,W

    dpi=4.d0*datan(1.d0)
    N=n_theta_max
    M=(N+1)/2
    XXM=0.d0
    XXL=1.d0
    EPS=3.D-14

    do I=1,M
        ZZ=dcos( dpi*( (dble(I)-0.25D0)/(dble(N)+0.5D0)) )

        ZZ1=0
        do while (dabs(ZZ-ZZ1).gt.EPS)
            P1=1.d0
            P2=0.d0
            do J=1,N
                P3=P2
                P2=P1
                P1=( dble(2*J-1)*ZZ*P2-dble(J-1)*P3 )/dble(J)
            enddo
            PP=dble(N)*(ZZ*P1-P2)/(ZZ*ZZ-1.D0)
            ZZ1=ZZ
            ZZ=ZZ1-P1/PP
        enddo

        XX(I)=dacos(XXM+XXL*ZZ)
        XX(N+1-I)=dacos(XXM-XXL*ZZ)
        W(I)=2.D0*XXL/((1.D0-ZZ*ZZ)*PP*PP)
        W(N+1-I)=W(I)
    enddo

    end subroutine
!---------------------------------------------------------------------
    subroutine get_blocks(l_max,lmP_max,n_m_max,minc, &
                        lm2l,lm2m,lmP2l,lmP2m,lmP2lm, &
                   lStartP,lStopP,lStart,lStop,lmOdd)

    implicit none

    integer:: l_max,lmP_max,minc,n_m_max
    integer, dimension(lmP_max) :: lmP2m,lmP2l,lmP2lm
    integer, dimension(l_max) :: lm2m,lm2l
    integer, dimension(n_m_max) :: lStart,lStop,lStartP,lStopP
    logical, dimension(n_m_max) :: lmOdd
    integer :: lm,mc,lmP,l,m

! lmP from s_blocking.f
    lm=0
    mc=0
    lmP=0
    do m=0,l_max,minc
        mc=mc+1
        do l=m,l_max
            lm         =lm+1
            lm2l(lm)   =l
            lm2m(lm)   =m
            lmP        =lmP+1
            lmP2l(lmP) =l
            lmP2m(lmP) =m
            lmP2lm(lmP)=lm
        enddo
        l=l_max+1    ! Extra l for lmP
        lmP=lmP+1
        lmP2l(lmP) =l
        lmP2m(lmP) =m
        lmP2lm(lmP)=-1
    enddo


    lStartP(1)=1
    lStopP(1) =l_max+2
    lStart(1) =1
    lStop(1)  =l_max+1
    IF ( MOD(l_max,2) == 0 ) THEN
        lmOdd(1) =.TRUE.
    ELSE
        lmOdd(1) =.FALSE.
    ENDIF
    DO mc=2,n_m_max
        m=(mc-1)*minc
        lStartP(mc)=lStopP(mc-1)+1
        lStopP(mc) =lStartP(mc) +l_max-m+1
        lStart(mc) =lStop(mc-1) +1
        lStop(mc)  =lStart(mc)  +l_max-m
        IF ( MOD(lStop(mc)-lStart(mc),2) == 0 ) THEN
            lmOdd(mc) =.TRUE.
        ELSE
            lmOdd(mc) =.FALSE.
        END IF
    END DO

    end subroutine
!--------------------------------------------------------
    subroutine plm_theta(theta,max_degree,max_order,m0, &
    plma,dtheta_plma,ndim_plma)
     
    implicit none
     
    real(kind=8) :: theta
    integer :: max_degree
    integer :: max_order
    integer :: m0
    integer :: ndim_plma

    real(kind=8) :: plma(ndim_plma)
    real(kind=8) :: dtheta_plma(ndim_plma)

    real(kind=8) :: sq2,dnorm,fac,plm,plm1,plm2
    integer :: l,m,j,pos
     
    dnorm=1.d0
    sq2=dsqrt(2.d0)

    pos=0
    do m=0,max_order,m0
        if (m == 0) then
            dnorm=1.d0
        else
            dnorm=sq2
        endif
         
        fac=1.d0
        do j=3,2*m+1,2
            fac=fac*dble(j)/dble(j-1)
        end do

        plm=dsqrt(fac)
        if( dsin(theta) /= 0.d0 ) then
            plm=plm*(-dsin(theta))**m
        elseif( m /= 0 ) then
            plm=0.d0
        endif
         
        l=m

        pos=pos+1
        plma(pos) = dnorm*plm
         
        plm1=0.d0
         
        do l=m+1,max_degree
            plm2=plm1
            plm1=plm
            plm= dcos(theta)* dsqrt( dble( (2*l-1)*(2*l+1) ) / &
            dble( (l-m)*(l+m) )  ) * plm1 - &
            dsqrt( dble( (2*l+1)*(l+m-1)*(l-m-1) ) / &
            dble( (2*l-3)*(l-m)*(l+m) ) ) * plm2
             
            pos=pos+1
            IF ( pos > ndim_plma ) THEN
                WRITE(*,*) '! Dimension ndim_plma too small'
                WRITE(*,*) '! in subroutine plm_theta!'
                STOP
            END IF
            plma(pos) = dnorm*plm
             
        end do

        l=max_degree+1
        plm2=plm1
        plm1=plm
        plm= dcos(theta)* dsqrt( dble( (2*l-1)*(2*l+1) ) / &
        dble( (l-m)*(l+m) )  ) * plm1 - &
        dsqrt( dble( (2*l+1)*(l+m-1)*(l-m-1) ) / &
        dble( (2*l-3)*(l-m)*(l+m) ) ) * plm2
        dtheta_plma(pos)=dnorm*plm
         
    end do    ! loop over order !
     
    pos=0
    do m=0,max_order,m0
        l=m
        pos=pos+1
        IF ( pos > ndim_plma ) THEN
            WRITE(*,*) '! Dimension ndim_plma too small'
            WRITE(*,*) '! in subroutine plm_theta!'
            STOP
        END IF
        if ( m < max_degree ) then
            dtheta_plma(pos)= l/dsqrt(dble(2*l+3)) * &
            plma(pos+1)
        else
            dtheta_plma(pos)= l/dsqrt(dble(2*l+3)) * &
            dtheta_plma(pos)
        end if
             
        do l=m+1,max_degree-1

            pos=pos+1
            IF ( pos > ndim_plma ) THEN
                WRITE(*,*) '! Dimension ndim_plma too small'
                WRITE(*,*) '! in subroutine plm_theta!'
                STOP
            END IF
            dtheta_plma(pos)= &
            l*dsqrt( dble((l+m+1)*(l-m+1)) / &
            dble((2*l+1)*(2*l+3)) &
            ) * plma(pos+1)  - &
            (l+1)*dsqrt( dble((l+m)*(l-m)) / &
            dble((2*l-1)*(2*l+1)) &
            ) * plma(pos-1)
        end do ! loop over degree

        IF ( m < max_degree ) THEN
            l=max_degree
            pos=pos+1
            IF ( pos > ndim_plma ) THEN
                WRITE(*,*) '! Dimension ndim_plma too small'
                WRITE(*,*) '! in subroutine plm_theta!'
                STOP
            END IF
            dtheta_plma(pos)= &
            l*dsqrt( dble((l+m+1)*(l-m+1)) / &
            dble((2*l+1)*(2*l+3)) &
            ) * dtheta_plma(pos)  - &
            (l+1)*dsqrt( dble((l+m)*(l-m)) / &
            dble((2*l-1)*(2*l+1)) &
            ) * plma(pos-1)
        END IF

    end do ! loop over order
     
    return
    end subroutine plm_theta
!---------------------------------------------------------------
    subroutine spatSpec(input,inputLM,wPlm,lStartP,lStopP,lStart, &
    lStop,np,nt,lmP_max,n_m_max)

    implicit none

! input variables
    integer :: np,nt,lmP_max,n_m_max
    real(kind=8), dimension(1:lmP_max,1:nt/2) :: wPlm
    complex(kind=8), dimension(np,nt) :: input
    integer, dimension(1:n_m_max) :: lStart,lStop,lStartP,lStopP

! output variables
    complex(kind=8), dimension(1:lmP_max) :: flm1,inputLM

! local variables
    complex(kind=8), dimension(np,nt) :: a1plus,a1minus,work2
    complex(kind=8)  :: ap_1,ap_2,am_1,am_2
    integer :: n_theta,n_theta_hs,n_theta_n,n_theta_s,mca
    integer :: n_theta_1,n_theta_2,n_theta_rel_1,n_theta_rel_2
    integer :: n_theta_max,lms,lm1,lmp

    n_theta_max=nt

    do n_theta=1,n_theta_max/2
        work2(:,2*n_theta-1)=input(:,n_theta)
        work2(:,2*n_theta)=input(:,n_theta_max-n_theta+1)
    end do

    n_theta_hs=0
    do n_theta_n=1,n_theta_max,2
        n_theta_s=n_theta_n+1
        n_theta_hs=n_theta_hs+1
        do mca=1,n_m_max
            a1plus(mca,n_theta_hs)=work2(mca,n_theta_n)+ &
            work2(mca,n_theta_s)
            a1minus(mca,n_theta_hs)=work2(mca,n_theta_n)- &
            work2(mca,n_theta_s)
        enddo
    enddo

    do mca=1,n_m_max
        lms=lStopP(mca)
        lm1=lms-1
        ap_1=a1plus(mca,1)
        ap_2=a1plus(mca,2)
        am_1=a1minus(mca,1)
        am_2=a1minus(mca,2)
        do lmp=lStartP(mca),lms-1,2
            lm1=lmp+1
            flm1(lmp)=ap_1*wPlm(lmp,1)+ap_2*wPlm(lmp,2)
            flm1(lm1)=am_1*wPlm(lm1,1)+am_2*wPlm(lm1,2)
        enddo

        if ( lm1 < lms) then
            flm1(lms)=ap_1*wPlm(lms,1)+ap_2*wPlm(lms,2)
        end if
    enddo

    n_theta_1=1
    do n_theta_rel_1=3,n_theta_max/2,2
        n_theta_rel_2=n_theta_rel_1+1
        n_theta_1=n_theta_1+2
        n_theta_2=n_theta_1+1
        do mca=1,n_m_max
            lms=lStopP(mca)
            lm1=lms-1
            ap_1=a1plus(mca,n_theta_rel_1)
            ap_2=a1plus(mca,n_theta_rel_2)
            am_1=a1minus(mca,n_theta_rel_1)
            am_2=a1minus(mca,n_theta_rel_2)
            do lmp=lStartP(mca),lms-1,2
                lm1=lmp+1
                flm1(lmp)=flm1(lmp)+ap_1*wPlm(lmp,n_theta_1)+ &
                ap_2*wPlm(lmp,n_theta_2)
                flm1(lm1)=flm1(lm1)+am_1*wPlm(lm1,n_theta_1)+ &
                am_2*wPlm(lm1,n_theta_2)
            enddo
                         
            if ( lm1 < lms) then
                flm1(lms)=flm1(lms)+ap_1*wPlm(lms,n_theta_1)+ &
                ap_2*wPlm(lms,n_theta_2)
            end if
        enddo
    enddo

    do mca=1,n_m_max
        inputLM(lStart(mca):lStop(mca))= &
        flm1(lStartP(mca):lStopP(mca)-1)
    enddo

    end subroutine
!-------------------------------------------------------------------
    subroutine specSpat(inputLM,cs2LM,br,bt,bp,Plm,dPlm,dPhi, &
                                    lStart,lStop,lmOdd,np,nt, &
                                      lm_max,lmP_max,n_m_max)

    implicit none
            
! input variables
    integer :: n_m_max,np,nt,lmP_max,lm_max
    integer, dimension(1:n_m_max) :: lStart,lStop
    logical, dimension(1:n_m_max) :: lmOdd
    real(kind=8), dimension(1:lm_max,1:nt/2) :: Plm,dPlm
    complex(kind=8), dimension(1:lm_max,1:nt/2) :: dPhi
    complex(kind=8), dimension(1:lmP_max) :: inputLM,cs2LM

! output variables
    complex(kind=8), dimension(np,nt) :: br,bt,bp

! local variables
    integer :: icn,icc,ic1,n_m,lm,lms,n_theta_max,n_phi_max
    complex(kind=8) :: s12,z12,t12,x12,u12,v12
    complex(kind=8), dimension(1:n_m_max,1:nt) :: work,work1,work2

    n_theta_max=nt
    n_phi_max=np

    icn=0
    do icc=1,n_theta_max,2
        ic1=icc+1
        icn=icn+1
        do n_m=1,n_m_max
            lms=lStop(n_m)
            s12=0.d0
            z12=0.d0
            t12=0.d0
            x12=0.d0
            u12=0.d0
            v12=0.d0
            do lm=lStart(n_m),lms-1,2
                s12=s12+inputLM(lm)*Plm(lm,icn)
                z12=z12+inputLM(lm+1)*Plm(lm+1,icn)
                t12=t12+cs2LM(lm)*dPlm(lm,icn)
                x12=x12+cs2LM(lm+1)*dPlm(lm+1,icn)
                u12=u12+cs2LM(lm)*dPhi(lm,icn)*Plm(lm,icn)
                v12=v12+cs2LM(lm+1)*dPhi(lm+1,icn)*Plm(lm+1,icn)
            enddo
            if ( lmOdd(n_m) ) then
                s12=s12+inputLM(lms)*Plm(lms,icn)
                t12=t12+cs2LM(lms)*dPlm(lms,icn)
                u12=u12+cs2LM(lms)*dPhi(lms,icn)*Plm(lms,icn)
            endif
            work(n_m,icc)=s12+z12
            work(n_m,ic1)=s12-z12
            work1(n_m,icc)=t12+x12
            work1(n_m,ic1)=-t12+x12
            work2(n_m,icc)=(u12+v12)
            work2(n_m,ic1)=(u12-v12)
        enddo
    enddo

    br(1:n_m_max,1:n_theta_max/2)=work(:,1:n_theta_max:2)
    br(1:n_m_max,n_theta_max/2+1:n_theta_max)=work(:,n_theta_max:2:-2)

    br(2:n_phi_max,:)=br(2:n_phi_max,:)/2.d0
    br(n_m_max+1:n_phi_max/2+1,:)=0.d0
    br(n_phi_max/2+2:n_phi_max,:)=conjg(br(n_phi_max/2:2:-1,:))

    bt(1:n_m_max,1:n_theta_max/2)=work1(:,1:n_theta_max:2)
    bt(1:n_m_max,n_theta_max/2+1:n_theta_max)=work1(:,n_theta_max:2:-2)

    bt(2:n_phi_max,:)=bt(2:n_phi_max,:)/2.d0
    bt(n_m_max+1:n_phi_max/2+1,:)=0.d0
    bt(n_phi_max/2+2:n_phi_max,:)=conjg(bt(n_phi_max/2:2:-1,:))

    bp(1:n_m_max,1:n_theta_max/2)=work2(:,1:n_theta_max:2)
    bp(1:n_m_max,n_theta_max/2+1:n_theta_max)=work2(:,n_theta_max:2:-2)

    bp(2:n_phi_max,:)=bp(2:n_phi_max,:)/2.d0
    bp(n_m_max+1:n_phi_max/2+1,:)=0.d0
    bp(n_phi_max/2+2:n_phi_max,:)=conjg(bp(n_phi_max/2:2:-1,:))

    end subroutine
