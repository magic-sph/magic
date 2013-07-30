!$Id$
!---------------------------------------------------------------------------
    SUBROUTINE fft_fac(a,b,c,d,trigs,nv,l1,l2,n,ifac,la,sinCos)
!---------------------------------------------------------------------------

!    !------------ This is release 2 level 1  --------------!
!    !------------ Created on 1/17/02  by JW. --------------!

!---------------------------------------------------------------------------

!     main part of Fourier / Chebychev transform

!     called in costf1, costf2

!---------------------------------------------------------------------------

    IMPLICIT NONE

!-- input:
    INTEGER :: nv,l1,l2,n,ifac,la
    REAL(kind=8) :: a(*),b(*),c(*),d(*)
    REAL(kind=8) :: trigs(2*n)
    REAL(kind=8) :: sinCos(5)

!-- output: a(*),b(*),c(*),d(*)

!-- local:
    INTEGER :: i,ia,ib,ic,id,ie
    INTEGER :: j,ja,jb,jc,jd,je
    INTEGER :: istart,istop,iink
    INTEGER :: jump,jadd,jink
    INTEGER :: k,kb,kc,kd,ke
    INTEGER :: nv2
    INTEGER :: l,m,lm1,lm2,ll,la1

    REAL(kind=8) :: c1,c2,c3,c4
    REAL(kind=8) :: s1,s2,s3,s4
    REAL(kind=8) :: sin36,cos36,sin72,cos72,sin60

!-- end of declaration
!---------------------------------------------------------------------------

    sin36=sinCos(1)
    cos36=sinCos(2)
    sin72=sinCos(3)
    cos72=sinCos(4)
    sin60=sinCos(5)

    m    =n/ifac
    iink =m*2*nv
    jink =la*2*nv
    jump =(ifac-1)*jink
    nv2  =2*nv
    lm1  =l1-1
    lm2  =l2-l1

    IF ( ifac == 2 ) THEN

        ia  =1
        ja  =1
        ib  =ia+iink
        jb  =ja+jink
        jadd=0
        DO l=1,la
            istart=(l-1)*nv2+lm1
            istop=istart+lm2
            DO i=istart,istop
                c(ja+i)=a(ia+i)+a(ib+i)
                d(ja+i)=b(ia+i)+b(ib+i)
                c(jb+i)=a(ia+i)-a(ib+i)
                d(jb+i)=b(ia+i)-b(ib+i)
            END DO
        END DO

        IF ( la == m ) RETURN

        la1=la+1
        DO k=la1,m,la
            jadd=jadd+jump
            kb  =k+k-2
            c1  =trigs(kb+1)
            s1  =trigs(kb+2)
            DO l=1,la
                ll    =k+l-1
                istart=(ll-1)*nv2+lm1
                istop =istart+lm2
                DO i=istart,istop
                    j=i+jadd
                    c(ja+j)=a(ia+i)+a(ib+i)
                    d(ja+j)=b(ia+i)+b(ib+i)
                    c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
                    d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
                END DO
            END DO
        END DO

    ELSE IF ( ifac == 3 ) THEN

        ia  =1
        ja  =1
        ib  =ia+iink
        jb  =ja+jink
        ic  =ib+iink
        jc  =jb+jink
        jadd=0
        DO l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            DO i=istart,istop
                c(ja+i)=a(ia+i)+(a(ib+i)+a(ic+i))
                d(ja+i)=b(ia+i)+(b(ib+i)+b(ic+i))
                c(jb+i)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i))) &
                             -(sin60*(b(ib+i)-b(ic+i)))
                c(jc+i)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i))) &
                             +(sin60*(b(ib+i)-b(ic+i)))
                d(jb+i)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i))) &
                             +(sin60*(a(ib+i)-a(ic+i)))
                d(jc+i)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i))) &
                             -(sin60*(a(ib+i)-a(ic+i)))
            END DO
        END DO

        IF ( la == m ) RETURN

        la1=la+1
        DO k=la1,m,la
            jadd=jadd+jump
            kb  =k+k-2
            kc  =kb+kb
            c1  =trigs(kb+1)
            s1  =trigs(kb+2)
            c2  =trigs(kc+1)
            s2  =trigs(kc+2)
            DO l=1,la
                ll    =k+l-1
                istart=(ll-1)*nv2+lm1
                istop =istart+lm2
                DO i=istart,istop
                    j=i+jadd
                    c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
                    d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
                    c(jb+j)= &
                        c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i))) &
                                 -(sin60*(b(ib+i)-b(ic+i))) ) &
                       -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i))) &
                                 +(sin60*(a(ib+i)-a(ic+i))) )
                    d(jb+j)= &
                        s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i))) &
                                 -(sin60*(b(ib+i)-b(ic+i))) ) &
                       +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i))) &
                                 +(sin60*(a(ib+i)-a(ic+i))) )
                    c(jc+j)= &
                        c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i))) &
                                 +(sin60*(b(ib+i)-b(ic+i))) ) &
                       -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i))) &
                                 -(sin60*(a(ib+i)-a(ic+i))) )
                    d(jc+j)= &
                        s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i))) &
                                 +(sin60*(b(ib+i)-b(ic+i))) ) &
                       +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i))) &
                                 -(sin60*(a(ib+i)-a(ic+i))) )
                END DO
            END DO
        END DO

    ELSE IF ( ifac == 4 ) THEN

        ia  =1
        ja  =1
        ib  =ia+iink
        jb  =ja+jink
        ic  =ib+iink
        jc  =jb+jink
        id  =ic+iink
        jd  =jc+jink
        jadd=0
        DO l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            DO i=istart,istop
                c(ja+i)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                c(jc+i)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
                d(ja+i)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                d(jc+i)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
                c(jb+i)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
                c(jd+i)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
                d(jb+i)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
                d(jd+i)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
            END DO
        END DO

        IF ( la == m ) RETURN

        la1=la+1
        DO k=la1,m,la
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
            DO l=1,la
                ll    =k+l-1
                istart=(ll-1)*nv2+lm1
                istop =istart+lm2
                DO i=istart,istop
                    j=i+jadd
                    c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
                    d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
                    c(jc+j)= &
                        c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                       -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                    d(jc+j)= &
                        s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))) &
                       +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
                    c(jb+j)= &
                        c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                       -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                    d(jb+j)= &
                        s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))) &
                       +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
                    c(jd+j)= &
                        c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                       -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                    d(jd+j)= &
                        s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))) &
                       +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
                END DO
            END DO
        END DO

    ELSE IF ( ifac == 5 ) THEN

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
        DO l=1,la
            istart=(l-1)*nv2+lm1
            istop =istart+lm2
            DO i=istart,istop
                c(ja+i)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                d(ja+i)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                c(jb+i)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i)) &
                                -cos36*(a(ic+i)+a(id+i)) ) &
                              -( sin72*(b(ib+i)-b(ie+i)) + &
                                 sin36*(b(ic+i)-b(id+i)) )
                c(je+i)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i)) &
                                -cos36*(a(ic+i)+a(id+i)) ) &
                              +( sin72*(b(ib+i)-b(ie+i)) + &
                                 sin36*(b(ic+i)-b(id+i)) )
                d(jb+i)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i)) &
                                -cos36*(b(ic+i)+b(id+i)) ) &
                              +( sin72*(a(ib+i)-a(ie+i)) + &
                                 sin36*(a(ic+i)-a(id+i)) )
                d(je+i)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i)) &
                                -cos36*(b(ic+i)+b(id+i)) ) &
                              -( sin72*(a(ib+i)-a(ie+i)) + &
                                 sin36*(a(ic+i)-a(id+i)) )
                c(jc+i)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i)) &
                                +cos72*(a(ic+i)+a(id+i)) ) &
                              -( sin36*(b(ib+i)-b(ie+i)) - &
                                 sin72*(b(ic+i)-b(id+i)) )
                c(jd+i)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i)) &
                                +cos72*(a(ic+i)+a(id+i)) ) &
                              +( sin36*(b(ib+i)-b(ie+i)) - &
                                 sin72*(b(ic+i)-b(id+i)) )
                d(jc+i)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i)) &
                                +cos72*(b(ic+i)+b(id+i)) ) &
                              +( sin36*(a(ib+i)-a(ie+i)) - &
                                 sin72*(a(ic+i)-a(id+i)) )
                d(jd+i)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i)) &
                                +cos72*(b(ic+i)+b(id+i)) ) &
                              -( sin36*(a(ib+i)-a(ie+i)) - &
                                 sin72*(a(ic+i)-a(id+i)) )
            END DO
        END DO

        IF ( la == m ) RETURN

        la1=la+1
        DO k=la1,m,la
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
            DO l=1,la
                ll=k+l-1
                istart=(ll-1)*nv2+lm1
                istop=istart+lm2
                DO i=istart,istop
                    j=i+jadd
                    c(ja+j)=a(ia+i) + &
                           (a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
                    d(ja+j)=b(ia+i) + &
                           (b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
                    c(jb+j)=c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i)) &
                                        -cos36*(a(ic+i)+a(id+i)) ) &
                                      -( sin72*(b(ib+i)-b(ie+i)) &
                                       + sin36*(b(ic+i)-b(id+i)) ) ) &
                           -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i)) &
                                        -cos36*(b(ic+i)+b(id+i)) ) &
                                      +( sin72*(a(ib+i)-a(ie+i)) &
                                        +sin36*(a(ic+i)-a(id+i)) ) )
                    d(jb+j)=s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i)) &
                                        -cos36*(a(ic+i)+a(id+i)) ) &
                                      -( sin72*(b(ib+i)-b(ie+i)) &
                                        +sin36*(b(ic+i)-b(id+i)) ) ) &
                           +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i)) &
                                        -cos36*(b(ic+i)+b(id+i)) ) &
                                      +( sin72*(a(ib+i)-a(ie+i)) &
                                        +sin36*(a(ic+i)-a(id+i)) ) )
                    c(je+j)=c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i)) &
                                        -cos36*(a(ic+i)+a(id+i)) ) &
                                      +( sin72*(b(ib+i)-b(ie+i)) &
                                        +sin36*(b(ic+i)-b(id+i)) ) ) &
                           -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i)) &
                                        -cos36*(b(ic+i)+b(id+i)) ) &
                                      -( sin72*(a(ib+i)-a(ie+i)) &
                                        +sin36*(a(ic+i)-a(id+i)) ) )
                    d(je+j)=s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i)) &
                                        -cos36*(a(ic+i)+a(id+i)) ) &
                                      +( sin72*(b(ib+i)-b(ie+i)) &
                                        +sin36*(b(ic+i)-b(id+i)) ) ) &
                           +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i)) &
                                        -cos36*(b(ic+i)+b(id+i)) ) &
                                      -( sin72*(a(ib+i)-a(ie+i)) &
                                        +sin36*(a(ic+i)-a(id+i)) ) )
                    c(jc+j)=c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i)) &
                                        +cos72*(a(ic+i)+a(id+i)) ) &
                                      -( sin36*(b(ib+i)-b(ie+i)) &
                                        -sin72*(b(ic+i)-b(id+i)) ) ) &
                           -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i)) &
                                        +cos72*(b(ic+i)+b(id+i)) ) &
                                      +( sin36*(a(ib+i)-a(ie+i)) &
                                        -sin72*(a(ic+i)-a(id+i)) ) )
                    d(jc+j)=s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i)) &
                                        +cos72*(a(ic+i)+a(id+i)) ) &
                                      -( sin36*(b(ib+i)-b(ie+i)) &
                                        -sin72*(b(ic+i)-b(id+i)) ) ) &
                           +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i)) &
                                        +cos72*(b(ic+i)+b(id+i)) ) &
                                      +( sin36*(a(ib+i)-a(ie+i)) &
                                        -sin72*(a(ic+i)-a(id+i)) ) )
                    c(jd+j)=c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i)) &
                                        +cos72*(a(ic+i)+a(id+i)) ) &
                                      +( sin36*(b(ib+i)-b(ie+i)) &
                                        -sin72*(b(ic+i)-b(id+i)) ) ) &
                           -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i)) &
                                        +cos72*(b(ic+i)+b(id+i)) ) &
                                      -( sin36*(a(ib+i)-a(ie+i)) &
                                        -sin72*(a(ic+i)-a(id+i)) ) )
                    d(jd+j)=s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i)) &
                                        +cos72*(a(ic+i)+a(id+i)) ) &
                                      +( sin36*(b(ib+i)-b(ie+i)) &
                                        -sin72*(b(ic+i)-b(id+i)) ) ) &
                           +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i)) &
                                        +cos72*(b(ic+i)+b(id+i)) ) &
                                      -( sin36*(a(ib+i)-a(ie+i)) &
                                        -sin72*(a(ic+i)-a(id+i)) ) )
                END DO
            END DO
        END DO

    END IF

    RETURN
    end SUBROUTINE fft_fac

!-----------------------------------------------------------------------------------
