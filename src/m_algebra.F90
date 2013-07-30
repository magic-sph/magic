!$Id$
module algebra

  implicit none

  contains
!***********************************************************************
    subroutine cgesl(a,ia,n,ip,bc1)
!***********************************************************************
!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  does the backward substitution into a lu-decomposed real         |
!  |  matrix a (to solve a * x = bc1) were bc1 is the right hand side  |
!  |  vector. On return x is stored in bc1.                            |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

!-----------------------------------------------------------------------

    implicit none

!-- input:
    integer :: n         ! dimension of problem
    integer :: ia        ! first dim of a
    integer :: ip(*)     ! pivot pointer of legth n
    real(kind=8) :: a(ia,*)    ! real n X n matrix
    complex(kind=8) :: bc1(*) ! on input RHS of problem

!-- output: solution stored in bc1(*)

!-- local:
    integer :: nm1,nodd,i,m
    integer :: k,k1
    complex(kind=8) :: c1

!-- end of declaration
!---------------------------------------------------------------------------
     
    nm1=n-1
    nodd=mod(n,2)

!     permute vectors b1

    do k=1,nm1
        m=ip(k)
        c1=bc1(m)
        bc1(m)=bc1(k)
        bc1(k)=c1
    end do

!     solve  l * y = b

    do k=1,n-2,2
        k1=k+1
        bc1(k1)=bc1(k1)-bc1(k)*a(k1,k)
        do i=k+2,n
            bc1(i)=bc1(i)-(bc1(k)*a(i,k)+bc1(k1)*a(i,k1))
        end do
    end do

    if ( nodd == 0 ) then
        bc1(n)=bc1(n)-bc1(nm1)*a(n,nm1)
    end if

!     solve  u * x = y

    do k=n,3,-2
        k1=k-1
        bc1(k)=bc1(k)*a(k,k)
        bc1(k1)=(bc1(k1)-bc1(k)*a(k1,k))*a(k1,k1)
        do i=1,k-2
            bc1(i)=bc1(i)-bc1(k)*a(i,k)-bc1(k1)*a(i,k1)
        end do
    end do

    if ( nodd == 0 ) then
        bc1(2)=bc1(2)*a(2,2)
        bc1(1)=(bc1(1)-bc1(2)*a(1,2))*a(1,1)
    else
        bc1(1)=bc1(1)*a(1,1)
    end if
     
    return
    end subroutine cgesl
!--------------------------------------------------------------------
    SUBROUTINE cgeslML(a,ia,n,ip,bc,ldBc,nRHSs)
!  +-------------+----------------+------------------------------------+
!  |                                                                   |
!  |  does the backward substitution into a lu-decomposed real         |
!  |  matrix a (to solve a * x = bc ) simultaneously for nRHSs complex |
!  |  vectors bc. On return the results are stored in                  |
!  |  the bc.                                                          |
!  |                                                                   |
!  +-------------------------------------------------------------------+
!  |  ruler                                                            |
!  |5 7 10   15   20   25   30   35   40   45   50   55   60   65   70 |
!--++-+--+----+----+----+----+----+----+----+----+----+----+----+----+-+

!---------------------------------------------------------------------------

    IMPLICIT NONE

!-- input:
    INTEGER :: n           ! dimension of problem
    INTEGER :: ia          ! leading dimension of a
    INTEGER :: ip(*)       ! pivot pointer of legth n
    INTEGER :: ldBc        ! leading dimension of bc
    REAL(kind=8) :: a(ia,*)      ! real n X n matrix
    COMPLEX(kind=8) :: bc(ldBc,*) ! on input RHS of problem
    INTEGER :: nRHSs       ! number of right-hand sides

!-- local:
    INTEGER :: nm1,nodd,i,m
    INTEGER :: k,k1,nRHS,nRHS2,noddRHS
    COMPLEX(kind=8) :: help

!-- end of declaration
!---------------------------------------------------------------------------
     
    nm1=n-1
    nodd   =MOD(n,2)
    noddRHS=MOD(nRHSs,2)

!     permute vectors bc
    DO nRHS=1,nRHSs
        DO k=1,nm1
            m=ip(k)
            help       =bc(m,nRHS)
            bc(m,nRHS) =bc(k,nRHS)
            bc(k,nRHS) =help
        END DO
    END DO


!     solve  l * y = b

    DO nRHS=1,nRHSs-1,2
        nRHS2=nRHS+1

        DO k=1,n-2,2
            k1=k+1
            bc(k1,nRHS) =bc(k1,nRHS)-bc(k,nRHS)*a(k1,k)
            bc(k1,nRHS2)=bc(k1,nRHS2)-bc(k,nRHS2)*a(k1,k)
            DO i=k+2,n
                bc(i,nRHS) =bc(i,nRHS) - &
                    (bc(k,nRHS)*a(i,k)+bc(k1,nRHS)*a(i,k1))
                bc(i,nRHS2)=bc(i,nRHS2) - &
                    (bc(k,nRHS2)*a(i,k)+bc(k1,nRHS2)*a(i,k1))
            END DO
        END DO
        IF ( nodd == 0 ) THEN
            bc(n,nRHS) =bc(n,nRHS) -bc(nm1,nRHS)*a(n,nm1)
            bc(n,nRHS2)=bc(n,nRHS2)-bc(nm1,nRHS2)*a(n,nm1)
        END IF
    
    !     solve  u * x = y
    
        DO k=n,3,-2
            k1=k-1
            bc(k,nRHS)  =bc(k,nRHS)*a(k,k)
            bc(k1,nRHS) =(bc(k1,nRHS)-bc(k,nRHS)*a(k1,k))*a(k1,k1)
            bc(k,nRHS2) =bc(k,nRHS2)*a(k,k)
            bc(k1,nRHS2)=(bc(k1,nRHS2)-bc(k,nRHS2)*a(k1,k))*a(k1,k1)
            DO i=1,k-2
                bc(i,nRHS)=bc(i,nRHS) - &
                    bc(k,nRHS)*a(i,k)-bc(k1,nRHS)*a(i,k1)
                bc(i,nRHS2)=bc(i,nRHS2) - &
                    bc(k,nRHS2)*a(i,k)-bc(k1,nRHS2)*a(i,k1)
            END DO
        END DO
        IF ( nodd == 0 ) THEN
            bc(2,nRHS)=bc(2,nRHS)*a(2,2)
            bc(1,nRHS)=(bc(1,nRHS)-bc(2,nRHS)*a(1,2))*a(1,1)
            bc(2,nRHS2)=bc(2,nRHS2)*a(2,2)
            bc(1,nRHS2)=(bc(1,nRHS2)-bc(2,nRHS2)*a(1,2))*a(1,1)
        ELSE
            bc(1,nRHS)=bc(1,nRHS)*a(1,1)
            bc(1,nRHS2)=bc(1,nRHS2)*a(1,1)
        END IF

    END DO

    IF ( noddRHS == 1 ) THEN
        nRHS=nRHSs

        DO k=1,n-2,2
            k1=k+1
            bc(k1,nRHS)=bc(k1,nRHS)-bc(k,nRHS)*a(k1,k)
            DO i=k+2,n
                bc(i,nRHS)=bc(i,nRHS) - &
                    (bc(k,nRHS)*a(i,k)+bc(k1,nRHS)*a(i,k1))
            END DO
        END DO
        IF ( nodd == 0 ) &
            bc(n,nRHS)=bc(n,nRHS)-bc(nm1,nRHS)*a(n,nm1)
        DO k=n,3,-2
            k1=k-1
            bc(k,nRHS) =bc(k,nRHS)*a(k,k)
            bc(k1,nRHS)=(bc(k1,nRHS)-bc(k,nRHS)*a(k1,k))*a(k1,k1)
            DO i=1,k-2
                bc(i,nRHS)=bc(i,nRHS) - &
                    bc(k,nRHS)*a(i,k)-bc(k1,nRHS)*a(i,k1)
            END DO
        END DO
        IF ( nodd == 0 ) THEN
            bc(2,nRHS)=bc(2,nRHS)*a(2,2)
            bc(1,nRHS)=(bc(1,nRHS)-bc(2,nRHS)*a(1,2))*a(1,1)
        ELSE
            bc(1,nRHS)=bc(1,nRHS)*a(1,1)
        END IF

    END IF


    RETURN
    end SUBROUTINE cgeslML
!--------------------------------------------------------------------
    SUBROUTINE sgesl(a,ia,n,ip,b)
!---------------------------------------------------------------------------

!     like the linpack routine
!     backward substitution of vector b into lu-decomposed matrix a
!     to solve  a * x = b for a single real vector b

!     sub sgefa must be called once first to initialize a and ip

!     a: (input)  nxn real matrix
!     n: (input)  size of a and b
!     ip: (input) pivot pointer array of length n
!     b: (in/output) rhs-vector on input, solution on output

!     called in supdate, wpupdate, tcond, jcond

!---------------------------------------------------------------------------

    IMPLICIT NONE

!-- input:
    INTEGER :: n      ! dim of problem
    INTEGER :: ia     ! first dim of a
    INTEGER :: ip(*)  ! pivot information
    REAL(kind=8) :: a(ia,*),b(*)

!-- output: solution stored in b(n)

!-- local:
    INTEGER :: nm1,i
    INTEGER :: k,k1,m,nodd
    REAL(kind=8) :: help

!-- end of declaration
!---------------------------------------------------------------------------

    nm1 =n-1
    nodd=MOD(n,2)

!     permute vector b

    DO k=1,nm1
        m   =ip(k)
        help=b(m)
        b(m)=b(k)
        b(k)=help
    END DO

!     solve  l * y = b

    DO k=1,n-2,2
        k1=k+1
        b(k1)=b(k1)-b(k)*a(k1,k)
        DO i=k+2,n
            b(i)=b(i)-(b(k)*a(i,k)+b(k1)*a(i,k1))
        END DO
    END DO
    IF ( nodd == 0 ) b(n)=b(n)-b(nm1)*a(n,nm1)

!     solve  u * x = y

    DO k=n,3,-2
        k1=k-1
        b(k) =b(k)*a(k,k)
        b(k1)=(b(k1)-b(k)*a(k1,k))*a(k1,k1)
        DO i=1,k-2
            b(i)=b(i)-(b(k)*a(i,k)+b(k1)*a(i,k1))
        END DO
    END DO
    IF ( nodd == 0 ) THEN
        b(2)=b(2)*a(2,2)
        b(1)=(b(1)-b(2)*a(1,2))*a(1,1)
    ELSE
        b(1)=b(1)*a(1,1)
    END IF

     
    RETURN
    end SUBROUTINE sgesl
!----------------------------------------------------------------------------
    SUBROUTINE sgefa(a,ia,n,ip,info)
!-------------------------------------------------------------------------------
!     like the linpack routine

!     lu decomposes the real matrix a(n,n) via gaussian elimination

!     a: (in/output) real nxn matrix on input, lu-decomposed matrix on output
!     ia: (input) first dimension of a (must be >= n)
!     n: (input) 2nd dimension and rank of a
!     ip: (output) pivot pointer array
!     info: (output) error message when .ne. 0

!-------------------------------------------------------------------------------

    IMPLICIT NONE

!-- input:
    INTEGER :: ia,n
    REAL(kind=8) :: a(ia,*)

!-- output:
!       REAL(kind=8) a(ia,*)  ! modified
    INTEGER :: ip(*)   ! pivoting information
    INTEGER :: info

!-- local:
    INTEGER :: nm1,k,kp1,l,i,j
    REAL(kind=8) :: help

!-- end of declaration
!--------------------------------------------------------------------------


    IF ( n <= 1 ) STOP '45'

    info=0
    nm1 =n-1
     
    DO k=1,nm1
        kp1=k+1
        l  =k

        DO i=kp1,n
            IF ( DABS(a(i,k)) > DABS(a(l,k)) ) l=i
        END DO

        ip(k)=l

        IF( a(l,k) /= 0.D0 ) THEN

            IF ( l /= k ) THEN
                DO i=1,n
                    help  =a(k,i)
                    a(k,i)=a(l,i)
                    a(l,i)=help
                END DO
            END IF

            help=1.D0/a(k,k)
            DO i=kp1,n
                a(i,k)=help*a(i,k)
            END DO

            DO j=kp1,n
                DO i=kp1,n
                    a(i,j)=a(i,j)-a(k,j)*a(i,k)
                END DO
            END DO

        ELSE
            info=k
        END IF
         
    END DO

    ip(n)=n
    IF( a(n,n) == 0.D0 ) info=n
    IF( info > 0 ) RETURN
    DO i=1,n
        a(i,i)=1.D0/a(i,i)
    END DO
     
    RETURN
    end SUBROUTINE sgefa
!-----------------------------------------------------------------------------
end module algebra
