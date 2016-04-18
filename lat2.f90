


      TYPE latt
         REAL(q) :: SCALE
         REAL(q) :: A(3,3),B(3,3)
         REAL(q) :: ANORM(3),BNORM(3)
         REAL(q) :: OMEGA

      END TYPE
      
      
      DO I=1,3
        READ(15,*,ERR=147,END=147) LATT_CUR%A(1,I),LATT_CUR%A(2,I),LATT_CUR%A(3,I)
      ENDDO
      
      
      SUBROUTINE LATTIC(Mylatt)
      USE prec

      IMPLICIT NONE

      TYPE(LATT) Mylatt
      REAL(q) Omega
      INTEGER I,J
      INTRINSIC SUM

      CALL EXPRO(Mylatt%B(1:3,1),Mylatt%A(1:3,2),Mylatt%A(1:3,3))
      CALL EXPRO(Mylatt%B(1:3,2),Mylatt%A(1:3,3),Mylatt%A(1:3,1))
      CALL EXPRO(Mylatt%B(1:3,3),Mylatt%A(1:3,1),Mylatt%A(1:3,2))

      Omega =Mylatt%B(1,1)*Mylatt%A(1,1)+Mylatt%B(2,1)*Mylatt%A(2,1) &
     &      +Mylatt%B(3,1)*Mylatt%A(3,1)

      DO I=1,3
      DO J=1,3
        Mylatt%B(I,J)=Mylatt%B(I,J)/Omega
      ENDDO
      ENDDO

      DO I=1,3
        Mylatt%ANORM(I)=SQRT(SUM(Mylatt%A(:,I)*Mylatt%A(:,I)))
        Mylatt%BNORM(I)=SQRT(SUM(Mylatt%B(:,I)*Mylatt%B(:,I)))
      ENDDO
      Mylatt%Omega=Omega
      RETURN
      END SUBROUTINE
      
      


! caclulates the x-product of two vectors

      SUBROUTINE EXPRO(H,U1,U2)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION H(3),U1(3),U2(3)

      H(1)=U1(2)*U2(3)-U1(3)*U2(2)
      H(2)=U1(3)*U2(1)-U1(1)*U2(3)
      H(3)=U1(1)*U2(2)-U1(2)*U2(1)

      RETURN
      END SUBROUTINE
      


! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)

      SUBROUTINE KARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION V(3,NMAX),BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE



! transform a set of vectors from
! ) direct lattice      (BASIS must be equal to A direct lattice)
! ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
! to cartesian coordinates

      SUBROUTINE DIRKAR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION V(3,NMAX),BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(1,2)+V(3,N)*BASIS(1,3)
        V2=V(1,N)*BASIS(2,1)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(2,3)
        V3=V(1,N)*BASIS(3,1)+V(2,N)*BASIS(3,2)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE


! bring all ions into the primitive cell

      SUBROUTINE TOPRIM(NIONS,POSION)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      DIMENSION POSION(3,NIONS)

      DO I=1,NIONS
      POSION(1,I)=MOD(POSION(1,I)+60,1._q)
      POSION(2,I)=MOD(POSION(2,I)+60,1._q)
      POSION(3,I)=MOD(POSION(3,I)+60,1._q)
      ENDDO
      END SUBROUTINE

      END MODULE



! transform a set of vectors from cartesian coordinates to
! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
! ) reciprocal lattice  (BASIS must be equal to A direct lattice)

      SUBROUTINE CKARDIR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(3,NMAX),V1,V2,V3
      DIMENSION  BASIS(3,3)

      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(2,1)+V(3,N)*BASIS(3,1)
        V2=V(1,N)*BASIS(1,2)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(3,2)
        V3=V(1,N)*BASIS(1,3)+V(2,N)*BASIS(2,3)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE


! transform a set of vectors from
! ) direct lattice      (BASIS must be equal to A direct lattice)
! ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
! to cartesian coordinates

      SUBROUTINE CDIRKAR(NMAX,V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      COMPLEX(q) V(3,NMAX),V1,V2,V3
      DIMENSION  BASIS(3,3)
      
      DO N=1,NMAX
        V1=V(1,N)*BASIS(1,1)+V(2,N)*BASIS(1,2)+V(3,N)*BASIS(1,3)
        V2=V(1,N)*BASIS(2,1)+V(2,N)*BASIS(2,2)+V(3,N)*BASIS(2,3)
        V3=V(1,N)*BASIS(3,1)+V(2,N)*BASIS(3,2)+V(3,N)*BASIS(3,3)
        V(1,N)=V1
        V(2,N)=V2
        V(3,N)=V3
      ENDDO

      RETURN
      END SUBROUTINE
      
      
