MODULE myModule
IMPLICIT NONE
SAVE
    INTEGER, PARAMETER :: PREC = KIND( 0.0D01 )
    REAL(KIND=PREC), PARAMETER :: EPS = EPSILON( 0.0_PREC )

CONTAINS
    FUNCTION expm( a )
        REAL(KIND=PREC), INTENT(IN), DIMENSION(:,:) :: a
        REAL(KIND=PREC), ALLOCATABLE, DIMENSION(:,:) :: expm
        REAL(KIND=PREC), ALLOCATABLE, DIMENSION(:,:) :: tmp
        REAL(KIND=PREC) :: prod
        INTEGER :: n
        INTEGER :: i

        n = SIZE( a, 1 )

        ALLOCATE ( expm(n,n) )
        expm = 0.0_PREC
        FORALL ( i = 1:n )
            expm( i, i ) = 1.0_PREC
        END FORALL


        tmp = a
        prod = 1.0_PREC
        DO
            IF ( SUM( tmp ** 2 ) < EPS ) EXIT
            expm = expm + tmp

            prod = prod + 1.0_PREC
            tmp = MATMUL( a, tmp ) / prod
        END DO

        DEALLOCATE ( tmp )
    END FUNCTION expm

    SUBROUTINE init_random_seed()
        INTEGER, ALLOCATABLE, DIMENSION(:) :: seed
        INTEGER :: clock
        INTEGER :: n
        INTEGER :: i

        CALL RANDOM_SEED( SIZE = n )
        ALLOCATE ( seed(n) )

        CALL SYSTEM_CLOCK( COUNT = clock )

        seed = clock + 37 * (/ ( i - 1, i = 1, n ) /)
        CALL RANDOM_SEED( PUT = seed )

        DEALLOCATE ( seed )
    END SUBROUTINE init_random_seed
END MODULE myModule
