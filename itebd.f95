! iTEBD, 1D quantum Ising chain with transverse magnetic field

PROGRAM itebd
USE myModule
IMPLICIT NONE
! parameters need to run itebd
    REAL(KIND=PREC), PARAMETER :: DELTA = 0.005_PREC
    INTEGER, PARAMETER :: NSTEPS = 5000
    INTEGER, PARAMETER :: D = 5
    REAL(KIND=PREC) :: hz = 0.5
    REAL(KIND=PREC), DIMENSION(4,4) :: hamiltonian
    REAL(KIND=PREC), DIMENSION(2,2,2,2) :: evolution
    REAL(KIND=PREC), DIMENSION(D,2,D,2) :: groundState
    REAL(KIND=PREC), DIMENSION(D,2) :: l
    REAL(KIND=PREC), DIMENSION(D,2,2,D) :: theta
! variables need to run SVD
    REAL(KIND=PREC), DIMENSION(2*D,2*D) :: u
    REAL(KIND=PREC), DIMENSION(2*D) :: sva
    REAL(KIND=PREC), DIMENSION(2*D,2*D) :: v
    REAL(KIND=PREC), DIMENSION(4*D) :: work
    INTEGER :: info
! some auxiliary variables
    INTEGER :: step
    INTEGER :: a, b
    INTEGER :: h, i, j, k

    CALL init_random_seed()
    CALL RANDOM_NUMBER( groundState )
    CALL RANDOM_NUMBER( l )


    hamiltonian = RESHAPE( (/ hz, 0.0_PREC, 0.0_PREC, 1.0_PREC, &
                              0.0_PREC, hz, 1.0_PREC, 0.0_PREC, &
                              0.0_PREC, 1.0_PREC, -hz, 0.0_PREC, &
                              1.0_PREC, 0.0_PREC, 0.0_PREC, -hz /), (/ 4, 4 /) )

    evolution = RESHAPE( expm( -DELTA * hamiltonian ), (/ 2, 2, 2, 2 /) )

    DO step = 1, NSTEPS
        a = MOD( step, 2 ) + 1
        b = MOD( step + 1, 2 ) + 1

        FORALL ( i = 1:2, j = 1:D )
            groundState( :, i, j, a ) = l( :, b ) * groundState( :, i, j, a )
            groundState( j, i, :, a ) = l( :, a ) * groundState( j, i, :, a )
            groundState( j, i, :, b ) = l( :, b ) * groundState( j, i, :, b )
        END FORALL

        FORALL ( i = 1:2, j = 1:2, h = 1:D, k = 1:D )
            theta( h, i, j, k ) = DOT_PRODUCT( groundState( h , i, :, a ), groundState( :, j, k, b ) )
            theta( h, i, j, k ) = SUM( theta( h, :, :, k ) * evolution( i, j, :, : ) )
        END FORALL

        u = RESHAPE( theta, (/ 2 * D, 2 * D /) )

        CALL DGESVJ( 'G', 'U', 'V', 2 * D, 2 * D, u, 2 * D, sva, 2 * D, v, 2 * D, work, 4 * D, info )

        l( :, a ) = sva( 1:D ) / SQRT( SUM( sva( 1:D ) ** 2 ) )
        groundState( :, :, :, a ) = RESHAPE( u( 1:, 1:D ), (/ D, 2, D /) )
        groundState( :, :, :, b ) = RESHAPE( TRANSPOSE( v( 1:, 1:D ) ), (/ D, 2, D /) )

        FORALL ( i = 1:2, j = 1:D )
            groundState( :, i, j, a ) = ( l( :, b ) ** ( -1 ) ) * groundState( :, i, j, a )
            groundState( j, i, :, b ) = ( l( :, b ) ** ( -1 ) ) * groundState( j, i, :, b )
        END FORALL
    END DO


    WRITE (*, '(F15.10)' ) -LOG( SUM( theta ** 2 ) ) / DELTA / 2.0_PREC
END PROGRAM itebd
