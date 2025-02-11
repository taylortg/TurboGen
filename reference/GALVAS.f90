PROGRAM GALVAS
    ! Declaration of variables and array used in the model
    integer :: SPLT, MIN, NVOVCR, K
    real :: LAMX, LOD, MUO, NONDES, MWRMS, INC, N
    real :: GAM, RGAS, POP, TOP, J, L, SW, THETA, DELTA
    real :: G1, G2, CP, FLFUNC, ROP, VCR, U1T
    real :: B2X, AL1MF, AL2, AL3, CHIH, CHIT, B1MFB
    real :: AMT(4) = [0.2, 0.4, 0.6, 0.8]
    real :: BARR(6) = [0.02, 0.04, 0.06, 0.08, .1, .12]

    ! These are all inputs from the paper itself
    SPLT = 1    ! 1 = splitter, 0 = no splitters
    NVOVCR = 15 ! Number of V/Vcr ratios to run
    ! Operating condition
    GAM     = 1.4
    RGAS    = 287.05
    POP     = 101325
    TOP     = 288.15
    N       = 72000.
    NONDES  = 1.
    ! Geometrical parameters
    B2X     = 25.       ! Backsweep impeller angle
    AL1MF   = 0.        ! Inducer meridional flow angle
    AL3     = 78.       ! Vaned diffuser inlet angle
    CHIH    = 7.        ! Inducer hub wall slope angle from axial
    CHIT    = 0.        ! Inducer tip wall slope angle from axial
    D1T     = 0.0813    ! Impeller tip diameter at inlet
    B1MFB   = 49.       ! RMS blade angle at impeller inlet


    ! Start of Galvas model
1   G1 = GAM+1.
    G2 = GAM-1.
    CP = GAM*RGAS/G2
    B2X = 0.01745*B2X
    AL1MF = 0.01745*AL1MF
    AL3 = 0.01745*AL3
    CHIH = 0.01745*CHIH
    CHIT = 0.01745*CHIT
    FLFUNC = SQRT(GAM/RGAS*(2./G1)**(G1/G2))
    ROP = POP/RGAS/TOP
    VCR = SQRT(2.*GAM/G1*RGAS*TOP)
    U1T = 3.14159*N*NONDES*D1T/60.
    DO 18 L=1,2
    MIN = 1
2   DO 17 J=MIN,NVOVCR
    K = 0


17  CONTINUE
    IF(L.EQ.1) GO TO 18
18  WCHK = SW*SQRT(THETA)/DELTA
    GO TO 1

    write(*,*) "GAMMA = ", GAM
    write(*,*) "G1    = ", G1
    write(*,*) "G2    =  ", G2
    write(*,*) "CP    = ", CP
    write(*,*) "FLFUNC= ", FLFUNC
    write(*,*) "ROP   = ", ROP
    write(*,*) "VCR   = ", VCR
    write(*,*) "U1T   = ", U1T

END PROGRAM GALVAS

