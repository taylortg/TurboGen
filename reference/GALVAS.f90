PROGRAM GALVAS
    ! Declaration of variables and array used in the model
    integer :: SPLT
    real :: LAMX, LOD, MUO, NONDES, MWRMS, INC, N
    real :: GAM, RGAS, POP, TOP
    real :: G1, G2, CP, FLFUNC, ROP, VCR
    real :: B2X, AL1MF, AL2, AL3, CHIH, CHIT
    real :: AMT(4) = [0.2, 0.4, 0.6, 0.8]
    real :: BARR(6) = [0.02, 0.04, 0.06, 0.08, .1, .12]

    ! These are all inputs from the paper itself
    SPLT = 1    ! 1 = splitter, 0 = no splitters
    ! Operating condition
    GAM     = 1.4
    RGAS    = 287.05
    POP     = 101325
    TOP     = 288.15
    ! Geometrical parameters
    B2X     = 25.0  ! Backsweep impeller angle
    AL1MF   = 0     ! Inducer meridional flow angle
    AL3     = 78.0  ! Vaned diffuser inlet angle
    CHIH    = 7.0   ! Inducer hub wall slope angle from axial
    CHIT    = 0     ! Inducer tip wall slope angle from axial


    ! Start of Galvas model
    G1 = GAM+1.
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

    write(*,*) "GAMMA = ", GAM
    write(*,*) "G1    = ", G1
    write(*,*) "G2    = ", G2
    write(*,*) "CP    = ", CP
    write(*,*) "FLFUNC= ", FLFUNC
    write(*,*) "ROP   = ", ROP
    write(*,*) "VOR   = ", VOR

END PROGRAM GALVAS