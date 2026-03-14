#P512 Structure
.set P512_M, 0
.set P512_DT, 64
.set P512_GR_PREFAC, 128
.set P512_GR_PREFAC2, 192
.set P512_m, 320
.set P512_X, 384
.set P512_Y, 448
.set P512_Z, 512
.set P512_VX, 576
.set P512_VY, 640
.set P512_VZ, 704
.set P512_MAT8_INERTIAL_TO_JACOBI, 768
.set P512_MAT8_JACOBI_TO_HELIOCENTRIC, 1280
.set P512_M0, 2304
.set P512_MASK, 2368

# Register use
.set X, %zmm10
.set Y, %zmm11
.set Z, %zmm12
.set VX, %zmm13
.set VY, %zmm14
.set VZ, %zmm15
.set ONE, %zmm16
.set DT, %zmm17
.set HALF_MASK, %zmm18
.set M, %zmm19

# Only used in Interaction step:
.set HVX, %zmm20
.set HVY, %zmm21
.set HVZ, %zmm22
.set HVXC, %zmm23
.set HVYC, %zmm24
.set HVZC, %zmm25
.set HX, %zmm26
.set HY, %zmm27
.set HZ, %zmm28

# Only used in Kepler step:
.set R, %zmm20
.set RI, %zmm21
.set ZETA, %zmm22
.set ETA, %zmm23
.set XX, %zmm24
.set GS0, %zmm25
.set GS1, %zmm26
.set GS2, %zmm27
.set GS3, %zmm28
.set BETA, %zmm29

