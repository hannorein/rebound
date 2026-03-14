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
.set X, %zmm29
.set Y, %zmm30
.set Z, %zmm31
.set VX, %zmm10
.set VY, %zmm11
.set VZ, %zmm12
.set ONE, %zmm26
.set DT, %zmm27
.set HALF_MASK, %zmm28
.set M, %zmm16

# Only used in Interaction step:
.set HVX, %zmm13
.set HVY, %zmm14
.set HVZ, %zmm15
.set HX, %zmm17
.set HY, %zmm18
.set HZ, %zmm19

# Onlt used in Kepler step:
.set R, %zmm9
.set RI, %zmm13
.set ZETA, %zmm14
.set ETA, %zmm15
.set XX, %zmm17
.set GS0, %zmm20
.set GS1, %zmm21
.set GS2, %zmm22
.set GS3, %zmm23
.set BETA, %zmm24

