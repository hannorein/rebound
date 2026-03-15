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
.set P512_MAT8_INERTIAL_TO_JACOBI_T, 2304
.set P512_M0, 2816
.set P512_MASK, 2880


# Only used in Interaction step:
.set HVX, %zmm13
.set HVY, %zmm14
.set HVZ, %zmm15
.set HVXC, %zmm16
.set HVYC, %zmm17
.set HVZC, %zmm18
.set HX, %zmm19
.set HY, %zmm20
.set HZ, %zmm21

# Only used in Kepler step:
.set R, %zmm13
.set RI, %zmm14
.set ZETA, %zmm15
.set ETA, %zmm16
.set XX, %zmm17
.set GS0, %zmm18
.set GS1, %zmm0     # Note: reusing register zmm0
.set GS2, %zmm19
.set GS3, %zmm20
.set BETA, %zmm21

# Common register use
.set X, %zmm22
.set Y, %zmm23
.set Z, %zmm24
.set VX, %zmm25
.set VY, %zmm26
.set VZ, %zmm27
.set ONE, %zmm28
.set DT, %zmm29
.set HALF_MASK, %zmm30
.set M, %zmm31
