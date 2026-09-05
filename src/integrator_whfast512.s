# file: integrator_whfast512.s
# Assembler code for WHFast512.
# This uses the GNU assembly syntax.
# 
# Copyright (c) 2026 Rishit Dagli, Hanno Rein
# 
# This file is part of rebound.
# 
# rebound is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rebound is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rebound.  If not, see <http://www.gnu.org/licenses/>.
#
#
.section .text
# Enable debug counter?
.equ DEBUG_AVX512, 0

#P512 Structure offsets
.set P512_M, 0
.set P512_DT, 64
.set P512_GR_PREFAC, 128
.set P512_m, 192
.set P512_X, 256
.set P512_Y, 320
.set P512_Z, 384
.set P512_VX, 448
.set P512_VY, 512
.set P512_VZ, 576
.set P512_MAT8_INERTIAL_TO_JACOBI, 640
.set P512_MAT8_JACOBI_TO_HELIOCENTRIC, 1152
.set P512_M0, 1664
.set P512_MASK, 1728
.set P512_EXIT_MAX_DISTANCE, 1792
.set P512_EXIT_MIN_DISTANCE, 1856
.set P512_EXIT_MIN_DISTANCE_R, 1920
.set P512_COUNTER, 2496

#####################################
# Register use
#                           Kepler   Interaction   Other
# 64    32    16   8      
# rax   eax   ax   ah,al                           return value 
# rbx   ebx   bx   bh,bl  
# rcx   ecx   cx   ch,cl                           Interrupt pointer
# rdx   edx   dx   dh,dl                           skip_first_kepler, corrector
# rsi   esi   si   sil                             number_of_steps 
# rdi   edi   di   dil      ------------pointer to simd_data--------- 
# rbp   ebp   bp   bpl      ------------frame pointer----------------
# rsp   esp   sp   spl      ------------stack pointer----------------
# r8    r8d   r8w  r8b                             corrector, status flag for encounters, ejections 
# r9    r9d   r9w  r9b      Netwon   matmul, exceptions        
# r10   r10d  r10w r10b                            corrector, step counter


#####################################
# Register alias
#####################################

# Only used in Interaction step:
.set .LHVX, %zmm13
.set .LHVY, %zmm14
.set .LHVZ, %zmm15
.set .LHVXC, %zmm16
.set .LHVYC, %zmm17
.set .LHVZC, %zmm18
.set .LHX, %zmm19
.set .LHY, %zmm20
.set .LHZ, %zmm21

# Only used in Kepler step:
.set .LR, %zmm13
.set .LIR, %zmm14
.set .LZETA, %zmm15
.set .LETA, %zmm16
.set .LXX, %zmm17
.set .LGS0, %zmm18
.set .LGS1, %zmm0     # Note: reusing register zmm0
.set .LGS2, %zmm19
.set .LGS3, %zmm20
.set .LBETA, %zmm21

# Common register use
.set .LX, %zmm22
.set .LY, %zmm23
.set .LZ, %zmm24
.set .LVX, %zmm25
.set .LVY, %zmm26
.set .LVZ, %zmm27
.set .LONE, %zmm28
.set .LDT, %zmm29
.set .LHALF, %zmm30
.set .LM, %zmm31
.set .LM_DT, %zmm12           # Only used once per step.
.set .LEPS, %zmm11
.set .LSIGN_ABS_MASK, %zmm10  # Only used once per step.
.set .LMM0_DT, %zmm9          # -dt*M0 Only used once per step.

#####################################
# Stack initialization 
#####################################
.macro alloc_stack64 bytes
    # Safe old base pointer
    pushq   %rbp
    movq    %rsp, %rbp
    # realign stack to nearest multiple of 64 bytes
    andq    $-64, %rsp
    subq    $\bytes, %rsp
.endm

.macro free_stack64
    # restore old base pointer
    movq    %rbp, %rsp
    popq    %rbp
.endm

#####################################
# Register initialization 
#####################################
.macro reb_whfast512_init_registers
    # Ignore exceptions
    subq $8, %rsp
    stmxcsr (%rsp)
    orl $0x8040, (%rsp)    # Set Global FTZ and DAZ
    ldmxcsr (%rsp)
    addq $8, %rsp

    # Load data
    kmovw           P512_MASK(%rdi), %k1
    vmovapd         P512_DT(%rdi), .LDT
    vmovapd         P512_M(%rdi), .LM
    vmulpd          .LDT, .LM, .LM_DT
    vmovapd         P512_M0(%rdi), .LMM0_DT
    vmulpd          .LDT, .LMM0_DT, .LMM0_DT
    vxorpd          .SIGN_FLIP_MASK(%rip){1to8}, .LMM0_DT, .LMM0_DT
    vbroadcastsd    .DOUBLE_ONE(%rip), .LONE
    vbroadcastsd    ..LHALF(%rip), .LHALF
    vbroadcastsd    .EPS(%rip), .LEPS
    vbroadcastsd    .SIGN_ABS_MASK(%rip), .LSIGN_ABS_MASK
    
    vmovapd     P512_X(%rdi), .LX
    vmovapd     P512_Y(%rdi), .LY
    vmovapd     P512_Z(%rdi), .LZ
    
    vmovapd     P512_VX(%rdi), .LVX
    vmovapd     P512_VY(%rdi), .LVY
    vmovapd     P512_VZ(%rdi), .LVZ
.endm  

.macro reb_whfast512_store_results
    vmovapd    .LVX, P512_VX(%rdi)
    vmovapd    .LVY, P512_VY(%rdi)
    vmovapd    .LVZ, P512_VZ(%rdi)
    vmovapd    .LX, P512_X(%rdi)
    vmovapd    .LY, P512_Y(%rdi)
    vmovapd    .LZ, P512_Z(%rdi)
.endm


#####################################
# Macros for Kepler step
#####################################
.macro vfnmadd_auto_inc reg1 reg2
    vfnmadd213pd    .IF0+(IF_offset*8)(%rip){1to8}, \reg1, \reg2
    .set IF_offset, IF_offset - 1
.endm

# High accuracy: (Gs1, Gs2, Gs3)
# Output: .LGS1==%zmm0, .LGS2, .LGS3
# numTerms must be an odd number
.macro mm_stiefel_Gs13_avx512 numTerms=19
    .set IF_offset, \numTerms
    vmulpd          .LXX, .LXX, %zmm2     # X^2
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm3
    .set IF_offset, IF_offset - 1
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm4
    .set IF_offset, IF_offset - 1
    vmulpd          %zmm2, .LBETA, %zmm0
    .set GS_iterations, (IF_offset -1)/2
    .rept GS_iterations
    vfnmadd_auto_inc %zmm0, %zmm3
    vfnmadd_auto_inc %zmm0, %zmm4
    .endr
    vmulpd          %zmm4, %zmm2, .LGS2
    vmulpd          %zmm3, .LXX, %zmm3
    vmulpd          %zmm3, %zmm2, .LGS3
    vfnmadd132pd    %zmm3, .LXX, %zmm0 # = .LGS1
.endm

.macro comp_horner_step C, CLO, ifoff, has_err=1
    vmulpd          %zmm0, \C, %zmm7                 # P = z*C
    vmovapd         %zmm7, %zmm8                     # zmm8 = P
    vfmsub231pd     %zmm0, \C, %zmm8                 # pe = z*C - P
    vbroadcastsd    .IF0+(\ifoff*8)(%rip), %zmm1     # IF[k]
    vsubpd          %zmm7, %zmm1, \C                 # Cn = IF[k] - P
    vsubpd          \C, %zmm1, %zmm1                 # IF[k] - Cn
    vsubpd          %zmm7, %zmm1, %zmm1              # se = (IF[k]-Cn) - P
    .if \has_err
    vbroadcastsd    .IF0_err+(\ifoff*8)(%rip), %zmm7 # err[k] = rounded-true
    vsubpd          %zmm7, %zmm1, %zmm1             # se - err
    .endif
    vsubpd          %zmm8, %zmm1, %zmm7             # (se[-err]) - pe
    vfnmadd213pd    %zmm7, %zmm0, \CLO             # \CLO = -z*\CLO + (se[-err]) - pe
.endm

# similar to mm_stiefel_Gs13_avx512
.macro mm_stiefel_Gs13_comp numTerms=19
    .set IF_offset, \numTerms
    vmulpd          .LXX, .LXX, %zmm2     # X^2
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm3
    .set IF_offset, IF_offset - 1
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm4
    .set IF_offset, IF_offset - 1
    vmulpd          %zmm2, .LBETA, %zmm0
    .set GS_iterations, (IF_offset -1)/2 - 2     # leave last 2 terms
    .rept GS_iterations
    vfnmadd_auto_inc %zmm0, %zmm3
    vfnmadd_auto_inc %zmm0, %zmm4
    .endr
    vxorpd          %zmm5, %zmm5, %zmm5          # cs3 lo = 0
    vxorpd          %zmm6, %zmm6, %zmm6          # cs2 lo = 0
    comp_horner_step %zmm3, %zmm5, 5
    comp_horner_step %zmm4, %zmm6, 4
    comp_horner_step %zmm3, %zmm5, 3
    comp_horner_step %zmm4, %zmm6, 2, 0          # .IF0_err[2]==0
    vaddpd          %zmm5, %zmm3, %zmm3
    vaddpd          %zmm6, %zmm4, %zmm4
    vmulpd          %zmm4, %zmm2, .LGS2
    vmulpd          %zmm3, .LXX, %zmm3
    vmulpd          %zmm3, %zmm2, .LGS3
    vfnmadd132pd    %zmm3, .LXX, %zmm0 # = .LGS1
.endm

# Low accuracy: (Gs0, Gs1, Gs2, Gs3)
# Output: .LGS0, .LGS1==%zmm0, .LGS2, .LGS3
# numTerms must be an odd number
.macro mm_stiefel_Gs03_avx512 numTerms=11
    .set IF_offset, \numTerms
    vmulpd          .LXX, .LXX, %zmm2     # X^2
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm3
    .set IF_offset, IF_offset - 1
    vbroadcastsd    .IF0+(IF_offset*8)(%rip), %zmm4
    .set IF_offset, IF_offset - 1
    vmulpd          %zmm2, .LBETA, %zmm0
    .set GS_iterations, (IF_offset -1)/2 -1
    .rept GS_iterations
    vfnmadd_auto_inc %zmm0, %zmm3
    vfnmadd_auto_inc %zmm0, %zmm4
    .endr
    vfnmadd213pd    .IF0+(3*8)(%rip){1to8}, %zmm0, %zmm3
    vmovapd         %zmm4, .LGS0
    vfnmadd213pd    .IF0+(2*8)(%rip){1to8}, %zmm0, %zmm4
    vfnmadd213pd    .IF0+(1*8)(%rip){1to8}, %zmm0, .LGS0
    vmulpd          %zmm4, %zmm2, .LGS2
    vmulpd          %zmm3, .LXX, %zmm3
    vmulpd          %zmm3, %zmm2, .LGS3
    vfnmadd132pd    %zmm3, .LXX, %zmm0 # = .LGS1
.endm

.macro halley
    # In: .LGS0,.LGS1,.LGS2,.LGS3
    # Out: .LXX
    # No other registers used. Destroys input.
    vfmsub213pd     .LDT, .LZETA, .LGS3
    vfmadd231pd     .LGS2, .LETA, .LGS3
    vfmadd231pd     .LXX, .LR, .LGS3              # f

    vfmadd132pd     .LZETA, .LR, .LGS2
    vfmadd231pd     .LETA, .LGS1, .LGS2           # fp

    vmulpd          .LGS0, .LETA, .LGS0
    vfmadd132pd     .LZETA, .LGS0, .LGS1          # fpp

    vmulpd          .LGS1, .LGS3, .LGS1           # f*fpp
    # 0.5*f*fpp via integer subtract of 1 from the exponent (4-cycle vmulpd -> 1-cycle vpsubq)
    vpsubq          .HALF_EXP_DECR(%rip){1to8}, .LGS1, .LGS1
    vfmsub231pd     .LGS2, .LGS2, .LGS1           # fp*fp-0.5*f*fpp
    vmulpd          .LGS3, .LGS2, .LGS3           # f*fp
    vdivpd          .LGS1, .LGS3, .LGS3
    vsubpd          .LGS3, .LXX, .LXX
.endm

.macro newton
    # In: .LGS1,.LGS2,.LGS3
    # Out: .LXX
    # No other registers used. Destroys input.
    vmulpd          .LGS1, .LETA, .LGS1
    vfmadd231pd     .LGS2, .LZETA, .LGS1
    vmulpd          .LGS1, .LXX, .LXX{%k4}
    vfnmadd132pd    .LETA, .LXX, .LGS2
    vaddpd          .LR, .LGS1, .LXX{%k4}
    vdivpd          .LXX, .LONE, .LXX{%k4}    # TODO: Hot spot
    vfnmadd231pd    .LGS3, .LZETA, .LGS2
    vaddpd          .LGS2, .LDT, .LGS2
    vmulpd          .LGS2, .LXX, .LXX{%k4}
.endm

###############################################################################
# Kepler Step
###############################################################################

.macro kepler_step 
    vmulpd          .LX, .LX, %zmm0
    vmulpd          .LVX, .LVX, %zmm1
    vfmadd231pd     .LY, .LY, %zmm0
    vfmadd231pd     .LVY, .LVY, %zmm1
    vfmadd231pd     .LZ, .LZ, %zmm0                 # r^2
    vfmadd231pd     .LVZ, .LVZ, %zmm1               # v^2
    vsqrtpd         %zmm0, .LR                    # r
    vdivpd          .LR, .LONE, .LIR                  # 1/r
    vaddpd          .LM, .LM, .LBETA                  # 2*M
    vfmsub132pd     .LIR, %zmm1, .LBETA             # beta
    vmulpd          .LVX, .LX, .LETA
    vfmadd231pd     .LVY, .LY, .LETA
    vfmadd231pd     .LVZ, .LZ, .LETA                  # eta
    vmovapd         .LBETA, .LZETA
    vfnmadd132pd    .LR, .LM, .LZETA                  # zeta
    vmulpd          .LIR, .LDT, .LXX                  # dt/r  = first order guess for .LXX
    # Second order alternative
    #    vmulpd          .LETA, %zmm5, %zmm4      # eta*dt/r
    #    vmulpd          .LHALF, %zmm4, %zmm4     # 0.5*eta*dt/r
    #    vfnmadd132pd    .LIR, .LONE, %zmm4        
    #    vmulpd          %zmm5, %zmm4, .LXX       # .LXX (second order initial guess)
    
    # Iterations to improve X
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mm_stiefel_Gs03_avx512 9
    halley
    mm_stiefel_Gs03_avx512 11
    halley
  
    movq            $0, %r9                     # Newton loop counter
    kxnorw          %k4, %k4, %k4               # k4 = all lanes active
.NewtonLoop\@:
    vmovapd         .LXX,     %zmm7               # Store old .LXX
    mm_stiefel_Gs13_avx512
    newton                                      # only updates .LXX for lanes still in k4

    vsubpd          .LXX, %zmm7, %zmm7            # Delta .LXX
    vpandq          .LSIGN_ABS_MASK, %zmm7, %zmm7 # abs(Delta .LXX)

    # Required precision reached? abs(Delta .LXX) < eps
    vcmppd          $0x11, %zmm7, .LEPS, %k4      # $11 = less than, ordered (nans fail), quiet, k4=1 for failed particles
    #vcmppd         $25, %zmm7, .LEPS, %k 4       # $25 = Not greater or equal, unordered (nans pass), quiet
.if DEBUG_AVX512 == 1 
    vmovdqa64       P512_COUNTER(%rdi), %zmm4
    vpaddq          .ONE_QUAD(%rip){1to8}, %zmm4, %zmm4{%k4}
    vmovdqa64       %zmm4, P512_COUNTER(%rdi)
.endif

    kortestw        %k4, %k4
    jz              .NewtonLoopDone\@

    # Maximum iterations reached?
    incq            %r9 
    cmpq            $5, %r9                     # max Newton iterations
    jne             .NewtonLoop\@

    # If not converged yet, fall back to bisection
    movq            $0, %r9 
    vxorpd          %zmm5, %zmm5, %zmm5         # X_MIN = 0

    # select lanes that require fallback
    vxorpd          %zmm0, %zmm0, %zmm0
    vcmppd          $0x12, %zmm0, .LBETA, %k5    # beta <= 0
    kandw           %k4, %k5, %k5
    knotw           %k5, %k6
    kandw           %k4, %k6, %k6

    # Elliptic bounds.
    vsqrtpd         .LBETA, %zmm1{%k6}{z}
    vmulpd          %zmm1, .LBETA, %zmm2{%k6}{z}  # sqrt(.LBETA)*.LBETA
    vmulpd          .TWOPI(%rip){1to8}, .LM, %zmm3
    vdivpd          %zmm3, %zmm2, %zmm2{%k6}{z}  # invperiod
    vmulpd          %zmm2, .LDT, %zmm2{%k6}{z}
    vrndscalepd     $0x1, %zmm2, %zmm2{%k6}{z}   # floor(dt*invperiod)

    vbroadcastsd    .TWOPI(%rip), %zmm3
    vdivpd          %zmm1, %zmm3, %zmm1{%k6}{z}  # X_per_period = 2*pi/sqrt(.LBETA)
    vmulpd          %zmm1, %zmm2, %zmm5{%k6}     # X_MIN = X_per_period*floor(dt_invperiod)
    vaddpd          %zmm1, %zmm5, %zmm1{%k6}     # X_MAX = X_MIN + X_per_period

    # Hyperbolic bounds.
    kortestw        %k5, %k5
    jz              .FallbackBoundsDone\@

    vmulpd          .LY, .LVZ, %zmm6{%k5}{z}
    vfnmadd231pd    .LZ, .LVY, %zmm6{%k5}        # h_x = y*vz - z*vy
    vmulpd          .LZ, .LVX, %zmm7{%k5}{z}
    vfnmadd231pd    .LX, .LVZ, %zmm7{%k5}        # h_y = z*vx - x*vz
    vmulpd          .LX, .LVY, %zmm8{%k5}{z}
    vfnmadd231pd    .LY, .LVX, %zmm8{%k5}        # h_z = x*vy - y*vx
    vmulpd          %zmm6, %zmm6, %zmm6{%k5}
    vfmadd231pd     %zmm7, %zmm7, %zmm6{%k5}
    vfmadd231pd     %zmm8, %zmm8, %zmm6{%k5}     # h^2
    vmovapd         %zmm6, %zmm7{%k5}{z}

    vmulpd          .LM, .LM, %zmm8{%k5}{z}
    vmulpd          .LBETA, %zmm7, %zmm2{%k5}{z}
    vdivpd          %zmm8, %zmm2, %zmm2{%k5}{z}
    vsubpd          %zmm2, .LONE, %zmm2{%k5}{z}  # 1 - h^2*beta/M^2
    vsqrtpd         %zmm2, %zmm2{%k5}{z}
    vaddpd          .LONE, %zmm2, %zmm2{%k5}{z}  # 1 + e

    vdivpd          .LM, %zmm7, %zmm6{%k5}{z}
    vdivpd          %zmm2, %zmm6, %zmm6{%k5}{z}  # q = h^2/M/(1+e)
    vsqrtpd         %zmm7, %zmm7{%k5}{z}
    vdivpd          %zmm6, %zmm7, %zmm7{%k5}{z}  # vq = sqrt(h^2)/q

    vmulpd          .LDT, %zmm7, %zmm8{%k5}{z}
    vpandq          .LSIGN_ABS_MASK, %zmm8, %zmm8
    vaddpd          .LR, %zmm8, %zmm8{%k5}{z}
    vdivpd          %zmm8, .LDT, %zmm5{%k5}      # X_MIN = dt/(abs(vq*dt)+r0)
    vdivpd          %zmm6, .LDT, %zmm1{%k5}      # X_MAX = dt/q

    vxorpd          %zmm0, %zmm0, %zmm0
    vcmppd          $0x11, %zmm0, .LDT, %k7      # dt < 0
    kandw           %k5, %k7, %k7
    vmovapd         %zmm5, %zmm2
    vmovapd         %zmm1, %zmm5{%k7}
    vmovapd         %zmm2, %zmm1{%k7}

.FallbackBoundsDone\@:

    vaddpd          %zmm5, %zmm1, .LXX{%k4}       # X_MIN + X_MAX
    vmulpd          .LHALF, .LXX, .LXX{%k4}           # X = (X_MIN + X_MAX)/2
.FallbackBisectionLoop\@:
.if DEBUG_AVX512 == 1 
    vmovdqa64       P512_COUNTER(%rdi), %zmm4
    vpaddq          .ONE_QUAD(%rip){1to8}, %zmm4, %zmm4{%k4}
    vmovdqa64       %zmm4, P512_COUNTER(%rdi)
.endif
    mm_stiefel_Gs13_avx512
    vmulpd          .LR, .LXX, %zmm2                # r0*X
    vfmadd231pd     .LGS2, .LETA, %zmm2
    vfmadd231pd     .LGS3, .LZETA, %zmm2            # r0*X + eta0*Gs2 + zeta0*Gs3

    vcmppd          $30, .LDT, %zmm2, %k2         # $30 = Greater than, ordered, quiet
    knotb           %k2, %k3
    vmovapd         .LXX, %zmm1{%k2}
    vmovapd         .LXX, %zmm5{%k3}
    vaddpd          %zmm5, %zmm1, .LXX{%k4}       # X_MIN + X_MAX
    vmulpd          .LHALF, .LXX, .LXX{%k4}           # X
    
    incq %r9 
    cmpq $52, %r9                               # Elliptic lanes use 52 bisection iterations.
    jl .FallbackBisectionLoop\@

    # for hyperbolic lanes, continue until stopping criterion.
    vsubpd          %zmm5, %zmm1, %zmm6
    vpandq          .LSIGN_ABS_MASK, %zmm6, %zmm6
    vaddpd          %zmm5, %zmm1, %zmm7
    vpandq          .LSIGN_ABS_MASK, %zmm7, %zmm7
    vmulpd          .BISECTION_EPS(%rip){1to8}, %zmm7, %zmm7
    vcmppd          $0x11, %zmm6, %zmm7, %k4    # tolerance < abs(X_MAX-X_MIN)
    kandw           %k5, %k4, %k4
    kortestw        %k4, %k4
    jz              .NewtonLoopDone\@

    cmpq            $2200, %r9
    jl              .FallbackBisectionLoop\@

.NewtonLoopDone\@:
    mm_stiefel_Gs13_comp
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    # Calculate r_new = .LR + .LGS1*.LETA + .LGS2*.LZETA, then 1/r.
    vfmadd231pd     .LGS1, .LETA, .LR
    vfmadd231pd     .LGS2, .LZETA, .LR
    vdivpd          .LR, .LONE, %zmm4          # 1/r
    
    # Calculate f and g functions
    vmulpd          .LGS2, .LM, %zmm5
    vmulpd          .LIR, %zmm5, %zmm3        # negative f
    vmulpd          %zmm5, %zmm4, %zmm2     # negative gd
    vmovapd         .LDT, %zmm1
    vfnmadd231pd    .LGS3, .LM, %zmm1           # g 
    vmulpd          .LGS1, .LM, %zmm0
    vmulpd          .LIR, %zmm0, %zmm0
    vmulpd          %zmm4, %zmm0, %zmm0     # negative fd

    vmovapd         %zmm3, %zmm4
    vmovapd         %zmm3, %zmm5
    // Calculate new x y z
    vfnmadd132pd    .LX, .LX, %zmm3
    vfnmadd132pd    .LY, .LY, %zmm4
    vfnmadd132pd    .LZ, .LZ, %zmm5
    vfmadd231pd     .LVX, %zmm1, %zmm3{%k1}{z}
    vfmadd231pd     .LVY, %zmm1, %zmm4{%k1}{z}
    vfmadd231pd     .LVZ, %zmm1, %zmm5{%k1}{z}
    // Calculate new vx vy vz
    vfnmadd132pd    %zmm2, .LVX, .LVX
    vfnmadd132pd    %zmm2, .LVY, .LVY
    vfnmadd132pd    %zmm2, .LVZ, .LVZ
    vfnmadd231pd    %zmm0, .LX, .LVX{%k1}{z}
    vfnmadd231pd    %zmm0, .LY, .LVY{%k1}{z}
    vfnmadd231pd    %zmm0, .LZ, .LVZ{%k1}{z}
    vmovapd    %zmm3, .LX 
    vmovapd    %zmm4, .LY
    vmovapd    %zmm5, .LZ
.endm

#####################################
# Macros for interaction step
#####################################
.macro gravity_prefactor multiplier encounterflag
    # Input:  zmm0=dx, zmm1=dy, zmm2=dz
    # Output: zmm6 = multiplier / r^3
    vmulpd      %zmm0, %zmm0, %zmm6
    vfmadd231pd %zmm1, %zmm1, %zmm6
    vfmadd231pd %zmm2, %zmm2, %zmm6     # zmm6 = r^2 = a

    vrsqrt14pd  %zmm6, %zmm7            # zmm7 = y_0 ~ 1/sqrt(a), ~14 bits
    vmulpd      .LHALF, %zmm6, %zmm8      # zmm8 = 0.5 * a  (constant across iters)

    # Newton iter 1
    vmulpd      %zmm7, %zmm7, %zmm5     # zmm5 = y_0^2
    vfnmadd213pd .ONE_AND_A_HALF(%rip){1to8}, %zmm8, %zmm5
                                        # zmm5 = -zmm5*zmm8 + 1.5 = 1.5 - 0.5*a*y_0^2
    vmulpd      %zmm7, %zmm5, %zmm7     # y_1 = y_0 * (1.5 - 0.5*a*y_0^2)

    # Newton iter 2
    vmulpd      %zmm7, %zmm7, %zmm5     # zmm5 = y_1^2
    vfnmadd213pd .ONE_AND_A_HALF(%rip){1to8}, %zmm8, %zmm5
                                        # zmm5 = 1.5 - 0.5*a*y_1^2
    vmulpd      %zmm7, %zmm5, %zmm7     # y_2 ~ 1/sqrt(a) to ~56 bits
    
    # Check for close encounters
    .if \encounterflag == 1
    vcmppd          $0x1E, P512_EXIT_MIN_DISTANCE_R(%rdi), %zmm7, %k4{%k1} # $1E = greater than, ordered (nans fail), quiet, k4=1 if distance<exit_min_distance
    kmovw           %k4, %r9d
    negq            %r9                     # If r9 =0, Carry Flag will be set
    movq            $3, %r9
    cmovnzq         %r9, %rax               # Set return value to REB_STATUS_ENCOUNTER
    .endif

    # y2*y2 -> *y2 -> *mult implemented as (y2*y2) || (mult*y2) -> mul.
    vmulpd      %zmm7, %zmm7, %zmm5         # zmm5 = y_2^2
    .ifc \multiplier,.LONE
    vmulpd      %zmm5, %zmm7, %zmm6         # zmm6 = y_2^3 ~ mult / r^3
    .else
    vmulpd      \multiplier, %zmm7, %zmm6   # zmm6 = mult * y_2  (parallel)
    vmulpd      %zmm5, %zmm6, %zmm6         # zmm6 = mult * y_2^3 ~ mult / r^3
    .endif
.endm


.macro mat8_mul3 in0, in1, in2, out0, out1, out2
    # 8x8 matrix multiplied with 3 different 8 vectors
    # in: r9  = vector to 64 matrix elements
    # Does not alter inputs
    # uses: zmm3-zmm7
    # The idea is to use embedded broadcast loads
    # Note: matrix needs to be transposed.
    vmovapd \in0,   0(%rsp)
    vmovapd \in1,  64(%rsp)
    vmovapd \in2, 128(%rsp)

    # Keeping six independent FMA chains going
    vmovapd        (%r9 ), %zmm4
    vmovapd      64(%r9 ), %zmm3
    vmulpd        0(%rsp){1to8}, %zmm4, \out0
    vmulpd       64(%rsp){1to8}, %zmm4, \out1
    vmulpd      128(%rsp){1to8}, %zmm4, \out2

    vmulpd        8(%rsp){1to8}, %zmm3, %zmm5
    vmulpd       72(%rsp){1to8}, %zmm3, %zmm6
    vmulpd      136(%rsp){1to8}, %zmm3, %zmm7

    vmovapd     128(%r9 ), %zmm4
    vmovapd     192(%r9 ), %zmm3
    vfmadd231pd  16(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd  80(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 144(%rsp){1to8}, %zmm4, \out2
    
    vfmadd231pd  24(%rsp){1to8}, %zmm3, %zmm5
    vfmadd231pd  88(%rsp){1to8}, %zmm3, %zmm6
    vfmadd231pd 152(%rsp){1to8}, %zmm3, %zmm7
    
    vmovapd     256(%r9 ), %zmm4
    vmovapd     320(%r9 ), %zmm3
    vfmadd231pd  32(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd  96(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 160(%rsp){1to8}, %zmm4, \out2
    
    vfmadd231pd  40(%rsp){1to8}, %zmm3, %zmm5
    vfmadd231pd 104(%rsp){1to8}, %zmm3, %zmm6
    vfmadd231pd 168(%rsp){1to8}, %zmm3, %zmm7
    
    vmovapd     384(%r9 ), %zmm4
    vmovapd     448(%r9 ), %zmm3
    vfmadd231pd  48(%rsp){1to8}, %zmm4, \out0
    vfmadd231pd 112(%rsp){1to8}, %zmm4, \out1
    vfmadd231pd 176(%rsp){1to8}, %zmm4, \out2
    
    vfmadd231pd  56(%rsp){1to8}, %zmm3, %zmm5
    vfmadd231pd 120(%rsp){1to8}, %zmm3, %zmm6
    vfmadd231pd 184(%rsp){1to8}, %zmm3, %zmm7

    # Using two accumulators, adding at end
    vaddpd  %zmm5, \out0, \out0
    vaddpd  %zmm6, \out1, \out1
    vaddpd  %zmm7, \out2, \out2
   .endm

###############################################################################
# Interaction Step
###############################################################################
.macro interaction_step grflag nsys encounterflag escapeflag
    # TODO: Floating point error accumulation might be less if Jacobi and GR are added after P-P perturbations
    # Add Jacobi term in Jacobi coordinates
    vmulpd      .LX, .LX, %zmm4     
    vfmadd231pd .LY, .LY, %zmm4      
    vfmadd231pd .LZ, .LZ, %zmm4             # r^2
    vsqrtpd     %zmm4, %zmm5            # r 
    vmulpd      %zmm4, %zmm5, %zmm4     # r^3
  
    vdivpd      %zmm4, .LM_DT, %zmm6      # M*dt/r^3 (where M=(m0, m0+m1, m0+m1+m2,...)
    
    vfmadd231pd     .LX, %zmm6, .LVX{%k1}{z} 
    vfmadd231pd     .LY, %zmm6, .LVY{%k1}{z} 
    vfmadd231pd     .LZ, %zmm6, .LVZ{%k1}{z} 
    
    leaq P512_MAT8_JACOBI_TO_HELIOCENTRIC(%rdi), %r9   # mat8_inertial_to_jacobi
    mat8_mul3 .LX, .LY, .LZ, .LHX, .LHY, .LHZ
    
    # Calculating r, r^2, r^3 for Jacobi term and GR
    vmulpd      .LHX, .LHX, %zmm6
    vfmadd231pd .LHY, .LHY, %zmm6
    vfmadd231pd .LHZ, .LHZ, %zmm6               # r^2
    vsqrtpd     %zmm6, %zmm7                # r
    
    # Check for escapes 
    .if \escapeflag == 1
    vcmppd          $0x1E, P512_EXIT_MAX_DISTANCE(%rdi), %zmm7, %k4{%k1} # $1E = greater than, ordered (nans fail), quiet, k4=1 if distance>exit_max_distance
    kmovw           %k4, %r9d
    negq            %r9                     # If r9 =0, Carry Flag will be set
    movq            $4, %r9
    cmovnzq         %r9, %rax               # Set return value to REB_STATUS_EJECTION
    .endif
    .if \encounterflag == 1
    vcmppd          $0x11, P512_EXIT_MIN_DISTANCE(%rdi), %zmm7, %k4{%k1} # $11 = less than, ordered (nans fail), quiet, k4=1 if distance<exit_min_distance
    kmovw           %k4, %r9d
    negq            %r9                     # If r9 =0, Carry Flag will be set
    movq            $3, %r9
    cmovnzq         %r9, %rax               # Set return value to REB_STATUS_ENCOUNTER
    .endif

    # Jacobi term
    vmulpd    %zmm6, %zmm7, %zmm7           # r^3    
    vdivpd    %zmm7, .LMM0_DT, %zmm8{%k1}{z}  # -m0*dt/r^3 (jacobi term)
        
    vmulpd    %zmm8, .LHX, .LHVX                # delta v_x due to Jacobi term, -x_j*m0*dt/r^3
    vmulpd    %zmm8, .LHY, .LHVY
    vmulpd    %zmm8, .LHZ, .LHVZ

    # GR term
    .if \grflag == 1
        vmulpd    P512_GR_PREFAC(%rdi), .LDT, %zmm3

        vmulpd    %zmm6, %zmm6, %zmm5           # r^4
        vdivpd    %zmm5, %zmm3, %zmm7{%k1}{z}   # -dt*6*m0*m0/(c*c) /r^4

        vfmadd231pd  %zmm7, .LHX, .LHVX{%k1}{z}     # -x_j*dt*6*m0*m0/(c*c) /r^4
        vfmadd231pd  %zmm7, .LHY, .LHVY{%k1}{z}
        vfmadd231pd  %zmm7, .LHZ, .LHVZ{%k1}{z}
    .endif

    #################################################################
    #// 0123 4567
    #// 3201 7645

    vmulpd  P512_m(%rdi), .LDT, %zmm3         # dt*m

  .if \nsys < 4                             # skip for 2-planet systems
    vpermpd $0x4B, .LHX, %zmm0                # 01234567 -> 32017645
    vpermpd $0x4B, .LHY, %zmm1
    vpermpd $0x4B, .LHZ, %zmm2
    vpermpd $0x4B, %zmm3, %zmm4

    vsubpd  %zmm0, .LHX, %zmm0                # d_x
    vsubpd  %zmm1, .LHY, %zmm1
    vsubpd  %zmm2, .LHZ, %zmm2

    gravity_prefactor .LONE \encounterflag    # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5         # dt*m/r^3

    vfnmadd231pd %zmm5, %zmm0,  .LHVX
    vfnmadd231pd %zmm5, %zmm1,  .LHVY
    vfnmadd231pd %zmm5, %zmm2,  .LHVZ

    vmulpd      %zmm6, %zmm3, %zmm5         # dt*m/r^3
    vpermpd $0x1E, %zmm0, %zmm0             # 32017645 -> 01234567
    vpermpd $0x1E, %zmm1, %zmm1
    vpermpd $0x1E, %zmm2, %zmm2
    vpermpd $0x1E, %zmm5, %zmm5

    #// 0123 4567
    #// 2310 6754

    vfmadd231pd %zmm5, %zmm0,  .LHVX
    vfmadd231pd %zmm5, %zmm1,  .LHVY
    vfmadd231pd %zmm5, %zmm2,  .LHVZ
  .endif

    #################################################################
    #// 0123 4567
    #// 1032 5476
    
    vshufpd $0x55, .LHX, .LHX, %zmm0                # 01234567 -> 10325476
    vshufpd $0x55, .LHY, .LHY, %zmm1                # Using vshufpd (1 cycle) rather than vpermpd (3 cycles) 
    vshufpd $0x55, .LHZ, .LHZ, %zmm2
    vshufpd $0x55, %zmm3, %zmm3, %zmm4 

    vsubpd  %zmm0, .LHX, %zmm0                    # d_x
    vsubpd  %zmm1, .LHY, %zmm1
    vsubpd  %zmm2, .LHZ, %zmm2
    
    gravity_prefactor %zmm4 \encounterflag      # zmm6 is 1/r^3
    
    vfnmadd231pd %zmm6, %zmm0,  .LHVX
    vfnmadd231pd %zmm6, %zmm1,  .LHVY
    vfnmadd231pd %zmm6, %zmm2,  .LHVZ

  .if \nsys == 1                                # only for a single 8-planet system
    #################################################################
    #// 0123 4567
    #// 4567 1230

    vmovdqa64 b3idx(%rip), %zmm7

    vpermpd .LHX, %zmm7, %zmm0                    # 01234567 -> 45671230 
    vpermpd .LHY, %zmm7, %zmm1
    vpermpd .LHZ, %zmm7, %zmm2
    vpermpd %zmm3, %zmm7, %zmm4 

    vsubpd  %zmm0, .LHX, %zmm0                    # d_x
    vsubpd  %zmm1, .LHY, %zmm1
    vsubpd  %zmm2, .LHZ, %zmm2
    
    gravity_prefactor .LONE \encounterflag        # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5             # m/r^3
  
    vfnmadd231pd %zmm5, %zmm0,  .LHVX
    vfnmadd231pd %zmm5, %zmm1,  .LHVY
    vfnmadd231pd %zmm5, %zmm2,  .LHVZ

    vmulpd      %zmm6, %zmm3, %zmm5             # m/r^3
    
    #// 4567 1230
    #// 0123 4567
    vmulpd %zmm5, %zmm0,  .LHVXC
    vmulpd %zmm5, %zmm1,  .LHVYC
    vmulpd %zmm5, %zmm2,  .LHVZC

    #################################################################
    #// 0123 4567
    #// 5674 2301
    
    vmovdqa64 b4idx(%rip), %zmm7

    vpermpd .LHX, %zmm7, %zmm0                    # 01234567 -> 56742301  
    vpermpd .LHY, %zmm7, %zmm1                    # TODO: Make this an in-line shuffle by reusing block3 data
    vpermpd .LHZ, %zmm7, %zmm2
    vpermpd %zmm3, %zmm7, %zmm4 

    vsubpd  %zmm0, .LHX, %zmm0                    # d_x
    vsubpd  %zmm1, .LHY, %zmm1
    vsubpd  %zmm2, .LHZ, %zmm2
    
    gravity_prefactor .LONE \encounterflag        # zmm6 is 1/r^3
    vmulpd      %zmm6, %zmm4, %zmm5             # m/r^3
  
    vfnmadd231pd %zmm5, %zmm0,  .LHVX
    vfnmadd231pd %zmm5, %zmm1,  .LHVY
    vfnmadd231pd %zmm5, %zmm2,  .LHVZ

    vmulpd      %zmm6, %zmm3, %zmm5             # m/r^3
    vpermpd $0x93, %zmm0, %zmm0                 # 5674 2301 -> 4567 1230
    vpermpd $0x93, %zmm1, %zmm1
    vpermpd $0x93, %zmm2, %zmm2
    vpermpd $0x93, %zmm5, %zmm5
    
    #// 4567 1230
    #// 3012 7456
    
    vfmadd231pd %zmm5, %zmm0,  .LHVXC
    vfmadd231pd %zmm5, %zmm1,  .LHVYC
    vfmadd231pd %zmm5, %zmm2,  .LHVZC
    
    #################################################################
    ## Final 256 bit lane crossing and add
    vmovdqa64 b34mergeidx(%rip), %zmm7

    vpermpd .LHVXC, %zmm7, %zmm0
    vpermpd .LHVYC, %zmm7, %zmm1
    vpermpd .LHVZC, %zmm7, %zmm2

    vaddpd %zmm0, .LHVX, %zmm0{%k1}{z}
    vaddpd %zmm1, .LHVY, %zmm1{%k1}{z}
    vaddpd %zmm2, .LHVZ, %zmm2{%k1}{z}
  .else                                        # for nsys 2 or 4
    vmovapd .LHVX, %zmm0{%k1}{z}
    vmovapd .LHVY, %zmm1{%k1}{z}
    vmovapd .LHVZ, %zmm2{%k1}{z}
  .endif

    # Convert accelerations (delta v) from heliocentric to Jacobi.
    leaq P512_MAT8_INERTIAL_TO_JACOBI(%rdi), %r9   # mat8_inertial_to_jacobi
   
    mat8_mul3 %zmm0, %zmm1, %zmm2, %zmm0, %zmm1, %zmm2

    # Update velocities
    # This could be combined with mat8_mul3.
    # However, that would increase floating point errors because sum(DVX) << .LVX
    vaddpd    .LVX, %zmm0, .LVX        
    vaddpd    .LVY, %zmm1, .LVY
    vaddpd    .LVZ, %zmm2, .LVZ
.endm 



###############################################################################
# Global functions
###############################################################################

.macro corrector_step grflag nsys
    reb_whfast512_init_registers                   # does not overwrite xmm0
    alloc_stack64   256                         # space for matricies (192) and direction (8), rounded up to nearest 64bytes
    movsd           %xmm0, 192(%rsp)            # store direction (1 or -1)
    leaq            .CORRECTOR17_AB(%rip), %r8
    leaq            504(%r8), %rdx              # end of array 63*8
    jmp             .L_CorrectorLoopK\@         # start with Kepler step

.L_CorrectorLoopI\@:
    vbroadcastsd    (%r8), %zmm0
    vmulpd          192(%rsp){1to8}, %zmm0, %zmm0
    # Interaction step uses .LDT, .LM_DT, .LMM0_DT
    vmulpd          P512_DT(%rdi), %zmm0, .LDT
    vmulpd          .LDT, .LM, .LM_DT
    vmovapd         P512_M0(%rdi), .LMM0_DT
    vmulpd          .LDT, .LMM0_DT, .LMM0_DT
    vxorpd          .SIGN_FLIP_MASK(%rip){1to8}, .LMM0_DT, .LMM0_DT
    interaction_step \grflag \nsys 0 0
    addq            $8, %r8

.L_CorrectorLoopK\@:
    vbroadcastsd    (%r8), %zmm0
    vmulpd          P512_DT(%rdi), %zmm0, .LDT
    vmulpd          .LDT, .LHALF, .LDT                # Reduce timestep for better convergence
    vmulpd          .LDT, .LHALF, .LDT
    movq            $4, %r10                     # Counter number of Kepler steps
.L_CorrectorLoopInnerKepler\@:
    kepler_step
    decq            %r10
    jnz             .L_CorrectorLoopInnerKepler\@
    addq            $8, %r8
    
    cmpq            %rdx, %r8
    jne             .L_CorrectorLoopI\@

    free_stack64
    reb_whfast512_store_results
    ret
.endm

# Also make the Kepler step available (for synchronization)
.globl reb_whfast512_kepler_step
reb_whfast512_kepler_step:
    reb_whfast512_init_registers
    kepler_step
    reb_whfast512_store_results
    ret

# Helper functions to load and move data
.globl reb_whfast512_set1_pd
reb_whfast512_set1_pd:
    # Input:
    #           rdi = pointer to 512bit memory
    #           xmm0 = double
    vpbroadcastq    %xmm0, %zmm0
    vmovapd         %zmm0, (%rdi)
    ret 

.globl reb_whfast512_movu_pd
reb_whfast512_movu_pd:
    # Input:
    #           rdi = destination
    #           rsi = source
    vmovupd     (%rsi), %zmm0
    vmovupd     %zmm0, (%rdi)
    ret


# Generate actual integration functions using macros.
# We do this to avoid branching during the inner loops.
# There is a GNU as bug which limits the number of nested irp loops, so we need to refactor this into macros.
# The basic idea is that we programatically create functions with all possible combinations of gr, nsys, encounter, and escape.
.macro full_steps grflag nsys encounterflag escapeflag
    # Input:
    #           rdi = p512
    #           rsi = pointer to number of steps
    #           rdx = skip_first_kepler_step
    #           rcx = pointer to reb_sigint (to check for interrupt)
    # Output: 
    #           rsi = pointer contains number of steps actually done
    
    .set ExceptionsCanOccur, 0
    .if \encounterflag == 1
    .set ExceptionsCanOccur, 1
    .endif
    .if \escapeflag == 1
    .set ExceptionsCanOccur, 1
    .endif
    
    movq        (%rsi), %r10    # counting down number of steps
    .if ExceptionsCanOccur == 1
    movq        $0, %rax        # flag to set return value
    .endif

    # Load constants
    reb_whfast512_init_registers
    # Allocate space on stack for matrix multiplications
    alloc_stack64 192
    # Ignore first Kepler step (if half timestep done manually)
    cmpq    $1, %rdx
    je      .LSkipFirstKeplerStep\@

    # Main loop
.LMainLoop\@:    
    kepler_step
.LSkipFirstKeplerStep\@:
    interaction_step \grflag \nsys \encounterflag \escapeflag
    cmpq    $0, (%rcx)
    jnz     .LInterruptOccured\@
    .if ExceptionsCanOccur == 1
    testq   %rax, %rax
    jnz     .LExceptionOccured\@
    .endif
    subq    $1, %r10
    jg      .LMainLoop\@
    movq    $-1, %rax       # Successfull completion of timestep status = REB_STATUS_RUNNING
    jmp     .LSuccess\@

.LInterruptOccured\@:
    movq    $6, %rax        # status = REB_STATUS_SIGINT
.LExceptionOccured\@:       # close mindistance or ejection
    subq    $1, %r10        # number of steps remaining (could be zero, but can't be negative)
    subq    %r10, (%rsi)    # steps done
.LSuccess\@:
    # Store final data in P512 structure
    reb_whfast512_store_results
    free_stack64
    ret
.endm

.macro reb_whfast512_full_steps_macro3 gr, nsys, encounterflag
.irp escapeflag,0,1
.globl reb_whfast512_full_steps_gr\gr\()_n\nsys\()_encounter\encounterflag\()_escape\escapeflag
reb_whfast512_full_steps_gr\gr\()_n\nsys\()_encounter\encounterflag\()_escape\escapeflag: full_steps \gr \nsys \encounterflag \escapeflag
.endr
.endm

.macro reb_whfast512_full_steps_macro2 gr, nsys
.irp encounterflag,0,1
reb_whfast512_full_steps_macro3 \gr \nsys \encounterflag
.endr
.endm

.macro reb_whfast512_full_steps_macro1 gr
.irp nsys,1,2,4
reb_whfast512_full_steps_macro2 \gr, \nsys
.global reb_whfast512_corrector_step_gr\gr\()_n\nsys
reb_whfast512_corrector_step_gr\gr\()_n\nsys: corrector_step \gr \nsys
.endr
.endm

.irp gr,0,1
reb_whfast512_full_steps_macro1 \gr
.endr


.section    .rodata


.if DEBUG_AVX512 == 1
.align 64
.ONE_QUAD:
    .quad 1
.endif

# Shuffle Indicies
# Each is eight 64-bit integers
.align 64
b3idx:
    .quad 4,5,6,7,1,2,3,0   
b4idx:
    .quad 5,6,7,4,2,3,0,1
b34mergeidx:
    .quad 7,4,5,6,0,1,2,3
# These two MASK are the inverse of each other and could be combined.
.align 64
.SIGN_ABS_MASK:
    .quad 0x7FFFFFFFFFFFFFFF
.align 64
.SIGN_FLIP_MASK:
    .quad 0x8000000000000000
.align 64
.EPS:
    .double 1e-11
.align 64
.BISECTION_EPS:
    .double 1e-15
.align 64
.TWOPI:
    .quad 0x401921fb54442d18
.align 64
..LHALF:
    .quad 0x3fe0000000000000
.align 64
# Exponent decrement: subtracting this from a normal double divides it by 2
.HALF_EXP_DECR:
    .quad 0x0010000000000000
.align 64
# 1.5 — Newton-Raphson constant for rsqrt iteration: y' = y*(1.5 - 0.5*a*y*y)
.ONE_AND_A_HALF:
    .quad 0x3ff8000000000000
.align 64
# Inverse factorial table
.IF0:
.DOUBLE_ONE:
    .quad 0x3ff0000000000000  # = 1.0000e+00 = 1/1
    .quad 0x3ff0000000000000  # = 1.0000e+00 = 1/1
    .quad 0x3fe0000000000000  # = 5.0000e-01 = 1/2
    .quad 0x3fc5555555555555  # = 1.6667e-01 = 1/6
    .quad 0x3fa5555555555555  # = 4.1667e-02 = 1/24
    .quad 0x3f81111111111111  # = 8.3333e-03 = 1/120
    .quad 0x3f56c16c16c16c17  # = 1.3889e-03 = 1/720
    .quad 0x3f2a01a01a01a01a  # = 1.9841e-04 = 1/5040
    .quad 0x3efa01a01a01a01a  # = 2.4802e-05 = 1/40320
    .quad 0x3ec71de3a556c734  # = 2.7557e-06 = 1/362880
    .quad 0x3e927e4fb7789f5c  # = 2.7557e-07 = 1/3628800
    .quad 0x3e5ae64567f544e4  # = 2.5052e-08 = 1/39916800
    .quad 0x3e21eed8eff8d898  # = 2.0877e-09 = 1/479001600
    .quad 0x3de6124613a86d09  # = 1.6059e-10 = 1/6227020800
    .quad 0x3da93974a8c07c9d  # = 1.1471e-11 = 1/87178291200
    .quad 0x3d6ae7f3e733b81f  # = 7.6472e-13 = 1/1307674368000
    .quad 0x3d2ae7f3e733b81f  # = 4.7795e-14 = 1/20922789888000
    .quad 0x3ce952c77030ad4a  # = 2.8115e-15 = 1/355687428096000
    .quad 0x3ca6827863b97d97  # = 1.5619e-16 = 1/6402373705728000
    .quad 0x3c62f49b46814157  # = 8.2206e-18 = 1/121645100408832000
    .quad 0x3c1e542ba4020225  # = 4.1103e-19 = 1/2432902008176640000
    .quad 0x3bd71b8ef6dcf572  # = 1.9573e-20 = 1/51090942171709440000
    .quad 0x3b90ce396db7f853  # = 8.8968e-22 = 1/1124000727777607680000
    .quad 0x3b4761b41316381a  # = 3.8682e-23 = 1/25852016738884976640000
    .quad 0x3aff2cf01972f578  # = 1.6117e-24 = 1/620448401733239439360000
    .quad 0x3ab3f3ccdd165fa9  # = 6.4470e-26 = 1/15511210043330985984000000
    .quad 0x3a688e85fc6a4e5a  # = 2.4796e-27 = 1/403291461126605635584000000
    .quad 0x3a1d1ab1c2dccea3  # = 9.1837e-29 = 1/10888869450418352160768000000
    .quad 0x39d0a18a2635085d  # = 3.2799e-30 = 1/304888344611713860501504000000
    .quad 0x398259f98b4358ad  # = 1.1310e-31 = 1/8841761993739701954543616000000
    .quad 0x3933932c5047d60e  # = 3.7700e-33 = 1/265252859812191058636308480000000
    .quad 0x38e434d2e783f5bc  # = 1.2161e-34 = 1/8222838654177922817725562880000000
    .quad 0x389434d2e783f5bc  # = 3.8004e-36 = 1/263130836933693530167218012160000000
    .quad 0x3843981254dd0d52  # = 1.1516e-37 = 1/8683317618811886495518194401280000000
    .quad 0x37f2710231c0fd7a  # = 3.3872e-39 = 1/295232799039604140847618609643520000000

.align 64
# rounding error of each .IF0 entry (true 1/k! - rounded double)
.IF0_err:
    .quad 0x0000000000000000  # = +0.0000e+00
    .quad 0x0000000000000000  # = +0.0000e+00
    .quad 0x0000000000000000  # = +0.0000e+00
    .quad 0xbc65555555555555  # = -9.2519e-18
    .quad 0xbc45555555555555  # = -2.3130e-18
    .quad 0xbc01111111111111  # = -1.1565e-19
    .quad 0x3bef49f49f49f49f  # = +5.3005e-20
    .quad 0xbb6a01a01a01a01a  # = -1.7210e-22
    .quad 0xbb3a01a01a01a01a  # = -2.1512e-23
    .quad 0x3b6c154f8ddc6c00  # = +1.8584e-22
    .quad 0xbb3cbbc05b4fa99a  # = -2.3768e-23
    .quad 0x3afc062e06d1f209  # = +1.4488e-24
    .quad 0x3ac2aec959e14c06  # = +1.2073e-25
    .quad 0xba8f28e0cc748ebe  # = -1.2585e-26
    .quad 0xba305d6f8a2efd1f  # = -2.0656e-28
    .quad 0xb9e1d8656b0ee8cb  # = -7.0387e-30
    .quad 0xb9a1d8656b0ee8cb  # = -4.3992e-31
    .quad 0xb98ac981465ddc6c  # = -1.6509e-31
    .quad 0xb94eec01221a8b0b  # = -1.1911e-32
    .quad 0xb8f2650f61dbdcb4  # = -2.2142e-34
    .quad 0xb87ea72b4afe3c2f  # = -1.4413e-36
    .quad 0x387d043ae40c4647  # = +1.3644e-36
    .quad 0x383aebcdbd20331c  # = +7.9114e-38
    .quad 0x37d3423c7d91404f  # = +8.8432e-40
    .quad 0x3789ada5fcc1ab14  # = +3.6847e-41
    .quad 0x37458ddadf344487  # = +1.9330e-42
    .quad 0x37071c37ebd16540  # = +1.2954e-43
    .quad 0xb6a054d0c78aea14  # = -1.4303e-45
    .quad 0xb66b9e2e28e1aa54  # = -1.5118e-46
    .quad 0xb62eaf8c39dd9bc5  # = -1.0498e-47
    .quad 0xb5d832b7b530a627  # = -2.5870e-49
    .quad 0xb580b87b91be9aff  # = -5.5863e-51
    .quad 0xb530b87b91be9aff  # = -1.7457e-52
    .quad 0x34e2b1f4c8015a2f  # = +6.0996e-54
    .quad 0xb473f8a2b4af9d6b  # = -5.0906e-56

.align 8
.CORRECTOR17_AB:        # 63 coefficients for 17th order corrector
    .quad 0xC00AC5EB3F7AB2F8    # Kepler
    .quad 0xBED22E64AF0557FF    # Interaction
    .quad 0x401AC5EB3F7AB2F8    # Kepler
    .quad 0x3ED22E64AF0557FF    # Interaction
    .quad 0xC019198C8B8307C8    # Kepler
    .quad 0x3F14098E956E85C7    # Interaction
    .quad 0x40176D2DD78B5C99    # Kepler
    .quad 0xBF14098E956E85C7    # Interaction
    .quad 0xC015C0CF2393B16A    # Kepler
    .quad 0xBF44D7273C9751C8    # Interaction
    .quad 0x401414706F9C063A    # Kepler
    .quad 0x3F44D7273C9751C8    # Interaction
    .quad 0xC0126811BBA45B0A    # Kepler
    .quad 0x3F6B2467AFD31964    # Interaction
    .quad 0x4010BBB307ACAFDB    # Kepler
    .quad 0xBF6B2467AFD31964    # Interaction
    .quad 0xC00E1EA8A76A0957    # Kepler
    .quad 0xBF88B9144F7F2B3F    # Interaction
    .quad 0x400AC5EB3F7AB2F8    # Kepler
    .quad 0x3F88B9144F7F2B3F    # Interaction
    .quad 0xC0076D2DD78B5C99    # Kepler
    .quad 0x3FA099A47793A306    # Interaction
    .quad 0x400414706F9C063A    # Kepler
    .quad 0xBFA099A47793A306    # Interaction
    .quad 0xC000BBB307ACAFDB    # Kepler
    .quad 0xBFB0B07AC0FE3DF1    # Interaction
    .quad 0x3FFAC5EB3F7AB2F8    # Kepler
    .quad 0x3FB0B07AC0FE3DF1    # Interaction
    .quad 0xBFF414706F9C063A    # Kepler
    .quad 0x3FB7D2865A643682    # Interaction
    .quad 0x3FEAC5EB3F7AB2F8    # Kepler
    .quad 0xBFC7D2865A643682    # Interaction (combined next three)
#    .quad 0xBFB7D2865A643682    # Interaction
#    .quad 0x0000000000000000    # Kepler
#    .quad 0xBFB7D2865A643682    # Interaction
    .quad 0xBFEAC5EB3F7AB2F8    # Kepler
    .quad 0x3FB7D2865A643682    # Interaction
    .quad 0x3FF414706F9C063A    # Kepler
    .quad 0x3FB0B07AC0FE3DF1    # Interaction
    .quad 0xBFFAC5EB3F7AB2F8    # Kepler
    .quad 0xBFB0B07AC0FE3DF1    # Interaction
    .quad 0x4000BBB307ACAFDB    # Kepler
    .quad 0xBFA099A47793A306    # Interaction
    .quad 0xC00414706F9C063A    # Kepler
    .quad 0x3FA099A47793A306    # Interaction
    .quad 0x40076D2DD78B5C99    # Kepler
    .quad 0x3F88B9144F7F2B3F    # Interaction
    .quad 0xC00AC5EB3F7AB2F8    # Kepler
    .quad 0xBF88B9144F7F2B3F    # Interaction
    .quad 0x400E1EA8A76A0957    # Kepler
    .quad 0xBF6B2467AFD31964    # Interaction
    .quad 0xC010BBB307ACAFDB    # Kepler
    .quad 0x3F6B2467AFD31964    # Interaction
    .quad 0x40126811BBA45B0A    # Kepler
    .quad 0x3F44D7273C9751C8    # Interaction
    .quad 0xC01414706F9C063A    # Kepler
    .quad 0xBF44D7273C9751C8    # Interaction
    .quad 0x4015C0CF2393B16A    # Kepler
    .quad 0xBF14098E956E85C7    # Interaction
    .quad 0xC0176D2DD78B5C99    # Kepler
    .quad 0x3F14098E956E85C7    # Interaction
    .quad 0x4019198C8B8307C8    # Kepler
    .quad 0x3ED22E64AF0557FF    # Interaction
    .quad 0xC01AC5EB3F7AB2F8    # Kepler
    .quad 0xBED22E64AF0557FF    # Interaction
    .quad 0x400AC5EB3F7AB2F8    # Kepler

# Now passed via --noexecstack to as
#.section .note.GNU-stack,"",@progbits
