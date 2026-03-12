	.file	"tmp.c"
	.text
	.section	.rodata
	.align 32
	.type	invfactorial, @object
	.size	invfactorial, 280
invfactorial:
	.long	0
	.long	1072693248
	.long	0
	.long	1072693248
	.long	0
	.long	1071644672
	.long	1431655765
	.long	1069897045
	.long	1431655765
	.long	1067799893
	.long	286331153
	.long	1065423121
	.long	381774871
	.long	1062650220
	.long	436314138
	.long	1059717536
	.long	436314138
	.long	1056571808
	.long	-1521039564
	.long	1053236707
	.long	-1216831652
	.long	1049787983
	.long	1744127204
	.long	1046144581
	.long	-268904296
	.long	1042411224
	.long	329805065
	.long	1038488134
	.long	-1463780195
	.long	1034500468
	.long	-416040929
	.long	1030416371
	.long	-416040929
	.long	1026222067
	.long	1882238282
	.long	1021924039
	.long	1673100695
	.long	1017545336
	.long	1182875991
	.long	1013118107
	.long	-1543372251
	.long	1008620587
	.long	-153291406
	.long	1003953038
	.long	1840773203
	.long	999345721
	.long	320223257
	.long	994533812
	.long	426964344
	.long	989801712
	.long	-585736279
	.long	984871884
	.long	-60141991
	.long	979930757
	.long	-1025716573
	.long	974985905
	.long	641009757
	.long	969974154
	.long	-1958520658
	.long	964844025
	.long	1346885134
	.long	959681324
	.long	-410782275
	.long	954479826
	.long	-410782275
	.long	949236946
	.long	1423773010
	.long	943953938
	.long	834731386
	.long	938635522
	.text
	.type	mm_stiefel_Gs13_avx512, @function
mm_stiefel_Gs13_avx512:
.LFB4831:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-64, %rsp
	subq	$1672, %rsp
	movq	%rdi, 64(%rsp)
	movq	%rsi, 56(%rsp)
	movq	%rdx, 48(%rsp)
	vmovapd	%zmm0, -56(%rsp)
	vmovapd	%zmm1, -120(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 1224(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 1160(%rsp)
	vmovapd	1224(%rsp), %zmm0
	vmulpd	1160(%rsp), %zmm0, %zmm0
	vmovapd	%zmm0, 1544(%rsp)
	vmovapd	1544(%rsp), %zmm0
	vmovapd	%zmm0, 1352(%rsp)
	vmovapd	-56(%rsp), %zmm0
	vmovapd	%zmm0, 1288(%rsp)
	vmovapd	1352(%rsp), %zmm0
	vmulpd	1288(%rsp), %zmm0, %zmm0
	vmovapd	%zmm0, 1480(%rsp)
	movl	$19, 1476(%rsp)
	movl	1476(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 1456(%rsp)
	vbroadcastsd	1456(%rsp), %zmm0
	movq	48(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movl	1476(%rsp), %eax
	decl	%eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 1464(%rsp)
	vbroadcastsd	1464(%rsp), %zmm0
	movq	56(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movl	1476(%rsp), %eax
	subl	$2, %eax
	movl	%eax, 1668(%rsp)
	jmp	.L6
.L11:
	movl	1668(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 704(%rsp)
	vbroadcastsd	704(%rsp), %zmm0
	vmovapd	%zmm0, %zmm2
	movq	48(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1480(%rsp), %zmm1
	vmovapd	%zmm1, 840(%rsp)
	vmovapd	%zmm0, 776(%rsp)
	vmovapd	%zmm2, 712(%rsp)
	vmovapd	840(%rsp), %zmm0
	vmovapd	776(%rsp), %zmm2
	vmovapd	712(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	movq	48(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movl	1668(%rsp), %eax
	decl	%eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 960(%rsp)
	vbroadcastsd	960(%rsp), %zmm0
	vmovapd	%zmm0, %zmm2
	movq	56(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1480(%rsp), %zmm1
	vmovapd	%zmm1, 1096(%rsp)
	vmovapd	%zmm0, 1032(%rsp)
	vmovapd	%zmm2, 968(%rsp)
	vmovapd	1096(%rsp), %zmm0
	vmovapd	1032(%rsp), %zmm2
	vmovapd	968(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	movq	56(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	subl	$2, 1668(%rsp)
.L6:
	cmpl	$2, 1668(%rsp)
	jg	.L11
	movq	48(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	%zmm0, 136(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 72(%rsp)
	vmovapd	136(%rsp), %zmm0
	vmulpd	72(%rsp), %zmm0, %zmm0
	movq	48(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	48(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1480(%rsp), %zmm1
	vmovapd	%zmm1, 328(%rsp)
	vmovapd	%zmm0, 264(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 200(%rsp)
	vmovapd	328(%rsp), %zmm0
	vmovapd	264(%rsp), %zmm2
	vmovapd	200(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k3}
	nop
	movq	64(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	48(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	%zmm0, 456(%rsp)
	vmovapd	1544(%rsp), %zmm0
	vmovapd	%zmm0, 392(%rsp)
	vmovapd	456(%rsp), %zmm0
	vmulpd	392(%rsp), %zmm0, %zmm0
	movq	48(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	56(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	%zmm0, 584(%rsp)
	vmovapd	1544(%rsp), %zmm0
	vmovapd	%zmm0, 520(%rsp)
	vmovapd	584(%rsp), %zmm0
	vmulpd	520(%rsp), %zmm0, %zmm0
	movq	56(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	nop
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4831:
	.size	mm_stiefel_Gs13_avx512, .-mm_stiefel_Gs13_avx512
	.type	mm_stiefel_Gs03_avx512, @function
mm_stiefel_Gs03_avx512:
.LFB4832:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	andq	$-64, %rsp
	subq	$1864, %rsp
	movq	%rdi, 48(%rsp)
	movq	%rsi, 40(%rsp)
	movq	%rdx, 32(%rsp)
	movq	%rcx, 24(%rsp)
	vmovapd	%zmm0, -56(%rsp)
	vmovapd	%zmm1, -120(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 1416(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 1352(%rsp)
	vmovapd	1416(%rsp), %zmm0
	vmulpd	1352(%rsp), %zmm0, %zmm0
	vmovapd	%zmm0, 1736(%rsp)
	vmovapd	1736(%rsp), %zmm0
	vmovapd	%zmm0, 1544(%rsp)
	vmovapd	-56(%rsp), %zmm0
	vmovapd	%zmm0, 1480(%rsp)
	vmovapd	1544(%rsp), %zmm0
	vmulpd	1480(%rsp), %zmm0, %zmm0
	vmovapd	%zmm0, 1672(%rsp)
	movl	$11, 1668(%rsp)
	movl	1668(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 1648(%rsp)
	vbroadcastsd	1648(%rsp), %zmm0
	movq	24(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movl	1668(%rsp), %eax
	decl	%eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 1656(%rsp)
	vbroadcastsd	1656(%rsp), %zmm0
	movq	32(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movl	1668(%rsp), %eax
	subl	$2, %eax
	movl	%eax, 1860(%rsp)
	jmp	.L21
.L26:
	movl	1860(%rsp), %eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 896(%rsp)
	vbroadcastsd	896(%rsp), %zmm0
	vmovapd	%zmm0, %zmm2
	movq	24(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1672(%rsp), %zmm1
	vmovapd	%zmm1, 1032(%rsp)
	vmovapd	%zmm0, 968(%rsp)
	vmovapd	%zmm2, 904(%rsp)
	vmovapd	1032(%rsp), %zmm0
	vmovapd	968(%rsp), %zmm2
	vmovapd	904(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	movq	24(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movl	1860(%rsp), %eax
	decl	%eax
	cltq
	leaq	0(,%rax,8), %rdx
	leaq	invfactorial(%rip), %rax
	vmovsd	(%rdx,%rax), %xmm0
	vmovsd	%xmm0, 1152(%rsp)
	vbroadcastsd	1152(%rsp), %zmm0
	vmovapd	%zmm0, %zmm2
	movq	32(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1672(%rsp), %zmm1
	vmovapd	%zmm1, 1288(%rsp)
	vmovapd	%zmm0, 1224(%rsp)
	vmovapd	%zmm2, 1160(%rsp)
	vmovapd	1288(%rsp), %zmm0
	vmovapd	1224(%rsp), %zmm2
	vmovapd	1160(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	movq	32(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	subl	$2, 1860(%rsp)
.L21:
	cmpl	$2, 1860(%rsp)
	jg	.L26
	vmovsd	.LC0(%rip), %xmm0
	vmovsd	%xmm0, 64(%rsp)
	vbroadcastsd	64(%rsp), %zmm0
	vmovapd	%zmm0, %zmm2
	movq	32(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1672(%rsp), %zmm1
	vmovapd	%zmm1, 200(%rsp)
	vmovapd	%zmm0, 136(%rsp)
	vmovapd	%zmm2, 72(%rsp)
	vmovapd	200(%rsp), %zmm0
	vmovapd	136(%rsp), %zmm2
	vmovapd	72(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k3}
	nop
	movq	48(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	24(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	%zmm0, 328(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 264(%rsp)
	vmovapd	328(%rsp), %zmm0
	vmulpd	264(%rsp), %zmm0, %zmm0
	movq	24(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	24(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	1672(%rsp), %zmm1
	vmovapd	%zmm1, 520(%rsp)
	vmovapd	%zmm0, 456(%rsp)
	vmovapd	-120(%rsp), %zmm0
	vmovapd	%zmm0, 392(%rsp)
	vmovapd	520(%rsp), %zmm0
	vmovapd	456(%rsp), %zmm2
	vmovapd	392(%rsp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k4
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k4}
	nop
	movq	40(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	24(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	%zmm0, 648(%rsp)
	vmovapd	1736(%rsp), %zmm0
	vmovapd	%zmm0, 584(%rsp)
	vmovapd	648(%rsp), %zmm0
	vmulpd	584(%rsp), %zmm0, %zmm0
	movq	24(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	movq	32(%rsp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	%zmm0, 776(%rsp)
	vmovapd	1736(%rsp), %zmm0
	vmovapd	%zmm0, 712(%rsp)
	vmovapd	776(%rsp), %zmm0
	vmulpd	712(%rsp), %zmm0, %zmm0
	movq	32(%rsp), %rax
	vmovapd	%zmm0, (%rax)
	nop
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4832:
	.size	mm_stiefel_Gs03_avx512, .-mm_stiefel_Gs03_avx512
	.globl	reb_whfast512_kepler_step
	.type	reb_whfast512_kepler_step, @function
reb_whfast512_kepler_step:
.LFB4833:
	.cfi_startproc
	leaq	8(%rsp), %r10
	.cfi_def_cfa 10, 0
	andq	$-64, %rsp
	pushq	-8(%r10)
	pushq	%rbp
	movq	%rsp, %rbp
	.cfi_escape 0x10,0x6,0x2,0x76,0
	pushq	%r10
	.cfi_escape 0xf,0x3,0x76,0x78,0x6
	subq	$16808, %rsp
	movq	%rdi, -16760(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	64(%rax), %zmm0
	vmovapd	%zmm0, -112(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	384(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	384(%rax), %zmm1
	vmovapd	%zmm1, -16432(%rbp)
	vmovapd	%zmm0, -16496(%rbp)
	vmovapd	-16432(%rbp), %zmm0
	vmulpd	-16496(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -176(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	448(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	448(%rax), %zmm1
	vmovapd	%zmm1, -16240(%rbp)
	vmovapd	%zmm0, -16304(%rbp)
	vmovapd	-176(%rbp), %zmm0
	vmovapd	%zmm0, -16368(%rbp)
	vmovapd	-16240(%rbp), %zmm0
	vmovapd	-16304(%rbp), %zmm2
	vmovapd	-16368(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	vmovapd	%zmm0, -176(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	512(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	512(%rax), %zmm1
	vmovapd	%zmm1, -16048(%rbp)
	vmovapd	%zmm0, -16112(%rbp)
	vmovapd	-176(%rbp), %zmm0
	vmovapd	%zmm0, -16176(%rbp)
	vmovapd	-16048(%rbp), %zmm0
	vmovapd	-16112(%rbp), %zmm2
	vmovapd	-16176(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	vmovapd	%zmm0, -176(%rbp)
	vmovapd	-176(%rbp), %zmm0
	vmovapd	%zmm0, -15920(%rbp)
	vmovapd	-15984(%rbp), %zmm0
	vmovapd	-15920(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vsqrtpd	%zmm1, %zmm0{%k3}
	nop
	vmovapd	%zmm0, -240(%rbp)
	vmovsd	.LC0(%rip), %xmm0
	vmovsd	%xmm0, -15800(%rbp)
	vbroadcastsd	-15800(%rbp), %zmm0
	vmovapd	%zmm0, -15728(%rbp)
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -15792(%rbp)
	vmovapd	-15728(%rbp), %zmm0
	vdivpd	-15792(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -304(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm1
	vmovapd	%zmm1, -15600(%rbp)
	vmovapd	%zmm0, -15664(%rbp)
	vmovapd	-15600(%rbp), %zmm0
	vmulpd	-15664(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -368(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm1
	vmovapd	%zmm1, -15408(%rbp)
	vmovapd	%zmm0, -15472(%rbp)
	vmovapd	-368(%rbp), %zmm0
	vmovapd	%zmm0, -15536(%rbp)
	vmovapd	-15408(%rbp), %zmm0
	vmovapd	-15472(%rbp), %zmm2
	vmovapd	-15536(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k4
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k4}
	nop
	vmovapd	%zmm0, -368(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm1
	vmovapd	%zmm1, -15216(%rbp)
	vmovapd	%zmm0, -15280(%rbp)
	vmovapd	-368(%rbp), %zmm0
	vmovapd	%zmm0, -15344(%rbp)
	vmovapd	-15216(%rbp), %zmm0
	vmovapd	-15280(%rbp), %zmm2
	vmovapd	-15344(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k5
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k5}
	nop
	vmovapd	%zmm0, -368(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	(%rax), %zmm1
	vmovsd	.LC1(%rip), %xmm0
	vmovsd	%xmm0, -15096(%rbp)
	vbroadcastsd	-15096(%rbp), %zmm0
	vmovapd	%zmm0, -15024(%rbp)
	vmovapd	%zmm1, -15088(%rbp)
	vmovapd	-15024(%rbp), %zmm0
	vmulpd	-15088(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -432(%rbp)
	vmovapd	-432(%rbp), %zmm0
	vmovapd	%zmm0, -14832(%rbp)
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm0, -14896(%rbp)
	vmovapd	-368(%rbp), %zmm0
	vmovapd	%zmm0, -14960(%rbp)
	vmovapd	-14832(%rbp), %zmm0
	vmovapd	-14896(%rbp), %zmm2
	vmovapd	-14960(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k6
	vfmsub132pd	%zmm2, %zmm1, %zmm0{%k6}
	nop
	vmovapd	%zmm0, -432(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	384(%rax), %zmm1
	vmovapd	%zmm1, -14704(%rbp)
	vmovapd	%zmm0, -14768(%rbp)
	vmovapd	-14704(%rbp), %zmm0
	vmulpd	-14768(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -496(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	448(%rax), %zmm1
	vmovapd	%zmm1, -14512(%rbp)
	vmovapd	%zmm0, -14576(%rbp)
	vmovapd	-496(%rbp), %zmm0
	vmovapd	%zmm0, -14640(%rbp)
	vmovapd	-14512(%rbp), %zmm0
	vmovapd	-14576(%rbp), %zmm2
	vmovapd	-14640(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k7
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k7}
	nop
	vmovapd	%zmm0, -496(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	512(%rax), %zmm1
	vmovapd	%zmm1, -14320(%rbp)
	vmovapd	%zmm0, -14384(%rbp)
	vmovapd	-496(%rbp), %zmm0
	vmovapd	%zmm0, -14448(%rbp)
	vmovapd	-14320(%rbp), %zmm0
	vmovapd	-14384(%rbp), %zmm2
	vmovapd	-14448(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	vmovapd	%zmm0, -496(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	(%rax), %zmm0
	vmovapd	-432(%rbp), %zmm1
	vmovapd	%zmm1, -14128(%rbp)
	vmovapd	-240(%rbp), %zmm1
	vmovapd	%zmm1, -14192(%rbp)
	vmovapd	%zmm0, -14256(%rbp)
	vmovapd	-14128(%rbp), %zmm0
	vmovapd	-14192(%rbp), %zmm2
	vmovapd	-14256(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	vmovapd	%zmm0, -560(%rbp)
	vmovapd	-112(%rbp), %zmm0
	vmovapd	%zmm0, -14000(%rbp)
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm0, -14064(%rbp)
	vmovapd	-14000(%rbp), %zmm0
	vmulpd	-14064(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -624(%rbp)
	vmovapd	-624(%rbp), %zmm0
	vmovapd	%zmm0, -13872(%rbp)
	vmovapd	-496(%rbp), %zmm0
	vmovapd	%zmm0, -13936(%rbp)
	vmovapd	-13872(%rbp), %zmm0
	vmulpd	-13936(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovsd	.LC2(%rip), %xmm0
	vmovsd	%xmm0, -13752(%rbp)
	vbroadcastsd	-13752(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -13680(%rbp)
	vmovapd	%zmm1, -13744(%rbp)
	vmovapd	-13680(%rbp), %zmm0
	vmulpd	-13744(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovsd	.LC0(%rip), %xmm0
	vmovsd	%xmm0, -13560(%rbp)
	vbroadcastsd	-13560(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -13424(%rbp)
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm0, -13488(%rbp)
	vmovapd	%zmm1, -13552(%rbp)
	vmovapd	-13424(%rbp), %zmm0
	vmovapd	-13488(%rbp), %zmm2
	vmovapd	-13552(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k3}
	nop
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-624(%rbp), %zmm0
	vmovapd	%zmm0, -13296(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -13360(%rbp)
	vmovapd	-13296(%rbp), %zmm0
	vmulpd	-13360(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-688(%rbp), %zmm1
	vmovapd	-432(%rbp), %zmm0
	leaq	-16688(%rbp), %rcx
	leaq	-16624(%rbp), %rdx
	leaq	-16560(%rbp), %rsi
	leaq	-16752(%rbp), %rax
	movq	%rax, %rdi
	call	mm_stiefel_Gs03_avx512
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -13104(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -13168(%rbp)
	vmovapd	-112(%rbp), %zmm0
	vmovapd	%zmm0, -13232(%rbp)
	vmovapd	-13104(%rbp), %zmm0
	vmovapd	-13168(%rbp), %zmm2
	vmovapd	-13232(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k4
	vfmsub132pd	%zmm2, %zmm1, %zmm0{%k4}
	nop
	vmovapd	%zmm0, -752(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -12912(%rbp)
	vmovapd	%zmm0, -12976(%rbp)
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -13040(%rbp)
	vmovapd	-12912(%rbp), %zmm0
	vmovapd	-12976(%rbp), %zmm2
	vmovapd	-13040(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k5
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k5}
	nop
	vmovapd	%zmm0, -752(%rbp)
	vmovapd	-16688(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -12720(%rbp)
	vmovapd	%zmm0, -12784(%rbp)
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -12848(%rbp)
	vmovapd	-12720(%rbp), %zmm0
	vmovapd	-12784(%rbp), %zmm2
	vmovapd	-12848(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k6
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k6}
	nop
	vmovapd	%zmm0, -752(%rbp)
	vmovapd	-16560(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -12528(%rbp)
	vmovapd	%zmm0, -12592(%rbp)
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -12656(%rbp)
	vmovapd	-12528(%rbp), %zmm0
	vmovapd	-12592(%rbp), %zmm2
	vmovapd	-12656(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k7
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k7}
	nop
	vmovapd	%zmm0, -816(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -12336(%rbp)
	vmovapd	%zmm0, -12400(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -12464(%rbp)
	vmovapd	-12336(%rbp), %zmm0
	vmovapd	-12400(%rbp), %zmm2
	vmovapd	-12464(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	vmovapd	%zmm0, -816(%rbp)
	vmovapd	-16752(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -12208(%rbp)
	vmovapd	%zmm0, -12272(%rbp)
	vmovapd	-12208(%rbp), %zmm0
	vmulpd	-12272(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -880(%rbp)
	vmovapd	-16560(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -12016(%rbp)
	vmovapd	%zmm0, -12080(%rbp)
	vmovapd	-880(%rbp), %zmm0
	vmovapd	%zmm0, -12144(%rbp)
	vmovapd	-12016(%rbp), %zmm0
	vmovapd	-12080(%rbp), %zmm2
	vmovapd	-12144(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	vmovapd	%zmm0, -880(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -11888(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -11952(%rbp)
	vmovapd	-11888(%rbp), %zmm0
	vmulpd	-11952(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -944(%rbp)
	vmovsd	.LC3(%rip), %xmm0
	vmovsd	%xmm0, -11768(%rbp)
	vbroadcastsd	-11768(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -11696(%rbp)
	vmovapd	%zmm1, -11760(%rbp)
	vmovapd	-11696(%rbp), %zmm0
	vmulpd	-11760(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -944(%rbp)
	vmovsd	.LC4(%rip), %xmm0
	vmovsd	%xmm0, -11576(%rbp)
	vbroadcastsd	-11576(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -11504(%rbp)
	vmovapd	-880(%rbp), %zmm0
	vmovapd	%zmm0, -11568(%rbp)
	vmovapd	-11504(%rbp), %zmm0
	vmulpd	-11568(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -11312(%rbp)
	vmovapd	%zmm1, -11376(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -11440(%rbp)
	vmovapd	-11312(%rbp), %zmm0
	vmovapd	-11376(%rbp), %zmm2
	vmovapd	-11440(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k3}
	nop
	vmovapd	%zmm0, -944(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -11184(%rbp)
	vmovapd	-11248(%rbp), %zmm0
	vmovapd	-11184(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k4
	vsqrtpd	%zmm1, %zmm0{%k4}
	nop
	vmovapd	%zmm0, -944(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -11056(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -11120(%rbp)
	vmovapd	-11056(%rbp), %zmm0
	vaddpd	-11120(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -944(%rbp)
	vmovsd	.LC5(%rip), %xmm0
	vmovsd	%xmm0, -10936(%rbp)
	vbroadcastsd	-10936(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -10864(%rbp)
	vmovapd	%zmm1, -10928(%rbp)
	vmovapd	-10864(%rbp), %zmm0
	vmulpd	-10928(%rbp), %zmm0, %zmm0
	vmovapd	-688(%rbp), %zmm1
	vmovapd	%zmm1, -10672(%rbp)
	vmovapd	-944(%rbp), %zmm1
	vmovapd	%zmm1, -10736(%rbp)
	vmovapd	%zmm0, -10800(%rbp)
	vmovapd	-10672(%rbp), %zmm0
	vmovapd	-10736(%rbp), %zmm2
	vmovapd	-10800(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k5
	vfmsub132pd	%zmm2, %zmm1, %zmm0{%k5}
	nop
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -10544(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -10608(%rbp)
	vmovapd	-10544(%rbp), %zmm0
	vdivpd	-10608(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-688(%rbp), %zmm1
	vmovapd	-432(%rbp), %zmm0
	leaq	-16688(%rbp), %rcx
	leaq	-16624(%rbp), %rdx
	leaq	-16560(%rbp), %rsi
	leaq	-16752(%rbp), %rax
	movq	%rax, %rdi
	call	mm_stiefel_Gs03_avx512
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -10352(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -10416(%rbp)
	vmovapd	-112(%rbp), %zmm0
	vmovapd	%zmm0, -10480(%rbp)
	vmovapd	-10352(%rbp), %zmm0
	vmovapd	-10416(%rbp), %zmm2
	vmovapd	-10480(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k6
	vfmsub132pd	%zmm2, %zmm1, %zmm0{%k6}
	nop
	vmovapd	%zmm0, -752(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -10160(%rbp)
	vmovapd	%zmm0, -10224(%rbp)
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -10288(%rbp)
	vmovapd	-10160(%rbp), %zmm0
	vmovapd	-10224(%rbp), %zmm2
	vmovapd	-10288(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k7
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k7}
	nop
	vmovapd	%zmm0, -752(%rbp)
	vmovapd	-16688(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -9968(%rbp)
	vmovapd	%zmm0, -10032(%rbp)
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -10096(%rbp)
	vmovapd	-9968(%rbp), %zmm0
	vmovapd	-10032(%rbp), %zmm2
	vmovapd	-10096(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	vmovapd	%zmm0, -752(%rbp)
	vmovapd	-16560(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -9776(%rbp)
	vmovapd	%zmm0, -9840(%rbp)
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -9904(%rbp)
	vmovapd	-9776(%rbp), %zmm0
	vmovapd	-9840(%rbp), %zmm2
	vmovapd	-9904(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	vmovapd	%zmm0, -816(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -9584(%rbp)
	vmovapd	%zmm0, -9648(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -9712(%rbp)
	vmovapd	-9584(%rbp), %zmm0
	vmovapd	-9648(%rbp), %zmm2
	vmovapd	-9712(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k3}
	nop
	vmovapd	%zmm0, -816(%rbp)
	vmovapd	-16752(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -9456(%rbp)
	vmovapd	%zmm0, -9520(%rbp)
	vmovapd	-9456(%rbp), %zmm0
	vmulpd	-9520(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -880(%rbp)
	vmovapd	-16560(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -9264(%rbp)
	vmovapd	%zmm0, -9328(%rbp)
	vmovapd	-880(%rbp), %zmm0
	vmovapd	%zmm0, -9392(%rbp)
	vmovapd	-9264(%rbp), %zmm0
	vmovapd	-9328(%rbp), %zmm2
	vmovapd	-9392(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k4
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k4}
	nop
	vmovapd	%zmm0, -880(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -9136(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -9200(%rbp)
	vmovapd	-9136(%rbp), %zmm0
	vmulpd	-9200(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -944(%rbp)
	vmovsd	.LC3(%rip), %xmm0
	vmovsd	%xmm0, -9016(%rbp)
	vbroadcastsd	-9016(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -8944(%rbp)
	vmovapd	%zmm1, -9008(%rbp)
	vmovapd	-8944(%rbp), %zmm0
	vmulpd	-9008(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -944(%rbp)
	vmovsd	.LC4(%rip), %xmm0
	vmovsd	%xmm0, -8824(%rbp)
	vbroadcastsd	-8824(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -8752(%rbp)
	vmovapd	-880(%rbp), %zmm0
	vmovapd	%zmm0, -8816(%rbp)
	vmovapd	-8752(%rbp), %zmm0
	vmulpd	-8816(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -8560(%rbp)
	vmovapd	%zmm1, -8624(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -8688(%rbp)
	vmovapd	-8560(%rbp), %zmm0
	vmovapd	-8624(%rbp), %zmm2
	vmovapd	-8688(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k5
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k5}
	nop
	vmovapd	%zmm0, -944(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -8432(%rbp)
	vmovapd	-8496(%rbp), %zmm0
	vmovapd	-8432(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k6
	vsqrtpd	%zmm1, %zmm0{%k6}
	nop
	vmovapd	%zmm0, -944(%rbp)
	vmovapd	-816(%rbp), %zmm0
	vmovapd	%zmm0, -8304(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -8368(%rbp)
	vmovapd	-8304(%rbp), %zmm0
	vaddpd	-8368(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -944(%rbp)
	vmovsd	.LC5(%rip), %xmm0
	vmovsd	%xmm0, -8184(%rbp)
	vbroadcastsd	-8184(%rbp), %zmm0
	vmovapd	%zmm0, %zmm1
	vmovapd	-752(%rbp), %zmm0
	vmovapd	%zmm0, -8112(%rbp)
	vmovapd	%zmm1, -8176(%rbp)
	vmovapd	-8112(%rbp), %zmm0
	vmulpd	-8176(%rbp), %zmm0, %zmm0
	vmovapd	-688(%rbp), %zmm1
	vmovapd	%zmm1, -7920(%rbp)
	vmovapd	-944(%rbp), %zmm1
	vmovapd	%zmm1, -7984(%rbp)
	vmovapd	%zmm0, -8048(%rbp)
	vmovapd	-7920(%rbp), %zmm0
	vmovapd	-7984(%rbp), %zmm2
	vmovapd	-8048(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k7
	vfmsub132pd	%zmm2, %zmm1, %zmm0{%k7}
	nop
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -7792(%rbp)
	vmovapd	-944(%rbp), %zmm0
	vmovapd	%zmm0, -7856(%rbp)
	vmovapd	-7792(%rbp), %zmm0
	vdivpd	-7856(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-688(%rbp), %zmm1
	vmovapd	-432(%rbp), %zmm0
	leaq	-16688(%rbp), %rdx
	leaq	-16624(%rbp), %rcx
	leaq	-16560(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	mm_stiefel_Gs13_avx512
	vmovapd	-16560(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -7664(%rbp)
	vmovapd	%zmm0, -7728(%rbp)
	vmovapd	-7664(%rbp), %zmm0
	vmulpd	-7728(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1008(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -7472(%rbp)
	vmovapd	%zmm0, -7536(%rbp)
	vmovapd	-1008(%rbp), %zmm0
	vmovapd	%zmm0, -7600(%rbp)
	vmovapd	-7472(%rbp), %zmm0
	vmovapd	-7536(%rbp), %zmm2
	vmovapd	-7600(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k1
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k1}
	nop
	vmovapd	%zmm0, -1008(%rbp)
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -7344(%rbp)
	vmovapd	-1008(%rbp), %zmm0
	vmovapd	%zmm0, -7408(%rbp)
	vmovapd	-7344(%rbp), %zmm0
	vaddpd	-7408(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1072(%rbp)
	vmovsd	.LC0(%rip), %xmm0
	vmovsd	%xmm0, -7224(%rbp)
	vbroadcastsd	-7224(%rbp), %zmm0
	vmovapd	%zmm0, -7152(%rbp)
	vmovapd	-1072(%rbp), %zmm0
	vmovapd	%zmm0, -7216(%rbp)
	vmovapd	-7152(%rbp), %zmm0
	vdivpd	-7216(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1072(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -7024(%rbp)
	vmovapd	-1008(%rbp), %zmm0
	vmovapd	%zmm0, -7088(%rbp)
	vmovapd	-7024(%rbp), %zmm0
	vmulpd	-7088(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -6832(%rbp)
	vmovapd	%zmm0, -6896(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -6960(%rbp)
	vmovapd	-6832(%rbp), %zmm0
	vmovapd	-6896(%rbp), %zmm2
	vmovapd	-6960(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k2
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k2}
	nop
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-16688(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -6640(%rbp)
	vmovapd	%zmm0, -6704(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -6768(%rbp)
	vmovapd	-6640(%rbp), %zmm0
	vmovapd	-6704(%rbp), %zmm2
	vmovapd	-6768(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k3
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k3}
	nop
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-112(%rbp), %zmm0
	vmovapd	%zmm0, -6512(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -6576(%rbp)
	vmovapd	-6512(%rbp), %zmm0
	vaddpd	-6576(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-1072(%rbp), %zmm0
	vmovapd	%zmm0, -6384(%rbp)
	vmovapd	-688(%rbp), %zmm0
	vmovapd	%zmm0, -6448(%rbp)
	vmovapd	-6384(%rbp), %zmm0
	vmulpd	-6448(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -688(%rbp)
	vmovapd	-688(%rbp), %zmm1
	vmovapd	-432(%rbp), %zmm0
	leaq	-16688(%rbp), %rdx
	leaq	-16624(%rbp), %rcx
	leaq	-16560(%rbp), %rax
	movq	%rcx, %rsi
	movq	%rax, %rdi
	call	mm_stiefel_Gs13_avx512
	vmovapd	-16560(%rbp), %zmm0
	vmovapd	-496(%rbp), %zmm1
	vmovapd	%zmm1, -6256(%rbp)
	vmovapd	%zmm0, -6320(%rbp)
	vmovapd	-6256(%rbp), %zmm0
	vmulpd	-6320(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1008(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	vmovapd	-560(%rbp), %zmm1
	vmovapd	%zmm1, -6064(%rbp)
	vmovapd	%zmm0, -6128(%rbp)
	vmovapd	-1008(%rbp), %zmm0
	vmovapd	%zmm0, -6192(%rbp)
	vmovapd	-6064(%rbp), %zmm0
	vmovapd	-6128(%rbp), %zmm2
	vmovapd	-6192(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k4
	vfmadd132pd	%zmm2, %zmm1, %zmm0{%k4}
	nop
	vmovapd	%zmm0, -1008(%rbp)
	vmovapd	-240(%rbp), %zmm0
	vmovapd	%zmm0, -5936(%rbp)
	vmovapd	-1008(%rbp), %zmm0
	vmovapd	%zmm0, -6000(%rbp)
	vmovapd	-5936(%rbp), %zmm0
	vaddpd	-6000(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1072(%rbp)
	vmovsd	.LC0(%rip), %xmm0
	vmovsd	%xmm0, -5816(%rbp)
	vbroadcastsd	-5816(%rbp), %zmm0
	vmovapd	%zmm0, -5744(%rbp)
	vmovapd	-1072(%rbp), %zmm0
	vmovapd	%zmm0, -5808(%rbp)
	vmovapd	-5744(%rbp), %zmm0
	vdivpd	-5808(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1072(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	(%rax), %zmm1
	vmovapd	%zmm1, -5616(%rbp)
	vmovapd	%zmm0, -5680(%rbp)
	vmovapd	-5616(%rbp), %zmm0
	vmulpd	-5680(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1136(%rbp)
	vmovapd	-1136(%rbp), %zmm0
	vmovapd	%zmm0, -5488(%rbp)
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm0, -5552(%rbp)
	vmovapd	-5488(%rbp), %zmm0
	vmulpd	-5552(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1136(%rbp)
	vmovapd	-16688(%rbp), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	(%rax), %zmm1
	vmovapd	%zmm1, -5296(%rbp)
	vmovapd	%zmm0, -5360(%rbp)
	vmovapd	-112(%rbp), %zmm0
	vmovapd	%zmm0, -5424(%rbp)
	vmovapd	-5296(%rbp), %zmm0
	vmovapd	-5360(%rbp), %zmm2
	vmovapd	-5424(%rbp), %zmm1
	movl	$-1, %eax
	kmovb	%eax, %k5
	vfnmadd132pd	%zmm2, %zmm1, %zmm0{%k5}
	nop
	vmovapd	%zmm0, -1200(%rbp)
	vmovapd	-16560(%rbp), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	(%rax), %zmm1
	vmovapd	%zmm1, -5168(%rbp)
	vmovapd	%zmm0, -5232(%rbp)
	vmovapd	-5168(%rbp), %zmm0
	vmulpd	-5232(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1264(%rbp)
	vmovapd	-1264(%rbp), %zmm0
	vmovapd	%zmm0, -5040(%rbp)
	vmovapd	-304(%rbp), %zmm0
	vmovapd	%zmm0, -5104(%rbp)
	vmovapd	-5040(%rbp), %zmm0
	vmulpd	-5104(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1264(%rbp)
	vmovapd	-1264(%rbp), %zmm0
	vmovapd	%zmm0, -4912(%rbp)
	vmovapd	-1072(%rbp), %zmm0
	vmovapd	%zmm0, -4976(%rbp)
	vmovapd	-4912(%rbp), %zmm0
	vmulpd	-4976(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1264(%rbp)
	vmovapd	-16624(%rbp), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	(%rax), %zmm1
	vmovapd	%zmm1, -4784(%rbp)
	vmovapd	%zmm0, -4848(%rbp)
	vmovapd	-4784(%rbp), %zmm0
	vmulpd	-4848(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1328(%rbp)
	vmovapd	-1328(%rbp), %zmm0
	vmovapd	%zmm0, -4656(%rbp)
	vmovapd	-1072(%rbp), %zmm0
	vmovapd	%zmm0, -4720(%rbp)
	vmovapd	-4656(%rbp), %zmm0
	vmulpd	-4720(%rbp), %zmm0, %zmm0
	vmovapd	%zmm0, -1328(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	384(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	384(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -4337(%rbp)
	vmovapd	-1136(%rbp), %zmm2
	vmovapd	%zmm2, -4464(%rbp)
	vmovapd	%zmm1, -4528(%rbp)
	vmovapd	%zmm0, -4592(%rbp)
	movzbl	-4337(%rbp), %eax
	vmovapd	-4464(%rbp), %zmm2
	vmovapd	-4528(%rbp), %zmm1
	vmovapd	-4592(%rbp), %zmm0
	kmovb	%eax, %k6
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k6}{z}
	nop
	vmovapd	%zmm0, -1392(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm0
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -4081(%rbp)
	vmovapd	-1200(%rbp), %zmm1
	vmovapd	%zmm1, -4208(%rbp)
	vmovapd	%zmm0, -4272(%rbp)
	vmovapd	-1392(%rbp), %zmm0
	vmovapd	%zmm0, -4336(%rbp)
	movzbl	-4081(%rbp), %eax
	vmovapd	-4208(%rbp), %zmm2
	vmovapd	-4272(%rbp), %zmm1
	vmovapd	-4336(%rbp), %zmm0
	kmovb	%eax, %k7
	vfmadd231pd	%zmm1, %zmm2, %zmm0{%k7}{z}
	nop
	vmovapd	%zmm0, -1392(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	448(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	448(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -3825(%rbp)
	vmovapd	-1136(%rbp), %zmm2
	vmovapd	%zmm2, -3952(%rbp)
	vmovapd	%zmm1, -4016(%rbp)
	vmovapd	%zmm0, -4080(%rbp)
	movzbl	-3825(%rbp), %eax
	vmovapd	-3952(%rbp), %zmm2
	vmovapd	-4016(%rbp), %zmm1
	vmovapd	-4080(%rbp), %zmm0
	kmovb	%eax, %k1
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k1}{z}
	nop
	vmovapd	%zmm0, -1456(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm0
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -3569(%rbp)
	vmovapd	-1200(%rbp), %zmm1
	vmovapd	%zmm1, -3696(%rbp)
	vmovapd	%zmm0, -3760(%rbp)
	vmovapd	-1456(%rbp), %zmm0
	vmovapd	%zmm0, -3824(%rbp)
	movzbl	-3569(%rbp), %eax
	vmovapd	-3696(%rbp), %zmm2
	vmovapd	-3760(%rbp), %zmm1
	vmovapd	-3824(%rbp), %zmm0
	kmovb	%eax, %k2
	vfmadd231pd	%zmm1, %zmm2, %zmm0{%k2}{z}
	nop
	vmovapd	%zmm0, -1456(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	512(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	512(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -3313(%rbp)
	vmovapd	-1136(%rbp), %zmm2
	vmovapd	%zmm2, -3440(%rbp)
	vmovapd	%zmm1, -3504(%rbp)
	vmovapd	%zmm0, -3568(%rbp)
	movzbl	-3313(%rbp), %eax
	vmovapd	-3440(%rbp), %zmm2
	vmovapd	-3504(%rbp), %zmm1
	vmovapd	-3568(%rbp), %zmm0
	kmovb	%eax, %k3
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k3}{z}
	nop
	vmovapd	%zmm0, -1520(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm0
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -3057(%rbp)
	vmovapd	-1200(%rbp), %zmm1
	vmovapd	%zmm1, -3184(%rbp)
	vmovapd	%zmm0, -3248(%rbp)
	vmovapd	-1520(%rbp), %zmm0
	vmovapd	%zmm0, -3312(%rbp)
	movzbl	-3057(%rbp), %eax
	vmovapd	-3184(%rbp), %zmm2
	vmovapd	-3248(%rbp), %zmm1
	vmovapd	-3312(%rbp), %zmm0
	kmovb	%eax, %k4
	vfmadd231pd	%zmm1, %zmm2, %zmm0{%k4}{z}
	nop
	vmovapd	%zmm0, -1520(%rbp)
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -2801(%rbp)
	vmovapd	-1328(%rbp), %zmm2
	vmovapd	%zmm2, -2928(%rbp)
	vmovapd	%zmm1, -2992(%rbp)
	vmovapd	%zmm0, -3056(%rbp)
	movzbl	-2801(%rbp), %eax
	vmovapd	-2928(%rbp), %zmm2
	vmovapd	-2992(%rbp), %zmm1
	vmovapd	-3056(%rbp), %zmm0
	kmovb	%eax, %k5
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k5}{z}
	nop
	movq	-16760(%rbp), %rax
	vmovapd	%zmm0, 576(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	576(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	384(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -2545(%rbp)
	vmovapd	-1264(%rbp), %zmm2
	vmovapd	%zmm2, -2672(%rbp)
	vmovapd	%zmm1, -2736(%rbp)
	vmovapd	%zmm0, -2800(%rbp)
	movzbl	-2545(%rbp), %eax
	vmovapd	-2672(%rbp), %zmm2
	vmovapd	-2736(%rbp), %zmm1
	vmovapd	-2800(%rbp), %zmm0
	kmovb	%eax, %k6
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k6}{z}
	nop
	movq	-16760(%rbp), %rax
	vmovapd	%zmm0, 576(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -2289(%rbp)
	vmovapd	-1328(%rbp), %zmm2
	vmovapd	%zmm2, -2416(%rbp)
	vmovapd	%zmm1, -2480(%rbp)
	vmovapd	%zmm0, -2544(%rbp)
	movzbl	-2289(%rbp), %eax
	vmovapd	-2416(%rbp), %zmm2
	vmovapd	-2480(%rbp), %zmm1
	vmovapd	-2544(%rbp), %zmm0
	kmovb	%eax, %k7
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k7}{z}
	nop
	movq	-16760(%rbp), %rax
	vmovapd	%zmm0, 640(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	640(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	448(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -2033(%rbp)
	vmovapd	-1264(%rbp), %zmm2
	vmovapd	%zmm2, -2160(%rbp)
	vmovapd	%zmm1, -2224(%rbp)
	vmovapd	%zmm0, -2288(%rbp)
	movzbl	-2033(%rbp), %eax
	vmovapd	-2160(%rbp), %zmm2
	vmovapd	-2224(%rbp), %zmm1
	vmovapd	-2288(%rbp), %zmm0
	kmovb	%eax, %k1
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k1}{z}
	nop
	movq	-16760(%rbp), %rax
	vmovapd	%zmm0, 640(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -1777(%rbp)
	vmovapd	-1328(%rbp), %zmm2
	vmovapd	%zmm2, -1904(%rbp)
	vmovapd	%zmm1, -1968(%rbp)
	vmovapd	%zmm0, -2032(%rbp)
	movzbl	-1777(%rbp), %eax
	vmovapd	-1904(%rbp), %zmm2
	vmovapd	-1968(%rbp), %zmm1
	vmovapd	-2032(%rbp), %zmm0
	kmovb	%eax, %k2
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k2}{z}
	nop
	movq	-16760(%rbp), %rax
	vmovapd	%zmm0, 704(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	704(%rax), %zmm0
	movq	-16760(%rbp), %rax
	vmovapd	512(%rax), %zmm1
	movq	-16760(%rbp), %rax
	movzbl	2368(%rax), %eax
	movzbl	%al, %eax
	movb	%al, -1521(%rbp)
	vmovapd	-1264(%rbp), %zmm2
	vmovapd	%zmm2, -1648(%rbp)
	vmovapd	%zmm1, -1712(%rbp)
	vmovapd	%zmm0, -1776(%rbp)
	movzbl	-1521(%rbp), %eax
	vmovapd	-1648(%rbp), %zmm2
	vmovapd	-1712(%rbp), %zmm1
	vmovapd	-1776(%rbp), %zmm0
	kmovb	%eax, %k3
	vfnmadd231pd	%zmm1, %zmm2, %zmm0{%k3}{z}
	nop
	movq	-16760(%rbp), %rax
	vmovapd	%zmm0, 704(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	-1392(%rbp), %zmm0
	vmovapd	%zmm0, 384(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	-1456(%rbp), %zmm0
	vmovapd	%zmm0, 448(%rax)
	movq	-16760(%rbp), %rax
	vmovapd	-1520(%rbp), %zmm0
	vmovapd	%zmm0, 512(%rax)
	nop
	movq	-8(%rbp), %r10
	.cfi_def_cfa 10, 0
	leave
	leaq	-8(%r10), %rsp
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE4833:
	.size	reb_whfast512_kepler_step, .-reb_whfast512_kepler_step
	.section	.rodata
	.align 8
.LC0:
	.long	0
	.long	1072693248
	.align 8
.LC1:
	.long	0
	.long	1073741824
	.align 8
.LC2:
	.long	0
	.long	1071644672
	.align 8
.LC3:
	.long	0
	.long	1076887552
	.align 8
.LC4:
	.long	0
	.long	1077149696
	.align 8
.LC5:
	.long	0
	.long	1075052544
	.ident	"GCC: (Debian 12.2.0-14+deb12u1) 12.2.0"
	.section	.note.GNU-stack,"",@progbits
