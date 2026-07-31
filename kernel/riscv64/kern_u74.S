/*
 * dgemm 4x4 scalar micro-kernel for SiFive U74 (RV64GC, lp64d).  v3:
 * software-pipelined, full operand double-buffering (P/Q ping-pong), one-
 * iteration lookahead, LOAD-BEFORE-FMA issue ordering (the U74 dual-issue
 * pairs a load/int in slot0 with an FP op in slot1, so each fld is emitted
 * immediately before the fmadd it co-issues with).  16 FP accumulators + 16 FP
 * operand regs = all 32 FP regs; alpha spilled for the k-loop.  Handles the
 * bm%4==0 && bn%4==0 grid with EVEN bk (bk>=2); ragged edges + odd bk later.
 *
 * int kern_asm(long bm,long bn,long bk,double alpha,double* A,double* B,double* C,long ldc)
 *   a0=bm a1=bn a2=bk fa0=alpha a3=A a4=B a5=C a6=ldc
 * acc: col0 ft0..ft3  col1 ft4..ft7  col2 ft8..ft11  col3 fa2..fa5
 * P: A=fa0,fa1,fa6,fa7  B=fs0,fs1,fs2,fs3    Q: A=fs4,fs5,fs6,fs7  B=fs8,fs9,fs10,fs11
 * int: t0=ldc*8 t1=bm/4 t2=j t3=ptrA t4=ptrB t5=i t6=pair-cnt  C0..3=a0,a1,a6,a7
 */
	.text
	.p2align 4
	.globl	kern_asm
	.type	kern_asm,@function
kern_asm:
	addi	sp, sp, -112
	fsd	fs0, 0(sp)
	fsd	fs1, 8(sp)
	fsd	fs2, 16(sp)
	fsd	fs3, 24(sp)
	fsd	fs4, 32(sp)
	fsd	fs5, 40(sp)
	fsd	fs6, 48(sp)
	fsd	fs7, 56(sp)
	fsd	fs8, 64(sp)
	fsd	fs9, 72(sp)
	fsd	fs10, 80(sp)
	fsd	fs11, 88(sp)
	fsd	fa0, 96(sp)		# spill alpha

	slli	t0, a6, 3		# ldc*8
	srli	t1, a0, 2		# bm/4
	srli	t2, a1, 2		# bn/4 (j)
	beqz	t2, .Ldone

.Lj:
	mv	a0, a5			# C0 = Crun
	add	a1, a0, t0		# C1
	add	a6, a1, t0		# C2
	add	a7, a6, t0		# C3
	mv	t3, a3			# ptrA = ba
	mv	t5, t1			# i = bm/4

.Li:
	mv	t4, a4			# ptrB = bb(run)
	fmv.d.x	ft0, zero
	fmv.d.x	ft1, zero
	fmv.d.x	ft2, zero
	fmv.d.x	ft3, zero
	fmv.d.x	ft4, zero
	fmv.d.x	ft5, zero
	fmv.d.x	ft6, zero
	fmv.d.x	ft7, zero
	fmv.d.x	ft8, zero
	fmv.d.x	ft9, zero
	fmv.d.x	ft10, zero
	fmv.d.x	ft11, zero
	fmv.d.x	fa2, zero
	fmv.d.x	fa3, zero
	fmv.d.x	fa4, zero
	fmv.d.x	fa5, zero

	# preload P from k=0, advance to k=1
	fld	fa0, 0(t3)
	fld	fa1, 8(t3)
	fld	fa6, 16(t3)
	fld	fa7, 24(t3)
	fld	fs0, 0(t4)
	fld	fs1, 8(t4)
	fld	fs2, 16(t4)
	fld	fs3, 24(t4)
	addi	t3, t3, 32
	addi	t4, t4, 32
	srli	t6, a2, 1		# bk/2
	addi	t6, t6, -1		# pair-cnt = bk/2 - 1
	beqz	t6, .Lepi

	.p2align 4
.Lbody:				# --- P current, load Q (load before paired fmadd) ---
	fld	fs4, 0(t3)
	fmadd.d	ft0, fa0, fs0, ft0
	fld	fs5, 8(t3)
	fmadd.d	ft1, fa1, fs0, ft1
	fld	fs6, 16(t3)
	fmadd.d	ft2, fa6, fs0, ft2
	fld	fs7, 24(t3)
	fmadd.d	ft3, fa7, fs0, ft3
	fld	fs8, 0(t4)
	fmadd.d	ft4, fa0, fs1, ft4
	fld	fs9, 8(t4)
	fmadd.d	ft5, fa1, fs1, ft5
	fld	fs10, 16(t4)
	fmadd.d	ft6, fa6, fs1, ft6
	fld	fs11, 24(t4)
	fmadd.d	ft7, fa7, fs1, ft7
	addi	t3, t3, 32
	fmadd.d	ft8, fa0, fs2, ft8
	addi	t4, t4, 32
	fmadd.d	ft9, fa1, fs2, ft9
	fmadd.d	ft10, fa6, fs2, ft10
	fmadd.d	ft11, fa7, fs2, ft11
	fmadd.d	fa2, fa0, fs3, fa2
	fmadd.d	fa3, fa1, fs3, fa3
	fmadd.d	fa4, fa6, fs3, fa4
	fmadd.d	fa5, fa7, fs3, fa5
				# --- Q current, load P ---
	fld	fa0, 0(t3)
	fmadd.d	ft0, fs4, fs8, ft0
	fld	fa1, 8(t3)
	fmadd.d	ft1, fs5, fs8, ft1
	fld	fa6, 16(t3)
	fmadd.d	ft2, fs6, fs8, ft2
	fld	fa7, 24(t3)
	fmadd.d	ft3, fs7, fs8, ft3
	fld	fs0, 0(t4)
	fmadd.d	ft4, fs4, fs9, ft4
	fld	fs1, 8(t4)
	fmadd.d	ft5, fs5, fs9, ft5
	fld	fs2, 16(t4)
	fmadd.d	ft6, fs6, fs9, ft6
	fld	fs3, 24(t4)
	fmadd.d	ft7, fs7, fs9, ft7
	addi	t3, t3, 32
	fmadd.d	ft8, fs4, fs10, ft8
	addi	t4, t4, 32
	fmadd.d	ft9, fs5, fs10, ft9
	fmadd.d	ft10, fs6, fs10, ft10
	fmadd.d	ft11, fs7, fs10, ft11
	fmadd.d	fa2, fs4, fs11, fa2
	fmadd.d	fa3, fs5, fs11, fa3
	fmadd.d	fa4, fs6, fs11, fa4
	fmadd.d	fa5, fs7, fs11, fa5
	addi	t6, t6, -1
	bnez	t6, .Lbody

.Lepi:				# --- P current, load last Q ---
	fld	fs4, 0(t3)
	fmadd.d	ft0, fa0, fs0, ft0
	fld	fs5, 8(t3)
	fmadd.d	ft1, fa1, fs0, ft1
	fld	fs6, 16(t3)
	fmadd.d	ft2, fa6, fs0, ft2
	fld	fs7, 24(t3)
	fmadd.d	ft3, fa7, fs0, ft3
	fld	fs8, 0(t4)
	fmadd.d	ft4, fa0, fs1, ft4
	fld	fs9, 8(t4)
	fmadd.d	ft5, fa1, fs1, ft5
	fld	fs10, 16(t4)
	fmadd.d	ft6, fa6, fs1, ft6
	fld	fs11, 24(t4)
	fmadd.d	ft7, fa7, fs1, ft7
	addi	t3, t3, 32
	fmadd.d	ft8, fa0, fs2, ft8
	addi	t4, t4, 32
	fmadd.d	ft9, fa1, fs2, ft9
	fmadd.d	ft10, fa6, fs2, ft10
	fmadd.d	ft11, fa7, fs2, ft11
	fmadd.d	fa2, fa0, fs3, fa2
	fmadd.d	fa3, fa1, fs3, fa3
	fmadd.d	fa4, fa6, fs3, fa4
	fmadd.d	fa5, fa7, fs3, fa5
				# --- Q current, compute only ---
	fmadd.d	ft0, fs4, fs8, ft0
	fmadd.d	ft1, fs5, fs8, ft1
	fmadd.d	ft2, fs6, fs8, ft2
	fmadd.d	ft3, fs7, fs8, ft3
	fmadd.d	ft4, fs4, fs9, ft4
	fmadd.d	ft5, fs5, fs9, ft5
	fmadd.d	ft6, fs6, fs9, ft6
	fmadd.d	ft7, fs7, fs9, ft7
	fmadd.d	ft8, fs4, fs10, ft8
	fmadd.d	ft9, fs5, fs10, ft9
	fmadd.d	ft10, fs6, fs10, ft10
	fmadd.d	ft11, fs7, fs10, ft11
	fmadd.d	fa2, fs4, fs11, fa2
	fmadd.d	fa3, fs5, fs11, fa3
	fmadd.d	fa4, fs6, fs11, fa4
	fmadd.d	fa5, fs7, fs11, fa5

	# ---- C += acc*alpha ; alpha->fs4, temps fs0/fs1 ----
	fld	fs4, 96(sp)
	fld	fs0, 0(a0)
	fmadd.d	fs0, ft0, fs4, fs0
	fsd	fs0, 0(a0)
	fld	fs1, 8(a0)
	fmadd.d	fs1, ft1, fs4, fs1
	fsd	fs1, 8(a0)
	fld	fs0, 16(a0)
	fmadd.d	fs0, ft2, fs4, fs0
	fsd	fs0, 16(a0)
	fld	fs1, 24(a0)
	fmadd.d	fs1, ft3, fs4, fs1
	fsd	fs1, 24(a0)
	fld	fs0, 0(a1)
	fmadd.d	fs0, ft4, fs4, fs0
	fsd	fs0, 0(a1)
	fld	fs1, 8(a1)
	fmadd.d	fs1, ft5, fs4, fs1
	fsd	fs1, 8(a1)
	fld	fs0, 16(a1)
	fmadd.d	fs0, ft6, fs4, fs0
	fsd	fs0, 16(a1)
	fld	fs1, 24(a1)
	fmadd.d	fs1, ft7, fs4, fs1
	fsd	fs1, 24(a1)
	fld	fs0, 0(a6)
	fmadd.d	fs0, ft8, fs4, fs0
	fsd	fs0, 0(a6)
	fld	fs1, 8(a6)
	fmadd.d	fs1, ft9, fs4, fs1
	fsd	fs1, 8(a6)
	fld	fs0, 16(a6)
	fmadd.d	fs0, ft10, fs4, fs0
	fsd	fs0, 16(a6)
	fld	fs1, 24(a6)
	fmadd.d	fs1, ft11, fs4, fs1
	fsd	fs1, 24(a6)
	fld	fs0, 0(a7)
	fmadd.d	fs0, fa2, fs4, fs0
	fsd	fs0, 0(a7)
	fld	fs1, 8(a7)
	fmadd.d	fs1, fa3, fs4, fs1
	fsd	fs1, 8(a7)
	fld	fs0, 16(a7)
	fmadd.d	fs0, fa4, fs4, fs0
	fsd	fs0, 16(a7)
	fld	fs1, 24(a7)
	fmadd.d	fs1, fa5, fs4, fs1
	fsd	fs1, 24(a7)

	addi	a0, a0, 32
	addi	a1, a1, 32
	addi	a6, a6, 32
	addi	a7, a7, 32
	addi	t5, t5, -1
	bnez	t5, .Li

	slli	t6, a2, 5		# bb(run) += bk*4
	add	a4, a4, t6
	slli	t6, t0, 2		# C(run) += 4*ldc
	add	a5, a5, t6
	addi	t2, t2, -1
	bnez	t2, .Lj

.Ldone:
	fld	fs0, 0(sp)
	fld	fs1, 8(sp)
	fld	fs2, 16(sp)
	fld	fs3, 24(sp)
	fld	fs4, 32(sp)
	fld	fs5, 40(sp)
	fld	fs6, 48(sp)
	fld	fs7, 56(sp)
	fld	fs8, 64(sp)
	fld	fs9, 72(sp)
	fld	fs10, 80(sp)
	fld	fs11, 88(sp)
	addi	sp, sp, 112
	li	a0, 0
	ret
	.size	kern_asm, .-kern_asm
