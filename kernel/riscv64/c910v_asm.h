#ifndef C910V_ASM_H
#define C910V_ASM_H

/*
 * Xuantie GCC accepts the unprefixed RVV 0.7 mnemonics. Upstream GCC
 * and binutils spell the same instructions with a th. prefix when the
 * architecture is xtheadvector. The macros are empty on the Xuantie
 * toolchain, so the assembly text is unchanged there.
 */
#if defined(__riscv_xtheadvector)
#define C910V_ASM_ENTER \
	".macro vle.v args:vararg\n\t" \
	"th.vle.v \\args\n\t" \
	".endm\n\t" \
	".macro vse.v args:vararg\n\t" \
	"th.vse.v \\args\n\t" \
	".endm\n\t" \
	".macro vsetvli args:vararg\n\t" \
	"th.vsetvli \\args\n\t" \
	".endm\n\t" \
	".macro vfmv.v.f args:vararg\n\t" \
	"th.vfmv.v.f \\args\n\t" \
	".endm\n\t" \
	".macro vfmacc.vv args:vararg\n\t" \
	"th.vfmacc.vv \\args\n\t" \
	".endm\n\t" \
	".macro vfadd.vv args:vararg\n\t" \
	"th.vfadd.vv \\args\n\t" \
	".endm\n\t" \
	".macro vfmul.vv args:vararg\n\t" \
	"th.vfmul.vv \\args\n\t" \
	".endm\n\t" \
	".macro vrgather.vi args:vararg\n\t" \
	"th.vrgather.vi \\args\n\t" \
	".endm\n\t"
#define C910V_ASM_LEAVE \
	".purgem vle.v\n\t" \
	".purgem vse.v\n\t" \
	".purgem vsetvli\n\t" \
	".purgem vfmv.v.f\n\t" \
	".purgem vfmacc.vv\n\t" \
	".purgem vfadd.vv\n\t" \
	".purgem vfmul.vv\n\t" \
	".purgem vrgather.vi\n\t"
#else
#define C910V_ASM_ENTER
#define C910V_ASM_LEAVE
#endif

#endif
