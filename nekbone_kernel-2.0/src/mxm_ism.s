/*
Routine to make general matrix matrix multiplication. lxd below should be the
maximum value of n2 for the multiplications
Registers:
	r15: a + (n2-1)*n1*8
	r14: b + 8*n2*n3 - 8*j*n2
	rdi: c + 8*n3*n1 - 8*j*n1
	rsi: &n1 at the beginning, then r15 - 8*k*n1 in loop_mult
	rcx: Temporary variable and loop counter in loop_mult and load_bs0
	rdx: unused (Contains b at calling time)
	r8 : 8*n1 (Contains c at calling time)
	r9 : n2, then 8*n2 (Contains &n3 at calling time)
	r11: Decreasing i from 8*n1.
	r10: n3, then j
Unused:
	rax, rbx, rdx, rbp, rsp, r12, r13
*/
.equ lxd, 12

#uses R  %r15, %rdi, %r8, %r9, %r11, %ymm0
#uses  W %rsi, %rcx, %ymm2, %ymm1
.macro loop_mult jmp_position, mask=0
.if \mask
	#We know we are at the top row i = 0
	movq	%r15, %rsi
.else
	movl	%r11d, %esi
	addq	%r15, %rsi
.endif # address a + (n2-1)*n1*8 + i*8 in rsi

.if \mask
	vmaskmovpd	(%rsi), %ymm0, %ymm1
.else
	vmovupd		(%rsi), %ymm1
.endif  # a loaded in ymm1

	movl	%r9d, %ecx
	subl	$8, %ecx
	vmulpd	bs0(,%rcx,4), %ymm1, %ymm1
	jz	\jmp_position\()_end	# Result of the subl

\jmp_position:
	subq		%r8, %rsi	# address a + i*8 + k*n1*8 in rsi
.if \mask
	vmaskmovpd	(%rsi), %ymm0, %ymm2
.else
	vmovupd		(%rsi), %ymm2
.endif
	subl		$8, %ecx
	vfmadd231pd	bs0(,%rcx,4), %ymm2, %ymm1
	jnz		\jmp_position	# Result of the subl
\jmp_position\()_end:
.if \mask
	vmaskmovpd	%ymm1, %ymm0, (%rdi)
.else
	vmovupd		%ymm1, (%rdi,%r11)
.endif
.endm

#uses R  %r9, %r14
#uses  W %rcx, %ymm1
.macro load_bs0_array
	movl		%r9d, %ecx

load_bs0_array:
	subl		$8, %ecx
	vbroadcastsd	(%r14,%rcx), %ymm1
	vmovapd		%ymm1, bs0(,%rcx,4)
	jnz		load_bs0_array	# Result of the subl
.endm


	.file "mxm_ism_test_assembly.s"
.section ".data", "aw"
	.align 32
	bs0:
	.fill (32*lxd)
.section ".rodata", "a"
	.align 32
mask:
	.long	0xffffffff,0xffffffff
	.long	0xffffffff,0xffffffff
	.long	0xffffffff,0xffffffff
	.long	0xffffffff,0xffffffff

	.long	0xffffffff,0xffffffff
	.long	0x00000000,0x00000000
	.long	0x00000000,0x00000000
	.long	0x00000000,0x00000000

	.long	0xffffffff,0xffffffff
	.long	0xffffffff,0xffffffff
	.long	0x00000000,0x00000000
	.long	0x00000000,0x00000000

	.long	0xffffffff,0xffffffff
	.long	0xffffffff,0xffffffff
	.long	0xffffffff,0xffffffff
	.long	0x00000000,0x00000000

	.long	0x00000000,0x00000000
	.long	0x00000000,0x00000000
	.long	0x00000000,0x00000000
	.long	0x00000000,0x00000000

	.type	mask,@object
	.size	mask,160

.text
	.align	2,0x90
	.globl mxm_ism_
mxm_ism_:
# parameter 1: %rdi  a
# parameter 2: %rsi  n1
# parameter 3: %rdx  b
# parameter 4: %rcx  n2
# parameter 5: %r8   c
# parameter 6: %r9   n3
	#pushq	%rbx
	#pushq	%rbp
	#pushq	%r12
	#pushq	%r13
	pushq	%r14
	pushq	%r15

	movq	%rdi, %r15	#Save a address in r15
	movq	%rdx, %r14	#Save b address in r14
	movq	%r8, %rdi	#Save c address in rdi

	# Make       n1 = r8d/8,  n2 = r9d/8,  n3 = r10d (we use them a lot)
	# And later:  i = r11d/8,   j= r10d,  k = ecx   (n3 is not actually used)
	# j <-> n3, i <-> n1, k <-> n2. loops nested in that order
	# a : n1 rows, n2 columns
	# b : n2 rows, n3 columns
	# c : n1 rows, n3 columns
	movl	(%r9), %r10d
	movl	(%rsi), %r8d
	movl	(%rcx), %r9d

	shll	$3, %r8d	#We only use 8*n1

	# Address a + (n2-1)*n1*8 in %r15
	movl	%r9d, %ecx
	decl	%ecx
	imull	%r8d, %ecx
	addq	%rcx, %r15

	shll	$3, %r9d	#We only use 8*n2

	# Address c + n3*n1(*8) in rdi
	movl	%r10d, %ecx
	imull	%r8d, %ecx
	addq	%rcx, %rdi

	#Create a mask for nbr of rows non-multiple of 4
	movl	%r8d, %ecx
	andl	$24, %ecx
	vmovupd	mask(,%rcx,4), %ymm0

	# Address b + n2*n3(*8) in r14
	movl	%r10d, %ecx
	imull	%r9d, %ecx
	addq	%rcx, %r14

	#Initialisation of j-loop.
	decl	%r10d

for_j_loop:
	subq	%r9, %r14
	load_bs0_array

	#Initialisation of i-loop
	movl	%r8d, %r11d
	subq	%r8, %rdi
	subl	$32, %r11d
	jle	while_i_loop_end

while_i_loop:
	loop_mult	for_k_loop
	subl		$32, %r11d
	jg		while_i_loop

while_i_loop_end:
	# Same as the above loop, but with maskload (in ymm0)
	loop_mult	for_k_loop_in_n1, 1
	decl		%r10d
	jge		for_j_loop

end_of_function:
	popq	%r15
	popq	%r14
	#popq	%r13
	#popq	%r12
	#popq	%rbp
	#popq	%rbx
	ret
	.align	2,0x90
	.type	mxm_ism_,@function
	.size	mxm_ism_,.-mxm_ism_
