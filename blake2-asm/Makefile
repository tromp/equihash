all: example_avx1 example_avx2

example_avx1: example_avx1.c asm/zcblake2_avx1.o
	gcc -o example_avx1 example_avx1.c asm/zcblake2_avx1.o

example_avx2: example_avx2.c asm/zcblake2_avx2.o
	gcc -o example_avx2 example_avx2.c asm/zcblake2_avx2.o

asm/zcblake2_avx1.o:
	make -C asm

asm/zcblake2_avx2.o:
	make -C asm

.PHONY: clean
clean:
	make -C asm clean && rm -f example_avx1 example_avx2 example_avx1.o example_avx2.o out.bin
