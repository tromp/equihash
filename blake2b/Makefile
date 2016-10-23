all: example_avx1 example_avx2

example_avx1: example_avx1.c zcblake2_avx1.o
	gcc -o example_avx1 example_avx1.c zcblake2_avx1.o

example_avx2: example_avx2.c zcblake2_avx2.o
	gcc -o example_avx2 example_avx2.c zcblake2_avx2.o

.PHONY: clean
clean:
	rm example_avx1 example_avx2
