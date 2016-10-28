all: zcblake2_avx1.o zcblake2_avx2.o

%.o: %.asm
	fasm $<

clean:
	rm -f zcblake2_avx1.o zcblake2_avx2.o
