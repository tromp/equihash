OPT   = -O3
FLAGS = -Wall -Wno-deprecated-declarations -D_POSIX_C_SOURCE=200112L $(OPT) -pthread 
GPP   = g++ -march=native -m64 -std=c++11 $(FLAGS)

all:	equi equi1 verify test spark test1445

equi:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DATOMIC equi_miner.cpp blake/blake2b.cpp -o equi

equi1:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) equi_miner.cpp blake/blake2b.cpp -o equi1

equix4:	equi.h equi_miner.h equi_miner.cpp blake2-avx2/blake2bip.c Makefile
	$(GPP) -mavx2 -DNBLAKES=4 -DATOMIC equi_miner.cpp blake/blake2b.cpp blake2-avx2/blake2bip.c -o equix4

equix41:	equi.h equi_miner.h equi_miner.cpp blake2-avx2/blake2bip.c Makefile
	$(GPP) -mavx2 -DNBLAKES=4 equi_miner.cpp blake/blake2b.cpp blake2-avx2/blake2bip.c -o equix41

equix81:	equi.h equi_miner.h equi_miner.cpp blake2-avx2/blake2bip.c Makefile
	$(GPP) -mavx2 -DNBLAKES=8 equi_miner.cpp blake/blake2b.cpp blake2-avx2/blake2bip.c -o equix81

equi1g:	equi.h equi_miner.h equi_miner.cpp Makefile
	g++ -g -std=c++11 -DLOGSPARK -DSPARKSCALE=11 equi_miner.cpp blake/blake2b.cpp -pthread -o equi1g

eq485:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DATOMIC -DWN=48 -DWK=5 -DRESTBITS=4 equi_miner.cpp blake/blake2b.cpp -o eq485

eq4851:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DWN=48 -DWK=5 -DRESTBITS=4 equi_miner.cpp blake/blake2b.cpp -o eq4851

eq1445:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DATOMIC -DRESTBITS=4 -DWN=144 -DWK=5 equi_miner.cpp blake/blake2b.cpp -o eq1445

eq1686:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DATOMIC -DRESTBITS=4 -DWN=168 -DWK=6 equi_miner.cpp blake/blake2b.cpp -o eq1686

eq1927:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DATOMIC -DRESTBITS=4 -DWN=192 -DWK=7 equi_miner.cpp blake/blake2b.cpp -o eq1927

eq21091:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DRESTBITS=5 -DWN=210 -DWK=9 equi_miner.cpp blake/blake2b.cpp -o eq21091

eq14451:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DRESTBITS=4 -DWN=144 -DWK=5 equi_miner.cpp blake/blake2b.cpp -o eq14451

eq1445x4:	equi.h equi_miner.h equi_miner.cpp blake2-avx2/blake2bip.c Makefile
	$(GPP) -DATOMIC -mavx2 -DNBLAKES=4 -DRESTBITS=4 -DWN=144 -DWK=5 equi_miner.cpp blake/blake2b.cpp blake2-avx2/blake2bip.c -o eq1445x4

eq1445x41:	equi.h equi_miner.h equi_miner.cpp blake2-avx2/blake2bip.c Makefile
	$(GPP) -mavx2 -DNBLAKES=4 -DRESTBITS=4 -DWN=144 -DWK=5 equi_miner.cpp blake/blake2b.cpp blake2-avx2/blake2bip.c -o eq1445x41

eqasm:	equi.h equi_miner.h equi_miner.cpp blake2-asm/asm/zcblake2_avx2.o Makefile
	$(GPP) -mavx2 -DASM_BLAKE -DATOMIC equi_miner.cpp blake/blake2b.cpp blake2-asm/asm/zcblake2_avx2.o -o eqasm

eqasm1:	equi.h equi_miner.h equi_miner.cpp blake2-asm/asm/zcblake2_avx2.o Makefile
	$(GPP) -mavx2 -DASM_BLAKE equi_miner.cpp blake/blake2b.cpp blake2-asm/asm/zcblake2_avx2.o -o eqasm1

eqcuda:	equi_miner.cu equi.h blake2b.cu Makefile
	nvcc -DXINTREE -arch sm_35 equi_miner.cu blake/blake2b.cpp -o eqcuda

eqcudah:	equi_miner.cu equi.h blake2b.cu Makefile
	nvcc -DHIST -DXINTREE -arch sm_35 equi_miner.cu blake/blake2b.cpp -o eqcudah

eqcuda1445:	equi_miner.cu equi.h blake2b.cu Makefile
	nvcc -DWN=144 -DWK=5 -arch sm_35 equi_miner.cu blake/blake2b.cpp -o eqcuda1445

eqcuda1927:	equi_miner.cu equi.h blake2b.cu Makefile
	nvcc -DWN=192 -DWK=7 -arch sm_35 equi_miner.cu blake/blake2b.cpp -o eqcuda1927

verify:	equi.h equi.c Makefile
	g++ -g equi.c blake/blake2b.cpp -o verify

verify1445:	equi.h equi.c Makefile
	g++ -DRESTBITS=4 -DWN=144 -DWK=5 -g equi.c blake/blake2b.cpp -o verify1445

verify1686:	equi.h equi.c Makefile
	g++ -DRESTBITS=4 -DWN=168 -DWK=6 -g equi.c blake/blake2b.cpp -o verify1686

bench:	equi1
	time ./equi1 -r 10

test:	equi1 verify Makefile
	time ./equi1 -h "" -n 0 -t 1 -s | grep ^Sol | ./verify -h "" -n 0

test1445:	eq14451 verify1445 Makefile
	time ./eq14451 -h "" -n 0 -t 1 -s | grep ^Sol | ./verify1445 -h "" -n 0

spark:	equi1g
	time ./equi1g

blake2-asm/asm/zcblake2_avx1.o:
	make -C blake2-asm

blake2-asm/asm/zcblake2_avx2.o:
	make -C blake2-asm

clean:	
	make -C blake2b clean && rm -f eqasm eqasm1 equi equi1 equix4 equix41 equi1g eq485 eq4851 eq1445 eq14451 eq1445x4 eq1445x41 eqcuda eqcuda1445 verify
