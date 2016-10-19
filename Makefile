OPT   = -O3
FLAGS = -Wall -Wno-deprecated-declarations -D_POSIX_C_SOURCE=200112L $(OPT) -pthread 
GPP   = g++ -march=native -m64 -std=c++11 $(FLAGS)

all:	equi equi1 verify test spark

equi:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DATOMIC equi_miner.cpp blake/blake2b.cpp -o equi

equi1:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) equi_miner.cpp blake/blake2b.cpp -o equi1

equi1g:	equi.h equi_miner.h equi_miner.cpp Makefile
	g++ -g -DLOGSPARK -DSPARKSCALE=11 equi_miner.cpp blake/blake2b.cpp -pthread -o equi1g

equi1445:	equi.h equi_miner.h equi_miner.cpp Makefile
	$(GPP) -DRESTBITS=4 -DWN=144 -DWK=5 equi_miner.cpp blake/blake2b.cpp -o equi1445

dev1:	equi.h dev_miner.h dev_miner.cpp Makefile
	$(GPP) -DRESTBITS=8 dev_miner.cpp blake/blake2b.cpp -o dev1

eqcuda:	equi_miner.cu equi.h blake2b.cu Makefile
	nvcc -DXINTREE -DUNROLL -arch sm_35 equi_miner.cu blake/blake2b.cpp -o eqcuda

devcuda:	dev_miner.cu equi.h blake2b.cu Makefile
	nvcc -DXINTREE -DUNROLL -arch sm_35 dev_miner.cu blake/blake2b.cpp -o devcuda

eqcuda1445:	equi_miner.cu equi.h blake2b.cu Makefile
	nvcc -DWN=144 -DWK=5 -arch sm_35 equi_miner.cu blake/blake2b.cpp -o eqcuda1445

verify:	equi.h equi.c Makefile
	g++ -g equi.c blake/blake2b.cpp -o verify

bench:	equi1
	time ./equi1 -r 10

test:	equi verify Makefile
	time ./equi -h "" -n 0 -t 1 -s | grep ^Sol | ./verify -h "" -n 0

spark:	equi1g
	time ./equi1g

clean:	
	rm equi equi1 equi1g equi1445 eqcuda eqcuda1445 feqcuda verify
