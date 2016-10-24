# equihash
Equihash proof-of-work solvers

Build with "make all"

The executables ending in 1 are compiled without atomics and thus
only suitable for single-threaded use (where they get some speedup over the generic version).

Options -h HEADER -n NONCE are self explanatory.
Add option -r RANGESIZE to search a range of nonces.

Build the faster version courtesy of xenoncat's AVX2 4-way parallel blake2b assembly code with
"make dev1" and bench with "time ./dev1 -n 255 -r 100"
