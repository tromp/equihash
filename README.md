# equihash
Equihash proof-of-work solvers


The executables ending in 1 are compiled without atomics and thus
only suitable for single-threaded use (where they get some speedup over the generic version).

Options -h <HEADER> -n <NONCE> are self explanatory.
Add option -r <RANGESIZE> to search a range of nonces.
