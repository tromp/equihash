# equihash
Equihash proof-of-work solvers

Build with "make all"

The executables ending in 1 are compiled without atomics and thus
only suitable for single-threaded use (where they get some speedup over the generic version).

Options -h HEADER -n NONCE are self explanatory.
Add option -r RANGESIZE to search a range of nonces.
For benching, options -n 255 -r 100 are useful as it gets exactly 188 solutions from 100 runs.

My original submission was triggered by seeing how xenoncat's
"has much of the same ideas as mine" so that making my open sourcing conditional
on getting sufficient funding for the Cuckoo Cycle Bounty Fund no longer made sense.

https://forum.z.cash/t/tromps-solvers/2465/76

I noticed that we both use bucket sorting with tree nodes stored as a directed acyclic graph.
Upon original submission, I wrote: Compared to xenoncat, my solver differs in
- having way more buckets,
- wasting some memory,
- having simpler pair compression,
- being multi-threaded,
- and supporting (144,5).
- And of course in not using any assembly.
- Oh, and having some cool visualization of bucket size distribution...

David Jaenson gave me the idea to disable atomics for single threaded operation,
which gave a nice speed boost for that case.

Since then I reduced the number of buckets in the cpu solver from 2^16 to 2^12,
which allowed for reducing the bucket space overhead. I borrowed from xenoncat
the idea to allocate all memory statically, and found a way to improve upon his memory layout,
reducing waste by about 7%. My solver now needs only 144MB compared to xenoncat's 178MB.

Seeing that my solver was spending 45% of runtime on hashing, I asked xenoncat if (s)he
could make their assembly blake2b implementation available through a C binding, which s(he)
very generously did. My solver executables using this are called dev1/dev.

Zooko had earlier suggested looking at Samuel Neves' blake2bp implemention for faster hashing.
After initially rejecting this approach due to different blake2bp semantics, I came back to 
to it in search of a more portable way to gain speed. I managed to bridge the semantic gap
and modify Samuel's source to serve Equihash's purposes.

On the morning of the submission deadline day, discussion on sorting with judge Solardiz
made me realize that my 2nd stage bucket sort could benefit from linking rather than listing
xor-able slots, which gave me the final speed boost.

More detailed documentation is available in the equi_miner.h source code.

Performance summary (on 4GHz i7-4790K and NVidia GTX980):

- equi1:        4.6 Sol/s
- eqavx21:      5.9 Sol/s
- dev1:         6.5 Sol/s

- equi -t 8:   13.7 Sol/s
- eqavx2 -t 8: 16.7 Sol/s
- dev -t 8:    17.2 Sol/s

- 8 x eqavx21: 20.3 Sol/s
- 8 x dev1:    20.6 Sol/s

- eqcuda:      23.6 Sol/s

And now, for something completely different: (144,5) taking 2.6 GB of memory

- eq1445 -t 8:     1.0 Sol/s
- eq1445avx2 -t 8: 1.2 Sol/s

- eqcuda1445:      2.2 Sol/s
