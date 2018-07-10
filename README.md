# equihash
Equihash proof-of-work solvers

Build with "make all"

The executables ending in 1 are compiled without atomics and thus
only suitable for single-threaded use (where they get some speedup over the generic version).

Options -h HEADER -n NONCE are self explanatory.
Non-ascii headers can be provided in hexadecimal with option -x HEXHEADER.
Add option -r RANGESIZE to search a range of nonces.
Add option -p PREFIX to change intial characters of the personalization string (not implemented in assembly and CUDA solvers).
For benching, options -n 255 -r 100 are useful as it gets exactly 188 solutions from 100 runs.
Finally, option -s shows the solutions.

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

On Thursday Nov 3, user elbandi on Slack reported a bug in verify() where it allows a non-zero
final digit in the top-level xor. That is now fixed.

On November 11, I implemented an interleaved 8-way blake, but this turned out to provide no gain.

On November 17, I added Cantor coding for slot pairs, as found in xenoncat's and morpav's solvers.
This allows the use of 2^10 buckets for (200,9) which turns out to be a small gain,
so I made this the new default.

I implemented prefetching for memory writes, but found no gain, and left the code out.

More detailed documentation is available in the equi_miner.h source code.

Performance summary (on 4GHz i7-4790K and NVidia GTX980):

- equi1:        4.9 Sol/s
- equix41:      6.2 Sol/s
- eqasm1:       6.9 Sol/s

- 8 x equi1:   22.2 Sol/s
- 8 x equix41: 27.1 Sol/s
- 8 x eqasm1:  27.2 Sol/s

- eqcuda:      27.2 Sol/s

And now, for something completely different: (144,5) taking 2.6 GB of memory

- eq1445   -t 8: 1.0 Sol/s
- eq1445x4 -t 8: 1.2 Sol/s

- eqcuda1445:    2.2 Sol/s

Contest judges requested the following information:

1. A brief self-assessment of your submission per the published judging criteria.

- testibility is integrated into the submission by provision of a 
    int verify(proof indices, const char *headernonce, const u32 headerlen);
  routine, and standalone verifier equi.c. This is part of the default make targets
  together with tests for both the (200,9) and (144,5) parameters.
- despite lack of implementation of the suggested API, the implemented API of
    equi(const u32 n_threads);
  for solver construction, with methods
    void setheadernonce(const char *headernonce, const u32 len);
    void digit0(const u32 id);
    void digitodd(const u32 r, const u32 id);
    void digiteven(const u32 r, const u32 id);
    void digitK(const u32 id);
  and specialized unrolled versions
    void digit1(const u32 id);
    ...
    void digit8(const u32 id);
  have proved practical enough to support integration into zcashd and nicehash miners.
- the submission is written with portability in mind, with no dependencies beyond pthreads,
  no architectural assumptions like word size or endian-ness, and using a subset of C++ features
  (i.e. no templates) for ease of porting to plain C.
- SIMD support is available in two ways:
  1) through an included blake2b reference impolementation that's been modified to make compression
     rounds strict rather than lazy, allowing for computation of an actual midstate
  2) through a custom 4-way blake2b implementation using intrinsics based on Samuel Neves' blake2bp code
- the implementation supports (200,9) and (144,5) out of the box, and can easily adapt to other
  parameters by changing a few lines to select the appropriate bit segements from the hash.
- memory is already minimized to the point of losing a tiny fraction of solutions (much less than 1%),
  but can trivially be reduced further with a compile time define, at the cost of more discarding.
- file equi_miner.h contains both a problem description as well as a very rough algorithm overview,
  followed by a slightly more detailed overview in lines 243--277. beyond that many single line
  comments can be found throughout the code.
- a list of post-deadline improvements may be found above, as well as expected performance.
- the solution rate has been measured as 1.88 Sol/run,
  with a fraction of about 0.002 of solutions discarded.
- due to static allocation, average and peak memory conincide.
- runtime varies only slightly (a few %) from run ro run.

2. An explanation about what you think are the strengths of your submission.

- relatively portable, concise and straightforward code
- support for multiple parameter sets including (144,5)
- multi-threading support (crucial for 144,5)
- support for a wide range of buckets
- support for storing part of hash in treenode to further optimize space use
- support for CUDA devices
- minimal staticically allocated memory use
- optional visualization of bucket size distribution

3. An explanation about what you think are the weaknesses of your submission.

- single threaded x86 performance on (200,9) lags somewhat behind other solvers
- lacking a 2-way blake SIMD implementation
- CUDA solver treats GPU as a many core cpu and fails to take better advantage of architectural
  features.
