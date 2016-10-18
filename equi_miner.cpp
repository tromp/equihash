// Wagner's algorithm for Generalized Birthday Paradox, a memory-hard proof-of-work
// Copyright (c) 2016 John Tromp

#include "equi_miner.h"
#include <unistd.h>

int main(int argc, char **argv) {
  int nthreads = 1;
  int nonce = 0;
  int range = 1;
  bool showsol = false;
  const char *header = "";
  int c;
  while ((c = getopt (argc, argv, "h:n:r:t:s")) != -1) {
    switch (c) {
      case 'h':
        header = optarg;
        break;
      case 'n':
        nonce = atoi(optarg);
        break;
      case 'r':
        range = atoi(optarg);
        break;
      case 's':
        showsol = true;
        break;
      case 't':
        nthreads = atoi(optarg);
        break;
    }
  }
#ifndef XWITHASH
  if (sizeof(tree) > 4)
    printf("WARNING: please compile with -DXWITHASH to shrink tree!\n");
#endif
#ifdef ATOMIC
  if (nthreads==1)
    printf("WARNING: use of atomics hurts single threaded performance!\n");
#else
  assert(nthreads==1);
#endif
  printf("Looking for wagner-tree on (\"%s\",%d", header, nonce);
  if (range > 1)
    printf("-%d", nonce+range-1);
  printf(") with %d %d-bit digits and %d threads\n", NDIGITS, DIGITBITS, nthreads);
  thread_ctx *threads = (thread_ctx *)calloc(nthreads, sizeof(thread_ctx));
  assert(threads);
  equi eq(nthreads);
  printf("Using %dMB of memory\n", eq.hta.alloced >> 20);
  u32 sumnsols = 0;
  for (int r = 0; r < range; r++) {
    eq.setnonce(header, strlen(header), nonce+r);
    for (int t = 0; t < nthreads; t++) {
      threads[t].id = t;
      threads[t].eq = &eq;
      int err = pthread_create(&threads[t].thread, NULL, worker, (void *)&threads[t]);
      assert(err == 0);
    }
    for (int t = 0; t < nthreads; t++) {
      int err = pthread_join(threads[t].thread, NULL);
      assert(err == 0);
    }
    u32 nsols = 0;
    for (unsigned s = 0; s < eq.nsols; s++) {
      nsols++;
      if (showsol) {
        printf("Solution");
        for (u32 i = 0; i < PROOFSIZE; i++)
          printf(" %jx", (uintmax_t)eq.sols[s][i]);
        printf("\n");
      }
    }
    printf("%d solutions\n", nsols);
    sumnsols += nsols;
  }
  free(threads);
  printf("%d total solutions\n", sumnsols);
  return 0;
}
