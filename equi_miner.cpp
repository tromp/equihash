// Wagner's algorithm for Generalized Birthday Paradox, a memory-hard proof-of-work
// Copyright (c) 2016 John Tromp

#include "equi_miner.h"
#include <unistd.h>
#include "ctype.h"

int hextobyte(const char * x) {
  u32 b = 0;
  for (int i = 0; i < 2; i++) {
    uchar c = tolower(x[i]);
    assert(isxdigit(c));
    b = (b << 4) | (c - (c >= '0' && c <= '9' ? '0' : ('a' - 10)));
  }
  return b;
}

int main(int argc, char **argv) {
  int nthreads = 1;
  int nonce = 0;
  int range = 1;
  bool showsol = false;
  const char *header = "";
  const char *hex = "";
  int c;
  while ((c = getopt (argc, argv, "h:n:r:t:x:s")) != -1) {
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
      case 'x':
        hex = optarg;
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
  printf("Looking for wagner-tree on (\"%s\",%d", hex ? "0x..." : header, nonce);
  if (range > 1)
    printf("-%d", nonce+range-1);
  printf(") with %d %d-bit digits and %d threads\n", NDIGITS, DIGITBITS, nthreads);
  thread_ctx *threads = (thread_ctx *)calloc(nthreads, sizeof(thread_ctx));
  assert(threads);
  equi eq(nthreads);
  printf("Using %dMB of memory", 1 + eq.hta.alloced / 0x100000);
#ifdef USE_AVX2
  printf(" and AVX2 intrinsics to compute 4-way blake2b\n");
#else
  printf(" and no AVX2\n");
#endif
  u32 sumnsols = 0;
  char headernonce[HEADERNONCELEN];
  u32 hdrlen = strlen(header);
  if (*hex) {
    assert(strlen(hex) == 2 * HEADERNONCELEN);
    for (int i = 0; i < HEADERNONCELEN; i++)
      headernonce[i] = hextobyte(&hex[2*i]);
  } else {
    memcpy(headernonce, header, hdrlen);
    memset(headernonce+hdrlen, 0, sizeof(headernonce)-hdrlen);
  }
  for (int r = 0; r < range; r++) {
    ((u32 *)headernonce)[32] = htole32(nonce+r);
    eq.setheadernonce(headernonce, sizeof(headernonce));
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
    u32 nsols, maxsols = min(MAXSOLS, eq.nsols);
    for (nsols = 0; nsols < maxsols; nsols++) {
      if (showsol) {
        printf("\nSolution");
        for (u32 i = 0; i < PROOFSIZE; i++)
          printf(" %jx", (uintmax_t)eq.sols[nsols][i]);
      }
    }
    printf("\n%d solutions\n", nsols);
    sumnsols += nsols;
  }
  free(threads);
  printf("%d total solutions\n", sumnsols);
  return 0;
}
