// Wagner's algorithm for Generalized Birthday Paradox, a memory-hard proof-of-work
// Copyright (c) 2016 John Tromp

#include "equi_miner.h"
#include <unistd.h>
#include "ctype.h"

#define COMPRESSED_SOL_SIZE (PROOFSIZE * (DIGITBITS + 1) / 8)

int hextobyte(const char * x) {
  u32 b = 0;
  for (int i = 0; i < 2; i++) {
    uchar c = tolower(x[i]);
    assert(isxdigit(c));
    b = (b << 4) | (c - (c >= '0' && c <= '9' ? '0' : ('a' - 10)));
  }
  return b;
}

void compress_solution(const proof sol, uchar* csol) {
  uchar b;
  for (u32 i = 0, j = 0, bits_left = DIGITBITS + 1;
      j < COMPRESSED_SOL_SIZE; csol[j++] = b) {
    if (bits_left >=8) {
      // Read next 8 bits, stay at same sol index
      b = sol[i] >> (bits_left -= 8);
    } else { // less than 8 bits to read
      // Read remaining bits and shift left to make space for next sol index
      b = sol[i];
      b <<= (8 - bits_left); // may also set b=0 if bits_left was 0, which is fine
      // Go to next sol index and read remaining bits
      bits_left += DIGITBITS + 1 - 8;
      b |= sol[++i] >> bits_left;
    }
  }
}

int main(int argc, char **argv) {
  int nthreads = 1;
  int nonce = 0;
  int range = 1;
  bool showsol = false;
  bool compress_sol = false;
  const char *header = "";
  const char *hex = "";
  char personal[] = "ZcashPoW";
  int prefixlen, c;
  while ((c = getopt (argc, argv, "h:n:p:r:t:x:sc")) != -1) {
    switch (c) {
      case 'h':
        header = optarg;
        break;
      case 'n':
        nonce = atoi(optarg);
        break;
      case 'p':
        prefixlen = strlen(optarg);
        assert(prefixlen <= 8);
        memcpy((void *)personal, optarg, prefixlen);
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
      case 'c':
        compress_sol = true;
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
  printf("Looking for wagner-tree on %s(\"%s\",%d", personal, hex ? "0x..." : header, nonce);
  if (range > 1)
    printf("-%d", nonce+range-1);
  printf(") with %d %d-bit digits and %d threads\n", NDIGITS, DIGITBITS, nthreads);
  thread_ctx *threads = (thread_ctx *)calloc(nthreads, sizeof(thread_ctx));
  assert(threads);
  equi eq(nthreads);
  printf("Using 2^%d buckets, %dMB of memory, and %d-way blake2b\n", BUCKBITS, 1 + eq.hta.alloced / 0x100000, NBLAKES);
#ifdef ASM_BLAKE
  printf("Using xenoncat's assembly blake code\n");
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
    ((u32 *)headernonce)[27] = htole32(nonce+r);
    eq.setheadernonce(headernonce, sizeof(headernonce), personal);
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
        if (compress_sol) {
          printf(" ");
          uchar csol[COMPRESSED_SOL_SIZE];
          compress_solution(eq.sols[nsols], csol);
          for (u32 i = 0; i < COMPRESSED_SOL_SIZE; ++i) {
            printf("%02hhx", csol[i]);
          }
        } else {
          for (u32 i = 0; i < PROOFSIZE; i++)
            printf(" %jx", (uintmax_t)eq.sols[nsols][i]);
        }
      }
    }
    printf("\n%d solutions\n", nsols);
    sumnsols += nsols;
  }
  free(threads);
  printf("%d total solutions\n", sumnsols);
  return 0;
}
