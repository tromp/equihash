#include "equi.h"
#include <inttypes.h> // for SCNx64 macro
#include <stdio.h>    // printf/scanf
#include <stdlib.h>   // exit
#include <unistd.h>   // getopt
#include <string.h>   // strlen
#include <assert.h>   // d'uh

int main(int argc, char **argv) {
  const char *header = "";
  char personal[] = "ZcashPoW";
  int prefixlen, nonce = 0, c;
  while ((c = getopt (argc, argv, "h:n:p:")) != -1) {
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
    }
  }
  printf("Verifying size %d proof for %s(\"%s\",%d)\n",
               PROOFSIZE, personal, header, nonce);
  char headernonce[HEADERNONCELEN];
  u32 hdrlen = strlen(header);
  memcpy(headernonce, header, hdrlen);
  memset(headernonce+hdrlen, 0, sizeof(headernonce)-hdrlen);
  ((u32 *)headernonce)[32] = htole32(nonce);
  for (int nsols=0; scanf(" Solution") == 0; nsols++) {
    u32 indices[PROOFSIZE];
    for (int n = 0; n < PROOFSIZE; n++) {
      int nscan = scanf(" %x", &indices[n]);
      assert(nscan == 1);
    }
    int pow_rc = verify(indices, headernonce, sizeof(headernonce), personal);
    if (pow_rc == POW_OK)
      printf("Verified\n");
    else
      printf("FAILED due to %s\n", errstr[pow_rc]);
  }
  return 0;
}
