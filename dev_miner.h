// Equihash solver
// Copyright (c) 2016 John Tromp

// This development version uses xenoncat's highly optimized
// 4-way parallel blake2b implementation

// Fix N, K, such that n = N/(k+1) is integer
// Fix M = 2^{n+1} hashes each of length N bits,
// H_0, ... , H_{M-1}, generated fom (n+1)-bit indices.
// Problem: find binary tree on 2^K distinct indices,
// for which the exclusive-or of leaf hashes is all 0s.
// Additionally, it should satisfy the Wagner conditions:
// for each height i subtree, the exclusive-or
// of its 2^i corresponding hashes starts with i*n 0 bits,
// and for i>0 the leftmost leaf of its left subtree
// is less than the leftmost leaf of its right subtree

// The algorithm below solves this by maintaining the tree
// in a graph of K layers, each split into buckets
// with buckets indexed by the first n-RESTBITS bits following
// the i*n 0s, each bucket having 4 * 2^RESTBITS slots,
// twice the number of subtrees expected to land there.

#include "equi.h"
#include <stdio.h>
#include <pthread.h>
#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif
void Blake2PrepareMidstate4(void *midstate, uchar *input);
#ifdef __cplusplus
}
#endif
//midstate: 256 bytes of buffer for output midstate, aligned by 32
//input: 140 bytes header, preferably aligned by 8

#ifdef __cplusplus
extern "C" {
#endif
void Blake2Run4(uchar *hashout, void *midstate, u32 indexctr);
#ifdef __cplusplus
}
#endif
//hashout: hash output buffer: 4*64 bytes
//midstate: 256 bytes from Blake2PrepareMidstate4
//indexctr: For n=200, k=9: {0, 4, 8, ..., 1048572}

#if defined __builtin_bswap32 && defined __LITTLE_ENDIAN
#undef htobe32
#define htobe32(x) __builtin_bswap32(x)
#elif defined __APPLE__
#undef htobe32
#define htobe32(x) OSSwapHostToBigInt32(x)
#endif

typedef uint16_t u16;
typedef uint64_t u64;

#ifdef ATOMIC
#include <atomic>
typedef std::atomic<u32> au32;
#else
typedef u32 au32;
#endif

#ifndef RESTBITS
#define RESTBITS	8
#endif

// 2_log of number of buckets
#define BUCKBITS (DIGITBITS-RESTBITS)

#ifndef SAVEMEM
#if RESTBITS == 4
// can't save memory in such small buckets
#define SAVEMEM 1
#elif RESTBITS >= 8
// take advantage of law of large numbers (sum of 2^8 random numbers)
// this reduces (200,9) memory to under 144MB, with negligible discarding
#define SAVEMEM 9/14
#endif
#endif

// number of buckets
static const u32 NBUCKETS = 1<<BUCKBITS;
// bucket mask
static const u32 BUCKMASK = NBUCKETS-1;
// 2_log of number of slots per bucket
static const u32 SLOTBITS = RESTBITS+1+1;
static const u32 SLOTRANGE = 1<<SLOTBITS;
static const u32 SLOTMSB = 1<<(SLOTBITS-1);
// number of slots per bucket
static const u32 NSLOTS = SLOTRANGE * SAVEMEM;
// SLOTBITS mask
static const u32 SLOTMASK = SLOTRANGE-1;
// number of possible values of xhash (rest of n) bits
static const u32 NRESTS = 1<<RESTBITS;
// nothing larger found in 100000 runs
static const u32 MAXSOLS = 8;

// tree node identifying its children as two different slots in
// a bucket on previous layer with the same rest bits (x-tra hash)
struct tree {
  u32 bid_s0_s1; // manual bitfields

  tree(const u32 idx) {
    bid_s0_s1 = idx;
  }
  tree(const u32 bid, const u32 s0, const u32 s1) {
#ifdef SLOTDIFF
    u32 ds10 = (s1 - s0) & SLOTMASK;
    if (ds10 & SLOTMSB) {
      bid_s0_s1 = (((bid << SLOTBITS) | s1) << (SLOTBITS-1)) | (SLOTMASK & ~ds10);
    } else {
      bid_s0_s1 = (((bid << SLOTBITS) | s0) << (SLOTBITS-1)) | (ds10 - 1);
    }
#else
    bid_s0_s1 = (((bid << SLOTBITS) | s0) << SLOTBITS) | s1;
#endif
  }
  u32 getindex() const {
    return bid_s0_s1;
  }
  u32 bucketid() const {
#ifdef SLOTDIFF
    return bid_s0_s1 >> (2 * SLOTBITS - 1);
#else
    return bid_s0_s1 >> (2 * SLOTBITS);
#endif
  }
  u32 slotid0() const {
#ifdef SLOTDIFF
    return (bid_s0_s1 >> (SLOTBITS-1)) & SLOTMASK;
#else
    return (bid_s0_s1 >> SLOTBITS) & SLOTMASK;
#endif
  }
  u32 slotid1() const {
#ifdef SLOTDIFF
    return (slotid0() + 1 + (bid_s0_s1 & (SLOTMASK>>1))) & SLOTMASK;
#else
    return bid_s0_s1 & SLOTMASK;
#endif
  }
};

union htunit {
  tree tag;
  u32 word;
  uchar bytes[sizeof(u32)];
};

#define WORDS(bits)	((bits + 31) / 32)
#define HASHWORDS0 WORDS(WN - DIGITBITS + RESTBITS)
#define HASHWORDS1 WORDS(WN - 2*DIGITBITS + RESTBITS)

// A slot is up to HASHWORDS0 hash units followed by a tag
typedef htunit slot0[HASHWORDS0+1];
typedef htunit slot1[HASHWORDS1+1];
// a bucket is NSLOTS treenodes
typedef slot0 bucket0[NSLOTS];
typedef slot1 bucket1[NSLOTS];
// the N-bit hash consists of K+1 n-bit "digits"
// each of which corresponds to a layer of NBUCKETS buckets
typedef bucket0 digit0[NBUCKETS];
typedef bucket1 digit1[NBUCKETS];
typedef au32 bsizes[NBUCKETS];

u32 min(const u32 a, const u32 b) {
  return a < b ? a : b;
}

// size (in bytes) of hash in round 0 <= r < WK
u32 hashsize(const u32 r) {
  const u32 hashbits = WN - (r+1) * DIGITBITS + RESTBITS;
  return (hashbits + 7) / 8;
}

u32 hashwords(u32 bytes) {
  return (bytes + 3) / 4;
}

// manages hash and tree data
struct htalloc {
  bucket0 *heap0;
  bucket1 *heap1;
  u32 alloced;
  htalloc() {
    alloced = 0;
  }
  void alloctrees() {
// optimize xenoncat's fixed memory layout, avoiding any waste
// digit  hashes   tree   hashes tree
// 0      A A A A A A 0   . . . . . .
// 1      A A A A A A 0   B B B B B 1
// 2      C C C C C 2 0   B B B B B 1
// 3      C C C C C 2 0   D D D D 3 1
// 4      E E E E 4 2 0   D D D D 3 1
// 5      E E E E 4 2 0   F F F 5 3 1
// 6      G G 6 . 4 2 0   F F F 5 3 1
// 7      G G 6 . 4 2 0   H H 7 5 3 1
// 8      I 8 6 . 4 2 0   H H 7 5 3 1
    assert(DIGITBITS >= 16); // ensures hashes shorten by 1 unit every 2 digits
    heap0 = (bucket0 *)alloc(NBUCKETS, sizeof(bucket0));
    heap1 = (bucket1 *)alloc(NBUCKETS, sizeof(bucket1));
  }
  void dealloctrees() {
    free(heap0);
    free(heap1);
  }
  void *alloc(const u32 n, const u32 sz) {
    void *mem  = calloc(n, sz);
    assert(mem);
    alloced += n * sz;
    return mem;
  }
};

struct equi {
  alignas(32) uchar blake_ctx[256];
  htalloc hta;
  bsizes *nslots;
  proof *sols;
  au32 nsols;
  u32 nthreads;
  u32 bfull;
  u32 hfull;
  pthread_barrier_t barry;
  equi(const u32 n_threads) {
    assert(sizeof(htunit) == 4);
    assert(WK&1); // assumed in candidate() calling indices1()
    nthreads = n_threads;
    const int err = pthread_barrier_init(&barry, NULL, nthreads);
    assert(!err);
    hta.alloctrees();
    nslots = (bsizes *)hta.alloc(2 * NBUCKETS, sizeof(au32));
    sols   =  (proof *)hta.alloc(MAXSOLS, sizeof(proof));
  }
  ~equi() {
    hta.dealloctrees();
    free(nslots);
    free(sols);
  }
  void setheadernonce(const char *headernonce, const u32 len) {
    alignas(8) uchar alignheader[HEADERNONCELEN];
    memcpy(alignheader, headernonce, len);
    assert(len == HEADERNONCELEN);
    alignas(32) uchar midstate[256];
    Blake2PrepareMidstate4(midstate, alignheader);
    memcpy(&blake_ctx, midstate, 256);
    memset(nslots, 0, NBUCKETS * sizeof(au32)); // only nslots[0] needs zeroing
    nsols = bfull = hfull = 0;
  }
  u32 getslot0(const u32 bucketi) {
#ifdef ATOMIC
    return std::atomic_fetch_add_explicit(&nslots[0][bucketi], 1U, std::memory_order_relaxed);
#else
    return nslots[0][bucketi]++;
#endif
  }
  u32 getslot1(const u32 bucketi) {
#ifdef ATOMIC
    return std::atomic_fetch_add_explicit(&nslots[1][bucketi], 1U, std::memory_order_relaxed);
#else
    return nslots[1][bucketi]++;
#endif
  }
  u32 getnslots0(const u32 bid) {
    au32 &nslot = nslots[0][bid];
    const u32 n = min(nslot, NSLOTS);
    nslot = 0;
    return n;
  }
  u32 getnslots1(const u32 bid) {
    au32 &nslot = nslots[1][bid];
    const u32 n = min(nslot, NSLOTS);
    nslot = 0;
    return n;
  }
#ifdef MERGESORT
  // if merged != 0, mergesort indices and return true if dupe found
  // if merged == 0, order indices as in Wagner condition
  bool orderindices(u32 *indices, u32 size, u32 *merged) {
    if (merged) {
      u32 i = 0, j = 0, k;
      for (k = 0; i<size && j<size; k++) {
        if (indices[i] == indices[size+j]) return true;
        merged[k] = indices[i] < indices[size+j] ? indices[i++] : indices[size+j++];
      }
      memcpy(merged+k, indices+i, (size-i) * sizeof(u32));
      memcpy(indices,  merged,    (size+j) * sizeof(u32));
      return false;
    } else {
      if (indices[0] > indices[size]) {
        for (u32 i=0; i < size; i++) {
          const u32 tmp = indices[i];
          indices[i] = indices[size+i];
          indices[size+i] = tmp;
        }
      }
      return false;
    }
  }
  // return true if dupe found
  bool listindices0(u32 r, const tree t, u32 *indices, u32 *merged) {
    if (r == 0) {
      *indices = t.getindex();
      return false;
    }
    const slot1 *buck = hta.heap1[t.bucketid()];
    const u32 size = 1 << --r;
    u32 *indices1 = indices + size;
    u32 tagi = hashwords(hashsize(r));
    return listindices1(r, buck[t.slotid0()][tagi].tag, indices,  merged)
        || listindices1(r, buck[t.slotid1()][tagi].tag, indices1, merged)
        || orderindices(indices, size, merged);
  }
  bool listindices1(u32 r, const tree t, u32 *indices, u32 *merged) {
    const slot0 *buck = hta.heap0[t.bucketid()];
    const u32 size = 1 << --r;
    u32 *indices1 = indices + size;
    u32 tagi = hashwords(hashsize(r));
    return listindices0(r, buck[t.slotid0()][tagi].tag, indices,  merged)
        || listindices0(r, buck[t.slotid1()][tagi].tag, indices1, merged)
        || orderindices(indices, size, merged);
  }
  void candidate(const tree t) {
    proof prf, merged;
    if (listindices1(WK, t, prf, merged)) return;
#ifdef ATOMIC
    u32 soli = std::atomic_fetch_add_explicit(&nsols, 1U, std::memory_order_relaxed);
#else
    u32 soli = nsols++;
#endif
    if (soli < MAXSOLS) listindices1(WK, t, sols[soli], 0);
  }
#else
  bool orderindices(u32 *indices, u32 size) {
    if (indices[0] > indices[size]) {
      for (u32 i=0; i < size; i++) {
        const u32 tmp = indices[i];
        indices[i] = indices[size+i];
        indices[size+i] = tmp;
      }
    }
    return false;
  }
  // if dupes != 0, list indices in arbitrary order and return true if dupe found
  // if dupes == 0, order indices as in Wagner condition
  bool listindices0(u32 r, const tree t, u32 *indices, u32 *dupes) {
    if (r == 0) {
      u32 idx = t.getindex();
      if (dupes) {
        u32 bin = idx & (PROOFSIZE-1);
        if (idx == dupes[bin]) return true;
        dupes[bin] = idx;
      }
      *indices = idx;
      return false;
    }
    const slot1 *buck = hta.heap1[t.bucketid()];
    const u32 size = 1 << --r;
    u32 tagi = hashwords(hashsize(r));
    return listindices1(r, buck[t.slotid0()][tagi].tag, indices,      dupes)
        || listindices1(r, buck[t.slotid1()][tagi].tag, indices+size, dupes)
        || (!dupes && orderindices(indices, size));
  }
  bool listindices1(u32 r, const tree t, u32 *indices, u32 *dupes) {
    const slot0 *buck = hta.heap0[t.bucketid()];
    const u32 size = 1 << --r;
    u32 tagi = hashwords(hashsize(r));
    return listindices0(r, buck[t.slotid0()][tagi].tag, indices,      dupes)
        || listindices0(r, buck[t.slotid1()][tagi].tag, indices+size, dupes)
        || (!dupes && orderindices(indices, size));
  }
  void candidate(const tree t) {
    proof prf, dupes;
    memset(dupes, 0xffff, sizeof(proof));
    if (listindices1(WK, t, prf, dupes)) return; // assume WK odd
    qsort(prf, PROOFSIZE, sizeof(u32), &compu32);
    for (u32 i=1; i<PROOFSIZE; i++) if (prf[i] <= prf[i-1]) return;
#ifdef ATOMIC
    u32 soli = std::atomic_fetch_add_explicit(&nsols, 1U, std::memory_order_relaxed);
#else
    u32 soli = nsols++;
#endif
    if (soli < MAXSOLS) listindices1(WK, t, sols[soli], 0); // assume WK odd
  }
#endif
  void showbsizes(u32 r) {
    printf(" b%d h%d\n", bfull, hfull);
    bfull = hfull = 0;
#if defined(HIST) || defined(SPARK) || defined(LOGSPARK)
    u32 binsizes[65];
    memset(binsizes, 0, 65 * sizeof(u32));
    for (u32 bucketid = 0; bucketid < NBUCKETS; bucketid++) {
      u32 bsize = min(nslots[r&1][bucketid], NSLOTS) >> (SLOTBITS-6);
      binsizes[bsize]++;
    }
    for (u32 i=0; i < 65; i++) {
#ifdef HIST
      printf(" %d:%d", i, binsizes[i]);
#else
#ifdef SPARK
      u32 sparks = binsizes[i] / SPARKSCALE;
#else
      u32 sparks = 0;
      for (u32 bs = binsizes[i]; bs; bs >>= 1) sparks++;
      sparks = sparks * 7 / SPARKSCALE;
#endif
      printf("\342\226%c", '\201' + sparks);
#endif
    }
    printf("\n");
#endif
    printf("Digit %d", r+1);
  }

  struct htlayout {
    htalloc hta;
    u32 prevhtunits;
    u32 nexthtunits;
    u32 dunits;
    u32 prevbo;
  
    htlayout(equi *eq, u32 r): hta(eq->hta), prevhtunits(0), dunits(0) {
      u32 nexthashbytes = hashsize(r);
      nexthtunits = hashwords(nexthashbytes);
      prevbo = 0;
      if (r) {
        u32 prevhashbytes = hashsize(r-1);
        prevhtunits = hashwords(prevhashbytes);
        prevbo = prevhtunits * sizeof(htunit) - prevhashbytes; // 0-3
        dunits = prevhtunits - nexthtunits;
      }
    }
    u32 getxhash0(const htunit* slot) const {
#if WN == 200 && RESTBITS == 4
      return slot->bytes[prevbo] >> 4;
#elif WN == 200 && RESTBITS == 8
      return (slot->bytes[prevbo] & 0xf) << 4 | slot->bytes[prevbo+1] >> 4;
#elif WN == 144 && RESTBITS == 4
      return slot->bytes[prevbo] & 0xf;
#else
#error non implemented
#endif
    }
    u32 getxhash1(const htunit* slot) const {
#if WN == 200 && RESTBITS == 4
      return slot->bytes[prevbo] & 0xf;
#elif WN == 200 && RESTBITS == 8
      return slot->bytes[prevbo];
#elif WN == 144 && RESTBITS == 4
      return slot->bytes[prevbo] & 0xf;
#else
#error non implemented
#endif
    }
    bool equal(const htunit *hash0, const htunit *hash1) const {
      return hash0[prevhtunits-1].word == hash1[prevhtunits-1].word;
    }
  };

  struct collisiondata {
#ifdef XBITMAP
#if NSLOTS > 64
#error cant use XBITMAP with more than 64 slots
#endif
    u64 xhashmap[NRESTS];
    u64 xmap;
#else
#if RESTBITS <= 6
    typedef uchar xslot;
#else
    typedef u16 xslot;
#endif
    static const xslot xnil = ~0;
    xslot xhashslots[NRESTS];
    xslot nextxhashslot[NSLOTS];
    xslot nextslot;
#endif
    u32 s0;

    void clear() {
#ifdef XBITMAP
      memset(xhashmap, 0, NRESTS * sizeof(u64));
#else
      memset(xhashslots, xnil, NRESTS * sizeof(xslot));
      memset(nextxhashslot, xnil, NSLOTS * sizeof(xslot));
#endif
    }
    void addslot(u32 s1, u32 xh) {
#ifdef XBITMAP
      xmap = xhashmap[xh];
      xhashmap[xh] |= (u64)1 << s1;
      s0 = -1;
#else
      nextslot = xhashslots[xh];
      nextxhashslot[s1] = nextslot;
      xhashslots[xh] = s1;
#endif
    }
    bool nextcollision() const {
#ifdef XBITMAP
      return xmap != 0;
#else
      return nextslot != xnil;
#endif
    }
    u32 slot() {
#ifdef XBITMAP
      const u32 ffs = __builtin_ffsll(xmap);
      s0 += ffs; xmap >>= ffs;
#else
      nextslot = nextxhashslot[s0 = nextslot];
#endif
      return s0;
    }
  };

  static const u32 BLAKESINPARALLEL = 4;
  // number of hashes extracted from BLAKESINPARALLEL blake2b outputs
  static const u32 HASHESPERBLOCK = BLAKESINPARALLEL*HASHESPERBLAKE;
  // number of blocks of parallel blake2b calls
  static const u32 NBLOCKS = (NHASHES+HASHESPERBLOCK-1)/HASHESPERBLOCK;

  void digit0(const u32 id) {
    htlayout htl(this, 0);
#ifndef HASHONLY
    const u32 hashbytes = hashsize(0);
#endif
    alignas(32) uchar midstate[256], hashes[256];
    //aligned256 midstate, hashes;
    memcpy((void *)midstate, blake_ctx, 256);
    for (u32 block = id; block < NBLOCKS; block += nthreads) {
      Blake2Run4(hashes, (void *)midstate, block * BLAKESINPARALLEL);
#ifndef HASHONLY
      for (u32 i = 0; i<BLAKESINPARALLEL; i++) {
        for (u32 j = 0; j<HASHESPERBLAKE; j++) {
          const uchar *ph = hashes+ i * 64 + j * WN/8;
          const u32 bucketid = ((u32)ph[0] << 4) | ph[1] >> 4;
          const u32 slot = getslot0(bucketid);
          if (slot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *s = hta.heap0[bucketid][slot] + htl.nexthtunits;
          memcpy(s->bytes-hashbytes, ph+WN/8-hashbytes, hashbytes);
          s->tag = tree((block * BLAKESINPARALLEL + i) * HASHESPERBLAKE + j);
        }
      }
#endif
    }
  }
  
  void digitodd(const u32 r, const u32 id) {
    htlayout htl(this, r);
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot0 *buck = htl.hta.heap0[bucketid];
      u32 bsize   = getnslots0(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, htl.getxhash0(slot1));
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (htl.equal(slot0, slot1)) {
            hfull++;
            continue;
          }
          u32 xorbucketid;
          const uchar *bytes0 = slot0->bytes, *bytes1 = slot1->bytes;
#if WN == 200 && BUCKBITS == 12 && RESTBITS == 8
          xorbucketid = (((u32)(bytes0[htl.prevbo+1] ^ bytes1[htl.prevbo+1]) & 0xf) << 8)
                             | (bytes0[htl.prevbo+2] ^ bytes1[htl.prevbo+2]);
#elif WN == 144 && BUCKBITS == 20 && RESTBITS == 4
          xorbucketid = ((((u32)(bytes0[htl.prevbo+1] ^ bytes1[htl.prevbo+1]) << 8)
                              | (bytes0[htl.prevbo+2] ^ bytes1[htl.prevbo+2])) << 4)
                              | (bytes0[htl.prevbo+3] ^ bytes1[htl.prevbo+3]) >> 4;
#elif WN == 96 && BUCKBITS == 12 && RESTBITS == 4
          xorbucketid = ((u32)(bytes0[htl.prevbo+1] ^ bytes1[htl.prevbo+1]) << 4)
                            | (bytes0[htl.prevbo+2] ^ bytes1[htl.prevbo+2]) >> 4;
#else
#error not implemented
#endif
          const u32 xorslot = getslot1(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *xs = htl.hta.heap1[xorbucketid][xorslot];
          for (u32 i=htl.dunits; i < htl.prevhtunits; i++)
            xs++->word = slot0[i].word ^ slot1[i].word;
          xs->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  
  void digiteven(const u32 r, const u32 id) {
    htlayout htl(this, r);
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot1 *buck = htl.hta.heap1[bucketid];
      u32 bsize   = getnslots1(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, htl.getxhash1(slot1));
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (htl.equal(slot0, slot1)) {
            hfull++;
            continue;
          }
          u32 xorbucketid;
          const uchar *bytes0 = slot0->bytes, *bytes1 = slot1->bytes;
#if WN == 200 && BUCKBITS == 12 && RESTBITS == 8
          xorbucketid = ((u32)(bytes0[htl.prevbo+1] ^ bytes1[htl.prevbo+1]) << 4)
                            | (bytes0[htl.prevbo+2] ^ bytes1[htl.prevbo+2]) >> 4;
#elif WN == 144 && BUCKBITS == 20 && RESTBITS == 4
          xorbucketid = ((((u32)(bytes0[htl.prevbo+1] ^ bytes1[htl.prevbo+1]) << 8)
                              | (bytes0[htl.prevbo+2] ^ bytes1[htl.prevbo+2])) << 4)
                              | (bytes0[htl.prevbo+3] ^ bytes1[htl.prevbo+3]) >> 4;
#elif WN == 96 && BUCKBITS == 12 && RESTBITS == 4
          xorbucketid = ((u32)(bytes0[htl.prevbo+1] ^ bytes1[htl.prevbo+1]) << 4)
                            | (bytes0[htl.prevbo+2] ^ bytes1[htl.prevbo+2]) >> 4;
#else
#error not implemented
#endif
          const u32 xorslot = getslot0(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *xs = htl.hta.heap0[xorbucketid][xorslot];
          for (u32 i=htl.dunits; i < htl.prevhtunits; i++)
            xs++->word = slot0[i].word ^ slot1[i].word;
          xs->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  
  void digit1(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot0 *buck = heaps.heap0[bucketid];
      u32 bsize   = getnslots0(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, htobe32(slot1->word) >> 20 & 0xff);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[5].word == slot1[5].word) {
            hfull++;
            continue;
          }
          u32 xorbucketid = htobe32(slot0->word ^ slot1->word) >> 8 & BUCKMASK;
          const u32 xorslot = getslot1(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          u64 *x  = (u64 *)heaps.heap1[xorbucketid][xorslot];
          u64 *x0 = (u64 *)slot0, *x1 = (u64 *)slot1;
          *x++ = x0[0] ^ x1[0];
          *x++ = x0[1] ^ x1[1];
          *x++ = x0[2] ^ x1[2];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit2(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot1 *buck = heaps.heap1[bucketid];
      u32 bsize   = getnslots1(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, slot1->bytes[3]);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[5].word == slot1[5].word) {
            hfull++;
            continue;
          }
          u32 xorbucketid = htobe32(slot0[1].word ^ slot1[1].word) >> 20;
          const u32 xorslot = getslot0(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *xs = heaps.heap0[xorbucketid][xorslot];
          xs++->word = slot0[1].word ^ slot1[1].word;
          u64 *x = (u64 *)xs, *x0 = (u64 *)slot0, *x1 = (u64 *)slot1;
          *x++ = x0[1] ^ x1[1];
          *x++ = x0[2] ^ x1[2];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit3(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot0 *buck = heaps.heap0[bucketid];
      u32 bsize   = getnslots0(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, htobe32(slot1->word) >> 12 & 0xff);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[4].word == slot1[4].word) {
            hfull++;
            continue;
          }
          u32 xorbucketid = htobe32(slot0[0].word ^ slot1[0].word) & BUCKMASK;
          const u32 xorslot = getslot1(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          u64 *x  = (u64 *)heaps.heap1[xorbucketid][xorslot];
          u64 *x0 = (u64 *)(slot0+1), *x1 = (u64 *)(slot1+1);
          *x++ = x0[0] ^ x1[0];
          *x++ = x0[1] ^ x1[1];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit4(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot1 *buck = heaps.heap1[bucketid];
      u32 bsize   = getnslots1(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, slot1->bytes[0]);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[3].word == slot1[3].word) {
            hfull++;
            continue;
          }
          u32 xorbucketid = htobe32(slot0[0].word ^ slot1[0].word) >> 12 & BUCKMASK;
          const u32 xorslot = getslot0(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          u64 *x  = (u64 *)heaps.heap0[xorbucketid][xorslot];
          u64 *x0 = (u64 *)slot0, *x1 = (u64 *)slot1;
          *x++ = x0[0] ^ x1[0];
          *x++ = x0[1] ^ x1[1];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit5(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot0 *buck = heaps.heap0[bucketid];
      u32 bsize   = getnslots0(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, htobe32(slot1->word) >> 4 & 0xff);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[3].word == slot1[3].word) {
            hfull++;
            continue;
          }
          u32 xor1 = slot0[1].word ^ slot1[1].word;
          u32 xorbucketid = (((u32)(slot0->bytes[3] ^ slot1->bytes[3]) & 0xf)
                               << 8) | (xor1 & 0xff);
          const u32 xorslot = getslot1(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *xs = heaps.heap1[xorbucketid][xorslot];
          xs++->word = xor1;
          u64 *x = (u64 *)xs, *x0 = (u64 *)slot0, *x1 = (u64 *)slot1;
          *x++ = x0[1] ^ x1[1];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit6(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot1 *buck = heaps.heap1[bucketid];
      u32 bsize   = getnslots1(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, slot1->bytes[1]);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[2].word == slot1[2].word) {
            hfull++;
            continue;
          }
          u32 xorbucketid = htobe32(slot0[0].word ^ slot1[0].word) >> 4 & BUCKMASK;
          const u32 xorslot = getslot0(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *xs = heaps.heap0[xorbucketid][xorslot];
          xs++->word = slot0[0].word ^ slot1[0].word;
          u64 *x = (u64 *)xs, *x0 = (u64 *)(slot0+1), *x1 = (u64 *)(slot1+1);
          *x++ = x0[0] ^ x1[0];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit7(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot0 *buck = heaps.heap0[bucketid];
      u32 bsize   = getnslots0(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, (slot1->bytes[3] & 0xf) << 4 | slot1->bytes[4] >> 4);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          if (slot0[2].word == slot1[2].word) {
            hfull++;
            continue;
          }
          u32 xorbucketid = htobe32(slot0[1].word ^ slot1[1].word) >> 16 & BUCKMASK;
          const u32 xorslot = getslot1(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          u64 *x  = (u64 *)heaps.heap1[xorbucketid][xorslot];
          u64 *x0 = (u64 *)(slot0+1), *x1 = (u64 *)(slot1+1);
          *x++ = x0[0] ^ x1[0];
          ((htunit *)x)->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  void digit8(const u32 id) {
    htalloc heaps = hta;
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot1 *buck = heaps.heap1[bucketid];
      u32 bsize   = getnslots1(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, slot1->bytes[2]);
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *slot0 = buck[s0];
          u32 xor1 = slot0[1].word ^ slot1[1].word;
          if (!xor1) {
            hfull++;
            continue;
          }
          u32 xorbucketid = ((u32)(slot0->bytes[3] ^ slot1->bytes[3]) << 4)
                          | (xor1 >> 4 & 0xf);
          const u32 xorslot = getslot0(xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          htunit *xs = heaps.heap0[xorbucketid][xorslot];
          xs++->word = xor1;
          xs->tag = tree(bucketid, s0, s1);
        }
      }
    }
  }
  
  void digitK(const u32 id) {
    collisiondata cd;
    htlayout htl(this, WK);
    u32 nc = 0;
    for (u32 bucketid = id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      slot0 *buck = htl.hta.heap0[bucketid];
      u32 bsize   = getnslots0(bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *slot1 = buck[s1];
        cd.addslot(s1, htl.getxhash0(slot1)); // assume WK odd
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          if (htl.equal(buck[s0], slot1)) { // EASY OPTIMIZE
            candidate(tree(bucketid, s0, s1));
            nc++;
          }
        }
      }
    }
    // printf(" %d candidates ", nc);
  }
};

typedef struct {
  u32 id;
  pthread_t thread;
  equi *eq;
} thread_ctx;

void barrier(pthread_barrier_t *barry) {
  const int rc = pthread_barrier_wait(barry);
  if (rc != 0 && rc != PTHREAD_BARRIER_SERIAL_THREAD) {
    printf("Could not wait on barrier\n");
    pthread_exit(NULL);
  }
}

void *worker(void *vp) {
  thread_ctx *tp = (thread_ctx *)vp;
  equi *eq = tp->eq;

  if (tp->id == 0) printf("Digit 0");
  eq->digit0(tp->id);
#ifdef HASHONLY
  pthread_exit(NULL);
#endif
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(0);
  barrier(&eq->barry);
#if WN == 200 && WK == 9 && RESTBITS == 8
  eq->digit1(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(1);
  barrier(&eq->barry);
  eq->digit2(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(2);
  barrier(&eq->barry);
  eq->digit3(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(3);
  barrier(&eq->barry);
  eq->digit4(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(4);
  barrier(&eq->barry);
  eq->digit5(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(5);
  barrier(&eq->barry);
  eq->digit6(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(6);
  barrier(&eq->barry);
  eq->digit7(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(7);
  barrier(&eq->barry);
  eq->digit8(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) eq->showbsizes(8);
  barrier(&eq->barry);
#else
  for (u32 r = 1; r < WK; r++) {
    r&1 ? eq->digitodd(r, tp->id) : eq->digiteven(r, tp->id);
    barrier(&eq->barry);
    if (tp->id == 0) eq->showbsizes(r);
    barrier(&eq->barry);
  }
#endif
  eq->digitK(tp->id);
  pthread_exit(NULL);
  return 0;
}
