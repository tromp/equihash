// Equihash solver
// Copyright (c) 2016 John Tromp

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

// The algorithm below solves this by maintaining the trees
// in a graph of K layers, each split into buckets
// with buckets indexed by the first n-RESTBITS bits following
// the i*n 0s, each bucket having 4 * 2^RESTBITS slots,
// twice the number of subtrees expected to land there.

#include "equi.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <assert.h>

typedef uint64_t u64;

#ifdef ATOMIC
#include <atomic>
typedef std::atomic<u32> au32;
#else
typedef u32 au32;
#endif

#ifndef RESTBITS
#define RESTBITS	4
#endif

// 2_log of number of buckets
#define BUCKBITS (DIGITBITS-RESTBITS)

// number of buckets
static const u32 NBUCKETS = 1<<BUCKBITS;
// 2_log of number of slots per bucket
static const u32 SLOTBITS = RESTBITS+1+1;
// number of slots per bucket
static const u32 NSLOTS = 1<<SLOTBITS;
// number of per-xhash slots
static const u32 XFULL = NSLOTS/4;
// SLOTBITS mask
static const u32 SLOTMASK = NSLOTS-1;
// number of possible values of xhash (rest of n) bits
static const u32 NRESTS = 1<<RESTBITS;
// number of blocks of hashes extracted from single 512 bit blake2b output
static const u32 NBLOCKS = (NHASHES+HASHESPERBLAKE-1)/HASHESPERBLAKE;
// nothing larger found in 100000 runs
static const u32 MAXSOLS = 8;

// scaling factor for showing bucketsize histogra as sparkline
#ifndef SPARKSCALE
#define SPARKSCALE	(40 << (BUCKBITS-12))
#endif

// tree node identifying its children as two different slots in
// a bucket on previous layer with the same rest bits (x-tra hash)
struct tree {
  unsigned bucketid : BUCKBITS;
  unsigned slotid0  : SLOTBITS;
  unsigned slotid1  : SLOTBITS;
#ifndef XWITHASH
  unsigned xhash    : RESTBITS;
#endif

// layer 0 has no children bit needs to encode index
  u32 getindex() const {
    return (bucketid << SLOTBITS) | slotid0;
  }
  void setindex(const u32 idx) {
    slotid0 = idx & SLOTMASK;
    bucketid = idx >> SLOTBITS;
  }
};

union htunit {
  tree attr;
  u32 hash;
  uchar bytes[sizeof(u32)];
};

// a bucket is NSLOTS treenodes
typedef htunit bucket[NSLOTS];
// the N-bit hash consists of K+1 n-bit "digits"
// each of which corresponds to a layer of NBUCKETS buckets

// size (in bytes) of hash in round 0 <= r < WK
u32 hashsize(const u32 r) {
#ifdef XWITHASH
  const u32 hashbits = WN - (r+1) * DIGITBITS + RESTBITS;
#else
  const u32 hashbits = WN - (r+1) * DIGITBITS;
#endif
  return (hashbits + 7) / 8;
}


u32 htunits(u32 bytes) {
  return (bytes + sizeof(htunit) - 1) / sizeof(htunit);
}

#ifdef JOINHT
u32 slotsize(const u32 r) {
  return 1 + htunits(hashsize(r));
}
// size (in htunits) of bucket in round 0 <= r < WK
u32 bucketsize(const u32 r) {
  return NSLOTS * slotsize(r);
}
#else
u32 slotsize(const u32 r) {
  return 1;
}
#endif

// manages hash and tree data
struct htalloc {
#ifdef JOINHT
  htunit *trees[WK];
#else
  bucket *trees[WK];
  htunit *hashes[WK];
#endif
  u64 alloced;
  htalloc() {
    alloced = 0;
  }
  void alloctrees() {
#ifdef JOINHT
    for (int r=0; r<WK; r++)
      trees[r] = (htunit *)alloc(NBUCKETS * NSLOTS * (1 + htunits(hashsize(r))), sizeof(htunit));
#endif
  }
  void dealloctrees() {
#ifdef JOINHT
    for (int r=0; r<WK; r++)
      dealloc(trees[r], NBUCKETS * NSLOTS * (1 + htunits(hashsize(r))), sizeof(htunit));
#else
    for (int r=0; r<WK; r++)
      dealloc(trees[r], NBUCKETS, sizeof(bucket));
#endif
  }
  void alloc_ht(u32 r) {
#ifndef JOINHT
    hashes[r] = (htunit *)alloc(NBUCKETS * NSLOTS * htunits(hashsize(r)), sizeof(htunit));
    trees[r] = (bucket *)alloc(NBUCKETS, sizeof(bucket));
#endif
  }
  void dealloc_ht(u32 r) {
#ifndef JOINHT
    dealloc(hashes[r], NBUCKETS * NSLOTS * htunits(hashsize(r)), sizeof(htunit));
#endif
  }
  htunit *getbucket(u32 r, u32 bid) const {
#ifdef JOINHT
    return &trees[r][bid * bucketsize(r)];
#else
    return trees[r][bid];
#endif
  }
  void *alloc(const u32 n, const u32 sz) {
    void *mem  = calloc(n, sz);
    assert(mem);
    alloced += (u64)n * sz;
    return mem;
  }
  void dealloc(void *mem, const u32 n, const u32 sz) {
    free(mem);
    alloced -= (u64)n * sz;
  }
};

typedef au32 bsizes[NBUCKETS];

u32 min(const u32 a, const u32 b) {
  return a < b ? a : b;
}

struct equi {
  blake2b_state blake_ctx;
  htalloc hta;
  bsizes *nslots;
  proof *sols;
  au32 nsols;
  u32 nthreads;
  u32 xfull;
  u32 hfull;
  u32 bfull;
  pthread_barrier_t barry;
  equi(const u32 n_threads) {
    nthreads = n_threads;
    const int err = pthread_barrier_init(&barry, NULL, nthreads);
    assert(!err);
    hta.alloctrees();
    nslots    = (bsizes *)hta.alloc(2 * NBUCKETS, sizeof(au32));
    sols = (proof *)hta.alloc(MAXSOLS, sizeof(proof));
  }
  ~equi() {
    free(nslots);
    free(sols);
#ifdef JOINHT
    hta.dealloctrees();
#endif
  }
  void setnonce(const char *header, u32 nonce) {
    setheader(&blake_ctx, header, nonce);
    memset(nslots, 0, NBUCKETS * sizeof(au32)); // only nslots[0] needs zeroing
    nsols = 0;
  }
  u32 findslot(const u32 r, const u32 bucketi) {
#ifdef ATOMIC
    return std::atomic_fetch_add_explicit(&nslots[r&1][bucketi], 1U, std::memory_order_relaxed);
#else
    return nslots[r&1][bucketi]++;
#endif
  }
  u32 getnslots(const u32 r, const u32 bid) {
    au32 &nslot = nslots[r&1][bid];
    const u32 n = min(nslot, NSLOTS);
    nslot = 0;
    return n;
  }
  void listindices(u32 r, const tree t, u32 *indices) {
    if (r == 0) {
      *indices = t.getindex();
      return;
    }
    const htunit *bt = hta.getbucket(--r,t.bucketid);
    const u32 size = 1 << r;
    u32 *indices1 = indices + size;
    listindices(r, bt[t.slotid0 * slotsize(r)].attr, indices);
    listindices(r, bt[t.slotid1 * slotsize(r)].attr, indices1);
    if (*indices > *indices1) {
      for (u32 i=0; i < size; i++) {
        const u32 tmp = indices[i];
        indices[i] = indices1[i];
        indices1[i] = tmp;
      }
    }
  }
  void candidate(const tree t) {
    proof prf;
    listindices(WK, t, prf);
    qsort(prf, PROOFSIZE, sizeof(u32), &compu32);
    for (u32 i=1; i<PROOFSIZE; i++)
      if (prf[i] <= prf[i-1])
        return;
#ifdef ATOMIC
    u32 soli = std::atomic_fetch_add_explicit(&nsols, 1U, std::memory_order_relaxed);
#else
    u32 soli = nsols++;
#endif
    if (soli < MAXSOLS)
      listindices(WK, t, sols[soli]);
  }
  void showbsizes(u32 r) {
#if defined(HIST) || defined(SPARK)
    u32 bsizes[NSLOTS+1];
    memset(bsizes, 0, (NSLOTS+1) * sizeof(u32));
    for (u32 bucketid = 0; bucketid < NBUCKETS; bucketid++) {
      u32 bsize = nslots[r&1][bucketid];
      if (bsize < NSLOTS)
        bsizes[bsize]++;
      else
        bsizes[NSLOTS]++;
    }
    for (u32 i=0; i<=NSLOTS; i++) {
#ifdef HIST
      printf(" %d:%d", i, bsizes[i]);
#else
      printf("\342\226%c", '\201'+bsizes[i]/SPARKSCALE);
#endif
    }
    printf(" %ld MB\n", hta.alloced >> 20);
#endif
  }

  struct htlayout {
    htalloc hta;
    u32 prevhtunits;
    u32 nexthtunits;
    u32 dunits;
    u32 prevbo;
    u32 nextbo;
    htunit *buck;
    htunit *hashbase;
  
    htlayout(equi *eq, u32 r): hta(eq->hta), prevhtunits(0), dunits(0) {
      u32 nexthashbytes = hashsize(r);
      nexthtunits = htunits(nexthashbytes);
      prevbo = 0;
      nextbo = nexthtunits * sizeof(htunit) - nexthashbytes; // 0-3
      if (r) {
        u32 prevhashbytes = hashsize(r-1);
        prevhtunits = htunits(prevhashbytes);
        prevbo = prevhtunits * sizeof(htunit) - prevhashbytes; // 0-3
        dunits = prevhtunits - nexthtunits;
      }
#ifdef JOINHT
      nexthtunits++;
      prevhtunits++;
#endif
    }
    void setbucket(u32 r, u32 bid) {
      buck = hta.getbucket(r, bid);
#ifdef JOINHT
      hashbase = buck + 1;
#else
      hashbase = hta.hashes[r] + (bid * NSLOTS) * prevhtunits;
#endif
    }
    u32 getxhash(const u32 slot, const htunit *hash) const {
#ifdef XWITHASH
      return hash->bytes[prevbo] & 0xf;
#elif defined JOINHT
      return buck[slot * prevhtunits].attr.xhash;
#else
      return buck[slot].attr.xhash;
#endif
    }
    u32 prevhashunits() const {
#ifdef JOINHT
      return prevhtunits - 1;
#else
      return prevhtunits;
#endif
    }
    bool equal(const htunit *hash0, const htunit *hash1) const {
      return hash0[prevhashunits()-1].hash == hash1[prevhashunits()-1].hash;
    }
    htunit *addtree(u32 r, tree t, u32 bid, u32 slot) {
      htunit *buck = hta.getbucket(r,bid);
#ifdef JOINHT
      htunit *slotree = buck + slot * nexthtunits;
      slotree->attr = t;
      return slotree + 1;
#else
      buck[slot].attr = t;
      return hta.hashes[r] + (bid * NSLOTS + slot) * nexthtunits;
#endif
    }
  };

  struct collisiondata {
#ifdef XBITMAP
    u64 xhashmap[NRESTS];
    u64 xmap;
#else
    typedef uchar xslot;
    xslot nxhashslots[NRESTS];
    xslot xhashslots[NRESTS][XFULL];
    xslot *xx;
    u32 n0;
    u32 n1;
#endif
    u32 s0;

    void clear() {
#ifdef XBITMAP
      memset(xhashmap, 0, NRESTS * sizeof(u64));
#else
      memset(nxhashslots, 0, NRESTS * sizeof(xslot));
#endif
    }
    bool addslot(u32 s1, u32 xh) {
#ifdef XBITMAP
      xmap = xhashmap[xh];
      xhashmap[xh] |= (u64)1 << s1;
      s0 = -1;
      return true;
#else
      n1 = (u32)nxhashslots[xh]++;
      if (n1 >= XFULL)
        return false;
      xx = xhashslots[xh];
      xx[n1] = s1;
      n0 = 0;
      return true;
#endif
    }
    bool nextcollision() const {
#ifdef XBITMAP
      return xmap != 0;
#else
      return n0 < n1;
#endif
    }
    u32 slot() {
#ifdef XBITMAP
      const u32 ffs = __builtin_ffsll(xmap);
      s0 += ffs; xmap >>= ffs;
      return s0;
#else
      return (u32)xx[n0++];
#endif
    }
  };

  void digit0(const u32 id) {
    uchar hash[HASHOUT];
    blake2b_state state;
    htlayout htl(this, 0);
    const u32 hashbytes = hashsize(0);
    for (u32 block = id; block < NBLOCKS; block += nthreads) {
      state = blake_ctx;
      const u32 leb = htole32(block);
      blake2b_update(&state, (uchar *)&leb, sizeof(u32));
      blake2b_final(&state, hash, HASHOUT);
      for (u32 i = 0; i<HASHESPERBLAKE; i++) {
        const uchar *ph = hash + i * WN/8;
#if BUCKBITS == 16 && RESTBITS == 4
        const u32 bucketid = ((u32)ph[0] << 8) | ph[1];
        const u32 xhash = ph[2] >> 4;
#elif BUCKBITS == 20 && RESTBITS == 4
        const u32 bucketid = ((((u32)ph[0] << 8) | ph[1]) << 4) | ph[2] >> 4;
#ifndef XWITHASH
        const u32 xhash = ph[2] & 0xf;
#endif
#elif BUCKBITS == 12 && RESTBITS == 4
        const u32 bucketid = ((u32)ph[0] << 4) | ph[1] >> 4;
        const u32 xhash = ph[1] & 0xf;
#else
#error not implemented
#endif
        const u32 slot = findslot(0, bucketid);
        if (slot >= NSLOTS) {
          bfull++;
          continue;
        }
        tree leaf;
        leaf.setindex(block*HASHESPERBLAKE+i);
#ifndef XWITHASH
        leaf.xhash = xhash;
#endif
        htunit *dest = htl.addtree(0, leaf, bucketid, slot);
        memcpy(dest->bytes+htl.nextbo, ph+WN/8-hashbytes, hashbytes);
      }
    }
  }
  
  void digitr(const u32 r, const u32 id) {
    htlayout htl(this, r);
    collisiondata cd;
    for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      htl.setbucket(r-1, bucketid);
      u32 bsize = getnslots(r-1, bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
        if (!cd.addslot(s1, htl.getxhash(s1, hash1))) {
          xfull++;
          continue;
        }
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
          if (htl.equal(hash0, hash1)) {
            hfull++;
            continue;
          }
          u32 xorbucketid;
          u32 xhash;
#if BUCKBITS == 16 && RESTBITS == 4
          xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
            | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
          xhash = hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2];
          if (r&1) {
            xorbucketid = ((xorbucketid & 0xfff) << 4) | (xhash >> 4);
            xhash &= 0xf;
          } else xhash >>= 4;
#elif BUCKBITS == 20 && RESTBITS == 4 && defined XWITHASH
          xhash = hash0->bytes[htl.prevbo+3] ^ hash1->bytes[htl.prevbo+3];
          xorbucketid = ((((u32)(hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]) << 8) | (hash0->bytes[htl.prevbo+2]^hash1->bytes[htl.prevbo+2])) << 4) | xhash >> 4;
          xhash &= 0xf;
#elif BUCKBITS == 12 && RESTBITS == 4
          xhash = hash0->bytes[htl.prevbo+1] ^ hash1->bytes[htl.prevbo+1];
          xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 4) | xhash >> 4;
          xhash &= 0xf;
#else
#error not implemented
#endif
          const u32 xorslot = findslot(r, xorbucketid);
          if (xorslot >= NSLOTS) {
            bfull++;
            continue;
          }
          tree xort; xort.bucketid = bucketid;
          xort.slotid0 = s0; xort.slotid1 = s1;
#ifndef XWITHASH
          xort.xhash = xhash;
#endif
          htunit *xorhash = htl.addtree(r, xort, xorbucketid, xorslot);
          for (u32 i=htl.dunits; i < htl.prevhashunits(); i++)
            xorhash[i-htl.dunits].hash = hash0[i].hash ^ hash1[i].hash;
        }
      }
    }
  }
  
  void digitK(const u32 id) {
    collisiondata cd;
    htlayout htl(this, WK);
    for (u32 bucketid = id; bucketid < NBUCKETS; bucketid += nthreads) {
      cd.clear();
      htl.setbucket(WK-1, bucketid);
      u32 bsize = getnslots(WK-1, bucketid);
      for (u32 s1 = 0; s1 < bsize; s1++) {
        const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
        if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
          continue;
        for (; cd.nextcollision(); ) {
          const u32 s0 = cd.slot();
          const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
          if (htl.equal(hash0, hash1)) {
            tree xort; xort.bucketid = bucketid;
            xort.slotid0 = s0; xort.slotid1 = s1;
            candidate(xort);
          }
        }
      }
    }
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

  if (tp->id == 0) {
    printf("Digit 0\n");
    eq->hta.alloc_ht(0);
  }
  barrier(&eq->barry);
  eq->digit0(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) {
    eq->xfull = eq->bfull = eq->hfull = 0;
    eq->showbsizes(0);
  }
  barrier(&eq->barry);
  for (u32 r = 1; r < WK; r++) {
    if (tp->id == 0) {
      printf("Digit %d", r);
      eq->hta.alloc_ht(r);
    }
    barrier(&eq->barry);
    eq->digitr(r, tp->id);
    barrier(&eq->barry);
    if (tp->id == 0) {
      printf(" x%d b%d h%d\n", eq->xfull, eq->bfull, eq->hfull);
      eq->xfull = eq->bfull = eq->hfull = 0;
      eq->showbsizes(r);
      eq->hta.dealloc_ht(r-1);
    }
    barrier(&eq->barry);
  }
  if (tp->id == 0)
    printf("Digit %d\n", WK);
  eq->digitK(tp->id);
  barrier(&eq->barry);
  if (tp->id == 0) {
    eq->hta.dealloc_ht(WK-1);
#ifndef JOINHT
    eq->hta.dealloctrees();
#endif
  }
  pthread_exit(NULL);
  return 0;
}
