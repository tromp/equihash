// Equihash CUDA solver
// Copyright (c) 2016 John Tromp

#include "equi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "blake2b.cu"

typedef uint64_t u64;

#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
  if (code != cudaSuccess) {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

#ifndef RESTBITS
#define RESTBITS	4
#endif

// 2_log of number of buckets
#define BUCKBITS	(DIGITBITS-RESTBITS)
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
  __device__ u32 getindex() const {
    return (bucketid << SLOTBITS) | slotid0;
  }
  __device__ void setindex(const u32 idx) {
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
typedef bucket digit[NBUCKETS];

// size (in bytes) of hash in round 0 <= r < WK
u32 hhashsize(const u32 r) {
#ifdef XWITHASH
  const u32 hashbits = WN - (r+1) * DIGITBITS + RESTBITS;
#else
  const u32 hashbits = WN - (r+1) * DIGITBITS;
#endif
  return (hashbits + 7) / 8;
}
// size (in bytes) of hash in round 0 <= r < WK
__device__ u32 hashsize(const u32 r) {
#ifdef XWITHASH
  const u32 hashbits = WN - (r+1) * DIGITBITS + RESTBITS;
#else
  const u32 hashbits = WN - (r+1) * DIGITBITS;
#endif
  return (hashbits + 7) / 8;
}

u32 hhtunits(u32 bytes) {
  return (bytes + sizeof(htunit) - 1) / sizeof(htunit);
}

__device__ u32 htunits(u32 bytes) {
  return (bytes + sizeof(htunit) - 1) / sizeof(htunit);
}

#ifdef JOINHT
__device__ u32 slotsize(const u32 r) {
  return 1 + htunits(hashsize(r));
}
// size (in htunits) of bucket in round 0 <= r < WK
__device__ u32 bucketsize(const u32 r) {
  return NSLOTS * slotsize(r);
}
#else
__device__ u32 slotsize(const u32 r) {
  return 1;
}
#endif

// manages hash and tree data
struct htalloc {
#ifdef JOINHT
  htunit *trees[WK];
#else
  digit *trees;
  htunit *hashes[2];
#endif
  __device__ htunit *getbucket(u32 r, u32 bid) const {
#ifdef JOINHT
    return &trees[r][bid * bucketsize(r)];
#else
    return trees[r][bid];
#endif
  }
};

typedef u32 bsizes[NBUCKETS];

//u32 __device__ min(const u32 a, const u32 b) {
//  return a < b ? a : b;
//}

struct equi {
  blake2b_state blake_ctx;
  htalloc hta;
  bsizes *nslots;
  proof *sols;
  u32 nsols;
  u32 nthreads;
  equi(const u32 n_threads) {
    nthreads = n_threads;
  }
  void setnonce(const char *header, u32 nonce) {
    setheader(&blake_ctx, header, nonce);
    checkCudaErrors(cudaMemset(nslots, 0, NBUCKETS * sizeof(u32)));
    nsols = 0;
  }
  __device__ u32 getnslots(const u32 r, const u32 bid) {
    u32 &nslot = nslots[r&1][bid];
    const u32 n = min(nslot, NSLOTS);
    nslot = 0;
    return n;
  }
  __device__ void orderindices(u32 *indices, u32 size) {
    if (indices[0] > indices[size]) {
      for (u32 i=0; i < size; i++) {
        const u32 tmp = indices[i];
        indices[i] = indices[size+i];
        indices[size+i] = tmp;
      }
    }
  }
  __device__ void listindices(u32 r, const tree t, u32 *indices) {
    if (r == 0) {
      *indices = t.getindex();
      return;
    }
    const htunit *bt = hta.getbucket(--r,t.bucketid);
    const u32 size = 1 << r;
    u32 *indices1 = indices + size;
    listindices(r, bt[t.slotid0 * slotsize(r)].attr, indices);
    listindices(r, bt[t.slotid1 * slotsize(r)].attr, indices1);
    orderindices(indices, size);
  }
  __device__ void listindices1(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(0,t.bucketid);
    const u32 size = 1 << 0;
    indices[0]    = bt[t.slotid0 * slotsize(0)].attr.getindex();
    indices[size] = bt[t.slotid1 * slotsize(0)].attr.getindex();
    orderindices(indices, size);
  }
  __device__ void listindices2(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(1,t.bucketid);
    const u32 size = 1 << 1;
    listindices1(bt[t.slotid0 * slotsize(1)].attr, indices);
    listindices1(bt[t.slotid1 * slotsize(1)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices3(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(2,t.bucketid);
    const u32 size = 1 << 2;
    listindices2(bt[t.slotid0 * slotsize(2)].attr, indices);
    listindices2(bt[t.slotid1 * slotsize(2)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices4(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(3,t.bucketid);
    const u32 size = 1 << 3;
    listindices3(bt[t.slotid0 * slotsize(3)].attr, indices);
    listindices3(bt[t.slotid1 * slotsize(3)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices5(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(4,t.bucketid);
    const u32 size = 1 << 4;
    listindices4(bt[t.slotid0 * slotsize(4)].attr, indices);
    listindices4(bt[t.slotid1 * slotsize(4)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices6(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(5,t.bucketid);
    const u32 size = 1 << 5;
    listindices5(bt[t.slotid0 * slotsize(5)].attr, indices);
    listindices5(bt[t.slotid1 * slotsize(5)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices7(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(6,t.bucketid);
    const u32 size = 1 << 6;
    listindices6(bt[t.slotid0 * slotsize(6)].attr, indices);
    listindices6(bt[t.slotid1 * slotsize(6)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices8(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(7,t.bucketid);
    const u32 size = 1 << 7;
    listindices7(bt[t.slotid0 * slotsize(7)].attr, indices);
    listindices7(bt[t.slotid1 * slotsize(7)].attr, indices+size);
    orderindices(indices, size);
  }
  __device__ void listindices9(const tree t, u32 *indices) {
    const htunit *bt = hta.getbucket(8,t.bucketid);
    const u32 size = 1 << 8;
    listindices8(bt[t.slotid0 * slotsize(8)].attr, indices);
    listindices8(bt[t.slotid1 * slotsize(8)].attr, indices+size);
    orderindices(indices, size);
  }
  void showbsizes(u32 r) {
#if defined(HIST) || defined(SPARK)
    u32 ns[NBUCKETS];
    checkCudaErrors(cudaMemcpy(ns, nslots[r&1], NBUCKETS * sizeof(u32), cudaMemcpyDeviceToHost));
    u32 bsizes[NSLOTS+1];
    memset(bsizes, 0, (NSLOTS+1) * sizeof(u32));
    for (u32 bucketid = 0; bucketid < NBUCKETS; bucketid++) {
      u32 bsize = ns[bucketid];
      if (bsize > NSLOTS)
        bsize = NSLOTS;
      bsizes[bsize]++;
    }
    for (u32 i=0; i<=NSLOTS; i++) {
#ifdef HIST
      printf(" %d:%d", i, bsizes[i]);
#else
      printf("\342\226%c", (uchar)'\201'+bsizes[i]/SPARKSCALE);
#endif
    }
    printf("\n");
#endif
  }
  // proper dupe test is a little costly on GPU, so allow false negatives
  __device__ bool probdupe(u32 *prf) {
    unsigned short susp[PROOFSIZE];
    memset(susp, 0xffff, PROOFSIZE * sizeof(unsigned short));
    for (u32 i=0; i<PROOFSIZE; i++) {
      u32 bin = prf[i] & (PROOFSIZE-1);
      unsigned short msb = prf[i]>>WK;
      if (msb == susp[bin])
        return true;
      susp[bin] = msb;
    }
    return false;
  }
  __device__ void candidate(const tree t) {
    proof prf;
#if WK==9
    listindices9(t, prf);
#elif WK==5
    listindices5(t, prf);
#else
#error not implemented
#endif
    if (probdupe(prf))
      return;
    u32 soli = atomicAdd(&nsols, 1);
    if (soli < MAXSOLS)
#if WK==9
      listindices9(t, sols[soli]);
#elif WK==5
      listindices5(t, sols[soli]);
#else
#error not implemented
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
  
    __device__ htlayout(equi *eq, u32 r): hta(eq->hta), prevhtunits(0), dunits(0) {
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
    __device__ void setbucket(u32 r, u32 bid) {
      buck = hta.getbucket(r, bid);
#ifdef JOINHT
      hashbase = buck + 1;
#else
      hashbase = hta.hashes[r&1] + (bid * NSLOTS) * prevhtunits;
#endif
    }
    __device__ u32 getxhash(const u32 slot, const htunit *hash) const {
#ifdef XWITHASH
      return hash->bytes[prevbo] & 0xf;
#elif defined JOINHT
      return buck[slot * prevhtunits].attr.xhash;
#else
      return buck[slot].attr.xhash;
#endif
    }
    __device__ u32 prevhashunits() const {
#ifdef JOINHT
      return prevhtunits - 1;
#else
      return prevhtunits;
#endif
    }
    __device__ bool equal(const htunit *hash0, const htunit *hash1) const {
      return hash0[prevhashunits()-1].hash == hash1[prevhashunits()-1].hash;
    }
    __device__ htunit *addtree(u32 r, tree t, u32 bid, u32 slot) {
      htunit *buck = hta.getbucket(r,bid);
#ifdef JOINHT
      htunit *slotree = buck + slot * nexthtunits;
      slotree->attr = t;
      return slotree + 1;
#else
      buck[slot].attr = t;
      return hta.hashes[r&1] + (bid * NSLOTS + slot) * nexthtunits;
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

    __device__ void clear() {
#ifdef XBITMAP
      memset(xhashmap, 0, NRESTS * sizeof(u64));
#else
      memset(nxhashslots, 0, NRESTS * sizeof(xslot));
#endif
    }
    __device__ bool addslot(u32 s1, u32 xh) {
#ifdef XBITMAP
      xmap = xhashmap[xh];
      xhashmap[xh] |= (u64)1 << s1;
      s0 = ~0;
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
    __device__ bool nextcollision() const {
#ifdef XBITMAP
      return xmap != 0;
#else
      return n0 < n1;
#endif
    }
    __device__ u32 slot() {
#ifdef XBITMAP
      const u32 ffs = __ffsll(xmap);
      s0 += ffs; xmap >>= ffs;
      return s0;
#else
      return (u32)xx[n0++];
#endif
    }
  };
};

__global__ void digit0(equi *eq) {
  uchar hash[HASHOUT];
  blake2b_state state;
  equi::htlayout htl(eq, 0);
  const u32 hashbytes = hashsize(0);
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 block = id; block < NBLOCKS; block += eq->nthreads) {
    state = eq->blake_ctx;
    blake2b_gpu_hash(&state, block, hash, HASHOUT);
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
      const u32 slot = atomicAdd(&eq->nslots[0][bucketid], 1);
      if (slot >= NSLOTS)
        continue;
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

__global__ void digitr(equi *eq, const u32 r) {
  equi::htlayout htl(eq, r);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(r-1, bucketid);
    u32 bsize = eq->getnslots(r-1, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
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
        const u32 xorslot = atomicAdd(&eq->nslots[r&1][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
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

#ifdef UNROLL
__global__ void digit1(equi *eq) {
  equi::htlayout htl(eq, 1);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(0, bucketid);
    u32 bsize = eq->getnslots(0, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2];
        xorbucketid = ((xorbucketid & 0xfff) << 4) | (xhash >> 4);
        xhash &= 0xf;
        const u32 xorslot = atomicAdd(&eq->nslots[1][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(1, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[1].hash ^ hash1[1].hash;
        xorhash[1].hash = hash0[2].hash ^ hash1[2].hash;
        xorhash[2].hash = hash0[3].hash ^ hash1[3].hash;
        xorhash[3].hash = hash0[4].hash ^ hash1[4].hash;
        xorhash[4].hash = hash0[5].hash ^ hash1[5].hash;
      }
    }
  }
}

__global__ void digit2(equi *eq) {
  equi::htlayout htl(eq, 2);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(1, bucketid);
    u32 bsize = eq->getnslots(1, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = (hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2]) >> 4;
        const u32 xorslot = atomicAdd(&eq->nslots[0][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(2, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[0].hash ^ hash1[0].hash;
        xorhash[1].hash = hash0[1].hash ^ hash1[1].hash;
        xorhash[2].hash = hash0[2].hash ^ hash1[2].hash;
        xorhash[3].hash = hash0[3].hash ^ hash1[3].hash;
        xorhash[4].hash = hash0[4].hash ^ hash1[4].hash;
      }
    }
  }
}

__global__ void digit3(equi *eq) {
  equi::htlayout htl(eq, 3);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(2, bucketid);
    u32 bsize = eq->getnslots(2, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2];
        xorbucketid = ((xorbucketid & 0xfff) << 4) | (xhash >> 4);
        xhash &= 0xf;
        const u32 xorslot = atomicAdd(&eq->nslots[1][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(3, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[1].hash ^ hash1[1].hash;
        xorhash[1].hash = hash0[2].hash ^ hash1[2].hash;
        xorhash[2].hash = hash0[3].hash ^ hash1[3].hash;
        xorhash[3].hash = hash0[4].hash ^ hash1[4].hash;
      }
    }
  }
}

__global__ void digit4(equi *eq) {
  equi::htlayout htl(eq, 4);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(3, bucketid);
    u32 bsize = eq->getnslots(3, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = (hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2]) >> 4;
        const u32 xorslot = atomicAdd(&eq->nslots[0][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(4, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[0].hash ^ hash1[0].hash;
        xorhash[1].hash = hash0[1].hash ^ hash1[1].hash;
        xorhash[2].hash = hash0[2].hash ^ hash1[2].hash;
        xorhash[3].hash = hash0[3].hash ^ hash1[3].hash;
      }
    }
  }
}

__global__ void digit5(equi *eq) {
  equi::htlayout htl(eq, 5);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(4, bucketid);
    u32 bsize = eq->getnslots(4, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2];
        xorbucketid = ((xorbucketid & 0xfff) << 4) | (xhash >> 4);
        xhash &= 0xf;
        const u32 xorslot = atomicAdd(&eq->nslots[1][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(5, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[1].hash ^ hash1[1].hash;
        xorhash[1].hash = hash0[2].hash ^ hash1[2].hash;
        xorhash[2].hash = hash0[3].hash ^ hash1[3].hash;
      }
    }
  }
}

__global__ void digit6(equi *eq) {
  equi::htlayout htl(eq, 6);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(5, bucketid);
    u32 bsize = eq->getnslots(5, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = (hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2]) >> 4;
        const u32 xorslot = atomicAdd(&eq->nslots[0][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(6, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[1].hash ^ hash1[1].hash;
        xorhash[1].hash = hash0[2].hash ^ hash1[2].hash;
      }
    }
  }
}

__global__ void digit7(equi *eq) {
  equi::htlayout htl(eq, 7);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(6, bucketid);
    u32 bsize = eq->getnslots(6, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2];
        xorbucketid = ((xorbucketid & 0xfff) << 4) | (xhash >> 4);
        xhash &= 0xf;
        const u32 xorslot = atomicAdd(&eq->nslots[1][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(7, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[0].hash ^ hash1[0].hash;
        xorhash[1].hash = hash0[1].hash ^ hash1[1].hash;
      }
    }
  }
}

__global__ void digit8(equi *eq) {
  equi::htlayout htl(eq, 8);
  equi::collisiondata cd;
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid=id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(7, bucketid);
    u32 bsize = eq->getnslots(7, bucketid);
    for (u32 s1 = 0; s1 < bsize; s1++) {
      const htunit *hash1 = htl.hashbase + s1 * htl.prevhtunits;
      if (!cd.addslot(s1, htl.getxhash(s1, hash1)))
        continue;
      for (; cd.nextcollision(); ) {
        const u32 s0 = cd.slot();
        const htunit *hash0 = htl.hashbase + s0 * htl.prevhtunits;
        if (htl.equal(hash0, hash1))
          continue;
        u32 xorbucketid;
        u32 xhash;
        xorbucketid = ((u32)(hash0->bytes[htl.prevbo]^hash1->bytes[htl.prevbo]) << 8)
          | (hash0->bytes[htl.prevbo+1]^hash1->bytes[htl.prevbo+1]);
        xhash = (hash0->bytes[htl.prevbo+2] ^ hash1->bytes[htl.prevbo+2]) >> 4;
        const u32 xorslot = atomicAdd(&eq->nslots[0][xorbucketid], 1);
        if (xorslot >= NSLOTS)
          continue;
        tree xort; xort.bucketid = bucketid;
        xort.slotid0 = s0; xort.slotid1 = s1;
        xort.xhash = xhash;
        htunit *xorhash = htl.addtree(8, xort, xorbucketid, xorslot);
        xorhash[0].hash = hash0[1].hash ^ hash1[1].hash;
      }
    }
  }
}
#endif

__global__ void digitK(equi *eq) {
  equi::collisiondata cd;
  equi::htlayout htl(eq, WK);
  const u32 id = blockIdx.x * blockDim.x + threadIdx.x;
  for (u32 bucketid = id; bucketid < NBUCKETS; bucketid += eq->nthreads) {
    cd.clear();
    htl.setbucket(WK-1, bucketid);
    u32 bsize = eq->getnslots(WK-1, bucketid);
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
          eq->candidate(xort);
        }
      }
    }
  }
}

#include <unistd.h>

int main(int argc, char **argv) {
  int nthreads = 8192;
  int nonce = 0;
  int tpb = 0;
  int range = 1;
  bool showsol = false;
  const char *header = "";
  int c;
  while ((c = getopt (argc, argv, "h:n:r:t:p:s")) != -1) {
    switch (c) {
      case 'h':
        header = optarg;
        break;
      case 'n':
        nonce = atoi(optarg);
        break;
      case 't':
        nthreads = atoi(optarg);
        break;
      case 'p':
        tpb = atoi(optarg);
        break;
      case 'r':
        range = atoi(optarg);
        break;
      case 's':
        showsol = true;
        break;
    }
  }
  if (!tpb) // if not set, then default threads per block to roughly square root of threads
    for (tpb = 1; tpb*tpb < nthreads; tpb *= 2) ;

  printf("Looking for wagner-tree on (\"%s\",%d", header, nonce);
  if (range > 1)
    printf("-%d", nonce+range-1);
  printf(") with %d %d-bits digits and %d threads (%d per block)\n", NDIGITS, DIGITBITS, nthreads, tpb);
  equi eq(nthreads);
#ifdef JOINHT
  for (u32 r=0; r < WK; r++)
    checkCudaErrors(cudaMalloc((void**)&eq.hta.trees[r], NBUCKETS * NSLOTS * (1 + hhtunits(hhashsize(r))) * sizeof(htunit)));
#else
  checkCudaErrors(cudaMalloc((void**)&eq.hta.trees, WK * NBUCKETS * NSLOTS * sizeof(tree)));
  for (u32 r=0; r < 2; r++)
    checkCudaErrors(cudaMalloc((void**)&eq.hta.hashes[r], NBUCKETS * NSLOTS * hhtunits(hhashsize(r)) * sizeof(htunit)));
#endif
  checkCudaErrors(cudaMalloc((void**)&eq.nslots, 2 * NBUCKETS * sizeof(u32)));
  checkCudaErrors(cudaMalloc((void**)&eq.sols, MAXSOLS * sizeof(proof)));

  equi *device_eq;
  checkCudaErrors(cudaMalloc((void**)&device_eq, sizeof(equi)));

  cudaEvent_t start, stop;
  checkCudaErrors(cudaEventCreate(&start));
  checkCudaErrors(cudaEventCreate(&stop));

  proof sols[MAXSOLS];
  u32 sumnsols = 0;
  for (int r = 0; r < range; r++) {
    cudaEventRecord(start, NULL);
    eq.setnonce(header, nonce+r);
    checkCudaErrors(cudaMemcpy(device_eq, &eq, sizeof(equi), cudaMemcpyHostToDevice));
    printf("Digit 0\n");
    digit0<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(0);
#if BUCKBITS == 16 && RESTBITS == 4 && defined(UNROLL)
    printf("Digit %d\n", 1);
    digit1<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(1);
    printf("Digit %d\n", 2);
    digit2<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(2);
    printf("Digit %d\n", 3);
    digit3<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(3);
    printf("Digit %d\n", 4);
    digit4<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(4);
    printf("Digit %d\n", 5);
    digit5<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(5);
    printf("Digit %d\n", 6);
    digit6<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(6);
    printf("Digit %d\n", 7);
    digit7<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(7);
    printf("Digit %d\n", 8);
    digit8<<<nthreads/tpb,tpb >>>(device_eq);
    eq.showbsizes(8);
#else
    for (u32 r=1; r < WK; r++) {
      printf("Digit %d\n", r);
      digitr<<<nthreads/tpb,tpb >>>(device_eq, r);
      eq.showbsizes(r);
    }
#endif
    printf("Digit %d\n", WK);
    digitK<<<nthreads/tpb,tpb >>>(device_eq);

    checkCudaErrors(cudaMemcpy(&eq, device_eq, sizeof(equi), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(sols, eq.sols, MAXSOLS * sizeof(proof), cudaMemcpyDeviceToHost));
    cudaEventRecord(stop, NULL);
    cudaEventSynchronize(stop);
    float duration;
    cudaEventElapsedTime(&duration, start, stop);
      printf("%d rounds completed in %.3f seconds.\n", WK, duration / 1000.0f);

    u32 nsols = 0;
    for (unsigned s = 0; s < eq.nsols; s++) {
      if (duped(sols[s])) {
        printf("Duped!\n");
        continue;
      }
      nsols++;
      if (showsol) {
        printf("Solution");
        for (int i = 0; i < PROOFSIZE; i++)
          printf(" %jx", (uintmax_t)sols[s][i]);
        printf("\n");
      }
    }
    printf("%d solutions\n", nsols);
    sumnsols += nsols;
  }
  checkCudaErrors(cudaFree(eq.nslots));
  checkCudaErrors(cudaFree(eq.sols));
#ifdef JOINHT
  for (u32 r=0; r < WK; r++)
    checkCudaErrors(cudaFree(eq.hta.trees[r]));
#else
  checkCudaErrors(cudaFree(eq.hta.trees));
  for (u32 r=0; r < 2; r++)
    checkCudaErrors(cudaFree(eq.hta.hashes[r]));
#endif

  printf("%d total solutions\n", sumnsols);
  return 0;
}
