/**
 * @file
 *
 * The hash spelling (CAP-13, AD-21): FNV-1a 64 over little-endian fields,
 * pinned by known answers from tools/gen-hash-vectors.py, which shares no code
 * with the library.
 *
 * Copyright 2026 by Corey Pennycuff
 */

#include "test_helpers.h"

#include "../known_answers.h"

#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

static uint64_t finish(const GSHGEOM_Hash & h) {
  return gshgeom_hash_value(&h);
}

TEST(Hash, ConstantsAreFnv1a64) {
  EXPECT_EQ(GSHGEOM_FNV_OFFSET_BASIS, UINT64_C(0xcbf29ce484222325));
  EXPECT_EQ(GSHGEOM_FNV_PRIME, UINT64_C(0x100000001b3));
}

TEST(Hash, FnvPublishedVectors) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  EXPECT_EQ(finish(h), kFnvEmpty);
  EXPECT_EQ(finish(h), UINT64_C(0xcbf29ce484222325));
  gshgeom_hash_bytes(&h, "a", 1);
  EXPECT_EQ(finish(h), kFnvA);
  gshgeom_hash_init(&h);
  gshgeom_hash_bytes(&h, "foobar", 6);
  EXPECT_EQ(finish(h), kFnvFoobar);
}

TEST(Hash, EveryFieldIsLittleEndian) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u8(&h, 0xA5);
  EXPECT_EQ(finish(h), kHashU8);
  gshgeom_hash_init(&h);
  gshgeom_hash_u16(&h, 0x1234);
  EXPECT_EQ(finish(h), kHashU16);
  gshgeom_hash_init(&h);
  gshgeom_hash_u32(&h, 0x01020304u);
  EXPECT_EQ(finish(h), kHashU32);
  gshgeom_hash_init(&h);
  gshgeom_hash_u64(&h, UINT64_C(0x0102030405060708));
  EXPECT_EQ(finish(h), kHashU64);
}

TEST(Hash, SignedAreTwosComplementAndDoublesAreTheirBits) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_i32(&h, -2);
  EXPECT_EQ(finish(h), kHashI32);
  gshgeom_hash_init(&h);
  gshgeom_hash_i64(&h, -3);
  EXPECT_EQ(finish(h), kHashI64);
  gshgeom_hash_init(&h);
  gshgeom_hash_f64(&h, 1.0);
  EXPECT_EQ(finish(h), kHashF64One);
  gshgeom_hash_init(&h);
  gshgeom_hash_f64(&h, -0.0);
  EXPECT_EQ(finish(h), kHashF64NegZero);
  gshgeom_hash_init(&h);
  gshgeom_hash_f64(&h, 0.1);
  EXPECT_EQ(finish(h), kHashF64Tenth);
  // -0.0 and +0.0 are different numbers in this spelling.
  EXPECT_NE(kHashF64NegZero, UINT64_C(0xcbf29ce484222325));
  EXPECT_EQ(gshgeom_f64_bits(1.0), UINT64_C(0x3ff0000000000000));
  EXPECT_EQ(gshgeom_f64_bits(-0.0), UINT64_C(0x8000000000000000));
  EXPECT_EQ(gshgeom_f64_bits(0.1), UINT64_C(0x3fb999999999999a));
}

TEST(Hash, StringIsLengthThenBytes) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_string(&h, "abc");
  EXPECT_EQ(finish(h), kHashString);
  gshgeom_hash_init(&h);
  gshgeom_hash_string(&h, "");
  EXPECT_EQ(finish(h), kHashEmptyString);
}

TEST(Hash, MixedStreamEqualsOneByteAtATime) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_u8(&h, 1);
  gshgeom_hash_u16(&h, 2);
  gshgeom_hash_u32(&h, 3);
  gshgeom_hash_u64(&h, 4);
  gshgeom_hash_f64(&h, 5.5);
  gshgeom_hash_string(&h, "xy");
  gshgeom_hash_bytes(&h, "zw", 2);
  EXPECT_EQ(finish(h), kHashMixed);

  // The same stream spelled as bytes, by hand.
  GSHGEOM_Hash g;
  gshgeom_hash_init(&g);
  std::vector<unsigned char> stream;
  auto le = [&](uint64_t v, int n) {
    for (int i = 0; i < n; i++) {
      stream.push_back((unsigned char)((v >> (8 * i)) & 0xff));
    }
  };
  le(1, 1);
  le(2, 2);
  le(3, 4);
  le(4, 8);
  le(UINT64_C(0x4016000000000000), 8);  // 5.5
  le(2, 4);
  stream.push_back('x');
  stream.push_back('y');
  stream.push_back('z');
  stream.push_back('w');
  gshgeom_hash_bytes(&g, stream.data(), stream.size());
  EXPECT_EQ(finish(g), kHashMixed);
}

TEST(Hash, ABytesCallOfSizeZeroMixesNothing) {
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_bytes(&h, nullptr, 0);
  gshgeom_hash_bytes(&h, "x", 0);
  EXPECT_EQ(finish(h), GSHGEOM_FNV_OFFSET_BASIS);
}

TEST(Hash, NullIsIgnored) {
  gshgeom_hash_init(nullptr);
  gshgeom_hash_u8(nullptr, 1);
  gshgeom_hash_u16(nullptr, 1);
  gshgeom_hash_u32(nullptr, 1);
  gshgeom_hash_u64(nullptr, 1);
  gshgeom_hash_i32(nullptr, 1);
  gshgeom_hash_i64(nullptr, 1);
  gshgeom_hash_f64(nullptr, 1.0);
  gshgeom_hash_bytes(nullptr, "a", 1);
  gshgeom_hash_string(nullptr, "a");
  EXPECT_EQ(gshgeom_hash_value(nullptr), GSHGEOM_FNV_OFFSET_BASIS);
  GSHGEOM_Hash h;
  gshgeom_hash_init(&h);
  gshgeom_hash_string(&h, nullptr);
  gshgeom_hash_bytes(&h, nullptr, 3);
  EXPECT_EQ(finish(h), GSHGEOM_FNV_OFFSET_BASIS);
}

TEST(Hash, AHashIsAValueThatForks) {
  GSHGEOM_Hash a;
  gshgeom_hash_init(&a);
  gshgeom_hash_u32(&a, 7);
  GSHGEOM_Hash b = a;
  gshgeom_hash_u8(&b, 1);
  EXPECT_NE(finish(a), finish(b));
  gshgeom_hash_u8(&a, 1);
  EXPECT_EQ(finish(a), finish(b));
}

TEST(Hash, EveryByteOfEveryFieldMatters) {
  // A defect that drops, repeats or reorders a byte of a wide field changes
  // one of these.
  GSHGEOM_Hash h;
  uint64_t seen[64];
  for (int i = 0; i < 64; i++) {
    gshgeom_hash_init(&h);
    gshgeom_hash_u64(&h, UINT64_C(1) << i);
    seen[i] = finish(h);
    for (int j = 0; j < i; j++) {
      EXPECT_NE(seen[i], seen[j]) << i << " " << j;
    }
  }
}

int main(int argc, char ** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
