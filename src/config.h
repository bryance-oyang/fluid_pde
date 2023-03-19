#ifndef CONFIG_H
#define CONFIG_H

#define NTHREAD 8

#define NU 320
#define NV 320

#define NSCALAR 0

#define NGHOST 4

// should be 0 unless large discontinuities
#define PPM_ALWAYS_LIM 0
// maybe good for nan cleaning
#define PPM_STRICT_LIM 1

// probably should be 0?
#define WEIRD_PPM 0

#define BROADCAST_PREIMAGE_MIN (-1)
#define BROADCAST_PREIMAGE_MAX 1

#endif /* CONFIG_H */
