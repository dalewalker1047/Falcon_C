#include <stdlib.h>
#include <stdint.h>


/*
 * Structure for a PRNG. This includes a large buffer so that values
 * get generated in advance. The 'state' is used to keep the current
 * PRNG algorithm state (contents depend on the selected algorithm).
 *
 * The unions with 'dummy_u64' are there to ensure proper alignment for
 * 64-bit direct access.
 */
typedef struct {
    union {
        uint8_t d[512]; /* MUST be 512, exactly */
        uint64_t dummy_u64;
    } buf;
    size_t ptr;
    union {
        uint8_t d[256];
        uint64_t dummy_u64;
    } state;
    int type;
} prng;

typedef struct {
    union {
        uint64_t A[25];
        uint8_t dbuf[200];
    } st;
    uint64_t dptr;
} inner_shake256_context;