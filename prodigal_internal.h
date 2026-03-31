/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Internal header -- not part of the public API.
*******************************************************************************/

#ifndef PRODIGAL_INTERNAL_H
#define PRODIGAL_INTERNAL_H

#include <stdlib.h>
#include "prodigal.h"
#include "sequence.h"
#include "node.h"
#include "gene.h"
#include "training.h"
#include "metagenomic.h"
#include "dprog.h"
#include "bitmap.h"

/*******************************************************************************
    Internal context structure (opaque to public API consumers)
*******************************************************************************/

struct prodigal_ctx {
    prodigal_config_t config;         /* Snapshot of caller config */

    /* Training data */
    struct _training tinf;            /* ~558KB training model */
    int trained;                      /* Nonzero if training is loaded/complete */

    /* Sequence buffers */
    unsigned char *seq;               /* 2-bit encoded forward sequence */
    unsigned char *rseq;              /* 2-bit encoded reverse complement */
    unsigned char *useq;              /* Ambiguity bitmap (N bases) */
    int slen;                         /* Current sequence length in bp */
    int max_slen;                     /* Max sequence length seen (for realloc) */
    double gc;                        /* GC content of current sequence */
    char cur_header[MAX_LINE];        /* Current sequence header */

    /* Masking */
    mask mlist[MAX_MASKS];
    int nmask;

    /* Node and gene working arrays */
    struct _node *nodes;
    int nn;                           /* Current node count */
    struct _gene *genes;

    /* Metagenomic bins (lazily initialized) */
    struct _metagenomic_bin *meta;
    int meta_initialized;

    /* Error state */
    char error_msg[1024];
    int error_code;
};

/*******************************************************************************
    Internal helpers
*******************************************************************************/

static inline void *pdg_alloc(prodigal_ctx_t *ctx, size_t size) {
    if (ctx->config.alloc_fn)
        return ctx->config.alloc_fn(size, ctx->config.allocator_user_data);
    return malloc(size);
}

static inline void pdg_free(prodigal_ctx_t *ctx, void *ptr) {
    if (ptr == NULL) return;
    if (ctx->config.free_fn)
        ctx->config.free_fn(ptr, ctx->config.allocator_user_data);
    else
        free(ptr);
}

#define PDG_LOG(ctx, fmt, ...) do { \
    if ((ctx)->config.log_callback) { \
        char _pdg_buf[512]; \
        snprintf(_pdg_buf, sizeof(_pdg_buf), fmt, ##__VA_ARGS__); \
        (ctx)->config.log_callback(_pdg_buf, (ctx)->config.log_user_data); \
    } \
} while(0)

#define PDG_SET_ERROR(ctx, code, fmt, ...) do { \
    (ctx)->error_code = (code); \
    snprintf((ctx)->error_msg, sizeof((ctx)->error_msg), fmt, ##__VA_ARGS__); \
} while(0)

#endif /* PRODIGAL_INTERNAL_H */
