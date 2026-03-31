/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************/

#ifndef PRODIGAL_API_H
#define PRODIGAL_API_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Symbol visibility */
#if defined(_WIN32) || defined(__CYGWIN__)
  #ifdef PRODIGAL_BUILDING_DLL
    #define PRODIGAL_API __declspec(dllexport)
  #else
    #define PRODIGAL_API __declspec(dllimport)
  #endif
#elif defined(__GNUC__) && __GNUC__ >= 4
  #define PRODIGAL_API __attribute__((visibility("default")))
#else
  #define PRODIGAL_API
#endif

/* Version */
#define PRODIGAL_VERSION_MAJOR 2
#define PRODIGAL_VERSION_MINOR 6
#define PRODIGAL_VERSION_PATCH 3
#define PRODIGAL_VERSION_STRING "2.6.3"

/* Error codes: zero = success, negative = error */
#define PRODIGAL_OK                   0
#define PRODIGAL_ERR_NOMEM           -1
#define PRODIGAL_ERR_INVALID_CONFIG  -2
#define PRODIGAL_ERR_INVALID_INPUT   -3
#define PRODIGAL_ERR_INTERNAL        -4
#define PRODIGAL_ERR_SEQ_TOO_SHORT   -5
#define PRODIGAL_ERR_CANCELLED       -6

/* Opaque context type */
typedef struct prodigal_ctx prodigal_ctx_t;

/*******************************************************************************
    Configuration
*******************************************************************************/

typedef struct {
    size_t struct_size;         /* Must be first field. Set by config_init. */

    int    trans_table;         /* NCBI translation table (default: 11) */
    int    closed_ends;         /* Nonzero: don't allow genes to run off edges */
    int    mask_regions;        /* Nonzero: treat runs of N as masked */
    int    force_nonsd;         /* Nonzero: bypass Shine-Dalgarno, use motif scan */
    int    meta_mode;           /* Nonzero: metagenomic mode */

    double start_weight;        /* Start score weight (default: 4.35) */

    /* Custom allocator (NULL = use system malloc/free) */
    void  *(*alloc_fn)(size_t size, void *user_data);
    void   (*free_fn)(void *ptr, void *user_data);
    void   *allocator_user_data;  /* Shared user_data for alloc_fn and free_fn */

    /* Logging callback (NULL = discard log messages) */
    void   (*log_callback)(const char *msg, void *user_data);
    void   *log_user_data;

    /* Progress callback (NULL = no progress reporting)
       Return nonzero from callback to cancel computation. */
    int    (*progress_callback)(const char *stage, double frac_done,
                                void *user_data);
    void   *progress_user_data;
} prodigal_config_t;

/*******************************************************************************
    Sequence info
*******************************************************************************/

typedef struct {
    int32_t length;             /* Encoded sequence length in bp */
    double  gc_content;         /* GC fraction [0, 1] */
} prodigal_seq_info_t;

/*******************************************************************************
    Output: Structure of Arrays (SOA) — primary for Arrow/columnar consumers
*******************************************************************************/

typedef struct {
    int32_t  n_genes;

    int32_t *begin;             /* 1-based left coordinate */
    int32_t *end;               /* 1-based right coordinate */
    int32_t *strand;            /* +1 forward, -1 reverse */

    int32_t *partial_left;      /* 1 if gene runs off left edge */
    int32_t *partial_right;     /* 1 if gene runs off right edge */
    int32_t *start_type;        /* 0=ATG, 1=GTG, 2=TTG, 3=Edge */

    double  *cscore;            /* Coding score (6-mer log-odds) */
    double  *sscore;            /* Start score (tscore+rscore+uscore) */
    double  *rscore;            /* RBS motif score */
    double  *uscore;            /* Upstream composition score */
    double  *tscore;            /* Start type score */
    double  *confidence;        /* Confidence [50, 100] */
    double  *gc_cont;           /* Per-gene GC content */

    const char **rbs_motif;     /* RBS motif string (static, not freed) */
    const char **rbs_spacer;    /* RBS spacer string (static, not freed) */

    void    *_base;             /* Single backing allocation (16-byte aligned) */
} prodigal_genes_soa_t;

/*******************************************************************************
    Output: Array of Structures (AOS) — convenience for per-gene iteration
*******************************************************************************/

typedef struct {
    int32_t begin;
    int32_t end;
    int32_t strand;
    int32_t partial_left;
    int32_t partial_right;
    int32_t start_type;
    double  cscore;
    double  sscore;
    double  rscore;
    double  uscore;
    double  tscore;
    double  confidence;
    double  gc_cont;
    const char *rbs_motif;
    const char *rbs_spacer;
} prodigal_gene_t;

typedef struct {
    int32_t         n_genes;
    prodigal_gene_t *genes;     /* Array of n_genes entries */
    void            *_base;     /* Single backing allocation */
} prodigal_genes_t;

/*******************************************************************************
    Statistics (pointer-free, safe to memcpy)
*******************************************************************************/

typedef struct {
    int32_t n_genes;
    int32_t n_nodes;
    double  gc_content;
    int32_t translation_table;
    int32_t uses_sd;
    int32_t best_meta_bin;      /* -1 if not metagenomic */
    char    best_meta_desc[512];
} prodigal_stats_t;

/*******************************************************************************
    API Functions
*******************************************************************************/

/* Config */
PRODIGAL_API void prodigal_config_init(prodigal_config_t *config);

/* Context lifecycle.
   The context struct itself is allocated with system malloc; the custom
   allocator (if provided) is used for internal working buffers only.
   prodigal_destroy(NULL) is a no-op. */
PRODIGAL_API prodigal_ctx_t *prodigal_create(const prodigal_config_t *config);
PRODIGAL_API void prodigal_destroy(prodigal_ctx_t *ctx);

/* Sequence input */
PRODIGAL_API int prodigal_set_sequence(prodigal_ctx_t *ctx,
                                       const char *seq, int32_t len,
                                       const char *header);
PRODIGAL_API int prodigal_set_training_sequences(prodigal_ctx_t *ctx,
                                                  const char **seqs,
                                                  const char **headers,
                                                  const int32_t *lens,
                                                  int32_t n_seqs);
PRODIGAL_API int prodigal_get_seq_info(const prodigal_ctx_t *ctx,
                                       prodigal_seq_info_t *info);

/* Training */
PRODIGAL_API int prodigal_train(prodigal_ctx_t *ctx);
PRODIGAL_API int prodigal_load_training(prodigal_ctx_t *ctx,
                                        const void *data, size_t len);
PRODIGAL_API int prodigal_export_training(const prodigal_ctx_t *ctx,
                                          void **data_out, size_t *len_out);

/* Training parameter setters (fine-grained control) */
PRODIGAL_API int prodigal_set_translation_table(prodigal_ctx_t *ctx, int table);
PRODIGAL_API int prodigal_set_start_weight(prodigal_ctx_t *ctx, double weight);
PRODIGAL_API int prodigal_set_gc(prodigal_ctx_t *ctx, double gc);
PRODIGAL_API int prodigal_set_uses_sd(prodigal_ctx_t *ctx, int uses_sd);

/* Gene finding */
PRODIGAL_API int prodigal_find_genes(prodigal_ctx_t *ctx,
                                     prodigal_genes_soa_t **genes_out,
                                     prodigal_stats_t *stats_out);
PRODIGAL_API int prodigal_find_genes_aos(prodigal_ctx_t *ctx,
                                         prodigal_genes_t **genes_out,
                                         prodigal_stats_t *stats_out);

/* Output cleanup.
   Output structs are always allocated with system malloc (not the custom
   allocator), matching the FastTree convention: output outlives the context
   and must be freeable without a context reference. */
PRODIGAL_API void prodigal_genes_free(prodigal_genes_soa_t *genes);
PRODIGAL_API void prodigal_genes_aos_free(prodigal_genes_t *genes);

/* Error reporting */
PRODIGAL_API const char *prodigal_strerror(int error_code);
PRODIGAL_API const char *prodigal_last_error(const prodigal_ctx_t *ctx);

/* Runtime version query */
PRODIGAL_API const char *prodigal_version_string(void);

#ifdef __cplusplus
}
#endif

#endif /* PRODIGAL_API_H */
