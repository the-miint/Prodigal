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

/**
 * @file prodigal.h
 * @brief Reentrant C library API for Prodigal gene prediction.
 *
 * Overview
 * --------
 * Prodigal (PROkaryotic DynamIc Programming Genefinding ALgorithm) predicts
 * protein-coding genes in prokaryotic genomes. This header defines a library
 * API that allows embedding Prodigal in other applications without shelling
 * out to the command-line tool.
 *
 * Two usage patterns are supported:
 *
 *   Metagenomic mode (no training):
 *     config_init -> create -> set_sequence -> find_genes -> free -> destroy
 *
 *   Single-genome mode (train, or load a model with load_training):
 *     config_init -> create -> set_training_sequences -> train ->
 *     set_sequence -> find_genes -> free -> (repeat per sequence) -> destroy
 *
 * Lifecycle
 * ---------
 *   1. Call prodigal_config_init() to populate a config with safe defaults.
 *   2. Modify config fields as needed (trans_table, meta_mode, callbacks...).
 *   3. Call prodigal_create() to allocate an opaque context.
 *   4. Load input:
 *      - Metagenomic mode: call prodigal_set_sequence() per contig, then
 *        prodigal_find_genes().  Training data is built-in (50 pre-trained
 *        models, lazily initialized on first use), so the calls that train,
 *        load or adjust a model are refused (prodigal_config_t.meta_mode).
 *      - Single-genome mode: call prodigal_set_training_sequences() with all
 *        contigs, then prodigal_train().  Alternatively, load a previously
 *        exported training file with prodigal_load_training().  Then call
 *        prodigal_set_sequence() + prodigal_find_genes() per contig.
 *   5. Free each output with prodigal_genes_free() / prodigal_genes_aos_free().
 *   6. Call prodigal_destroy() when done.
 *
 * Memory model
 * ------------
 *   - The context struct itself is allocated with system malloc/free.
 *   - The context's working memory uses the custom allocator if provided,
 *     otherwise system malloc/free: the sequence, nodes, genes and masked
 *     regions, training's scratch and, in metagenomic mode, the 50 bins
 *     (sizes and lifetimes: prodigal_create()).
 *   - Output structs (prodigal_genes_soa_t, prodigal_genes_t) never use the
 *     custom allocator, so they can be freed without a context: with
 *     prodigal_genes_free() / prodigal_genes_aos_free(), never free() (the
 *     SOA's backing block is an aligned allocation, _aligned_malloc on
 *     Windows).
 *   - prodigal_destroy(NULL) and prodigal_genes_free(NULL) are safe no-ops.
 *
 * Output formats
 * --------------
 *   - SOA (prodigal_genes_soa_t): Structure-of-Arrays layout ideal for
 *     columnar / Apache Arrow consumers, in one backing allocation (see the
 *     struct).
 *   - AOS (prodigal_genes_t): Array-of-Structures layout for convenient
 *     per-gene iteration.  Also backed by a single allocation.
 *
 * Error handling
 * --------------
 *   - All fallible functions return int: PRODIGAL_OK (0) on success,
 *     negative error codes on failure.
 *   - prodigal_strerror() maps an error code to a static category string.
 *   - prodigal_last_error() returns a detailed message from the context;
 *     read it right after the call that failed (see prodigal_last_error()).
 *   - The context is reusable after an error -- call prodigal_set_sequence()
 *     again and retry.
 *   - The library never calls exit(), abort(), or assert().
 *   - The library never writes to stdout or stderr.  prodigal_train()
 *     reports its progress to log_callback.
 *
 * Thread safety
 * -------------
 *   - Each context is independent.  Multiple contexts can be used
 *     concurrently from different threads without synchronization.
 *   - A single context must not be used from multiple threads simultaneously.
 *   - No global mutable state exists in the library.
 *   - Callbacks and alloc_fn/free_fn run on the thread that made the call;
 *     a callback must not call the library on the same context.  Contexts
 *     that share an allocator call it concurrently.
 *   - prodigal_config_validate(), prodigal_training_info(),
 *     prodigal_strerror() and the version and ABI functions need no
 *     context.
 *
 * GPL-boundary integration
 * ------------------------
 *   - prodigal_config_init_v() refuses a caller with another ABI version
 *     or config size, and prodigal_abi_layout() lets bindings check every
 *     struct's field offsets in their tests.
 *   - Build with -DPRODIGAL_NO_MAIN to exclude the CLI entry point when
 *     compiling as a static library for embedding.
 *   - The Makefile builds libprodigal.a (static) and libprodigal.so
 *     (shared) when asked for them (plain make builds the CLI).  Neither
 *     needs zlib: only the CLI is built with gzip support.
 *   - Build with -ffp-contract=off, as the Makefile does: without it ARM64
 *     builds fuse multiply-adds and print some scores differently from
 *     x86-64 in the last digit.
 *
 * Quick example (metagenomic mode)
 * --------------------------------
 * @code
 *   prodigal_config_t config;
 *   prodigal_config_init(&config);
 *   config.meta_mode = 1;
 *
 *   prodigal_ctx_t *ctx = prodigal_create(&config);
 *   if (ctx == NULL) { handle error }
 *
 *   prodigal_set_sequence(ctx, seq_chars, seq_len, "contig_1");
 *
 *   prodigal_genes_soa_t *genes = NULL;
 *   prodigal_stats_t stats;
 *   int rc = prodigal_find_genes(ctx, &genes, &stats);
 *   if (rc != PRODIGAL_OK) { handle error }
 *
 *   for (int i = 0; i < genes->n_genes; i++) {
 *       printf("gene %d: %d..%d strand=%d score=%.1f\n", i,
 *              genes->begin[i], genes->end[i], genes->strand[i],
 *              genes->cscore[i] + genes->sscore[i]);
 *   }
 *
 *   prodigal_genes_free(genes);
 *   prodigal_destroy(ctx);
 * @endcode
 *
 * Quick example (single-genome with training)
 * --------------------------------------------
 * @code
 *   prodigal_config_t config;
 *   prodigal_config_init(&config);
 *
 *   prodigal_ctx_t *ctx = prodigal_create(&config);
 *
 *   // Load all contigs for training
 *   prodigal_set_training_sequences(ctx, seqs, headers, lens, n_seqs);
 *   int rc = prodigal_train(ctx);
 *   if (rc != PRODIGAL_OK) { handle error }
 *
 *   // Find genes per contig
 *   for (int s = 0; s < n_seqs; s++) {
 *       prodigal_set_sequence(ctx, seqs[s], lens[s], headers[s]);
 *       prodigal_genes_soa_t *genes = NULL;
 *       prodigal_find_genes(ctx, &genes, NULL);
 *       // ... process genes ...
 *       prodigal_genes_free(genes);
 *   }
 *
 *   prodigal_destroy(ctx);
 * @endcode
 *
 * Quick example (custom allocator)
 * --------------------------------
 * @code
 *   static void *my_alloc(size_t size, void *ud) { return my_pool_alloc(ud, size); }
 *   static void  my_free(void *ptr, void *ud)    { my_pool_free(ud, ptr); }
 *
 *   prodigal_config_t config;
 *   prodigal_config_init(&config);
 *   config.alloc_fn = my_alloc;
 *   config.free_fn  = my_free;
 *   config.allocator_user_data = my_pool;
 *   // ... create, use, destroy as usual ...
 * @endcode
 */

#ifndef PRODIGAL_API_H
#define PRODIGAL_API_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Symbol visibility.
   For static linking (the default), PRODIGAL_API is empty.
   For building a shared library (DLL), define PRODIGAL_BUILDING_DLL.
   For consuming a shared library on Windows, define PRODIGAL_DLL. */
#if defined(_WIN32) || defined(__CYGWIN__)
  #if defined(PRODIGAL_BUILDING_DLL)
    #define PRODIGAL_API __declspec(dllexport)
  #elif defined(PRODIGAL_DLL)
    #define PRODIGAL_API __declspec(dllimport)
  #else
    #define PRODIGAL_API
  #endif
#elif defined(__GNUC__) && __GNUC__ >= 4
  #define PRODIGAL_API __attribute__((visibility("default")))
#else
  #define PRODIGAL_API
#endif

/** @name Compile-time version constants */
/**@{*/
#define PRODIGAL_VERSION_MAJOR 2
#define PRODIGAL_VERSION_MINOR 6
#define PRODIGAL_VERSION_PATCH 5
#define PRODIGAL_VERSION_STRING "2.6.5"
/**@}*/

/**
 * @name Error codes
 * Zero is success; all errors are negative.  Use prodigal_strerror() for a
 * human-readable category and prodigal_last_error() for a detailed message.
 */
/**@{*/
#define PRODIGAL_OK                   0
#define PRODIGAL_ERR_NOMEM           -1  /**< Memory allocation failed */
#define PRODIGAL_ERR_INVALID_CONFIG  -2  /**< Bad config (NULL, wrong struct_size, invalid field) */
#define PRODIGAL_ERR_INVALID_INPUT   -3  /**< Bad argument, or a call the context's state does not allow (see each function) */
#define PRODIGAL_ERR_INTERNAL        -4  /**< Internal error (should not happen) */
#define PRODIGAL_ERR_SEQ_TOO_SHORT   -5  /**< Sequence < 20000 bp for single-genome training */
#define PRODIGAL_ERR_CANCELLED       -6  /**< Cancelled by progress_callback returning nonzero */
/**@}*/

/**
 * Opaque context.  All computation state lives here.
 * Created by prodigal_create(), destroyed by prodigal_destroy().
 */
typedef struct prodigal_ctx prodigal_ctx_t;

/*******************************************************************************
    Configuration
*******************************************************************************/

/**
 * Configuration struct passed to prodigal_create().
 *
 * Initialize with prodigal_config_init() which sets struct_size and all
 * defaults.  Modify fields as needed, then pass to prodigal_create().
 * The context takes a snapshot of the config; further changes to the config
 * struct after create() have no effect.
 *
 * New fields are always appended at the end.  Never reorder or remove fields.
 * The ABI is checked by prodigal_config_init(), and by prodigal_create() for
 * struct_size (the first field); see also prodigal_abi_layout().
 */
typedef struct {
    /** Must be first field.  Set by prodigal_config_init().
     *  prodigal_create() rejects configs where struct_size !=
     *  sizeof(prodigal_config_t); that misses an addition that fills
     *  padding, which prodigal_config_init() catches. */
    size_t struct_size;

    /** NCBI translation table (genetic code) number.  Default: 11
     *  (Bacterial, Archaeal and Plant Plastid).  Supported: 4, 5, 9, 11,
     *  13, 15, 21, 25, 26, 32, the NCBI codes Prodigal uses exactly as NCBI
     *  defines them (docs/genetic-codes.md).  An unsupported table
     *  makes prodigal_create() return NULL; prodigal_config_validate()
     *  says why.  With meta_mode only 11 (see meta_mode): a genome that
     *  uses another table needs single mode. */
    int    trans_table;

    /** Nonzero: closed ends -- do not allow genes to run off sequence edges.
     *  Default: 0 (allow edge genes). */
    int    closed_ends;

    /** Nonzero: treat runs of >= 50 N's as masked regions; don't build genes
     *  across them.  Default: 0. */
    int    mask_regions;

    /** Nonzero: bypass Shine-Dalgarno trainer and force a full upstream
     *  motif scan.  Default: 0 (auto-detect SD usage).  Refused with
     *  meta_mode. */
    int    force_nonsd;

    /** Nonzero: metagenomic mode (50 pre-trained bins, each with its own
     *  model, table 11 or 4).  Zero: single-genome mode (requires
     *  training).  Default: 0.  Metagenomic mode refuses the single-genome
     *  model's settings, which it would ignore: in the config, a
     *  trans_table other than 11, a nonzero force_nonsd, or a start_weight
     *  other than exactly 4.35 (prodigal_config_validate() says which); on
     *  the context, prodigal_train(), prodigal_load_training() and the
     *  setters, whatever their arguments, return PRODIGAL_ERR_INVALID_INPUT
     *  with the reason and change nothing, except
     *  prodigal_set_translation_table(11), which is accepted and changes
     *  nothing.  So it has no model to export.
     *  prodigal_set_training_sequences() is accepted: the sequences, spacers
     *  included, become the one sequence to search. */
    int    meta_mode;

    /** Start score weight, finite and > 0.  Affects the balance between
     *  coding potential and start signal strength.  Default: 4.35.  Rarely
     *  needs changing.  With meta_mode only exactly 4.35 is accepted. */
    double start_weight;

    /** Custom allocator for internal working buffers.  If NULL, system
     *  malloc is used.  Must return memory suitable for any alignment
     *  (16-byte aligned recommended).  What it serves: Memory model, at
     *  the top of this file.  Set both alloc_fn and free_fn, or neither: free_fn without alloc_fn
     *  is refused, and with alloc_fn but no free_fn the library frees
     *  alloc_fn's blocks with system free(). */
    void  *(*alloc_fn)(size_t size, void *user_data);

    /** Custom deallocator paired with alloc_fn.  If NULL, system free is
     *  used.  Called with pointers previously returned by alloc_fn. */
    void   (*free_fn)(void *ptr, void *user_data);

    /** User data pointer passed as the second argument to both alloc_fn
     *  and free_fn.  Typically a pool or arena handle. */
    void   *allocator_user_data;

    /** Logging callback.  Receives prodigal_train()'s progress messages,
     *  the only ones the library logs.  If NULL, they are discarded.  The
     *  msg string is valid only for the duration of the callback. */
    void   (*log_callback)(const char *msg, void *user_data);

    /** User data pointer passed to log_callback. */
    void   *log_user_data;

    /** Progress callback.  Called by prodigal_find_genes() in metagenomic
     *  mode only, before each of the 50 bins, with stage "metagenomic
     *  scoring" and the fraction of bins done, i/50 (0.0 to 0.98; never
     *  1.0).  Return 0 to continue, nonzero to cancel: prodigal_find_genes()
     *  then returns PRODIGAL_ERR_CANCELLED, with the sequence still loaded.
     *  If NULL, no progress reporting. */
    int    (*progress_callback)(const char *stage, double frac_done,
                                void *user_data);

    /** User data pointer passed to progress_callback. */
    void   *progress_user_data;

    /** Nonzero: prodigal_find_genes() also returns each gene's protein
     *  (prodigal_genes_soa_t's protein_offsets and protein_data;
     *  prodigal_gene_t's protein).  Default: 0 (no proteins). */
    int    want_proteins;
} prodigal_config_t;

/*******************************************************************************
    Sequence info
*******************************************************************************/

/** Basic information about the currently loaded sequence. */
typedef struct {
    /** Encoded sequence length in bp: the input's characters in 'A'..'z'
     *  (the rest are dropped), plus, after prodigal_set_training_sequences
     *  with several sequences, 12 for each spacer (one after each). */
    int32_t length;
    double  gc_content;         /**< GC fraction [0, 1] of length (ambiguous bases and spacers count as not GC).  0.0 if no sequence loaded. */
} prodigal_seq_info_t;

/*******************************************************************************
    Output: Structure of Arrays (SOA)
*******************************************************************************/

/**
 * Gene prediction results in Structure-of-Arrays layout.
 *
 * Preferred output format for columnar / Apache Arrow consumers.  All numeric
 * arrays are carved from a single 16-byte-aligned backing allocation pointed
 * to by _base.  Call prodigal_genes_free() to release.
 *
 * Gene i is the CLI's gene i+1 for the sequence (GFF ID <seqnum>_<i+1>).
 *
 * Every value is the one the CLI prints in its GFF attributes (partial=,
 * start_type=, rbs_motif=, ..., tscore=): the values that chose the gene.
 * In metagenomic mode the CLI's GFF column 6 differs from cscore + sscore
 * (an upstream quirk, docs/upstream-behavior.md).
 *
 * The strings (rbs_motif, rbs_spacer) and the proteins are stored in the
 * same allocation: they are valid until prodigal_genes_free() and are not
 * freed on their own.
 *
 * If n_genes == 0, every array pointer and _base are NULL.
 */
typedef struct {
    int32_t  n_genes;           /**< Number of predicted genes */

    /** 1-based, inclusive coordinates, begin < end on both strands, the stop
     *  codon included (unless the gene runs off the edge there); in the
     *  sequence as encoded (characters outside 'A'..'z' dropped: see
     *  prodigal_set_sequence()).  Length n_genes. */
    int32_t *begin;             /**< Left coordinate */
    int32_t *end;               /**< Right coordinate */
    int32_t *strand;            /**< +1 forward, -1 reverse, length n_genes */

    int32_t *partial_left;      /**< 1 if gene runs off left edge, length n_genes */
    int32_t *partial_right;     /**< 1 if gene runs off right edge, length n_genes */
    int32_t *start_type;        /**< 0=ATG, 1=GTG, 2=TTG, 3=Edge, length n_genes */

    double  *cscore;            /**< Coding score (6-mer log-odds), length n_genes */
    double  *sscore;            /**< Start score as the CLI reports it; in general not tscore + rscore + uscore (adjustments: docs/internals.md, "Scoring").  Length n_genes */
    double  *rscore;            /**< RBS motif score, length n_genes */
    double  *uscore;            /**< Upstream composition score, length n_genes */
    double  *tscore;            /**< Start type score, length n_genes */
    double  *confidence;        /**< Confidence in [50, 100], length n_genes */
    double  *gc_cont;           /**< GC fraction [0, 1] of the gene as Prodigal counts it: ambiguous bases count as GC, and a reverse-strand gene's window is 2 bases off (docs/upstream-behavior.md, quirks 3 and 14).  Length n_genes */

    const char **rbs_motif;     /**< RBS motif, as the CLI prints it (e.g. "AGGAG", "TAA", "None"), length n_genes */
    const char **rbs_spacer;    /**< RBS spacer, as the CLI prints it (e.g. "5-10bp", "8bp", "None"), length n_genes */

    /** Proteins, when the config set want_proteins; otherwise NULL.
     *  Arrow-style (int32 offsets, as Arrow's Utf8): gene i's protein is
     *  the protein_offsets[i+1] - protein_offsets[i] letters at
     *  protein_data + protein_offsets[i]; there are n_genes + 1 offsets,
     *  the first 0, and protein_data has no NUL bytes.  With n_genes == 0
     *  both are NULL (an Arrow consumer supplies the single 0 offset).
     *  A protein is the CLI's -a record without its line breaks: one
     *  letter per codon of the gene, from its start, with X for a codon
     *  holding an ambiguous base, M for the start codon (unless the gene
     *  runs off the edge there), * for the stop codon (none when the gene
     *  runs off the edge at its stop end), and between them the amino acids
     *  of the table the genes were called with
     *  (prodigal_stats_t.translation_table), as NCBI defines it
     *  (docs/genetic-codes.md). */
    int32_t *protein_offsets;
    char    *protein_data;

    void    *_base;             /**< Single backing allocation (16-byte aligned).  Internal. */
} prodigal_genes_soa_t;

/*******************************************************************************
    Output: Array of Structures (AOS)
*******************************************************************************/

/** Per-gene data in struct form.  See prodigal_genes_soa_t for field docs. */
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
    const char *rbs_motif;      /**< Valid until prodigal_genes_aos_free() */
    const char *rbs_spacer;     /**< Valid until prodigal_genes_aos_free() */
    /** The protein, NUL-terminated (see prodigal_genes_soa_t), or NULL
     *  unless the config set want_proteins.  Valid until
     *  prodigal_genes_aos_free(). */
    const char *protein;
} prodigal_gene_t;

/**
 * Gene prediction results in Array-of-Structures layout.
 * Call prodigal_genes_aos_free() to release.
 */
typedef struct {
    int32_t         n_genes;    /**< Number of predicted genes */
    prodigal_gene_t *genes;     /**< Array of n_genes entries; NULL if n_genes == 0 */
    void            *_base;     /**< Single backing allocation (NULL if n_genes == 0).  Internal. */
} prodigal_genes_t;

/*******************************************************************************
    Statistics
*******************************************************************************/

/**
 * Computation statistics.  Pointer-free, safe to memcpy.
 * Passed by pointer to prodigal_find_genes(); may be NULL if not needed.
 * Filled when it succeeds; not to be read after an error.
 */
typedef struct {
    int32_t n_genes;            /**< Number of genes found */
    int32_t n_nodes;            /**< Number of start/stop nodes evaluated */
    double  gc_content;         /**< Sequence GC fraction, as prodigal_get_seq_info() */
    int32_t translation_table;  /**< Table the genes were called with: the model's, or in metagenomic mode the bin's */
    int32_t uses_sd;            /**< Nonzero if that model or bin uses Shine-Dalgarno motifs */
    int32_t best_meta_bin;      /**< The bin that called the genes, 0..49 (-1 in single-genome mode) */
    char    best_meta_desc[512];/**< Its description ("" in single-genome mode) */
} prodigal_stats_t;

/*******************************************************************************
    Model settings
*******************************************************************************/

/** A model's settings, as prodigal_training_info() reads them from a blob.
 *  Pointer-free, safe to memcpy. */
typedef struct {
    int32_t translation_table;  /**< NCBI translation table */
    int32_t uses_sd;            /**< 1 if the model scores Shine-Dalgarno motifs, else 0 */
    double  gc_content;         /**< The model's GC fraction [0, 1]: the training sequence's, unless prodigal_set_gc changed it */
    double  start_weight;       /**< Start score weight (> 0) */
} prodigal_training_info_t;

/*******************************************************************************
    ABI checks
*******************************************************************************/

/** This header's ABI version.  Bumped whenever a public struct's layout, a
 *  function's signature or a public constant's value changes (the version
 *  macros aside: a release changes them).
 *  prodigal_config_init() checks it; bindings with their own mirrors (such
 *  as GPL-boundary's Rust #[repr(C)] structs) pass theirs to
 *  prodigal_config_init_v() or compare it with prodigal_abi_version(). */
#define PRODIGAL_ABI_VERSION 1

/** One row of prodigal_abi_layout(): a public struct, or one of its fields.
 *  This struct is in the table too, but a binding must mirror it to read
 *  the table, so its layout changes only with PRODIGAL_ABI_VERSION. */
typedef struct {
    const char *type;    /**< Struct name, e.g. "prodigal_config_t"; NULL ends the table */
    const char *field;   /**< Field name, or NULL for the struct's own row */
    size_t      offset;  /**< offsetof the field (0 for the struct's row) */
    size_t      size;    /**< sizeof the field, or of the struct */
} prodigal_abi_field_t;

/*******************************************************************************
    API Functions
*******************************************************************************/

/**
 * Initialize a config struct with safe defaults.
 *
 * Sets struct_size, trans_table=11, start_weight=4.35, all other fields to
 * zero/NULL.  A config initialized this way is ready for prodigal_create()
 * without further modification (defaults to single-genome mode, table 11).
 *
 * C callers use the macro prodigal_config_init(config), which passes this
 * header's PRODIGAL_ABI_VERSION and sizeof(*config).  A library with
 * another ABI then refuses and writes nothing, so the check can't be
 * forgotten.  Bindings in other languages pass the ABI version they were
 * written for and their struct's size.
 *
 * @param config       Config to initialize.
 * @param abi_version  The caller's PRODIGAL_ABI_VERSION.
 * @param struct_size  The caller's sizeof(prodigal_config_t).
 * @return PRODIGAL_OK, or PRODIGAL_ERR_INVALID_CONFIG (config unchanged)
 *         for a NULL config, another ABI version or another size.
 */
PRODIGAL_API int prodigal_config_init_v(prodigal_config_t *config,
                                        int abi_version, size_t struct_size);
#define prodigal_config_init(config) \
    prodigal_config_init_v((config), PRODIGAL_ABI_VERSION, sizeof(*(config)))

/**
 * Create a new Prodigal context.
 *
 * The context starts small (about 1.4 MB: its model, and its working
 * arrays at their floors).  The arrays grow as each sequence needs (in
 * prodigal_set_sequence, prodigal_set_training_sequences, prodigal_train
 * and prodigal_find_genes, any of which can then return
 * PRODIGAL_ERR_NOMEM).  The nodes and genes go back to their floors when
 * prodigal_train() or prodigal_find_genes() returns, and the sequence's
 * bitmaps (0.625-1.25 B a base) at the next prodigal_set_sequence() or
 * prodigal_set_training_sequences().  The masked-region list and, in
 * metagenomic mode, the bins (about 28 MB, from the first
 * prodigal_find_genes()) are kept until prodigal_destroy().  The config
 * is snapshotted; the caller may modify or discard it after this call.
 *
 * Returns NULL on failure: a config prodigal_config_validate() rejects, or
 * out of memory.  No context exists then to carry an error message; call
 * prodigal_config_validate() for the reason.
 *
 * @param config  Pointer to initialized config.  Must not be NULL.
 * @return New context, or NULL on failure.
 */
PRODIGAL_API prodigal_ctx_t *prodigal_create(const prodigal_config_t *config);

/**
 * Check a config the way prodigal_create() does, and say what is wrong.
 * Needs no context; thread-safe.
 *
 * Checked in this order, and the first failure is reported: NULL; a
 * struct_size other than this library's sizeof(prodigal_config_t); an
 * unsupported trans_table (docs/genetic-codes.md), or with meta_mode any
 * but 11; a start_weight that is not finite and > 0; with meta_mode, a
 * nonzero force_nonsd, then a start_weight other than exactly 4.35 (the
 * default); free_fn without alloc_fn.  alloc_fn without free_fn is
 * accepted (see alloc_fn); the ABI version is prodigal_config_init()'s
 * check.
 *
 * @param config   Config to check.  May be NULL (rejected).
 * @param msg      Receives the reason, or "" if the config is valid; may be
 *                 NULL.  Truncated to msg_len bytes, NUL included; nothing
 *                 is written when msg_len is 0.
 * @param msg_len  Size of msg in bytes; messages fit in 256.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_CONFIG.
 */
PRODIGAL_API int prodigal_config_validate(const prodigal_config_t *config,
                                          char *msg, size_t msg_len);

/**
 * Destroy a context and free all associated memory.
 *
 * Safe to call with NULL (no-op).  After this call, the context pointer is
 * invalid.  Any output structs (prodigal_genes_soa_t, prodigal_genes_t)
 * previously returned remain valid until explicitly freed.
 *
 * @param ctx  Context to destroy, or NULL.
 */
PRODIGAL_API void prodigal_destroy(prodigal_ctx_t *ctx);

/**
 * Load a single sequence for gene finding.
 *
 * Encodes the raw ASCII nucleotide string into Prodigal's internal 2-bit
 * representation, computes the reverse complement, and calculates GC content.
 * A, C, G and T, in either case, are bases.  Characters outside 'A'..'z'
 * (digits, '-', '.', '*', whitespace, ...) are dropped, as the CLI drops
 * them, which shifts the coordinates of everything after them:
 * prodigal_get_seq_info() reports how many characters were kept, so a
 * caller can compare that with len.  Every other character in 'A'..'z'
 * (N, IUPAC codes, '[', '_', ...) is an ambiguous base, encoded as C with
 * an ambiguity flag (how GC counts it: prodigal_seq_info_t.gc_content,
 * prodigal_genes_soa_t.gc_cont; training counts it as GC).
 *
 * At most 31,989,999 characters may be kept (MAX_SEQ - MAX_LINE - 1), the
 * CLI's limit; a longer sequence is an error (PRODIGAL_ERR_INVALID_INPUT).
 * PRODIGAL_ERR_NOMEM: no memory for the sequence or, with mask_regions,
 * its masked regions.  On an error no
 * sequence is loaded.
 *
 * Resets internal state from any previous sequence.  The context can be
 * reused across multiple set_sequence + find_genes cycles.
 *
 * @param ctx     Context.
 * @param seq     Raw nucleotide string (not null-terminated required; len used).
 * @param len     Number of characters in seq.  Must be > 0.
 * @param header  Sequence name/header.  May be NULL.  Nothing the library
 *                returns uses it.
 * @return PRODIGAL_OK; PRODIGAL_ERR_INVALID_INPUT for a NULL ctx or seq,
 *         len <= 0, no character in 'A'..'z', or more than 31,989,999 kept;
 *         PRODIGAL_ERR_NOMEM.
 */
PRODIGAL_API int prodigal_set_sequence(prodigal_ctx_t *ctx,
                                       const char *seq, int32_t len,
                                       const char *header);

/**
 * Load multiple sequences concatenated for single-genome training.
 *
 * Concatenates all sequences with TTAATTAATTAA stop-codon spacers (forcing
 * stops in all 6 reading frames), matching Prodigal's training behavior.
 * After this call, use prodigal_train() to build the model.  Characters are
 * kept or dropped as in prodigal_set_sequence(), and the same limit applies
 * to the whole, spacers included: a longer input is an error, where the
 * CLI warns and trains on the start of it.  With more than one sequence,
 * a 12-base spacer follows each, and a NULL entry, or one of length <= 0,
 * adds only its spacer (so spacers alone are accepted).  On an error no
 * sequence is loaded.
 *
 * @param ctx      Context.
 * @param seqs     Array of n_seqs raw nucleotide strings.
 * @param headers  Array of n_seqs header strings (may be NULL).  Only
 *                 headers[0] is kept, as prodigal_set_sequence()'s header.
 * @param lens     Array of n_seqs sequence lengths.
 * @param n_seqs   Number of sequences.  Must be > 0.
 * @return As prodigal_set_sequence(), and PRODIGAL_ERR_INVALID_INPUT for a
 *         NULL seqs or lens, or n_seqs <= 0.
 */
PRODIGAL_API int prodigal_set_training_sequences(prodigal_ctx_t *ctx,
                                                  const char **seqs,
                                                  const char **headers,
                                                  const int32_t *lens,
                                                  int32_t n_seqs);

/**
 * Get information about the currently loaded sequence.
 *
 * Returns length=0 and gc_content=0.0 if no sequence is loaded.
 *
 * @param ctx   Context.
 * @param info  Output struct to populate.
 * @return PRODIGAL_OK, or PRODIGAL_ERR_INVALID_INPUT for a NULL argument
 *         (no message: the context is const).
 */
PRODIGAL_API int prodigal_get_seq_info(const prodigal_ctx_t *ctx,
                                       prodigal_seq_info_t *info);

/**
 * Train a gene-finding model from the loaded training sequences.
 *
 * Trains on the loaded sequence: normally every contig, loaded with
 * prodigal_set_training_sequences(); one loaded with
 * prodigal_set_sequence() works too.  It needs at least 20,000 bases
 * (characters kept, spacers included).  Runs the full Prodigal training pipeline:
 * node creation, GC bias analysis, initial dynamic programming, dicodon
 * statistics, RBS scoring, and start codon weight training.
 *
 * Training starts from a fresh model: the translation table and start
 * weight are the context's (the config's, or the last value given to
 * prodigal_set_translation_table() or prodigal_set_start_weight()), the GC
 * is the training sequence's, and nothing else carries over from a model
 * trained, loaded or adjusted before.  If training fails after checking
 * the sequence, the context has no model until it is trained again or a
 * model is loaded.
 *
 * For the same input and options, the model equals the training file the
 * CLI writes with -t, provided the input is within the limits listed under
 * "Input limits" in docs/internals.md.
 *
 * After training, call prodigal_set_sequence() + prodigal_find_genes()
 * for each contig.  Alternatively, export the model with
 * prodigal_export_training() for later reuse.
 *
 * @param ctx  Context with training sequences loaded.
 * @return PRODIGAL_OK; PRODIGAL_ERR_INVALID_INPUT for a NULL ctx, a
 *         metagenomic context or no sequence loaded;
 *         PRODIGAL_ERR_SEQ_TOO_SHORT under 20,000 bases; PRODIGAL_ERR_NOMEM.
 */
PRODIGAL_API int prodigal_train(prodigal_ctx_t *ctx);

/**
 * Load a pre-trained model from a binary blob.
 *
 * The blob is a model as produced by prodigal_export_training() or the
 * CLI's -t flag: Prodigal's model struct as it lies in memory, 558,392
 * bytes where doubles are 8-byte aligned, in the machine's byte order
 * (docs/internals.md, "The model file").  So a blob moves between x86-64,
 * ARM64 and wasm32 builds, all little-endian, but not to a big-endian
 * machine or to 32-bit x86 Linux, which aligns doubles to 4 bytes.  It
 * need not be aligned.
 * Rejected, leaving the context unchanged: any model in a metagenomic
 * context (prodigal_config_t.meta_mode); a wrong size, an unsupported
 * translation table (docs/genetic-codes.md), a GC outside [0, 1], a start
 * weight that is not finite and > 0, and an SD flag other than 0 or 1.
 * Other values are not checked: they can be NaN
 * (docs/upstream-behavior.md, quirk 11).
 *
 * @param ctx   Context.
 * @param data  Pointer to training data blob.
 * @param len   Size of blob in bytes.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_load_training(prodigal_ctx_t *ctx,
                                        const void *data, size_t len);

/**
 * Export the current training model as a binary blob.
 *
 * The caller receives a malloc'd buffer that must be freed with free().
 * The blob is binary-compatible with Prodigal's -t training file format.
 *
 * @param ctx       Context with a model, from prodigal_train() or
 *                  prodigal_load_training() (a metagenomic context never
 *                  has one).
 * @param data_out  Receives pointer to malloc'd blob.
 * @param len_out   Receives size of blob in bytes.
 * @return PRODIGAL_OK, or, with no message (the context is const),
 *         PRODIGAL_ERR_NOMEM, or PRODIGAL_ERR_INVALID_INPUT when an
 *         argument is NULL or the context has no model.
 */
PRODIGAL_API int prodigal_export_training(const prodigal_ctx_t *ctx,
                                          void **data_out, size_t *len_out);

/**
 * Read a model's settings from a blob without loading it.
 *
 * Accepts exactly the blobs prodigal_load_training() accepts in a
 * single-genome context (the blob of
 * prodigal_export_training() or a CLI -t file; it need not be aligned) and
 * gives the same reasons for the others, so it also checks a model file
 * before any context exists.
 *
 * @param data     The blob.
 * @param len      Its size in bytes.
 * @param info     Receives the settings; unchanged when the blob is rejected.
 * @param msg      Receives the reason, or "" if the blob is accepted; may be
 *                 NULL.  Truncated to msg_len bytes, NUL included.
 * @param msg_len  Size of msg in bytes; messages fit in 256.
 * @return PRODIGAL_OK, or PRODIGAL_ERR_INVALID_INPUT for a NULL data or
 *         info, or a blob prodigal_load_training() rejects.
 */
PRODIGAL_API int prodigal_training_info(const void *data, size_t len,
                                        prodigal_training_info_t *info,
                                        char *msg, size_t msg_len);

/** @name Training parameter setters
 *  Fine-grained control over individual training parameters.  These modify
 *  the current model directly.  Useful for experimentation or for tweaking
 *  a loaded model.  The translation table and start weight also apply to
 *  later training; prodigal_train() replaces the GC and SD usage, so set
 *  those after training or loading.  A metagenomic context refuses them
 *  all, except prodigal_set_translation_table(11) (see
 *  prodigal_config_t.meta_mode).
 */
/**@{*/

/**
 * Set the NCBI translation table.
 * @param table  Supported: those listed at prodigal_config_t.trans_table
 *               (docs/genetic-codes.md), and in a metagenomic context 11
 *               only; the error message says why another is not.
 * @return PRODIGAL_OK or PRODIGAL_ERR_INVALID_INPUT.
 */
PRODIGAL_API int prodigal_set_translation_table(prodigal_ctx_t *ctx, int table);

/**
 * Set the start score weight.
 * @param weight  Must be finite and > 0.  Default: 4.35.
 * @return PRODIGAL_OK, or PRODIGAL_ERR_INVALID_INPUT (NULL ctx, a weight
 *         that is not finite and > 0, or a metagenomic context).
 */
PRODIGAL_API int prodigal_set_start_weight(prodigal_ctx_t *ctx, double weight);

/**
 * Set the GC content in the training model.
 * @param gc  Must be in [0.0, 1.0].
 * @return PRODIGAL_OK, or PRODIGAL_ERR_INVALID_INPUT (NULL ctx, gc outside
 *         [0, 1] or NaN, or a metagenomic context).
 */
PRODIGAL_API int prodigal_set_gc(prodigal_ctx_t *ctx, double gc);

/**
 * Set whether the model uses Shine-Dalgarno motifs.
 * Zero on a model trained with SD motifs scores upstream motifs with
 * weights its training never set (zero).
 * @param uses_sd  Nonzero: use SD.  Zero: use upstream motif scan.
 * @return PRODIGAL_OK, or PRODIGAL_ERR_INVALID_INPUT (NULL ctx, or a
 *         metagenomic context).
 */
PRODIGAL_API int prodigal_set_uses_sd(prodigal_ctx_t *ctx, int uses_sd);
/**@}*/

/**
 * Find genes in the currently loaded sequence (SOA output).
 *
 * In single-genome mode, requires prior training (prodigal_train() or
 * prodigal_load_training()).  In metagenomic mode the first call allocates
 * the 50 built-in models (see prodigal_create()).
 *
 * On success, *genes_out is set to a newly allocated SOA struct.  The caller
 * must free it with prodigal_genes_free().  If stats_out is non-NULL, it is
 * populated with computation statistics.
 *
 * On error, *genes_out is set to NULL (when genes_out isn't NULL), and the
 * sequence and model stay loaded, so the call can be repeated.
 *
 * @param ctx        Context with a sequence loaded.
 * @param genes_out  Receives pointer to SOA output (caller frees).
 * @param stats_out  Receives computation stats, or NULL to skip.
 * @return PRODIGAL_OK; PRODIGAL_ERR_INVALID_INPUT for a NULL ctx or
 *         genes_out, no model in single-genome mode, or no sequence loaded;
 *         PRODIGAL_ERR_NOMEM; PRODIGAL_ERR_CANCELLED (metagenomic mode,
 *         progress_callback); PRODIGAL_ERR_INTERNAL (not expected: protein
 *         letters past 32-bit offsets).
 */
PRODIGAL_API int prodigal_find_genes(prodigal_ctx_t *ctx,
                                     prodigal_genes_soa_t **genes_out,
                                     prodigal_stats_t *stats_out);

/**
 * Find genes in the currently loaded sequence (AOS output).
 *
 * Same as prodigal_find_genes() but returns an Array-of-Structures layout.
 * Free with prodigal_genes_aos_free().
 *
 * @param ctx        Context with a sequence loaded.
 * @param genes_out  Receives pointer to AOS output (caller frees).
 * @param stats_out  Receives computation stats, or NULL to skip.
 * @return As prodigal_find_genes(), and PRODIGAL_ERR_NOMEM for the copy.
 */
PRODIGAL_API int prodigal_find_genes_aos(prodigal_ctx_t *ctx,
                                         prodigal_genes_t **genes_out,
                                         prodigal_stats_t *stats_out);

/**
 * Free a SOA gene output struct.
 * Safe to call with NULL (no-op).  Never uses the custom allocator.
 * @param genes  SOA struct to free, or NULL.
 */
PRODIGAL_API void prodigal_genes_free(prodigal_genes_soa_t *genes);

/**
 * Free an AOS gene output struct.
 * Safe to call with NULL (no-op).  Never uses the custom allocator.
 * @param genes  AOS struct to free, or NULL.
 */
PRODIGAL_API void prodigal_genes_aos_free(prodigal_genes_t *genes);

/**
 * Get a static human-readable string for an error code.
 *
 * Thread-safe; returns a pointer to a string literal.
 *
 * @param error_code  An error code (PRODIGAL_OK, PRODIGAL_ERR_*, etc.).
 * @return Static string like "Success", "Out of memory", "Unknown error".
 */
PRODIGAL_API const char *prodigal_strerror(int error_code);

/**
 * Get a detailed error message from the last failed operation.
 *
 * Returns a pointer to an internal buffer in the context, which the next
 * failing call on it overwrites; valid until prodigal_destroy().  Read it
 * right after a call returns an error: a successful call does not clear
 * it (prodigal_set_sequence() and prodigal_set_training_sequences() do),
 * and the errors of prodigal_get_seq_info() and prodigal_export_training(),
 * whose context is const, set none, leaving an older message.
 * Returns "" on a fresh context with no errors.
 * Returns "NULL context" if ctx is NULL.
 *
 * @param ctx  Context, or NULL.
 * @return Error message string (never NULL).
 */
PRODIGAL_API const char *prodigal_last_error(const prodigal_ctx_t *ctx);

/**
 * The library's ABI version: PRODIGAL_ABI_VERSION of the header it was
 * built with.
 */
PRODIGAL_API int prodigal_abi_version(void);

/**
 * The layout of every public struct as the library was compiled: for each
 * struct its own row (field NULL), then a row per field in declaration
 * order, and a final row whose type is NULL.  A binding that mirrors the
 * structs can compare every row with its own layout (in Rust: offset_of!
 * and size_of), which catches reordered fields as well as added or removed
 * ones, and fields struct_size can't see: those an addition puts in
 * padding.  The table is static and constant.
 *
 * @param n  Receives the number of rows before the final one; may be NULL.
 * @return The table.
 */
PRODIGAL_API const prodigal_abi_field_t *prodigal_abi_layout(size_t *n);

/**
 * Get the library version string at runtime.
 *
 * Useful for FFI consumers that cannot use compile-time macros.
 * Returns a pointer to a string literal (e.g., "2.6.5").
 *
 * @return Static version string.
 */
PRODIGAL_API const char *prodigal_version_string(void);

#ifdef __cplusplus
}
#endif

#endif /* PRODIGAL_API_H */
