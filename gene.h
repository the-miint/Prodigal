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

#ifndef _GENE_H
#define _GENE_H

#include <stdio.h>
#include "training.h"
#include "node.h"
#include "dprog.h"

struct _gene {
  int begin;               /* Left end of the gene */
  int end;                 /* Right end of the gene */
  int start_ndx;           /* Index to the start node in the nodes file */
  int stop_ndx;            /* Index to the stop node in the nodes file */

  /* The gene as record_gene_data reports it, from the nodes that chose it.
     The library returns these; the CLI prints them, formatted by
     pdg_gene_data and pdg_score_data. */
  int strand;              /* 1 = forward, -1 = reverse */
  int partial_left;        /* Runs off the left edge */
  int partial_right;       /* Runs off the right edge */
  int start_type;          /* 0 = ATG, 1 = GTG, 2 = TTG, 3 = Edge */
  double cscore, sscore, rscore, uscore, tscore;
  double confidence;
  double gc_cont;
  struct _rbs_text rbs;
};

int add_genes(struct _gene *, int, struct _node *, int);
void record_gene_data(struct _gene *, int, struct _node *, struct _training *);

/* The bytes of pdg_gene_data's and pdg_score_data's buffers */
#define PDG_GENE_TEXT 500
void pdg_gene_data(char *, const struct _gene *, int, int);
void pdg_score_data(char *, const struct _gene *);
void tweak_final_starts(struct _gene *, int, struct _node *, int, struct
                       _training *);

void print_genes(FILE *, struct _gene *, int, struct _node *, int, int, int,
                 int, char *, struct _training *, char *, char *, char *);
/* Each codon's protein letter under one translation table, by the codon's
   mer_ndx(3, ...) index: aa[0] inside a gene, aa[1] as the first codon of a
   gene that is not edge-started (M at a start). */
struct _codon_letters {
  int table;                  /* the trans_table aa holds; 0: none yet */
  char aa[2][64];
};

int pdg_protein_length(const struct _gene *);
void pdg_codon_letters(struct _training *, struct _codon_letters *);
char pdg_protein_letter(const struct _gene *, int, unsigned char *,
                        unsigned char *, unsigned char *, int,
                        const struct _codon_letters *);
void write_translations(FILE *, struct _gene *, int, unsigned char *,
                        unsigned char *, unsigned char *, int,
                        struct _training *, struct _codon_letters *, int,
                        char *);
void write_nucleotide_seqs(FILE *, struct _gene *, int, struct _node *, 
                           unsigned char *, unsigned char *, unsigned char *,
                           int, struct _training *, int, char *);
double calculate_confidence(double, double);

#endif
