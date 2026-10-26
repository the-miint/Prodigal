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

#ifndef BITMAP_H_
#define BITMAP_H_

/* Inline (the fork's): every read of a base tests bits (docs/internals.md) */

/* Test a bit, 0 = not set, 1 = set */
static inline unsigned char test(const unsigned char *bm, int ndx) {
  return (bm[ndx>>3] >> (ndx&0x07)) & 1;
}

/* Clear a bit (set it to 0) */
static inline void clear(unsigned char *bm, int ndx) {
  bm[ndx>>3] &= ~(1 << (ndx&0x07));
}

/* Set a bit to 1 */
static inline void set(unsigned char *bm, int ndx) {
  bm[ndx>>3] |= (1 << (ndx&0x07));
}

/* Flip a bit's value 0->1 or 1->0 */
static inline void toggle(unsigned char *bm, int ndx) {
  bm[ndx>>3] ^= (1 << (ndx&0x07));
}

#endif
