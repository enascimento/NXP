/* ===================================================================== */
/* This file is part of Daredevil                                        */
/* Daredevil is a side-channel analysis tool                             */
/* Copyright (C) 2016                                                    */
/* Original author:   Paul Bottinelli <paulbottinelli@hotmail.com>       */
/* Contributors:      Joppe Bos <joppe_bos@hotmail.com>                  */
/*                                                                       */
/* This program is free software: you can redistribute it and/or modify  */
/* it under the terms of the GNU General Public License as published by  */
/* the Free Software Foundation, either version 3 of the License, or     */
/* any later version.                                                    */
/*                                                                       */
/* This program is distributed in the hope that it will be useful,       */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/* GNU General Public License for more details.                          */
/*                                                                       */
/* You should have received a copy of the GNU General Public License     */
/* along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/* ===================================================================== */
#ifndef CPA_H
#define CPA_H

#include <stdio.h>
#include <stdint.h>
#include <iostream>
#include "memUtils.h"
#include "aes.h"
#include "des.h"


#define ALG_AES                 0
#define ALG_DES                 1
/*
#define ALG_DES_AFTER           2
#define ALG_DES_BEFORE_SMALL    3
#define ALG_DES_AFTER_SMALL     4
*/

/* Given the messages stored in m, use the bytenum-th byte to construct
 * the guesses for round R at position pos for algorithm alg and store the guesses in guess.
 */
template <class TypeGuess> int construct_guess (TypeGuess ***guess, uint32_t alg, Matrix *m, uint32_t n_m, uint32_t bytenum, uint32_t R, uint32_t pos, uint16_t * sbox, uint32_t n_keys, int8_t bit);

/* Stucture to store a first order correlation element to be put in the
 * priority queue. Such an element is defined by its correlation, and
 * the two time sample and key that led to this correlation.
 */
template <typename Type>
struct CorrSecondOrder {

  Type corr;
  int time1;
  int time2;
  int key;

  CorrSecondOrder(Type c, int t1, int t2, int k) : corr(c), time1(t1), time2(t2), key(k) {
  }

  CorrSecondOrder() : corr(0), time1(0), time2(0), key(0) {
  }

  bool operator<(const struct CorrSecondOrder<Type> & other) const {
    return fabs(this->corr) < fabs(other.corr);
  }

  /* Not really correct in a logical PoV, but this is to get the rank
   * of the highest correct key in the PriorityQueue.
   */
  bool operator==(const int other_key) const {
    return this->key == other_key;
  }

  friend std::ostream& operator<<( std::ostream& out, const CorrSecondOrder& b ){
    return out << setw(16) << b.corr << setw(6) << "0x" << setw(4) << left << hex << b.key << right << setw(8) << dec << b.time1 <<setw(8) << dec << b.time2;
  }

  void corr2str(string sep){
    cout << time1 << sep << time2 << sep << "0x" << hex << key << dec << sep << corr << endl;
  }

};

/* Stucture to store a first order correlation element to be put in the
 * priority queue. Such an element is defined by its correlation, and
 * the time sample and key that led to this correlation.
 */
template <typename Type>
struct CorrFirstOrder {

  Type corr;
  int time;
  int key;

  CorrFirstOrder() : corr(0), time(0), key(0) {
  }

  CorrFirstOrder(Type c, int t, int k) : corr(c), time(t), key(k) {
  }

  bool operator<(const struct CorrFirstOrder & other) const {
    return fabs(this->corr) < fabs(other.corr);
  }

  bool operator==(const int other_key) const {
    return this->key == other_key;
  }

  friend std::ostream& operator<<( std::ostream& out, const CorrFirstOrder& b ){
    return out << setw(16) << b.corr << setw(6) << "0x" << setfill('0') << setw(2) << hex << b.key << setw(6) << setfill(' ') << right << setw(8) << dec << b.time;
  }

  void corr2str(string sep){
    cout << time << sep << "0x" << hex << key << dec << sep << corr << endl;
  }

};

template <typename Type>
struct SecondOrderQueues {
  PriorityQueue<CorrSecondOrder<Type> > * pqueue;
  CorrSecondOrder<Type> * top_corr;

  SecondOrderQueues(PriorityQueue<CorrSecondOrder<Type> > * q, CorrSecondOrder<Type> * t):
    pqueue(q), top_corr(t) {
    }
};


template <typename Type>
struct FirstOrderQueues {
  PriorityQueue<CorrFirstOrder<Type> > * pqueue;
  CorrFirstOrder<Type> * top_corr;

  FirstOrderQueues(PriorityQueue<CorrFirstOrder<Type> > * q, CorrFirstOrder<Type> * t):
    pqueue(q), top_corr(t) {
    }
};

template <typename TypeTrace>
struct PrecompTraces {

  int start;
  int end;
  int length;
  TypeTrace ** trace;

  PrecompTraces(int st, int en, int nt, TypeTrace ** tr):
    start(st), end(en), length(nt), trace(tr) {
  }
};

template <typename TypeTrace>
struct PrecompTracesNorm {

  int start;
  int end;
  int length;
  uint8_t exponent;
  TypeTrace ** trace;

  PrecompTracesNorm(int st, int en, int nt, int ex, TypeTrace ** tr):
    start(st), end(en), length(nt), exponent(ex), trace(tr) {
  }
};

template <typename TypeGuess, typename TypeReturn>
struct PrecompGuesses {

  int start;
  int end;
  int n_traces;
  TypeGuess ** guess;
  TypeReturn ** precomp_k;

  PrecompGuesses(int st, int en, int n_t, TypeGuess ** gu, TypeReturn ** pk):
    start(st), end(en), n_traces(n_t), guess(gu), precomp_k(pk) {
  }
};


/* This function precomputes the mean for the traces and subtract this mean
 * from every element of the traces. This is to be used by the newer v_5 of
 * SOCPA.
 */
  template <class TypeTrace, class TypeReturn>
void * precomp_traces_v_2(void * args_in);

/* This functions simply splits the total work (n_rows) into an equal number of
 * threads, creates this amount of threads and starts them to precompute the
 * distance of means for each row of the matrix trace.
 * ! We expect a matrix where the number of traces is n_rows
 */
  template <class TypeTrace, class TypeReturn>
int p_precomp_traces(TypeTrace ** trace, int n_rows, int n_columns, int n_threads, int offset=0);

template <class TypeTrace, class TypeReturn>
void * precomp_traces_v_2_norm(void * args_in);

template <class TypeTrace, class TypeReturn>
int p_precomp_traces_norm(TypeTrace ** trace, int n_rows, int n_columns, int n_threads,  uint8_t exponent, int offset=0);

#endif
