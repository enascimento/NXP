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

#ifndef MEMUTILS_H
#define MEMUTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include "utils.h"

#ifndef RESOURCES
#define RESOURCES "/usr/share/daredevil"
#endif //RESOURCES

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define RESET "\033[0m"

#define GIGA 1e9
#define MEGA 1e6

using namespace std;

/* Structure used to stores all the arrays.
 */
template <typename TypeTrace, typename TypeReturn, typename TypeGuess>
struct MatArgs {

  TypeTrace ** trace;
  TypeGuess ** guess;
  TypeReturn ** results;

  MatArgs(TypeTrace ** tr, TypeGuess ** gues, TypeReturn ** res):
    trace(tr), guess(gues), results(res) {
    }
};

/* Structures used to stores all the *common* information needed by the
 * threads to compute the correlation, like the pointers to the traces
 * and guesses arrays and their sizes.
 */
template <typename TypeTrace, typename TypeReturn, typename TypeGuess>
struct Args {

  TypeTrace ** trace;
  int n_samples;
  TypeGuess ** guess;
  int n_keys;
  TypeReturn ** results;
  int n_traces;
  int nsqr;

  Args(TypeTrace ** tr, int n_s, TypeGuess ** gues, int nk, TypeReturn ** res, int nt, int ns):
    trace(tr), n_samples(n_s), guess(gues), n_keys(nk), results(res), n_traces(nt), nsqr(ns) {
    }
};

/* Structure used by threads to access the common information, and the
 * indices in the array of traces where every individual thread will start
 * and stop computing correlations.
 */
template <typename TypeTrace, typename TypeReturn, typename TypeGuess>
struct ThreadArgs {

    Args<TypeTrace, TypeReturn, TypeGuess> * args;
    int start;
    int length;

    ThreadArgs(Args<TypeTrace, TypeReturn, TypeGuess> * a, int st, int len):
      args(a), start(st), length(len) {
    }
};

/* Homemade Priority Queue used to store the best correlations
 */
template <typename Type>
class PriorityQueue
{
  Type * array;
  int size, max_size, index_min, total;

  public:
  PriorityQueue(int s)
  {
    init(s);
  }
  PriorityQueue(){}

  void init(int s)
  {
    max_size = s;
    array = (Type *) malloc(max_size * sizeof(Type));
    size = 0;
    index_min = 0;
    total = 0;
  }
  void insert(const Type& elem)
  {
    if (size == max_size) {
      if(array[index_min] < elem) {
        array[index_min] = elem;
        update_smallest_ind();
      }
    }else {
      array[size] = elem;
      if (elem < array[index_min]) {
        index_min = size;
      }
      size += 1;
    }
    total++;
  }

  void print(int length = -1, int key = -1)
  {
    int i;
    uint8_t seen_key = 0;

    if (length == -1 || length > size)
      length = size;

    cout << "[INFO]\t" << total <<" correlations computed in total." << endl;
    cout << "[INFO]\tGlobal top " << length << " correlations." << endl;
    sort(array, array + size);
    update_smallest_ind();
    cout << endl;
    for (i = size - 1; i >= (size-length); i--) {
      if (array[i] == key){
        seen_key = 1;
        cout << KGRN << array[i] << RESET << endl;
      }else
        cout << array[i] << endl;
    }

    if(length != size && !seen_key){
      while (i > 0){
        if (array[i] == key){
          seen_key = 1;
          for(int j = 0; j < 3; j++)
            cout << setw(13) << '.' << setw(10) << '.' <<  setw(10) << '.' <<setw(8) << '.' << endl;
          cout << KGRN << array[i] << RESET << "\tat rank " << size-i << "." << endl;
          break;
        }
        i--;
      }
    }
    if (!seen_key && key != -1 && size != 0){
      cout << endl;
      cout << "Key 0x" << hex << key << " does not appear in the top " << dec << size << " correlations." << endl;
    }
    cout << endl;
  }
  private:
  void update_smallest_ind()
  {
    int i;
    for (i = 0; i < size; i++) {
      if (array[i] < array[index_min]) {
        index_min = i;
      }
    }
  }

};

/* Structure used to store ALL the general and common information
 */
template <typename TypeTrace, typename TypeReturn, typename TypeGuess>
struct FinalConfig {

  MatArgs<TypeTrace, TypeReturn, TypeGuess> * mat_args;
  Config * conf;
  void * queues;

  FinalConfig(MatArgs<TypeTrace, TypeReturn, TypeGuess> * m_a, Config * c, void * q):
    mat_args(m_a), conf(c), queues(q){
    }
};

//The files from socpa.h start here
template <typename TypeTrace, typename TypeReturn, typename TypeGuess>
struct General {

  int start;
  int length;
  int n_traces;
  int global_offset;
  /* The number of columns that we can store in memory. Needed to know when to
   * stop computing correlations in the last slice when computing on big files.
   */
  int n_samples;
  TypeReturn ** precomp_guesses;
  FinalConfig<TypeTrace, TypeReturn, TypeGuess> * fin_conf;

  General(int st, int len, int nt, int go, int nc, TypeReturn ** pg, FinalConfig<TypeTrace, TypeReturn, TypeGuess> * s):
    start(st), length(len), n_traces(nt), global_offset(go), n_samples(nc), precomp_guesses(pg), fin_conf(s){
  }
};

/* Frees a matrix
 */
template <class Type>
void free_matrix(Type *** matrix, int n_rows);

/* Allocates memory for a matrix
 */
template <class Type>
int allocate_matrix(Type *** matrix, int n_rows, int n_columns);

template <typename Type>
int get_ncol(long int memsize, int ntraces);

template <class TypeTrace, class TypeReturn, class TypeGuess>
int split_work(FinalConfig<TypeTrace, TypeReturn, TypeGuess> & fin_conf, void * (*fct)(void *), TypeReturn ** precomp_k, int total_work, int offset=0);

template <class TypeTrace, class TypeReturn, class TypeGuess>
void * precomp_guesses(void * args_in);

#endif
