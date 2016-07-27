/* ===================================================================== */
/* This file is part of Daredevil                                        */
/* Daredevil is a side-channel analysis tool                             */
/* Copyright (C) 2016                                                    */
/* Original author:   Paul Bottinelli <paulbottinelli@hotmail.com>       */
/* Contributors:      Joppe Bos <joppe_bos@hotmail.com>                  */
/*                    Philippe Teuwen <phil@teuwen.org>                  */
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


/* This file deals with general function that help with memory management
 * which are used throughout the rest of the code.
 */


#include <typeinfo>
#include <fstream>
#include <sstream>
#include <string>
#include "memUtils.h"

/* Returns the number of columns that can be loaded in memory.
 */
  template <typename Type>
int get_ncol(long int memsize, int ntraces)
{
  // We use 60% of the available memory. We never know what can happen :)
  return (0.6*memsize)/(sizeof(Type)*ntraces);
}


  template <class Type>
void free_matrix(Type *** matrix, int n_rows)
{
  for (int i=0; i < n_rows; i++) {
    free((*matrix)[i]);
  }
  free(*matrix);
}

/* Allocates the array matrix
 */
  template <class Type>
int allocate_matrix(Type *** matrix, int n_rows, int n_columns)
{
  *matrix = (Type **)malloc(n_rows * sizeof(Type *));
  if(*matrix == NULL)
    return -1;

  for (int i=0; i < n_rows; i++) {
    (*matrix)[i] = (Type *) malloc (n_columns * sizeof(Type));
    if((*matrix)[i] == NULL)
      return -1;
  }
  return 0;
}

/* This functions simply splits the total_work (usually represents the number
 * of columns of the matrix we're processing) into an equal number of threads,
 * creates this amount of threads and starts them with the function fct.
 */
  template <class TypeTrace, class TypeReturn, class TypeGuess>
int split_work(FinalConfig<TypeTrace, TypeReturn, TypeGuess> & fin_conf, void * (*fct)(void *), TypeReturn ** precomp_k, int total_work, int offset)
{
  int n, rc,
      workload = 0,
      n_threads = fin_conf.conf->n_threads,
      /* Can be changed later in order to compute on less traces.
       */
      n_traces = fin_conf.conf->total_n_traces;


  /* If the total work by thread is smaller than 1, only the last thread would
   * work, which is against the sole principle of multithreading. Thus, we
   * reduce the number of threads until the workload is larger than 1.
   *
   * This is quite a naive approach, and it would be better to look into more
   * efficient load balancing algorithms.
   *
   */
  workload = (total_work/n_threads);
  while (workload < 1) {
    n_threads -= 1;
    workload = (total_work/n_threads);
  }

  pthread_t threads[n_threads];
  General<TypeTrace, TypeReturn, TypeGuess> *ta = NULL;

  ta = (General<TypeTrace, TypeReturn, TypeGuess> *) malloc(n_threads * sizeof(General<TypeTrace, TypeReturn, TypeGuess>));
  if (ta == NULL) {
    fprintf (stderr, "[ERROR] Memory alloc failed.\n");
    return -1;
  }

  for (n = 0; n < n_threads; n++) {
    //printf(" Thread_%i [%i-%i]\n", n, n*workload + offset, offset + n*workload + workload + ((n + 1) / n_threads)*(total_work % n_threads));
    ta[n] = General<TypeTrace, TypeReturn, TypeGuess>(n*workload, workload + ((n + 1) / n_threads) * (total_work % n_threads), n_traces, offset, total_work, precomp_k, &fin_conf);
    rc = pthread_create(&threads[n], NULL, (*fct), (void *) &ta[n]);
    if (rc != 0) {
      fprintf(stderr, "[ERROR] Creating thread.\n");
      free (ta);
      return -1;
    }
  }

  for (n = 0; n < n_threads; n++) {
    rc = pthread_join(threads[n], NULL);
    if (rc != 0) {
      fprintf(stderr, "[ERROR] Joining thread.\n");
      free (ta);
      return -1;
    }
  }
  free (ta);
  return 0;
}

/* This function precomputes the sum and the sum of squares for all guesses
 * which will later be used in the correlation computation.
 */
  template <class TypeTrace, class TypeReturn, class TypeGuess>
void * precomp_guesses(void * args_in)
{
  int i, j;
  TypeReturn tmp;
  General<TypeTrace, TypeReturn, TypeGuess> * G = (General<TypeTrace, TypeReturn, TypeGuess> *) args_in;

  for (i = G->start; i < G->start + G->length; i++) {
    for (j = 0; j < G->n_traces; j++) {
      tmp = G->fin_conf->mat_args->guess[i][j];
      G->precomp_guesses[i][0] += tmp;
      G->precomp_guesses[i][1] += tmp*tmp;
    }
  }
  return NULL;
}

template int get_ncol<int8_t>(long int memsize, int ntraces);
template int get_ncol<float>(long int memsize, int ntraces);
template int get_ncol<double>(long int memsize, int ntraces);

template void free_matrix(float *** matrix, int n_rows);
template void free_matrix(double *** matrix, int n_rows);
template void free_matrix(uint8_t *** matrix, int n_rows);
template void free_matrix(int8_t *** matrix, int n_rows);
template void free_matrix(int *** matrix, int n_rows);

template int allocate_matrix(float *** matrix, int n_rows, int n_columns);
template int allocate_matrix(double *** matrix, int n_rows, int n_columns);
template int allocate_matrix(uint8_t *** matrix, int n_rows, int n_columns);
template int allocate_matrix(int8_t *** matrix, int n_rows, int n_columns);
template int allocate_matrix(int *** matrix, int n_rows, int n_columns);

template int split_work<double, double, uint8_t>(FinalConfig<double, double, uint8_t> & fin_conf, void * (*fct)(void *), double ** precomp_k, int total_work, int offset);
template int split_work<float, double, uint8_t>(FinalConfig<float, double, uint8_t> & fin_conf, void * (*fct)(void *), double ** precomp_k, int total_work, int offset);
template int split_work<int8_t, double, uint8_t>(FinalConfig<int8_t, double, uint8_t> & fin_conf, void * (*fct)(void *), double ** precomp_k, int total_work, int offset);
template int split_work<float, float, uint8_t>(FinalConfig<float, float, uint8_t> & fin_conf, void * (*fct)(void *), float ** precomp_k, int total_work, int offset);
template int split_work<int8_t, float, uint8_t>(FinalConfig<int8_t, float, uint8_t> & fin_conf, void * (*fct)(void *), float ** precomp_k, int total_work, int offset);

template void * precomp_guesses<int8_t, double, uint8_t>(void * args_in);
template void * precomp_guesses<float, double, uint8_t>(void * args_in);
template void * precomp_guesses<double, double, uint8_t>(void * args_in);
template void * precomp_guesses<int8_t, float, uint8_t>(void * args_in);
template void * precomp_guesses<float, float, uint8_t>(void * args_in);
