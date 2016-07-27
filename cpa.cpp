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

/* This file is used in every correlation power attack. It has functions
 * like split_work or precomp_traces that are used in focpa and socpa.
 * Also all the used structures in focpa and socpa can be found in cpa.h.
 */

#include "cpa.h"
#include "omp.h"

/* Given the messages stored in m, use the bytenum-th byte to construct
 * the guesses for round R for algorithm alg and store the guesses in guess.
 * des_switch is only used by DES
 */
template <class TypeGuess>
int construct_guess (TypeGuess ***guess, uint32_t alg, Matrix *m, uint32_t n_m, uint32_t bytenum, uint32_t R, uint32_t des_switch, uint16_t * sbox, uint32_t n_keys, int8_t bit) {
  int ret;

  switch (alg) {
    case ALG_AES:
      ret = construct_guess_AES (guess, m, n_m, bytenum, R, sbox, n_keys, bit); //construct_guess_AES can be found in aes.cpp
      if (ret < 0) return -1;
      break;
    case ALG_DES:
      ret = construct_guess_DES (guess, m, n_m, bytenum, R, des_switch, sbox, n_keys, bit); //construct_guess_DES can be found in des.cpp
      if (ret < 0) return -1;
      break;
    default:
      fprintf (stderr, "Algorithm is not supported (yet).\n");
      return -1;
  }
  return 1;

}

/* This functions simply splits the total work (n_rows) into an equal number of
 * threads, creates this amount of threads and starts them to precompute the
 * distance of means for each row of the matrix trace. If the offset value is
 * specified, we start splitting the work starting at offset.
 *
 * ! We expect a matrix where the number of traces is n_rows
 */
  template <class TypeTrace, class TypeReturn>
int p_precomp_traces(TypeTrace ** trace, int n_rows, int n_columns, int n_threads, int offset/*, int n_traces_from_offset*/)
{
  int n, rc,
      workload = 0,
      n_traces = n_columns;

  //printf("Offset: %i\n", offset);

  /* If the total work by thread is smaller than 1, only the last thread would
   * work, which is against the sole principle of multithreading. Thus, we
   * reduce the number of threads until the workload is larger than 1.
   */
  workload = ((n_rows-offset)/n_threads);
  while (workload < 1) {
    n_threads -= 1;
    workload = ((n_rows-offset)/n_threads);
  }

  pthread_t threads[n_threads];
  PrecompTraces<TypeTrace> *ta = NULL;

  ta = (PrecompTraces<TypeTrace>*) malloc(n_threads * sizeof(PrecompTraces<TypeTrace>));
  if (ta == NULL) {
    fprintf (stderr, "[ERROR] Memory alloc failed.\n");
    return -1;
  }

  for (n = 0; n < n_threads; n++) {
    //printf(" Thread_%i [%i-%i]\n",n , offset+ n*workload, offset+n*workload + workload + ((n + 1) / n_threads)*(n_rows % n_threads));
    ta[n] = PrecompTraces<TypeTrace>(offset + n*workload, workload + ((n + 1) / n_threads) * (n_rows % n_threads), n_traces, trace);
    rc = pthread_create(&threads[n], NULL, precomp_traces_v_2<TypeTrace, TypeReturn>, (void *) &ta[n]);
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

/* This function precomputes the mean for the traces and subtract this mean
 * from every element of the traces. This is to be used by the newer v_5 of
 * SOCPA.
 */
  template <class TypeTrace, class TypeReturn>
void * precomp_traces_v_2(void * args_in)
{

  int i, j;
  TypeReturn mean = 0.0;

  PrecompTraces<TypeTrace> * G = (PrecompTraces<TypeTrace> *) args_in;

  for (i = G->start; i < G->start + G->end; i++) {
    for (j = 0; j < G->length; j++) {
      mean += G->trace[i][j];
    }
    mean /= G->length;
    for (j = 0; j < G->length; j++) {
      G->trace[i][j] -= mean;
    }
  }
  return NULL;
}

template int construct_guess (uint8_t ***guess, uint32_t alg, Matrix *m, uint32_t n_m, uint32_t bytenum, uint32_t R, uint32_t des_switch, uint16_t * sbox, uint32_t n_keys, int8_t bit);

template int p_precomp_traces<int8_t, double>(int8_t ** trace, int n_rows, int n_columns, int n_threads, int offset);
template int p_precomp_traces<double, double>(double ** trace, int n_rows, int n_columns, int n_threads, int offset);
template int p_precomp_traces<float, float>(float ** trace, int n_rows, int n_columns, int n_threads, int offset);
