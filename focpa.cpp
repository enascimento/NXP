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

#include "cpa.h"
#include "string.h"
#include "focpa.h"
#include "memUtils.h"
#include "des.h"

pthread_mutex_t lock;

/* Implements first order CPA in a faster and multithreaded way on big files,
 * using the vertical partitioning approach.
 */
  template <class TypeTrace, class TypeReturn, class TypeGuess>
int first_order(Config & conf)
{

  long int memory = conf.memory;

  double start, end = 0;

  int res,
      n_keys = conf.total_n_keys,
      n_samples = conf.n_samples,
      nmat = conf.n_file_trace,
      nrows = conf.total_n_traces,
      ncol = min(get_ncol<TypeTrace>(memory-(nrows*n_keys*sizeof(TypeGuess)), nrows), n_samples), //get_ncol is defined in memUtils.cpp
      col_incr = ncol,
      col_offset = 0,
      row_offset = 0,
      sample_offset = 0,
      cur_n_rows, cur_n_cols,
      samples_loaded = 0,
      to_load = ncol;

  uint8_t is_last_iter = 0;
  unsigned int max_n_rows = 0;


  TypeTrace ** traces = NULL;
  TypeTrace ** tmp = NULL;
  TypeGuess ** guesses = NULL;
  TypeReturn ** precomp_k;

  if (col_incr <= 0) {
    fprintf(stderr, "[ERROR] Invalid parameters ncol(=%i).\n", ncol);
    return -1;
  }

  /* We determine the size of the file having the largest number of rows, to
   * allocate memory for tmp.
   */
  for (int i = 0; i < nmat; i++){
    if(conf.traces[i].n_rows > max_n_rows)
      max_n_rows = conf.traces[i].n_rows;
  }

  res = allocate_matrix(&tmp, max_n_rows, ncol); //allocate_matrix can be found in memUtils.cpp
  if (res != 0) {
    fprintf (stderr, "[ERROR] Allocating matrix in focpa vp.\n");
    return -1;
  }

  res = allocate_matrix(&traces, ncol, nrows); //allocate_matrix can be found in memUtils.cpp
  if (res != 0) {
    fprintf (stderr, "[ERROR] Allocating matrix in focpa vp.\n");
    return -1;
  }

  res = allocate_matrix(&precomp_k, n_keys, 2); //allocate_matrix can be found in memUtils.cpp
  if (res != 0){
    fprintf(stderr, "[ERROR] Memory allocation failed in focpa vp\n");
    return -1;
  }

  /* We initialize the priority queues to store the highest correlations.
   */
  PriorityQueue<CorrFirstOrder <TypeReturn> > * pqueue = new PriorityQueue<CorrFirstOrder <TypeReturn> >; //PriorityQueue can be found in memUtils.h
  (*pqueue).init(conf.top);

  CorrFirstOrder <TypeReturn> * top_r_by_key; //CorrFirstOrder can be found in cpa.h

  /* If we initialize with malloc, the default constructor is not called,
   * leading to possible issued when inserting/comparing elements.
   */
  top_r_by_key = new CorrFirstOrder <TypeReturn> [n_keys];
  if (top_r_by_key == NULL){
    fprintf(stderr, "[ERROR] Allocating memory for top correlations.\n");
    return -1;
  }


  MatArgs<TypeTrace, TypeReturn, TypeGuess> mat_args = MatArgs<TypeTrace, TypeReturn, TypeGuess> (traces, guesses, NULL); //MatArgs can be found in memUtils.h

  FirstOrderQueues<TypeReturn>* queues = new FirstOrderQueues<TypeReturn>(pqueue, top_r_by_key); //FirstOrderQueues can be found in cpa.h
  if(queues == NULL){
    fprintf(stderr, "[ERROR] Allocating memory for the priority queues.\n");
    return -1;
  }

  FinalConfig<TypeTrace, TypeReturn, TypeGuess> fin_conf = FinalConfig<TypeTrace, TypeReturn, TypeGuess>(&mat_args, &conf, (void*)queues); //FinalConfig can be found in memUtils.h
  pthread_mutex_init(&lock, NULL);

  /* We loop over all the key bytes.
   */

  for (int bn = 0; bn < conf.key_size; bn++){
    ostringstream best_out;
    int lowest_rank = 16;
    double sum_bit_cor[256] = {0};
    double peak_bit_cor[256] = {0};

    /* We keep time for each key byte individually;
     */
    start = omp_get_wtime();

    if (conf.key_size == 1){
      bn = conf.bytenum;
    }
    else if (conf.bytenum != -1 && conf.bytenum != bn){
      continue;
    }

    if (conf.sep == "") printf("[ATTACK] Key byte number %i\n\n", bn);
    else if (conf.key_size > 1) printf("%i%s", bn, conf.sep.c_str());

    /* Potentially attack each bit individually. */
    int bitsperbyte = 0;

    if (conf.algo == ALG_AES) bitsperbyte = 8;
    else if (conf.algo == ALG_DES) bitsperbyte = 4;

    for (int bit=0; bit >= 0 && bit < bitsperbyte; bit=(bit!=-1)?bit+1:bit) {
      if (conf.bitnum == -2) bit = -1;
      else if (conf.bitnum >= 0 && conf.bitnum != bit) continue;

      if (conf.bitnum != -2) {
        if (conf.sep == "") printf("[ATTACK] Target bit number %i\n\n", bit);
        else if (conf.key_size > 1) printf("%i%s", bit, conf.sep.c_str());
      }

      res = construct_guess (&fin_conf.mat_args->guess, conf.algo, conf.guesses, conf.n_file_guess, bn, conf.round, conf.des_switch, conf.sbox, conf.total_n_keys, bit); // construct_guess can be found in cpa.cpp
      if (res < 0) {
        fprintf (stderr, "[ERROR] Constructing guess.\n");
        return -1;
      }

      res = split_work(fin_conf, precomp_guesses<TypeTrace, TypeReturn, TypeGuess>, precomp_k, n_keys); //split_work can be found in memUtils.cpp
      if (res != 0) {
        fprintf(stderr, "[ERROR] Precomputing sum and sum of square for the guesses.\n");
        return -1;
      }


      /* We iterate over the all the files, loading ncol columns to memory at a
       * time.
       */
      while (!is_last_iter) {

        /* If the number of samples loaded so far + what we will load in this
         * iteration is larger than the number of samples, it's the last iter.
         */
        if (samples_loaded + to_load >= n_samples){
          is_last_iter = 1;
          to_load = n_samples - samples_loaded;
        }

        /* We iterate over all the files, loading to_load samples at a time and
         * starting at offset 'conf.index_sample + sample_offset + row_offset'
         * in the files. This offset depends on the iteration and the variable
         * to_load depends on whether it is the first iteration or not
         * (we have to load more in the first iteration)
         */
        for (int i = 0; i < nmat; i++){
          cur_n_rows = conf.traces[i].n_rows;
          cur_n_cols = conf.traces[i].n_columns;

          res = load_file_v_1(conf.traces[i].filename, &tmp, cur_n_rows, to_load, conf.index_sample + sample_offset + row_offset, cur_n_cols); //load_file_v_1 can be found in utils.cpp
          if (res != 0) {
            fprintf (stderr, "[ERROR] Loading file.\n");
            return -1;
          }

          /* We copy the array tmp in the array traces at the good offset, and we
           * transpose it at the same time.
           * row_offset is used to make the distinction between the first iteration
           * and the following.
           */
          for (int j = 0; j < cur_n_rows; j++){
            for (int k = 0; k < to_load; k++){
              fin_conf.mat_args->trace[k + row_offset][j + col_offset] = tmp[j][k];
            }
          }
          col_offset += conf.traces[i].n_rows;
        }

        samples_loaded += to_load;

        /* We set to_load to col_incr. So that only in the very first iteration
         * we load ncol.
         */
        to_load = col_incr;

        /* Same principle for row_offset
         */
        row_offset = 0;

        col_offset = 0;
        if(conf.attack_moment == 0)
          res = split_work(fin_conf, correlation_first_order<TypeTrace, TypeReturn, TypeGuess>, precomp_k, is_last_iter ? (n_samples - sample_offset) : col_incr, sample_offset); //split_work can be found in memUtils.cpp, correlation_first_order is defined below
        else
          res = split_work(fin_conf, higher_moments_correlation<TypeTrace, TypeReturn, TypeGuess>, precomp_k, is_last_iter ? (n_samples - sample_offset) : col_incr, sample_offset); //split_work can be found in memUtils.cpp, higher_moments_correlation is defined below
        if (res != 0) {
          fprintf(stderr, "[ERROR] Computing correlations.\n");
          return -1;
        }
        sample_offset += col_incr;

        /* If we are at the last iteration at that point, no need to do more
         * work.
         */
        if (is_last_iter)
          break;
      }

      /* Warning, when using DES, the correct key doesn't correspond to the actual
       * good key, as we are predicting the input state based on a round key.
       */
      int correct_key;
      if (conf.key_size == 1) {
        if (conf.des_switch == DES_4_BITS && conf.correct_key != -1) correct_key = get_4_middle_bits(conf.correct_key); //get_4_middle_bits can be found in des.cpp
        else correct_key = conf.correct_key;
        pqueue->print(conf.top, correct_key);
        print_top_r(top_r_by_key, n_keys, correct_key); //print_top_r can be found in utils.cpp
      } else {
        if (conf.des_switch == DES_4_BITS) correct_key = get_4_middle_bits(conf.complete_correct_key[bn]); //get_4_middle_bits can be found in des.cpp
        else correct_key = conf.complete_correct_key[bn];

        if (conf.bitnum == -1) {
          sort(top_r_by_key, top_r_by_key + n_keys);
          for (int i = n_keys - 1; i >= 0; i--) {
            if (top_r_by_key[i] == correct_key) {
              if (n_keys - i - 1 < lowest_rank) {
                lowest_rank = n_keys - i - 1;
                best_out.str(std::string());  /* Clear best_out. */
                best_out << "Best bit: " << bit << " rank: " << n_keys - i - 1 << "." << setw(-2) << top_r_by_key[i] << endl;
              }
            }
          }
        } else {
          print_top_r(top_r_by_key, n_keys, correct_key, conf.sep); //print_top_r can be found in utils.cpp
        }
      }

      int key_guess_used[256] = {0};
      for (int i = n_keys - 1; i >= 0; i--) {
        if (key_guess_used[top_r_by_key[i].key] == 0) {
          key_guess_used[top_r_by_key[i].key] = 1;
          sum_bit_cor[top_r_by_key[i].key] += abs(top_r_by_key[i].corr);
          if (abs(top_r_by_key[i].corr) > peak_bit_cor[top_r_by_key[i].key])
            peak_bit_cor[top_r_by_key[i].key] = abs(top_r_by_key[i].corr);
        }
      }

      if (bit == bitsperbyte-1) {
        int nbest=10; // TODO: make it a config parameter
        double sum_bit_cor_sort[256];
        memcpy (sum_bit_cor_sort, sum_bit_cor, n_keys*sizeof(double));
        sort (sum_bit_cor_sort, sum_bit_cor_sort + n_keys);
        cout << "Best " << nbest << " candidates for key byte #" << bn << " according to sum(abs(bit_correlations)):" << endl;
        for (int best=n_keys-1; best > n_keys -1 - nbest;) {
          for (int i=0; (i < n_keys) && (best > n_keys -1 - nbest); i++) {
            if (sum_bit_cor_sort[best] == sum_bit_cor[i]) {
              cout << setfill(' ') << setw(2) << n_keys -1 - best << ": 0x" << setfill('0') << setw(2) << hex << i;
              cout << setfill(' ') << dec << "  sum: " << setw(8) << left << sum_bit_cor_sort[best] << right;
              if (i == correct_key)
                cout << "  <==";
              cout << endl;
              best--;
            }
          }
        }
        cout << endl;
        double peak_bit_cor_sort[256];
        memcpy (peak_bit_cor_sort, peak_bit_cor, n_keys*sizeof(double));
        sort (peak_bit_cor_sort, peak_bit_cor_sort + n_keys);
        cout << "Best " << nbest << " candidates for key byte #" << bn << " according to highest abs(bit_correlations):" << endl;
        for (int best=n_keys-1; best > n_keys -1 - nbest;) {
          for (int i=0; (i < n_keys) && (best > n_keys -1 - nbest); i++) {
            if (peak_bit_cor_sort[best] == peak_bit_cor[i]) {
              cout << setfill(' ') << setw(2) << n_keys -1 - best << ": 0x" << setfill('0') << setw(2) << hex << i;
              cout << setfill(' ') << dec << " peak: " << setw(8) << left << peak_bit_cor_sort[best] << right;
              if (i == correct_key)
                cout << "  <==";
              cout << endl;
              best--;
            }
          }
        }
        cout << endl;
      }


      /* We reset the variables and arrays.
       */
      for (int k = 0; k < n_keys; k++){
        precomp_k[k][0] = 0;
        precomp_k[k][1] = 0;
        top_r_by_key[k].corr = 0.0;
      }

      end = omp_get_wtime();

      is_last_iter = 0;
      col_offset = 0;
      row_offset = 0;
      sample_offset = 0;
      samples_loaded = 0;
      to_load = ncol;
    }
    if (conf.sep == ""){
      printf("[INFO] Attack of byte number %i done in %lf seconds.\n", bn, end - start);
      fflush(stdout);
    }
    if (conf.bitnum == -1) {
      cout << best_out.str() << endl;
    }
  }
  delete[] top_r_by_key;
  delete pqueue;
  delete queues;
  free_matrix(&precomp_k, n_keys); //free_matrix can be found in utils.cpp
  free_matrix(&traces, ncol);
  free_matrix(&tmp, max_n_rows);
  free_matrix(&fin_conf.mat_args->guess, n_keys);
  pthread_mutex_destroy(&lock);
  return 0;
}


/* This function computes the first order correlation between a subset
 * of the traces defined in the structure passed as argument and all the
 * key guesses.
 */
  template <class TypeTrace, class TypeReturn, class TypeGuess>
void * correlation_first_order(void * args_in)
{

  General<TypeTrace, TypeReturn, TypeGuess> * G = (General<TypeTrace, TypeReturn, TypeGuess> *) args_in; //General can be found in memUtils.h
  FirstOrderQueues<TypeReturn> * queues = (FirstOrderQueues<TypeReturn> *)(G->fin_conf->queues); //FirstOrderQueues can be found in cpa.h
  int i, j, k,
      n_keys = G->fin_conf->conf->total_n_keys,
      n_traces = G->fin_conf->conf->n_traces,
      first_sample = G->fin_conf->conf->index_sample,
      offset = G->global_offset;

  TypeReturn corr,
    sum_trace,
    sum_sq_trace,
    tmp,
    std_dev_t;
  CorrFirstOrder<TypeReturn> * q = (CorrFirstOrder<TypeReturn> *) malloc(n_keys * sizeof(CorrFirstOrder<TypeReturn>)); //CorrFirstOrder can be found in cpa.h
  if (q == NULL){
    fprintf (stderr, "[ERROR] Allocating memory for q in correlation\n");
  }

  for (i = G->start; i < G->start + G->length; i++) {
    sum_trace = 0.0;
    sum_sq_trace = 0.0;
    // calculate the sum and sum of squares of the i-th trace
    for (j = 0; j < n_traces; j++){
      tmp = G->fin_conf->mat_args->trace[i][j];
      sum_trace += tmp;
      sum_sq_trace += tmp*tmp;
    }

    // with the sum and sum of squares we can calculate the standard deviation of the trace
    std_dev_t = sqrt(n_traces*sum_sq_trace - sum_trace*sum_trace);
    for (k = 0; k < n_keys; k++) {
      // calculate the standard deviation of the guesses, this is done by the sum and sum of squares via precomp_guesses
      tmp = sqrt(n_traces * G->precomp_guesses[k][1] - G->precomp_guesses[k][0] * G->precomp_guesses[k][0]); //precomp_guesses can be found in memUtils.cpp

      // calculate the correlation between the guess and the measured trace
      // arguments: x[], sum_x, std_dev_x, y[], sum_y, std_dev_y, length
      corr = pearson_v_2_2<TypeReturn, TypeTrace, TypeGuess>(G->fin_conf->mat_args->guess[k],\
        G->precomp_guesses[k][0], tmp, G->fin_conf->mat_args->trace[i], sum_trace, std_dev_t, n_traces); //pearson_v_2_2 can be found in pearson.h

      if (!isnormal(corr)) corr = (TypeReturn) 0;

      q[k].corr  = corr;
      q[k].time  = i + first_sample + offset;
      q[k].key   = k;
    }

    pthread_mutex_lock(&lock);
    for (int key=0; key < n_keys; key++) {
      if (G->fin_conf->conf->key_size == 1)
        queues->pqueue->insert(q[key]);
      if (queues->top_corr[key] < q[key]){
        queues->top_corr[key] = q[key];
      }
    }
    pthread_mutex_unlock(&lock);
  }
  free (q);
  return NULL;
}

/* This function computes the higher moments correlation between a subset
 * of the traces defined in the structure passed as argument and all the
 * key guesses.
 */
  template <class TypeTrace, class TypeReturn, class TypeGuess>
void * higher_moments_correlation(void * args_in)
{

  General<TypeTrace, TypeReturn, TypeGuess> * G = (General<TypeTrace, TypeReturn, TypeGuess> *) args_in; //General can be found in memUtils.h
  FirstOrderQueues<TypeReturn> * queues = (FirstOrderQueues<TypeReturn> *)(G->fin_conf->queues); //FirstOrderQueues can be found in cpa.h
  int i, k,
      n_keys = G->fin_conf->conf->total_n_keys,
      n_traces = G->fin_conf->conf->n_traces,
      first_sample = G->fin_conf->conf->index_sample,
      offset = G->global_offset,
      exponent = G->fin_conf->conf->attack_moment;

  TypeReturn corr,
    sum_trace,
    sum_sq_trace,
    tmp,
    std_dev_t,
    mean_t,
    sigma_n;
  TypeReturn * t = (TypeReturn *) malloc(n_traces * sizeof(TypeReturn));
  if (t == NULL){
    fprintf (stderr, "[ERROR] Allocating memory for t in correlation\n");
  }

  CorrFirstOrder<TypeReturn> * q = (CorrFirstOrder<TypeReturn> *) malloc(n_keys * sizeof(CorrFirstOrder<TypeReturn>)); //CorrFirstOrder can be found in cpa.h
  if (q == NULL){
    fprintf (stderr, "[ERROR] Allocating memory for q in correlation\n");
  }


  for (i = G->start; i < G->start + G->length; i++) {
      sum_trace = 0.0;
      sum_sq_trace = 0.0;
      mean_t = 0.0;
      sigma_n = 0.0;
      for (k = 0; k < n_traces; k++) {
        tmp = G->fin_conf->mat_args->trace[i][k];
        mean_t += tmp;
        sigma_n += tmp*tmp;
      }
      mean_t /= n_traces;
      sigma_n = pow(sqrt(sigma_n/n_traces - mean_t*mean_t), exponent);
      for (k = 0; k < n_traces; k++) {
        tmp = pow((G->fin_conf->mat_args->trace[i][k] - mean_t), exponent)/sigma_n;
        t[k] = tmp;
        sum_trace += tmp;
        sum_sq_trace += tmp*tmp;
      }
      std_dev_t = sqrt(n_traces*sum_sq_trace - sum_trace*sum_trace);
      for (k = 0; k < n_keys; k++) {
        tmp = sqrt(n_traces * G->precomp_guesses[k][1] - G->precomp_guesses[k][0] * G->precomp_guesses[k][0]);

        corr = pearson_v_2_2<TypeReturn, TypeReturn, TypeGuess>(G->fin_conf->mat_args->guess[k],\
          G->precomp_guesses[k][0], tmp, t, sum_trace, std_dev_t, n_traces); //pearson_v_2_2 can be found in pearson.h

        if (!isnormal(corr)) corr = (TypeReturn) 0;

        q[k].corr  = corr;
        q[k].time = i + first_sample + offset;
        q[k].key   = k;
      }
      pthread_mutex_lock(&lock);
      for (int key=0; key < n_keys; key++) {
        if (G->fin_conf->conf->key_size == 1)
          queues->pqueue->insert(q[key]);
        if (queues->top_corr[key] < q[key]){
          queues->top_corr[key] = q[key];
        }
      }
      pthread_mutex_unlock(&lock);


  }
  free (t);
  free (q);
  return NULL;
}


template int first_order<float, double, uint8_t>(Config & conf);
template int first_order<double, double, uint8_t>(Config & conf);
template int first_order<int8_t, double, uint8_t>(Config & conf);
template int first_order<int8_t, float, uint8_t>(Config & conf);

template void * correlation_first_order<int8_t, double, uint8_t> (void * args_in);
template void * correlation_first_order<int8_t, float, uint8_t> (void * args_in);
template void * correlation_first_order<float, double, uint8_t> (void * args_in);
template void * correlation_first_order<double, double, uint8_t> (void * args_in);
