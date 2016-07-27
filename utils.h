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

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <queue>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

#ifndef RESOURCES
#define RESOURCES "/usr/share/daredevil"
#endif //RESOURCES

#define QUEUE_INIT 1024
#define QUEUE_PRINT 100

#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define RESET "\033[0m"

#define GIGA 1e9
#define MEGA 1e6
/*
#define ALG_AES                 0
#define ALG_DES_BEFORE          1
#define ALG_DES_AFTER           2
#define ALG_DES_BEFORE_SMALL    3
#define ALG_DES_AFTER_SMALL     4
*/

using namespace std;

/* Used to store the filename and dimensions of the matrix for loading them
 * in files later on.
 */
struct Matrix {

  const char * filename;
  unsigned int n_rows, n_columns;

  Matrix(const char * f_name, unsigned int rows, unsigned int columns):
    filename(f_name), n_rows(rows), n_columns(columns) {
    }
};

/* Structure used to store all the configuration information, used by the
 * config file at the moment.
 */
struct Config {

  /* The number of threads
   */
  int n_threads;

  /* The index of the first sample to start computing from. Useful when we
   * want to target only a subset of the time samples.
   */
  int index_sample;

  /* The number of samples after index_sample we want to correlate.
   */
  int n_samples;

  /* The number of traces we want to analyze, in case we don't want to
   * compute correlation on all of them.
   */
  int n_traces;

  /* The total number of traces, might be useless. To be removed if this is
   * the case.
   */
  int total_n_traces;

  /* The total number of samples, might be useless. To be removed if this is
   * the case.
   */
  int total_n_samples;

  /* The total number of keys, might be useless. To be removed if this is
   * the case.
   */
  int total_n_keys;

  /* The total number of columns of the keys guesses
   */
  int n_col_keys;

  /* Whether we want to transpose the array of traces and guesses.
   */
  bool transpose_traces;
  bool transpose_guesses;

  /* The number of trace and guess files.
   */
  int n_file_trace;
  int n_file_guess;

  /* The type of the traces and guesses, represented by a char.
   * u: uint8_t
   * f: float
   * d: double
   * i: int8_t
   */
  char type_trace;
  char type_guess;
  char type_return;

  /* The matrices structures containing file informations.
   */
  Matrix * traces;
  Matrix * guesses;

  /* The order of the attack
   */
  uint8_t attack_order;

  /* The moment of the attack
   */
  uint8_t attack_moment;

  /* The algorithm to attack.
   * A: AES
   * D: DES
   */
  uint8_t algo;

  /* The round of the algorithm to attack.
   */
  uint32_t round;

  /* The position where to attack.
   */
  uint32_t position;

  /* The list of all position we want to attack.
   */
  vector<uint32_t> all_positions;

  /* The bytenumber to contruct the guesses.
   */
  int bytenum;

  /* The window size when computing higher order attacks.
   */
  int window;

  /* The correct key byte. If specified, the correct key will be highlighted when
   * displaying the results. Could also serve later on when doing known key
   * attack.
   */
  int correct_key;

  /* The complete correct key, in bytes.
   */
  uint8_t * complete_correct_key;

  /* The original correct key, in bytes. This is used for DES as we correlate
   * to the round key and not to the input key directly. However, this is only
   * when printing the configuration.
   */
  uint8_t * original_correct_key;

  /* The key size in bytes.
   */
  int key_size;

  /* The memory dedicated to the attack.
   */
  long int memory;

  /* The number of top element we keep track of globally.
   */
  int top;

  /* The SBOX is specified.
   */
  uint16_t * sbox;

  /* Array to store the multiple sboxes.
   */
  vector<string> all_sboxes;

  /* Switch to specify what des layout for the sboxes we want
   * des_switch = 0 => [8][64]
   * des_switch = 1 => [32][16]
   */
  uint8_t des_switch;

  /* Separator for printing
   */
  string sep;

  /* Do we want to target an individual bit?
   * If so, what bit?
   * -2 = none
   * -1 = all
   */
  int8_t bitnum;

};

/* Parse a file describing an SBOX
 */
int parse_sbox_file(const char * fname, uint16_t ** sbox);

/* Parse the command line arguments. For now, only supports the path to the
 * configuration file.
 */
int parse_args(int argc, char * argv[], char ** config_file);

/* Loads the configuration from a config file.
 */
int load_config(Config & conf, const char * conf_file);

/* Prints the current configuration
 */
void print_config(Config &conf);

  /* Latest version of load file. This function is used to load chunks in the
   * chunk partitioning approach.
   *
   * @param str           Path to the file to be loaded
   * @param mem           Pointer to the array in which to load the chunk
   * @param chunk_size    Number of rows to load
   * @param chunk_offset  Initial position in the rows from which we start loading
   * @param n_columns     Number of columns to load
   * @param col_offset    Initial position in the columns from which we start loading
   * @param tot_n_cols    Total number of columns in the file
   *
   * @return The number of lines read
   */
  template <class Type>
size_t fload(const char str[], Type *** mem, int chunk_size, long int chunk_offset, int n_columns, long int col_offset, int tot_n_cols);

  /* Like load_file but doens't allocate new memory each time.
   */
  template <class Type>
int load_file_v_1(const char str[], Type *** mem, int n_rows, int n_columns, long int offset, int total_n_columns);

/* Loads the file located at str in the 2D array mem, whose dimensions
 * are specified by n_rows and n_columns
 */
template <class Type>
int load_file(const char str[], Type *** mem, int n_rows, int n_columns, long int offset=0, int total_n_columns=0);

/*
 * Loads in the array mem the matrices contained in the array of Matrix
 * matrices (which represent files). A smaller subset can be selected by
 * setting the parameters first_sample and n_samples.
 * We assume that n_columns among all matrices is equal, or that n_rows
 * is equal (or both), but it makes no sense if they're both unequal.
 *
 * Warning: No check is done on the bounds if a smaller subset is selected!
 *
 * @param n_matrices: length of the array matrices
 * @param transpose: If set to true, the resulting array "mem"  will be transposed
 * @param first_sample: Index of the first time_sample we want
 * @param n_samples: number of time samples we want
 */
template <class Type>
int import_matrices(Type *** mem, Matrix * matrices,
    unsigned int n_matrices, bool transpose,
    int first_sample = 0, int n_samples = 0);

/* Prints the top correlations by key, ranked by the correlation value. If the
 * correct key is specified, colors it :).
 */
template <class Type>
void print_top_r(Type corrs[], int n_keys, int correct_key=-1, string csv = "");

#endif
