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
#ifndef PEARSON_H
#define PEARSON_H

#include <math.h>



/* Computes the correlation between the vectors x and y, given the
 * precomputed values sum_* and std_dev_*, using the single pass approach. The
 * precomputed values can be calculated by the function precomp_guesses.
 *
 * The formula can be viewed in its mathematical form on https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
 * Note that std_dev_x and std_dev_y are not the standard deviation but n times the standard deviation
 */
  template <class Type1, class Type2, class Type3>
Type1 pearson_v_2_2(Type3 x[], Type1 sum_x, Type1 std_dev_x, Type2 y[], Type1 sum_y, Type1 std_dev_y, int n)
{
  Type1 sum_xy = 0.0;

  for(int i = 0; i < n; i++) {
    sum_xy += (Type1) x[i] * (Type1) y[i];
  }

  return n * (( sum_xy - (sum_x * sum_y)/n ) /
   (std_dev_x * std_dev_y));

}

#endif
