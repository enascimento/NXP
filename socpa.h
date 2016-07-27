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
#ifndef SOCPA_H
#define SOCPA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include "utils.h"
#include "pearson.h"

/* Second order and higher order moments cpa for large files
 */
template <class TypeTrace, class TypeReturn, class TypeGuess>
int second_order(Config & conf);

/* Correlation function used when computing second order
 * for moment = 1
 */
template <class TypeTrace, class TypeReturn, class TypeGuess>
void * second_order_correlation(void * args_in);

#endif
