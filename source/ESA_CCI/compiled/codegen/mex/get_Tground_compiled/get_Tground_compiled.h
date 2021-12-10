/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * get_Tground_compiled.h
 *
 * Code generation for function 'get_Tground_compiled'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void get_Tground_compiled(const emlrtStack *sp, const real_T energy[15400],
  const real_T c_thawed[14350], const real_T c_frozen[14350], const real_T
  E_frozen[14350], const real_T T_end_freezing[14350], real_T T[14000]);

/* End of code generation (get_Tground_compiled.h) */
