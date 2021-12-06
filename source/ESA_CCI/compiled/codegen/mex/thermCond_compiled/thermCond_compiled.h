/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * thermCond_compiled.h
 *
 * Code generation for function 'thermCond_compiled'
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
void thermCond_compiled(const emlrtStack *sp, const real_T T[15750], const
  real_T T_end_freezing[14350], const real_T k_frozen[14350], const real_T
  k_freezing[14350], const real_T k_thawed[14350], real_T thermCond[14350]);

/* End of code generation (thermCond_compiled.h) */
