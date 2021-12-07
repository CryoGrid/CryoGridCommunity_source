/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * thermCond_eff_compiled.h
 *
 * Code generation for function 'thermCond_eff_compiled'
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
void thermCond_eff_compiled(const emlrtStack *sp, const real_T thermCond[15750],
  const real_T layerThick[15750], real_T thermCond_eff[14000]);

/* End of code generation (thermCond_eff_compiled.h) */
