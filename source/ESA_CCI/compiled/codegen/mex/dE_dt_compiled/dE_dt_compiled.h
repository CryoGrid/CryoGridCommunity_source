/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dE_dt_compiled.h
 *
 * Code generation for function 'dE_dt_compiled'
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
void dE_dt_compiled(const emlrtStack *sp, real_T d_energy[15400], const real_T
                    thermCond_eff[15400], const real_T T[15750], const real_T
                    layerDistance[15400]);

/* End of code generation (dE_dt_compiled.h) */
