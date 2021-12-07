/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * thermCond_eff_compiled.c
 *
 * Code generation for function 'thermCond_eff_compiled'
 *
 */

/* Include files */
#include "thermCond_eff_compiled.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void thermCond_eff_compiled(const emlrtStack *sp, const real_T thermCond[15750],
  const real_T layerThick[15750], real_T thermCond_eff[14000])
{
  real_T d;
  real_T d1;
  real_T d2;
  real_T d3;
  int32_T i;
  int32_T i1;
  int32_T i2;
  (void)sp;
  for (i = 0; i < 350; i++) {
    for (i1 = 0; i1 < 40; i1++) {
      i2 = i1 + 45 * i;
      d = thermCond[i2 + 4];
      d1 = thermCond[i2 + 5];
      d2 = layerThick[i2 + 4];
      d3 = layerThick[i2 + 5];
      thermCond_eff[i1 + 40 * i] = d * d1 * (d2 / 2.0 + d3 / 2.0) / (d * d3 /
        2.0 + d1 * d2 / 2.0);
    }
  }

  /* size N */
}

/* End of code generation (thermCond_eff_compiled.c) */
