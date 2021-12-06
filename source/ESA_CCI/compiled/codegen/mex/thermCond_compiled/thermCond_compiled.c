/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * thermCond_compiled.c
 *
 * Code generation for function 'thermCond_compiled'
 *
 */

/* Include files */
#include "thermCond_compiled.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void thermCond_compiled(const emlrtStack *sp, const real_T T[15750], const
  real_T T_end_freezing[14350], const real_T k_frozen[14350], const real_T
  k_freezing[14350], const real_T k_thawed[14350], real_T thermCond[14350])
{
  real_T d;
  real_T d1;
  int32_T i;
  int32_T i1;
  int32_T i2;
  (void)sp;
  for (i = 0; i < 350; i++) {
    for (i1 = 0; i1 < 41; i1++) {
      d = T[(i1 + 45 * i) + 4];
      i2 = i1 + 41 * i;
      d1 = T_end_freezing[i2];
      thermCond[i2] = ((real_T)(d < d1) * k_frozen[i2] + (real_T)(d > 0.0) *
                       k_thawed[i2]) + (real_T)((d >= d1) && (d <= 0.0)) *
        k_freezing[i2];
    }
  }
}

/* End of code generation (thermCond_compiled.c) */
