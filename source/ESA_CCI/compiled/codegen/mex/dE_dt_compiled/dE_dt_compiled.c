/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dE_dt_compiled.c
 *
 * Code generation for function 'dE_dt_compiled'
 *
 */

/* Include files */
#include "dE_dt_compiled.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void dE_dt_compiled(const emlrtStack *sp, real_T d_energy[15400], const real_T
                    thermCond_eff[15400], const real_T T[15750], const real_T
                    layerDistance[15400])
{
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  (void)sp;

  /*  function res   = dE_dt_compiled(k_eff, T, cT_delta) */
  /*   */
  /*           res = (k_eff(2:end,:).*(T(3:end,:)-T(2:end-1,:))./cT_delta(2:end,:) -... */
  /*                  k_eff(1:end-1,:).*(T(2:end-1,:)-T(1:end-2,:))./cT_delta(1:end-1,:));  %size N-1 */
  /* downwards flux */
  /* upwards flux, lower boundary already added */
  for (i = 0; i < 350; i++) {
    for (i1 = 0; i1 < 44; i1++) {
      i2 = i1 + 44 * i;
      i3 = i1 + 45 * i;
      d_energy[i2] -= thermCond_eff[i2] * (T[i3 + 1] - T[i3]) / layerDistance[i2];
    }

    for (i1 = 0; i1 < 43; i1++) {
      i2 = i1 + 44 * i;
      i3 = i1 + 45 * i;
      d_energy[i2] += thermCond_eff[i2 + 1] * (T[i3 + 2] - T[i3 + 1]) /
        layerDistance[i2 + 1];
    }
  }
}

/* End of code generation (dE_dt_compiled.c) */
