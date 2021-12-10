/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * get_Tground_compiled.c
 *
 * Code generation for function 'get_Tground_compiled'
 *
 */

/* Include files */
#include "get_Tground_compiled.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void get_Tground_compiled(const emlrtStack *sp, const real_T energy[15400],
  const real_T c_thawed[14350], const real_T c_frozen[14350], const real_T
  E_frozen[14350], const real_T T_end_freezing[14350], real_T T[14000])
{
  real_T d;
  real_T d1;
  real_T d2;
  int32_T i;
  int32_T i1;
  int32_T i2;
  (void)sp;

  /*  function res = get_T_compiled(E, c_thawed, c_frozen, E_frozen, T_end_freezing) */
  /*   */
  /*  res = double(E(5:end,:)>=0) .* E(5:end,:) ./ c_thawed(2:end,:) + ... */
  /*                  double(E(5:end,:) <= E_frozen(2:end,:)) .* ((E(5:end,:) - E_frozen(2:end,:)) ./ c_frozen(2:end,:) + T_end_freezing(2:end,:)) + ... */
  /*                  double(E(5:end,:) < 0 & E(5:end,:) > E_frozen(2:end,:)) .* E(5:end,:)./E_frozen(2:end,:) .*(T_end_freezing(2:end,:)); */
  for (i = 0; i < 350; i++) {
    for (i1 = 0; i1 < 40; i1++) {
      d = energy[(i1 + 44 * i) + 4];
      i2 = (i1 + 41 * i) + 1;
      d1 = E_frozen[i2];
      d2 = T_end_freezing[i2];
      T[i1 + 40 * i] = ((real_T)(d >= 0.0) * d / c_thawed[i2] + (real_T)(d <= d1)
                        * ((d - d1) / c_frozen[i2] + d2)) + (real_T)((d < 0.0) &&
        (d > d1)) * d / d1 * d2;
    }
  }
}

/* End of code generation (get_Tground_compiled.c) */
