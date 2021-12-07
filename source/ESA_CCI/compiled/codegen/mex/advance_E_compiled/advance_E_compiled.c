/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * advance_E_compiled.c
 *
 * Code generation for function 'advance_E_compiled'
 *
 */

/* Include files */
#include "advance_E_compiled.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void advance_E_compiled(real32_T energy[15400], const real32_T d_energy[15400],
  real32_T timestep)
{
  int32_T i;

  /*  function E = advance_E_compiled(E, timestep, dE_dt) */
  /*  E = E + timestep .* dE_dt; */
  for (i = 0; i < 15400; i++) {
    energy[i] += d_energy[i] * timestep;
  }
}

/* End of code generation (advance_E_compiled.c) */
