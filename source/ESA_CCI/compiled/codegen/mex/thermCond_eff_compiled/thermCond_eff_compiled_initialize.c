/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * thermCond_eff_compiled_initialize.c
 *
 * Code generation for function 'thermCond_eff_compiled_initialize'
 *
 */

/* Include files */
#include "thermCond_eff_compiled_initialize.h"
#include "_coder_thermCond_eff_compiled_mex.h"
#include "rt_nonfinite.h"
#include "thermCond_eff_compiled_data.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;

/* Function Definitions */
void thermCond_eff_compiled_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mex_InitInfAndNan();
  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (thermCond_eff_compiled_initialize.c) */
