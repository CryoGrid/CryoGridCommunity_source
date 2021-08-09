/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_thermCond_compiled_mex.c
 *
 * Code generation for function '_coder_thermCond_compiled_mex'
 *
 */

/* Include files */
#include "_coder_thermCond_compiled_mex.h"
#include "_coder_thermCond_compiled_api.h"
#include "rt_nonfinite.h"
#include "thermCond_compiled_data.h"
#include "thermCond_compiled_initialize.h"
#include "thermCond_compiled_terminate.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(&thermCond_compiled_atexit);

  /* Module initialization. */
  thermCond_compiled_initialize();

  /* Dispatch the entry-point. */
  thermCond_compiled_mexFunction(nlhs, plhs, nrhs, prhs);

  /* Module termination. */
  thermCond_compiled_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

void thermCond_compiled_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[5])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *outputs[1];
  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        18, "thermCond_compiled");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 18,
                        "thermCond_compiled");
  }

  /* Call the function. */
  thermCond_compiled_api(prhs, outputs);

  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, plhs, outputs);
}

/* End of code generation (_coder_thermCond_compiled_mex.c) */
