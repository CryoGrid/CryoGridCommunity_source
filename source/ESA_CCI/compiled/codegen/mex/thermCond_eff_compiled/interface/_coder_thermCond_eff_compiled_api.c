/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_thermCond_eff_compiled_api.c
 *
 * Code generation for function '_coder_thermCond_eff_compiled_api'
 *
 */

/* Include files */
#include "_coder_thermCond_eff_compiled_api.h"
#include "rt_nonfinite.h"
#include "thermCond_eff_compiled.h"
#include "thermCond_eff_compiled_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15750];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15750];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *thermCond,
  const char_T *identifier))[15750];
static const mxArray *emlrt_marshallOut(const real_T u[14000]);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15750]
{
  real_T (*y)[15750];
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15750]
{
  static const int32_T dims[2] = { 45, 350 };

  real_T (*ret)[15750];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[15750])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *thermCond,
  const char_T *identifier))[15750]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[15750];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(thermCond), &thisId);
  emlrtDestroyArray(&thermCond);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u[14000])
{
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 40, 350 };

  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

void thermCond_eff_compiled_api(const mxArray * const prhs[2], const mxArray
  *plhs[1])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  real_T (*layerThick)[15750];
  real_T (*thermCond)[15750];
  real_T (*thermCond_eff)[14000];
  st.tls = emlrtRootTLSGlobal;
  thermCond_eff = (real_T (*)[14000])mxMalloc(sizeof(real_T [14000]));

  /* Marshall function inputs */
  thermCond = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "thermCond");
  layerThick = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "layerThick");

  /* Invoke the target function */
  thermCond_eff_compiled(&st, *thermCond, *layerThick, *thermCond_eff);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*thermCond_eff);
}

/* End of code generation (_coder_thermCond_eff_compiled_api.c) */
