/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_dE_dt_compiled_api.c
 *
 * Code generation for function '_coder_dE_dt_compiled_api'
 *
 */

/* Include files */
#include "_coder_dE_dt_compiled_api.h"
#include "dE_dt_compiled.h"
#include "dE_dt_compiled_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15400];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *T, const
  char_T *identifier))[15750];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15750];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15400];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *d_energy,
  const char_T *identifier))[15400];
static void emlrt_marshallOut(const real_T u[15400], const mxArray *y);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15750];

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15400]
{
  real_T (*y)[15400];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *T,
  const char_T *identifier))[15750]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[15750];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(T), &thisId);
  emlrtDestroyArray(&T);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15750]
{
  real_T (*y)[15750];
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15400]
{
  static const int32_T dims[2] = { 44, 350 };

  real_T (*ret)[15400];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[15400])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *d_energy,
  const char_T *identifier))[15400]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[15400];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(d_energy), &thisId);
  emlrtDestroyArray(&d_energy);
  return y;
}
  static void emlrt_marshallOut(const real_T u[15400], const mxArray *y)
{
  static const int32_T iv[2] = { 44, 350 };

  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, iv, 2);
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15750]
{
  static const int32_T dims[2] = { 45, 350 };

  real_T (*ret)[15750];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[15750])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  void dE_dt_compiled_api(const mxArray * const prhs[4], const mxArray *plhs[1])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *prhs_copy_idx_0;
  real_T (*T)[15750];
  real_T (*d_energy)[15400];
  real_T (*layerDistance)[15400];
  real_T (*thermCond_eff)[15400];
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_0 = emlrtProtectR2012b(prhs[0], 0, true, -1);

  /* Marshall function inputs */
  d_energy = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_0), "d_energy");
  thermCond_eff = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "thermCond_eff");
  T = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "T");
  layerDistance = emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "layerDistance");

  /* Invoke the target function */
  dE_dt_compiled(&st, *d_energy, *thermCond_eff, *T, *layerDistance);

  /* Marshall function outputs */
  emlrt_marshallOut(*d_energy, prhs_copy_idx_0);
  plhs[0] = prhs_copy_idx_0;
}

/* End of code generation (_coder_dE_dt_compiled_api.c) */
