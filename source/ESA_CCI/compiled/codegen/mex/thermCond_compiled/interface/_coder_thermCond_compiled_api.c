/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_thermCond_compiled_api.c
 *
 * Code generation for function '_coder_thermCond_compiled_api'
 *
 */

/* Include files */
#include "_coder_thermCond_compiled_api.h"
#include "rt_nonfinite.h"
#include "thermCond_compiled.h"
#include "thermCond_compiled_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15750];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *T_end_freezing, const char_T *identifier))[14350];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[14350];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15750];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *T, const
  char_T *identifier))[15750];
static const mxArray *emlrt_marshallOut(const real_T u[14350]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14350];

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15750]
{
  real_T (*y)[15750];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *T_end_freezing, const char_T *identifier))[14350]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[14350];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(T_end_freezing), &thisId);
  emlrtDestroyArray(&T_end_freezing);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[14350]
{
  real_T (*y)[14350];
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15750]
{
  static const int32_T dims[2] = { 45, 350 };

  real_T (*ret)[15750];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[15750])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *T, const
  char_T *identifier))[15750]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[15750];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(T), &thisId);
  emlrtDestroyArray(&T);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u[14350])
{
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 41, 350 };

  const mxArray *m;
  const mxArray *y;
  y = NULL;
  m = emlrtCreateNumericArray(2, &iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14350]
{
  static const int32_T dims[2] = { 41, 350 };

  real_T (*ret)[14350];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[14350])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  void thermCond_compiled_api(const mxArray * const prhs[5], const mxArray *
  plhs[1])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  real_T (*T)[15750];
  real_T (*T_end_freezing)[14350];
  real_T (*k_freezing)[14350];
  real_T (*k_frozen)[14350];
  real_T (*k_thawed)[14350];
  real_T (*thermCond)[14350];
  st.tls = emlrtRootTLSGlobal;
  thermCond = (real_T (*)[14350])mxMalloc(sizeof(real_T [14350]));

  /* Marshall function inputs */
  T = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "T");
  T_end_freezing = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "T_end_freezing");
  k_frozen = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "k_frozen");
  k_freezing = c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "k_freezing");
  k_thawed = c_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "k_thawed");

  /* Invoke the target function */
  thermCond_compiled(&st, *T, *T_end_freezing, *k_frozen, *k_freezing, *k_thawed,
                     *thermCond);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*thermCond);
}

/* End of code generation (_coder_thermCond_compiled_api.c) */
