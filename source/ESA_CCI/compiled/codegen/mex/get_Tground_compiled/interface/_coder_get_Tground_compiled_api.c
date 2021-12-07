/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_get_Tground_compiled_api.c
 *
 * Code generation for function '_coder_get_Tground_compiled_api'
 *
 */

/* Include files */
#include "_coder_get_Tground_compiled_api.h"
#include "get_Tground_compiled.h"
#include "get_Tground_compiled_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15400];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c_thawed,
  const char_T *identifier))[14350];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[14350];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15400];
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *energy,
  const char_T *identifier))[15400];
static const mxArray *emlrt_marshallOut(const real_T u[14000]);
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14350];

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[15400]
{
  real_T (*y)[15400];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *c_thawed, const char_T *identifier))[14350]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[14350];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(c_thawed), &thisId);
  emlrtDestroyArray(&c_thawed);
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
  const emlrtMsgIdentifier *msgId))[15400]
{
  static const int32_T dims[2] = { 44, 350 };

  real_T (*ret)[15400];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[15400])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *energy,
  const char_T *identifier))[15400]
{
  emlrtMsgIdentifier thisId;
  real_T (*y)[15400];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(energy), &thisId);
  emlrtDestroyArray(&energy);
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
  void get_Tground_compiled_api(const mxArray * const prhs[5], const mxArray
  *plhs[1])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  real_T (*energy)[15400];
  real_T (*E_frozen)[14350];
  real_T (*T_end_freezing)[14350];
  real_T (*c_frozen)[14350];
  real_T (*c_thawed)[14350];
  real_T (*T)[14000];
  st.tls = emlrtRootTLSGlobal;
  T = (real_T (*)[14000])mxMalloc(sizeof(real_T [14000]));

  /* Marshall function inputs */
  energy = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "energy");
  c_thawed = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "c_thawed");
  c_frozen = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "c_frozen");
  E_frozen = c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "E_frozen");
  T_end_freezing = c_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "T_end_freezing");

  /* Invoke the target function */
  get_Tground_compiled(&st, *energy, *c_thawed, *c_frozen, *E_frozen,
                       *T_end_freezing, *T);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*T);
}

/* End of code generation (_coder_get_Tground_compiled_api.c) */
