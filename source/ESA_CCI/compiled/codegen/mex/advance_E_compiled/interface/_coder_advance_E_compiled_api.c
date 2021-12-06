/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_advance_E_compiled_api.c
 *
 * Code generation for function '_coder_advance_E_compiled_api'
 *
 */

/* Include files */
#include "_coder_advance_E_compiled_api.h"
#include "advance_E_compiled.h"
#include "advance_E_compiled_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[15400];
static real32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *timestep,
  const char_T *identifier);
static real32_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real32_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15400];
static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *energy,
  const char_T *identifier))[15400];
static void emlrt_marshallOut(const real32_T u[15400], const mxArray *y);
static real32_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real32_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[15400]
{
  real32_T (*y)[15400];
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real32_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *timestep, const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real32_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(timestep), &thisId);
  emlrtDestroyArray(&timestep);
  return y;
}

static real32_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real32_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real32_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[15400]
{
  static const int32_T dims[2] = { 44, 350 };

  real32_T (*ret)[15400];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 2U, dims);
  ret = (real32_T (*)[15400])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real32_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *energy,
  const char_T *identifier))[15400]
{
  emlrtMsgIdentifier thisId;
  real32_T (*y)[15400];
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(energy), &thisId);
  emlrtDestroyArray(&energy);
  return y;
}

static void emlrt_marshallOut(const real32_T u[15400], const mxArray *y)
{
  static const int32_T iv[2] = { 44, 350 };

  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, iv, 2);
}

static real32_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real32_T ret;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "single", false, 0U, &dims);
  ret = *(real32_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void advance_E_compiled_api(const mxArray * const prhs[3], const mxArray *plhs[1])
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *prhs_copy_idx_0;
  real32_T (*d_energy)[15400];
  real32_T (*energy)[15400];
  real32_T timestep;
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_0 = emlrtProtectR2012b(prhs[0], 0, true, -1);

  /* Marshall function inputs */
  energy = emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_0), "energy");
  d_energy = emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "d_energy");
  timestep = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "timestep");

  /* Invoke the target function */
  advance_E_compiled(*energy, *d_energy, timestep);

  /* Marshall function outputs */
  emlrt_marshallOut(*energy, prhs_copy_idx_0);
  plhs[0] = prhs_copy_idx_0;
}

/* End of code generation (_coder_advance_E_compiled_api.c) */
