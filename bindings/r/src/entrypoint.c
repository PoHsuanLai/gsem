// Minimal entrypoint — let extendr handle symbol registration dynamically.

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

extern SEXP wrap__make_gsemr_wrappers(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"wrap__make_gsemr_wrappers", (DL_FUNC) &wrap__make_gsemr_wrappers, 2},
    {NULL, NULL, 0}
};

void R_init_gsemr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}
