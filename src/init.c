#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

extern void getsizeratios(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void increasingreg(void *, void *);
extern void sigmashrink(void *, void *, void *, void *, void *, void *, void *);
extern void wipeout(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"getsizeratios", (DL_FUNC) &getsizeratios, 12},
    {"increasingreg", (DL_FUNC) &increasingreg,  2},
    {"sigmashrink",   (DL_FUNC) &sigmashrink,    7},
    {"wipeout",       (DL_FUNC) &wipeout,        4},
    {NULL, NULL, 0}
};

void R_init_ruv(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
