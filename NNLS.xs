#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"
#include "nnls.h"

MODULE = Math::NNLS		PACKAGE = Math::NNLS

AV *
nnls(AV * A, AV * b, int n, int p, int maxIter, IN_OUT double rNorm)
    CODE:
        double * _A, * _b, * _x, * _workW, * _workZ;
        int * _workIndex;
        int i, j, workMode;

        int A_len = n * p;

        Newxz(_A, A_len, double);
        Newxz(_b, n, double);
        Newxz(_x, p, double);
        Newxz(_workW, p, double);
        Newxz(_workZ, n, double);
        Newxz(_workIndex, p, int);

        for (i = 0; i < A_len; i++)
            _A[i] = SvNV( *av_fetch(A, i, 0) );

        for (i = 0; i < n; i++)
            _b[i] = SvNV( *av_fetch(b, i, 0) );

        _nnls(_A, n, n, p, _b, _x, &rNorm, _workW, _workZ, _workIndex, &workMode, maxIter);

        // printf("rNorm %f\n", rNorm);

        RETVAL = newAV();
        for (i = 0; i < p; i++)
            av_push( RETVAL, newSVnv(_x[i]) );

        Safefree(_A);
        Safefree(_b);
        Safefree(_x);
        Safefree(_workW);
        Safefree(_workZ);
        Safefree(_workIndex);

    OUTPUT:
        RETVAL
