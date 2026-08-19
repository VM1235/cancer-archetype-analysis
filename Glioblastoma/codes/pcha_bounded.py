"""py_pcha PCHA with a bounded line search.

Stock py_pcha can hang on 1-D data: S_update/C_update use `while True` and
halve the step forever if SSE never improves. Same algorithm otherwise.
Used only for hypothesis-driven GBM k=2.
"""

import time
from datetime import datetime as dt

import numpy as np
from py_pcha.furthest_sum import furthest_sum
from scipy.sparse import csr_matrix

MAX_LINESEARCH = 64


def PCHA(X, noc, I=None, U=None, delta=0, verbose=False, conv_crit=1e-6, maxiter=500):
    def S_update(S, XCtX, CtXtXC, muS, SST, SSE, niter):
        noc, J = S.shape
        e = np.ones((noc, 1))
        for _ in range(niter):
            SSE_old = SSE
            g = (np.dot(CtXtXC, S) - XCtX) / (SST / J)
            g = g - e * np.sum(g.A * S.A, axis=0)
            S_old = S
            accepted = False
            for _ in range(MAX_LINESEARCH):
                S = (S_old - g * muS).clip(min=0)
                colsum = np.dot(e, np.sum(S, axis=0))
                S = np.divide(S, colsum, out=np.ones_like(S) / noc, where=np.asarray(colsum) != 0)
                SSt = S * S.T
                SSE = SST - 2 * np.sum(XCtX.A * S.A) + np.sum(CtXtXC.A * SSt.A)
                if SSE <= SSE_old * (1 + 1e-9):
                    muS = muS * 1.2
                    accepted = True
                    break
                muS = muS / 2
                if muS < 1e-20:
                    break
            if not accepted:
                S = S_old
                SSt = S * S.T
                SSE = SSE_old
        return S, SSE, muS, SSt

    def C_update(X, XSt, XC, SSt, C, delta, muC, mualpha, SST, SSE, niter=1):
        J, nos = C.shape
        if delta != 0:
            alphaC = np.sum(C, axis=0).A[0]
            C = np.dot(C, np.diag(1 / alphaC))
        e = np.ones((J, 1))
        XtXSt = np.dot(X.T, XSt)
        for _ in range(niter):
            SSE_old = SSE
            g = (np.dot(X.T, np.dot(XC, SSt)) - XtXSt) / SST
            if delta != 0:
                g = np.dot(g, np.diag(alphaC))
            g = g.A - e * np.sum(g.A * C.A, axis=0)
            C_old = C
            accepted = False
            for _ in range(MAX_LINESEARCH):
                C = (C_old - muC * g).clip(min=0)
                nC = np.sum(C, axis=0) + np.finfo(float).eps
                C = np.dot(C, np.diag(1 / nC.A[0]))
                Ct = C * np.diag(alphaC) if delta != 0 else C
                XC = np.dot(X, Ct)
                CtXtXC = np.dot(XC.T, XC)
                SSE = SST - 2 * np.sum(XC.A * XSt.A) + np.sum(CtXtXC.A * SSt.A)
                if SSE <= SSE_old * (1 + 1e-9):
                    muC = muC * 1.2
                    accepted = True
                    break
                muC = muC / 2
                if muC < 1e-20:
                    break
            if not accepted:
                C = C_old
                SSE = SSE_old
            if delta != 0:
                SSE_old = SSE
                g = (np.diag(CtXtXC * SSt).T / alphaC - np.sum(C.A * XtXSt.A)) / (SST * J)
                alphaC_old = alphaC
                for _ in range(MAX_LINESEARCH):
                    alphaC = alphaC_old - mualpha * g
                    alphaC[alphaC < 1 - delta] = 1 - delta
                    alphaC[alphaC > 1 + delta] = 1 + delta
                    XCt = np.dot(XC, np.diag(alphaC / alphaC_old))
                    CtXtXC = np.dot(XCt.T, XCt)
                    SSE = SST - 2 * np.sum(XCt.A * XSt.A) + np.sum(CtXtXC.A * SSt.A)
                    if SSE <= SSE_old * (1 + 1e-9):
                        mualpha = mualpha * 1.2
                        XC = XCt
                        break
                    mualpha = mualpha / 2
                    if mualpha < 1e-20:
                        alphaC = alphaC_old
                        break
        if delta != 0:
            C = C * np.diag(alphaC)
        return C, SSE, muC, mualpha, CtXtXC, XC

    N, M = X.shape
    if I is None:
        I = range(M)
    if U is None:
        U = range(M)
    SST = np.sum(X[:, U] * X[:, U])
    try:
        i = furthest_sum(X[:, I], noc, [int(np.ceil(len(I) * np.random.rand()))])
    except IndexError as err:
        raise RuntimeError("Initialization does not converge. Too few examples in dataset.") from err

    j = range(noc)
    C = csr_matrix((np.ones(len(i)), (i, j)), shape=(len(I), noc)).todense()
    XC = np.dot(X[:, I], C)
    muS, muC, mualpha = 1, 1, 1
    XCtX = np.dot(XC.T, X[:, U])
    CtXtXC = np.dot(XC.T, XC)
    S = -np.log(np.random.random((noc, len(U))))
    S = S / np.dot(np.ones((noc, 1)), np.mat(np.sum(S, axis=0)))
    SSt = np.dot(S, S.T)
    SSE = SST - 2 * np.sum(XCtX.A * S.A) + np.sum(CtXtXC.A * SSt.A)
    S, SSE, muS, SSt = S_update(S, XCtX, CtXtXC, muS, SST, SSE, 25)

    iter_ = 0
    dSSE = np.inf
    t1 = dt.now()
    varexpl = (SST - SSE) / SST
    while np.abs(dSSE) >= conv_crit * np.abs(SSE) and iter_ < maxiter and varexpl < 0.9999:
        told = t1
        iter_ += 1
        SSE_old = SSE
        XSt = np.dot(X[:, U], S.T)
        C, SSE, muC, mualpha, CtXtXC, XC = C_update(
            X[:, I], XSt, XC, SSt, C, delta, muC, mualpha, SST, SSE, 10
        )
        XCtX = np.dot(XC.T, X[:, U])
        S, SSE, muS, SSt = S_update(S, XCtX, CtXtXC, muS, SST, SSE, 10)
        dSSE = SSE_old - SSE
        t1 = dt.now()
        time.sleep(0.000001)
        varexpl = (SST - SSE) / SST
        if verbose:
            print(iter_, varexpl, SSE, (t1 - told).seconds)

    varexpl = (SST - SSE) / SST
    ind, vals = zip(*sorted(enumerate(np.sum(S, axis=1)), key=lambda x: x[0], reverse=1))
    S = S[ind, :]
    C = C[:, ind]
    XC = XC[:, ind]
    return XC, S, C, SSE, varexpl
