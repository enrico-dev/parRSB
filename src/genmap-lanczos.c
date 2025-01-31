#include "genmap-impl.h"

#include <math.h>
#include <stdio.h>
//
// Orthogonalize by 1-vector (vector of all 1's)
//
int GenmapOrthogonalizebyOneVector(GenmapHandle h, GenmapComm c,
                                   GenmapVector q1, GenmapLong n) {
  GenmapInt i;
  GenmapScalar sum = 0.0;
  for(i = 0;  i < q1->size; i++) {
    sum += q1->data[i];
  }

  GenmapGop(c, &sum, 1, GENMAP_SCALAR, GENMAP_SUM);
  sum /= n;
  for(i = 0;  i < q1->size; i++) {
    q1->data[i] -= sum;
  }

  return 0;
}
//
// Lanczos version used in Paul's original genmap code.
//
int GenmapLanczosLegendary(GenmapHandle h, GenmapComm c, GenmapVector f,
                           GenmapInt niter, GenmapVector **rr, GenmapVector diag,
                           GenmapVector upper) {
  assert(diag->size == niter);
  assert(diag->size == upper->size + 1);
  assert(f->size == GenmapGetNLocalElements(h));

  if(GenmapGetNGlobalElements(h) < niter) {
    niter = GenmapGetNGlobalElements(h);
    diag->size = niter;
    upper->size = niter - 1;
  }

  GenmapScalar eps = 1.e-5;
  GenmapScalar alpha, beta;
  GenmapScalar rnorm, rtol, rni, rtr, rtz1, rtz2, pap = 0.0, pap_old;
  GenmapVector r, p, w, weights;

  rtz1 = 1.0;

  GenmapScalar tmp;
  GenmapInt lelt = GenmapGetNLocalElements(h);

  // Store Local Laplacian weights
  GenmapCreateVector(&weights, lelt);
  GenmapCreateZerosVector(&p, lelt);
  GenmapCreateVector(&w, lelt);
  GenmapInitLaplacian(h, c, weights);

  // Create vector r orthogonalizing init in 1-norm to (1,1,1...)
  GenmapCreateVector(&r, lelt);
  GenmapCopyVector(r, f);
  GenmapOrthogonalizebyOneVector(h, c, r, GenmapGetNGlobalElements(h));
  rtr = GenmapDotVector(r, r);
  GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
  rnorm = sqrt(rtr);
  rtol = rnorm * eps;
  rni = 1.0 / rnorm;

  // Allocate memory for q-vectors
  if(*rr == NULL) {
    GenmapMalloc((size_t)(niter + 1), rr);
    GenmapInt i;
    for(i = 0; i < niter + 1; ++i)(*rr)[i] = NULL;
  }
  GenmapCreateVector(&(*rr)[0], lelt);

  GenmapScaleVector((*rr)[0], r, rni);

  int iter;
  for(iter = 0; iter < niter; iter++) {
    rtz2 = rtz1;
    rtz1 = rtr;
    beta = rtz1 / rtz2;
    if(iter == 0) beta = 0.0;

    GenmapAxpbyVector(p, p, beta, r, 1.0);
    GenmapOrthogonalizebyOneVector(h, c, p, GenmapGetNGlobalElements(h));

    // Multiplication by the laplacian
    GenmapLaplacian(h, c, p, weights, w);
    GenmapScaleVector(w, w, -1.0);

    pap_old = pap;
    pap = GenmapDotVector(w, p);
    GenmapGop(c, &pap, 1, GENMAP_SCALAR, GENMAP_SUM);
    alpha = rtz1 / pap;
    GenmapAxpbyVector(r, r, 1.0, w, -1.0 * alpha);

    rtr = GenmapDotVector(r, r);
    GenmapGop(c, &rtr, 1, GENMAP_SCALAR, GENMAP_SUM);
    if(rtr < 0) {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
    rnorm = sqrt(rtr);
    rni = 1.0 / rnorm;
    GenmapCreateVector(&(*rr)[iter + 1], lelt);
    GenmapScaleVector((*rr)[iter + 1], r, rni);

    if(iter == 0) {
      diag->data[iter] = pap / rtz1;
    } else {
      diag->data[iter] = (beta * beta * pap_old + pap) / rtz1;
      upper->data[iter - 1] = -beta * pap_old / sqrt(rtz2 * rtz1);
    }
#if defined(GENMAP_DEBUG)
    if(GenmapCommRank(GenmapGetGlobalComm(h)) == 0) {
      printf("diag[%d]="GenmapScalarFormat"\n", iter + 1, diag->data[iter]);
    }
#endif

    if(rnorm < rtol)  {
      diag->size = iter + 1;
      upper->size = iter;
      iter = iter + 1;
      break;
    }
  }

  GenmapDestroyVector(weights);
  GenmapDestroyVector(p);
  GenmapDestroyVector(w);
  GenmapDestroyVector(r);

  return iter;
}
//
// Lanczos version in Introduction to Sci. Comp by Prof. Heath.
//
int GenmapLanczos(GenmapHandle h, GenmapComm c, GenmapVector init,
                  GenmapInt iter, GenmapVector **q, GenmapVector alpha,
                  GenmapVector beta) {
  assert(alpha->size == iter);
  assert(alpha->size == beta->size + 1);
  assert(init->size == GenmapGetNLocalElements(h));

  if(GenmapGetNGlobalElements(h) < iter) {
    iter = GenmapGetNGlobalElements(h);
    alpha->size = iter;
    beta->size = iter - 1;
  }

  GenmapVector q0, q1, u;
  GenmapScalar normq1 = 0., b = 0.;

  GenmapInt lelt = GenmapGetNLocalElements(h);

  // Create vector q1 orthogonalizing init in 1-norm to (1,1,1...)
  GenmapCreateVector(&q1, lelt);
  GenmapCopyVector(q1, init);
  GenmapOrthogonalizebyOneVector(h, c, q1, GenmapGetNGlobalElements(h));
  normq1 = GenmapDotVector(q1, q1);
  GenmapGop(c, &normq1, 1, GENMAP_SCALAR, GENMAP_SUM);
  normq1 = sqrt(normq1);
  GenmapScaleVector(q1, q1, 1. / normq1);

  // Create vector u
  GenmapCreateVector(&u, lelt);

  // Set q_0 and beta_0 to zero (both uses 0-indexing)
  GenmapCreateZerosVector(&q0, lelt);
  beta->data[0] = 0.;

  // Allocate memory for q-vectors
  if(*q == NULL) {
    GenmapMalloc((size_t)iter, q);
    GenmapInt i;
    for(i = 0; i < iter; ++i)
      (*q)[i] = NULL;
  }

  // Store Local Laplacian weights
  GenmapVector weights;
  GenmapCreateVector(&weights, lelt);
  GenmapInitLaplacian(h, c, weights);

  int k;
  for(k = 0; k < iter; k++) {
    // Store q1
    GenmapCreateVector(&(*q)[k], lelt);
    GenmapCopyVector((*q)[k], q1);

    // Multiplication by the laplacian
    GenmapLaplacian(h, c, q1, weights, u);

    alpha->data[k] = GenmapDotVector(q1, u);
    GenmapGop(c, &alpha->data[k], 1, GENMAP_SCALAR, GENMAP_SUM);

    GenmapAxpbyVector(u, u, 1., q0, -b);
    GenmapAxpbyVector(u, u, 1., q1, -alpha->data[k]);

    b = GenmapDotVector(u, u);
    GenmapGop(c, &b, 1, GENMAP_SCALAR, GENMAP_SUM);
    b = sqrt(b);

    if(k < iter - 1) {
      beta->data[k] = b;
      if(fabs(b) < normq1 * GENMAP_TOL) {
        beta->size = k;
        alpha->size = k + 1;
        k = k + 1;
        break;
      }

      GenmapCopyVector(q0, q1);
      GenmapScaleVector(q1, u, 1. / beta->data[k]);
    }
  }

  GenmapDestroyVector(q0);
  GenmapDestroyVector(q1);
  GenmapDestroyVector(u);
  GenmapDestroyVector(weights);

  return k;
}
