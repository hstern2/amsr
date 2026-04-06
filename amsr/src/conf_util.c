/*
 * conf_util.c — C utilities for AMSR conformer generation.
 *
 * Two entry points:
 *   cost_and_grad()  — cost function + analytical gradient for L-BFGS-B
 *   embed()          — distance-geometry embedding with AMSR dihedrals
 *
 * Build:
 *   macOS:  cc -O3 -shared -fPIC -o conf_util.dylib conf_util.c -lm
 *   Linux:  cc -O3 -shared -fPIC -o conf_util.so conf_util.c -lm
 *
 * Portable C99, no external dependencies.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================== */
/* 3-vector helpers                                                    */
/* ================================================================== */

static inline void cross3(const double *a, const double *b, double *out) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static inline double norm3(const double *v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline double dot3(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline void sub3(const double *a, const double *b, double *out) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    out[2] = a[2] - b[2];
}

/* ================================================================== */
/* Torsion angle (degrees) for four points                             */
/* ================================================================== */

static double measure_torsion(const double *p0, const double *p1,
                              const double *p2, const double *p3) {
    double b1[3], b2[3], b3[3], n1[3], n2[3], cn[3];
    sub3(p1, p0, b1);
    sub3(p2, p1, b2);
    sub3(p3, p2, b3);
    cross3(b1, b2, n1);
    cross3(b2, b3, n2);
    double n1n = norm3(n1), n2n = norm3(n2);
    if (n1n < 1e-10 || n2n < 1e-10) return 0.0;
    for (int k = 0; k < 3; k++) { n1[k] /= n1n; n2[k] /= n2n; }
    double b2n = norm3(b2);
    double b2h[3] = {b2[0]/b2n, b2[1]/b2n, b2[2]/b2n};
    cross3(n1, n2, cn);
    return atan2(dot3(cn, b2h), dot3(n1, n2)) * (180.0 / M_PI);
}

/* ================================================================== */
/* Cost function + gradient                                            */
/* ================================================================== */

static inline const double *coord(const double *x, const double *fixed,
                                  int n_free, int slot) {
    if (slot < n_free) return x + 3*slot;
    return fixed + 3*(slot - n_free);
}

static inline void scatter(double *grad, int slot, int n_free,
                           const double *contrib) {
    if (slot < n_free) {
        grad[3*slot]   += contrib[0];
        grad[3*slot+1] += contrib[1];
        grad[3*slot+2] += contrib[2];
    }
}

double cost_and_grad(
    const double *x, double *grad, int n_free,
    const double *fixed, int n_fixed,
    const int *bond_pairs, const double *ideal_lengths, int n_bonds,
    const int *angle_triples, const double *ideal_angles_deg, int n_angles,
    const int *planar_groups, int n_planar,
    const int *chiral_info, const double *chiral_target_vols, int n_chiral,
    const int *dih_quads, const double *dih_targets_deg, int n_dih,
    const int *ez_quads, const double *ez_targets_deg, int n_ez,
    const int *linear_triples, int n_linear,
    double w_bond, double w_angle, double w_planar, double w_chiral,
    double w_dih, double w_ez, double w_linear)
{
    double cost = 0.0;
    memset(grad, 0, 3 * n_free * sizeof(double));

    /* Bond terms */
    for (int ib = 0; ib < n_bonds; ib++) {
        int si = bond_pairs[2*ib], sj = bond_pairs[2*ib+1];
        const double *ri = coord(x, fixed, n_free, si);
        const double *rj = coord(x, fixed, n_free, sj);
        double d[3]; sub3(ri, rj, d);
        double dist = norm3(d);
        if (dist < 1e-10) continue;
        double r = w_bond * (dist - ideal_lengths[ib]);
        cost += r * r;
        double scale = 2.0 * w_bond * r / dist;
        double g[3] = {scale*d[0], scale*d[1], scale*d[2]};
        scatter(grad, si, n_free, g);
        double ng[3] = {-g[0], -g[1], -g[2]};
        scatter(grad, sj, n_free, ng);
    }

    /* Angle terms */
    for (int ia = 0; ia < n_angles; ia++) {
        int sa = angle_triples[3*ia];
        int sb = angle_triples[3*ia+1];
        int sc = angle_triples[3*ia+2];
        const double *ra = coord(x, fixed, n_free, sa);
        const double *rb = coord(x, fixed, n_free, sb);
        const double *rc = coord(x, fixed, n_free, sc);
        double v1[3], v2[3]; sub3(ra, rb, v1); sub3(rc, rb, v2);
        double n1 = norm3(v1), n2 = norm3(v2);
        if (n1 < 1e-10 || n2 < 1e-10) continue;
        double L = 0.5 * (n1 + n2);
        double cos_a = dot3(v1, v2) / (n1 * n2);
        if (cos_a > 1.0) cos_a = 1.0;
        if (cos_a < -1.0) cos_a = -1.0;
        double theta = acos(cos_a);
        double theta0 = ideal_angles_deg[ia] * (M_PI / 180.0);
        double dh = (theta - theta0) / 2.0;
        double r = w_angle * 2.0 * L * sin(dh);
        cost += r * r;
        double sin_th = sin(theta);
        if (sin_th < 1e-10) continue;
        double v1h[3] = {v1[0]/n1, v1[1]/n1, v1[2]/n1};
        double v2h[3] = {v2[0]/n2, v2[1]/n2, v2[2]/n2};
        double dr_dtheta = w_angle * L * cos(dh);
        double dr_dL = w_angle * 2.0 * sin(dh);
        double dth_dra[3], dth_drc[3], dth_drb[3];
        for (int k = 0; k < 3; k++) {
            dth_dra[k] = (cos_a * v1h[k] - v2h[k]) / (sin_th * n1);
            dth_drc[k] = (cos_a * v2h[k] - v1h[k]) / (sin_th * n2);
            dth_drb[k] = -(dth_dra[k] + dth_drc[k]);
        }
        double ga[3], gb[3], gc[3];
        double sc2r = 2.0 * r;
        for (int k = 0; k < 3; k++) {
            ga[k] = sc2r * (dr_dtheta * dth_dra[k] + dr_dL * 0.5 * v1h[k]);
            gc[k] = sc2r * (dr_dtheta * dth_drc[k] + dr_dL * 0.5 * v2h[k]);
            gb[k] = sc2r * (dr_dtheta * dth_drb[k] - dr_dL * 0.5 * (v1h[k] + v2h[k]));
        }
        scatter(grad, sa, n_free, ga);
        scatter(grad, sb, n_free, gb);
        scatter(grad, sc, n_free, gc);
    }

    /* Planarity terms */
    for (int ip = 0; ip < n_planar; ip++) {
        int sj = planar_groups[4*ip], sa = planar_groups[4*ip+1];
        int sb = planar_groups[4*ip+2], sc = planar_groups[4*ip+3];
        const double *rj = coord(x, fixed, n_free, sj);
        const double *ra = coord(x, fixed, n_free, sa);
        const double *rb = coord(x, fixed, n_free, sb);
        const double *rc = coord(x, fixed, n_free, sc);
        double v1[3], v2[3], v3[3];
        sub3(ra, rj, v1); sub3(rb, rj, v2); sub3(rc, rj, v3);
        double c23[3], c31[3], c12[3];
        cross3(v2, v3, c23);
        double vol = dot3(v1, c23);
        double nn1 = norm3(v1), nn2 = norm3(v2), nn3 = norm3(v3);
        double nrm = nn1 * nn2 * nn3 + 1e-10;
        double r = w_planar * vol / nrm;
        cost += r * r;
        cross3(v3, v1, c31); cross3(v1, v2, c12);
        double inv = 1.0 / nrm;
        double q = vol * inv;
        double dnorm_a[3] = {0,0,0}, dnorm_b[3] = {0,0,0}, dnorm_c[3] = {0,0,0};
        if (nn1 > 1e-10) { double f = nn2*nn3/nn1; for (int k=0;k<3;k++) dnorm_a[k]=f*v1[k]; }
        if (nn2 > 1e-10) { double f = nn1*nn3/nn2; for (int k=0;k<3;k++) dnorm_b[k]=f*v2[k]; }
        if (nn3 > 1e-10) { double f = nn1*nn2/nn3; for (int k=0;k<3;k++) dnorm_c[k]=f*v3[k]; }
        double ga[3], gb[3], gc[3], gj[3];
        double sc2r = 2.0 * r;
        for (int k = 0; k < 3; k++) {
            ga[k] = sc2r * w_planar * (c23[k]*inv - q*dnorm_a[k]*inv);
            gb[k] = sc2r * w_planar * (c31[k]*inv - q*dnorm_b[k]*inv);
            gc[k] = sc2r * w_planar * (c12[k]*inv - q*dnorm_c[k]*inv);
            gj[k] = -(ga[k]+gb[k]+gc[k]);
        }
        scatter(grad, sj, n_free, gj);
        scatter(grad, sa, n_free, ga);
        scatter(grad, sb, n_free, gb);
        scatter(grad, sc, n_free, gc);
    }

    /* Chirality terms */
    for (int ic = 0; ic < n_chiral; ic++) {
        int sj = chiral_info[5*ic], sa = chiral_info[5*ic+1];
        int sb = chiral_info[5*ic+2], sc = chiral_info[5*ic+3];
        double sign = (double)chiral_info[5*ic+4];
        double target = chiral_target_vols[ic];
        const double *rj = coord(x, fixed, n_free, sj);
        const double *ra = coord(x, fixed, n_free, sa);
        const double *rb = coord(x, fixed, n_free, sb);
        const double *rc = coord(x, fixed, n_free, sc);
        double v1[3], v2[3], v3[3], c23[3];
        sub3(ra, rj, v1); sub3(rb, rj, v2); sub3(rc, rj, v3);
        cross3(v2, v3, c23);
        double vol = dot3(v1, c23);
        double raw = sign * (target - vol);
        if (raw <= 0.0) continue;
        double r = w_chiral * raw;
        cost += r * r;
        double c31[3], c12[3];
        cross3(v3, v1, c31); cross3(v1, v2, c12);
        double dr_dvol = -w_chiral * sign;
        double sc2r = 2.0 * r;
        double ga[3], gb[3], gc[3], gj[3];
        for (int k = 0; k < 3; k++) {
            ga[k] = sc2r * dr_dvol * c23[k];
            gb[k] = sc2r * dr_dvol * c31[k];
            gc[k] = sc2r * dr_dvol * c12[k];
            gj[k] = -(ga[k]+gb[k]+gc[k]);
        }
        scatter(grad, sj, n_free, gj);
        scatter(grad, sa, n_free, ga);
        scatter(grad, sb, n_free, gb);
        scatter(grad, sc, n_free, gc);
    }

    /* Dihedral + EZ terms (shared loop) */
    struct { const int *quads; const double *targets; int n; double w; }
    dih_sets[2] = {
        {dih_quads, dih_targets_deg, n_dih, w_dih},
        {ez_quads, ez_targets_deg, n_ez, w_ez}
    };
    for (int iset = 0; iset < 2; iset++) {
        const int *quads = dih_sets[iset].quads;
        const double *targets = dih_sets[iset].targets;
        int n = dih_sets[iset].n;
        double w = dih_sets[iset].w;
        for (int id = 0; id < n; id++) {
            int s0=quads[4*id], s1=quads[4*id+1], s2=quads[4*id+2], s3=quads[4*id+3];
            const double *p0=coord(x,fixed,n_free,s0), *p1=coord(x,fixed,n_free,s1);
            const double *p2=coord(x,fixed,n_free,s2), *p3=coord(x,fixed,n_free,s3);
            double actual = measure_torsion(p0, p1, p2, p3);
            double diff = fmod(actual - targets[id] + 540.0, 360.0) - 180.0;
            double r = w * diff;
            cost += r * r;
            double b1[3], b2[3], b3[3], n1v[3], n2v[3];
            sub3(p1,p0,b1); sub3(p2,p1,b2); sub3(p3,p2,b3);
            cross3(b1,b2,n1v); cross3(b2,b3,n2v);
            double n1n=norm3(n1v), n2n=norm3(n2v), b2n=norm3(b2);
            if (n1n<1e-10 || n2n<1e-10 || b2n<1e-10) continue;
            double f0 = -b2n/(n1n*n1n), f3 = b2n/(n2n*n2n);
            double dt_dp0[3], dt_dp3[3];
            for (int k=0;k<3;k++) { dt_dp0[k]=f0*n1v[k]; dt_dp3[k]=f3*n2v[k]; }
            double b1db2=dot3(b1,b2), b3db2=dot3(b3,b2), b2sq=b2n*b2n;
            double c1=b1db2/b2sq, c2=b3db2/b2sq;
            double dt_dp1[3], dt_dp2[3];
            for (int k=0;k<3;k++) {
                dt_dp1[k] = -(c1+1.0)*dt_dp0[k] + c2*dt_dp3[k];
                dt_dp2[k] = c1*dt_dp0[k] - (c2+1.0)*dt_dp3[k];
            }
            double sc = 2.0 * w * r * (180.0 / M_PI);
            double g0[3], g1[3], g2[3], g3[3];
            for (int k=0;k<3;k++) {
                g0[k]=sc*dt_dp0[k]; g1[k]=sc*dt_dp1[k];
                g2[k]=sc*dt_dp2[k]; g3[k]=sc*dt_dp3[k];
            }
            scatter(grad,s0,n_free,g0); scatter(grad,s1,n_free,g1);
            scatter(grad,s2,n_free,g2); scatter(grad,s3,n_free,g3);
        }
    }

    /* Linearity terms */
    for (int il = 0; il < n_linear; il++) {
        int sa=linear_triples[3*il], sb=linear_triples[3*il+1], sc=linear_triples[3*il+2];
        const double *ra=coord(x,fixed,n_free,sa);
        const double *rb=coord(x,fixed,n_free,sb);
        const double *rc=coord(x,fixed,n_free,sc);
        double v1[3], v2[3]; sub3(ra,rb,v1); sub3(rc,rb,v2);
        double n1sq=dot3(v1,v1), n2sq=dot3(v2,v2), d12=dot3(v1,v2);
        double denom=n1sq*n2sq+1e-20;
        double cross_sq=n1sq*n2sq-d12*d12;
        double sin2=cross_sq/denom;
        double w2=w_linear*w_linear;
        cost += w2*sin2;
        double inv=1.0/denom;
        double gb[3];
        for (int k=0;k<3;k++) {
            double da = w2*2.0*((n2sq*v1[k]-d12*v2[k])*inv - sin2*v1[k]/n1sq);
            double dc = w2*2.0*((n1sq*v2[k]-d12*v1[k])*inv - sin2*v2[k]/n2sq);
            gb[k] = -(da+dc);
        }
        scatter(grad, sb, n_free, gb);
    }

    return cost;
}

/* ================================================================== */
/* L-BFGS optimizer (via liblbfgs)                                     */
/* ================================================================== */

#include "lbfgs.h"

/* Context passed to the liblbfgs evaluate callback. */
typedef struct {
    int n_free; const double *fixed; int n_fixed;
    const int *bond_pairs; const double *ideal_lengths; int n_bonds;
    const int *angle_triples; const double *ideal_angles_deg; int n_angles;
    const int *planar_groups; int n_planar;
    const int *chiral_info; const double *chiral_target_vols; int n_chiral;
    const int *dih_quads; const double *dih_targets_deg; int n_dih;
    const int *ez_quads; const double *ez_targets_deg; int n_ez;
    const int *linear_triples; int n_linear;
    double w_bond, w_angle, w_planar, w_chiral, w_dih, w_ez, w_linear;
} opt_ctx;

static lbfgsfloatval_t _evaluate(void *instance,
    const lbfgsfloatval_t *x, lbfgsfloatval_t *g,
    const int n, const lbfgsfloatval_t step)
{
    opt_ctx *c = (opt_ctx*)instance;
    return cost_and_grad(x, g, c->n_free, c->fixed, c->n_fixed,
        c->bond_pairs, c->ideal_lengths, c->n_bonds,
        c->angle_triples, c->ideal_angles_deg, c->n_angles,
        c->planar_groups, c->n_planar,
        c->chiral_info, c->chiral_target_vols, c->n_chiral,
        c->dih_quads, c->dih_targets_deg, c->n_dih,
        c->ez_quads, c->ez_targets_deg, c->n_ez,
        c->linear_triples, c->n_linear,
        c->w_bond, c->w_angle, c->w_planar, c->w_chiral,
        c->w_dih, c->w_ez, c->w_linear);
}

double lbfgs_optimize(
    double *x, int ndim,
    int n_free, const double *fixed, int n_fixed,
    const int *bond_pairs, const double *ideal_lengths, int n_bonds,
    const int *angle_triples, const double *ideal_angles_deg, int n_angles,
    const int *planar_groups, int n_planar,
    const int *chiral_info, const double *chiral_target_vols, int n_chiral,
    const int *dih_quads, const double *dih_targets_deg, int n_dih,
    const int *ez_quads, const double *ez_targets_deg, int n_ez,
    const int *linear_triples, int n_linear,
    double w_bond, double w_angle, double w_planar, double w_chiral,
    double w_dih, double w_ez, double w_linear,
    double ftol, double gtol, int max_iter)
{
    opt_ctx ctx = {
        n_free, fixed, n_fixed,
        bond_pairs, ideal_lengths, n_bonds,
        angle_triples, ideal_angles_deg, n_angles,
        planar_groups, n_planar,
        chiral_info, chiral_target_vols, n_chiral,
        dih_quads, dih_targets_deg, n_dih,
        ez_quads, ez_targets_deg, n_ez,
        linear_triples, n_linear,
        w_bond, w_angle, w_planar, w_chiral, w_dih, w_ez, w_linear
    };

    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    param.epsilon = gtol;
    param.delta = ftol;
    param.max_iterations = max_iter;
    param.m = 10;

    lbfgsfloatval_t fx;
    lbfgs(ndim, x, &fx, _evaluate, NULL, &ctx, &param);
    return fx;
}

/* ================================================================== */
/* Distance geometry embedding                                         */
/* ================================================================== */

static double dist_from_angle(double d_ab, double d_bc, double angle_deg) {
    double a = angle_deg * M_PI / 180.0;
    return sqrt(d_ab*d_ab + d_bc*d_bc - 2.0*d_ab*d_bc*cos(a));
}

static double dist_from_dihedral(double d_ab, double d_bc, double d_cd,
                                  double ang_abc, double ang_bcd,
                                  double dih_deg) {
    double a1 = ang_abc * M_PI / 180.0;
    double a2 = ang_bcd * M_PI / 180.0;
    double tau = dih_deg * M_PI / 180.0;
    double ax = -d_ab * cos(M_PI - a1);
    double ay =  d_ab * sin(M_PI - a1);
    double dx = d_bc + d_cd * cos(M_PI - a2);
    double dy = d_cd * sin(M_PI - a2) * cos(tau);
    double dz = d_cd * sin(M_PI - a2) * sin(tau);
    double ex = dx - ax, ey = dy - ay, ez = dz;
    return sqrt(ex*ex + ey*ey + ez*ez);
}

static void smooth_bounds(double *lower, double *upper, int n) {
    int changed = 1;
    for (int iter = 0; iter < 10 && changed; iter++) {
        changed = 0;
        for (int k = 0; k < n; k++)
            for (int i = 0; i < n; i++) {
                if (i == k) continue;
                for (int j = i+1; j < n; j++) {
                    if (j == k) continue;
                    int ij=i*n+j, ji=j*n+i, ik=i*n+k, kj=k*n+j;
                    double u_new = upper[ik] + upper[kj];
                    if (u_new < upper[ij]) { upper[ij]=upper[ji]=u_new; changed=1; }
                    double l1=lower[ik]-upper[kj], l2=lower[kj]-upper[ik];
                    double l_new = l1>l2 ? l1 : l2;
                    if (l_new > lower[ij]) { lower[ij]=lower[ji]=l_new; changed=1; }
                    if (lower[ij] > upper[ij]) {
                        double mid=0.5*(lower[ij]+upper[ij]);
                        lower[ij]=lower[ji]=mid; upper[ij]=upper[ji]=mid;
                    }
                }
            }
    }
}

static unsigned int _rng_state;
static void seed_rng(unsigned int s) { _rng_state = s; }
static double rand_uniform(void) {
    _rng_state = _rng_state * 1103515245u + 12345u;
    return (_rng_state >> 16) / 32768.0;
}

static void jacobi_rotate(double *A, double *V, int n, int p, int q) {
    double app=A[p*n+p], aqq=A[q*n+q], apq=A[p*n+q];
    if (fabs(apq) < 1e-15) return;
    double tau_j = (aqq-app)/(2.0*apq);
    double t = (tau_j>=0?1.0:-1.0)/(fabs(tau_j)+sqrt(1.0+tau_j*tau_j));
    double c = 1.0/sqrt(1.0+t*t), s = t*c;
    A[p*n+p] = app-t*apq; A[q*n+q] = aqq+t*apq; A[p*n+q]=A[q*n+p]=0.0;
    for (int r=0;r<n;r++) {
        if (r==p||r==q) continue;
        double arp=A[r*n+p], arq=A[r*n+q];
        A[r*n+p]=A[p*n+r]=c*arp-s*arq;
        A[r*n+q]=A[q*n+r]=s*arp+c*arq;
    }
    for (int r=0;r<n;r++) {
        double vrp=V[r*n+p], vrq=V[r*n+q];
        V[r*n+p]=c*vrp-s*vrq; V[r*n+q]=s*vrp+c*vrq;
    }
}

static void eigen_symmetric(double *A, double *evals, double *evecs, int n) {
    double *M = (double*)malloc(n*n*sizeof(double));
    memcpy(M, A, n*n*sizeof(double));
    memset(evecs, 0, n*n*sizeof(double));
    for (int i=0;i<n;i++) evecs[i*n+i]=1.0;
    for (int iter=0;iter<100;iter++) {
        double off=0;
        for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) off+=M[i*n+j]*M[i*n+j];
        if (off < 1e-20) break;
        for (int p=0;p<n;p++) for (int q=p+1;q<n;q++) jacobi_rotate(M,evecs,n,p,q);
    }
    for (int i=0;i<n;i++) evals[i]=M[i*n+i];
    free(M);
}

void embed(int n, double *coords_out,
           int n_bonds, const int *bond_pairs, const double *bond_lengths,
           int n_angles, const int *angle_triples, const double *angle_values,
           int n_dihedrals, const int *dihedral_quads, const double *dihedral_values,
           unsigned int seed) {
    if (n <= 0) return;

    double *lower = (double*)calloc(n*n, sizeof(double));
    double *upper = (double*)malloc(n*n*sizeof(double));

    /* Default upper bound: rough molecular diameter estimate. */
    double avg_bond = 0;
    for (int i = 0; i < n_bonds; i++) avg_bond += bond_lengths[i];
    avg_bond = n_bonds > 0 ? avg_bond / n_bonds : 1.5;
    double max_dist = sqrt((double)n) * avg_bond * 2.0;
    if (max_dist < 5.0) max_dist = 5.0;

    for (int i = 0; i < n*n; i++) upper[i] = max_dist;
    for (int i = 0; i < n; i++) { lower[i*n+i] = 0; upper[i*n+i] = 0; }

    double min_dist = 1.5;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++) lower[i*n+j]=lower[j*n+i]=min_dist;

    /* Build bond length lookup matrix for O(1) access. */
    double *blen = (double*)calloc(n*n, sizeof(double));
    for (int k=0;k<n_bonds;k++) {
        int i=bond_pairs[2*k], j=bond_pairs[2*k+1];
        blen[i*n+j]=blen[j*n+i]=bond_lengths[k];
        lower[i*n+j]=lower[j*n+i]=bond_lengths[k];
        upper[i*n+j]=upper[j*n+i]=bond_lengths[k];
    }

    /* 1-3 distances from angles */
    for (int k=0;k<n_angles;k++) {
        int a=angle_triples[3*k], b=angle_triples[3*k+1], c=angle_triples[3*k+2];
        double d_ab=blen[a*n+b], d_bc=blen[b*n+c];
        if (d_ab>0 && d_bc>0) {
            double d=dist_from_angle(d_ab, d_bc, angle_values[k]);
            lower[a*n+c]=lower[c*n+a]=d; upper[a*n+c]=upper[c*n+a]=d;
        }
    }

    /* 1-4 distances from dihedrals */
    for (int k=0;k<n_dihedrals;k++) {
        int a=dihedral_quads[4*k], b=dihedral_quads[4*k+1];
        int c=dihedral_quads[4*k+2], dd=dihedral_quads[4*k+3];
        double d_ab=blen[a*n+b], d_bc=blen[b*n+c], d_cd=blen[c*n+dd];
        if (d_ab<=0 || d_bc<=0 || d_cd<=0) continue;
        double ang_abc=-1, ang_bcd=-1;
        for (int m=0;m<n_angles;m++) {
            int ta=angle_triples[3*m], tb=angle_triples[3*m+1], tc=angle_triples[3*m+2];
            if (tb==b && ((ta==a&&tc==c)||(ta==c&&tc==a))) ang_abc=angle_values[m];
            if (tb==c && ((ta==b&&tc==dd)||(ta==dd&&tc==b))) ang_bcd=angle_values[m];
        }
        if (ang_abc>0 && ang_bcd>0) {
            double d=dist_from_dihedral(d_ab,d_bc,d_cd,ang_abc,ang_bcd,dihedral_values[k]);
            double tol=0.1;
            double lo = d-tol>min_dist ? d-tol : min_dist;
            double hi = d+tol;
            if (lo > lower[a*n+dd]) lower[a*n+dd]=lower[dd*n+a]=lo;
            if (hi < upper[a*n+dd]) upper[a*n+dd]=upper[dd*n+a]=hi;
        }
    }

    smooth_bounds(lower, upper, n);

    /* Sample distances.  Use lower bounds plus a small random
     * perturbation to produce a compact, nearly-realizable matrix.
     * Large random offsets create distance matrices that can't be
     * embedded in 3D, causing the MDS to collapse atoms together. */
    seed_rng(seed);
    double *D = (double*)malloc(n*n*sizeof(double));
    for (int i=0;i<n;i++) {
        D[i*n+i]=0;
        for (int j=i+1;j<n;j++) {
            double lo=lower[i*n+j], hi=upper[i*n+j];
            if (lo>hi) { double m=0.5*(lo+hi); lo=hi=m; }
            double d;
            if (hi-lo < 0.5) {
                d=0.5*(lo+hi);
            } else {
                /* Sample in the middle third of the range */
                double r = rand_uniform();
                d = lo + (0.33 + 0.34 * r) * (hi - lo);
            }
            D[i*n+j]=D[j*n+i]=d;
        }
    }

    /* Metric matrix embedding */
    double *D2 = (double*)malloc(n*n*sizeof(double));
    for (int i=0;i<n*n;i++) D2[i]=D[i]*D[i];
    double *row_mean = (double*)calloc(n, sizeof(double));
    double grand_mean = 0;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) row_mean[i]+=D2[i*n+j];
        row_mean[i]/=n; grand_mean+=row_mean[i];
    }
    grand_mean /= n;
    double *G = (double*)malloc(n*n*sizeof(double));
    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
            G[i*n+j] = -0.5*(D2[i*n+j]-row_mean[i]-row_mean[j]+grand_mean);

    double *evals = (double*)malloc(n*sizeof(double));
    double *evecs = (double*)malloc(n*n*sizeof(double));
    eigen_symmetric(G, evals, evecs, n);

    /* Sort eigenvalues descending (top 3 needed) */
    for (int i=0;i<3&&i<n;i++) {
        int best=i;
        for (int j=i;j<n;j++) if (evals[j]>evals[best]) best=j;
        if (best!=i) {
            double tmp=evals[i]; evals[i]=evals[best]; evals[best]=tmp;
            for (int r=0;r<n;r++) {
                double t=evecs[r*n+i]; evecs[r*n+i]=evecs[r*n+best]; evecs[r*n+best]=t;
            }
        }
    }

    int dim = 3<n ? 3 : n;
    for (int i=0;i<n;i++) {
        for (int k=0;k<dim;k++) {
            double ev = evals[k]>0 ? evals[k] : 0;
            coords_out[i*3+k] = evecs[i*n+k]*sqrt(ev);
        }
        for (int k=dim;k<3;k++) coords_out[i*3+k]=0;
    }

    /* Iterative projection onto bond + angle constraints */
    for (int iter=0;iter<50;iter++) {
        for (int k=0;k<n_bonds;k++) {
            int bi=bond_pairs[2*k], bj=bond_pairs[2*k+1];
            double dx=coords_out[bj*3]-coords_out[bi*3];
            double dy=coords_out[bj*3+1]-coords_out[bi*3+1];
            double dz=coords_out[bj*3+2]-coords_out[bi*3+2];
            double dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (dist<1e-10) { coords_out[bj*3]+=0.1; continue; }
            double sc=0.5*(bond_lengths[k]-dist)/dist;
            coords_out[bi*3]-=sc*dx; coords_out[bi*3+1]-=sc*dy; coords_out[bi*3+2]-=sc*dz;
            coords_out[bj*3]+=sc*dx; coords_out[bj*3+1]+=sc*dy; coords_out[bj*3+2]+=sc*dz;
        }
    }
    for (int iter=0;iter<20;iter++) {
        for (int k=0;k<n_angles;k++) {
            int a=angle_triples[3*k], b=angle_triples[3*k+1], c=angle_triples[3*k+2];
            double d_ab=blen[a*n+b], d_bc=blen[b*n+c];
            if (d_ab<=0||d_bc<=0) continue;
            double target_ac=dist_from_angle(d_ab,d_bc,angle_values[k]);
            double dx=coords_out[c*3]-coords_out[a*3];
            double dy=coords_out[c*3+1]-coords_out[a*3+1];
            double dz=coords_out[c*3+2]-coords_out[a*3+2];
            double dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (dist<1e-10) continue;
            double sc=0.25*(target_ac-dist)/dist;
            coords_out[a*3]-=sc*dx; coords_out[a*3+1]-=sc*dy; coords_out[a*3+2]-=sc*dz;
            coords_out[c*3]+=sc*dx; coords_out[c*3+1]+=sc*dy; coords_out[c*3+2]+=sc*dz;
        }
        for (int k=0;k<n_bonds;k++) {
            int bi=bond_pairs[2*k], bj=bond_pairs[2*k+1];
            double dx=coords_out[bj*3]-coords_out[bi*3];
            double dy=coords_out[bj*3+1]-coords_out[bi*3+1];
            double dz=coords_out[bj*3+2]-coords_out[bi*3+2];
            double dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (dist<1e-10) continue;
            double sc=0.5*(bond_lengths[k]-dist)/dist;
            coords_out[bi*3]-=sc*dx; coords_out[bi*3+1]-=sc*dy; coords_out[bi*3+2]-=sc*dz;
            coords_out[bj*3]+=sc*dx; coords_out[bj*3+1]+=sc*dy; coords_out[bj*3+2]+=sc*dz;
        }
    }

    free(lower); free(upper); free(D); free(D2);
    free(row_mean); free(G); free(evals); free(evecs);
    free(blen);
}
