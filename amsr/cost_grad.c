/*
 * cost_grad.c — C implementation of the ring geometry cost function and
 * analytical gradient for amsr/zmatrix.py.
 *
 * Replaces _cost_and_grad (Python/numpy) with a single C function call,
 * eliminating per-iteration Python→numpy dispatch overhead.
 *
 * Build:
 *   macOS:  cc -O3 -shared -fPIC -o cost_grad.dylib cost_grad.c -lm
 *   Linux:  cc -O3 -shared -fPIC -o cost_grad.so cost_grad.c -lm
 *
 * Portable C99, no external dependencies.  Structured for future CUDA porting
 * (each constraint loop maps to a GPU kernel).
 */

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* 3-vector helpers (inlined for performance)                         */
/* ------------------------------------------------------------------ */

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

/* ------------------------------------------------------------------ */
/* Torsion angle (degrees) for four points                            */
/* ------------------------------------------------------------------ */

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

/* ------------------------------------------------------------------ */
/* Lookup: slot -> coordinate pointer                                 */
/* ------------------------------------------------------------------ */

static inline const double *coord(const double *x, const double *fixed,
                                  int n_free, int slot) {
    if (slot < n_free) return x + 3*slot;
    return fixed + 3*(slot - n_free);
}

/* ------------------------------------------------------------------ */
/* Accumulate gradient contribution                                   */
/* ------------------------------------------------------------------ */

static inline void scatter(double *grad, int slot, int n_free,
                           const double *contrib) {
    if (slot < n_free) {
        grad[3*slot]   += contrib[0];
        grad[3*slot+1] += contrib[1];
        grad[3*slot+2] += contrib[2];
    }
}

/* ------------------------------------------------------------------ */
/* Main cost + gradient function                                      */
/* ------------------------------------------------------------------ */

double cost_and_grad(
    /* Free atom coordinates (n_free * 3, row-major) */
    const double *x,
    /* Output gradient (n_free * 3) */
    double *grad,
    int n_free,
    /* Fixed atom coordinates (n_fixed * 3) */
    const double *fixed,
    int n_fixed,
    /* Bond constraints: pairs[n_bonds*2], ideal_lengths[n_bonds] */
    const int *bond_pairs, const double *ideal_lengths, int n_bonds,
    /* Angle constraints: triples[n_angles*3], ideal_angles_deg[n_angles] */
    const int *angle_triples, const double *ideal_angles_deg, int n_angles,
    /* Planarity: groups[n_planar*4] */
    const int *planar_groups, int n_planar,
    /* Chirality: info[n_chiral*5] (center, a, b, c, sign) */
    const int *chiral_info, int n_chiral,
    /* Dihedral: quads[n_dih*4], targets_deg[n_dih] */
    const int *dih_quads, const double *dih_targets_deg, int n_dih,
    /* EZ: quads[n_ez*4], targets_deg[n_ez] */
    const int *ez_quads, const double *ez_targets_deg, int n_ez,
    /* Weights */
    double w_bond, double w_angle, double w_planar, double w_chiral,
    double w_dih, double w_ez)
{
    double cost = 0.0;
    memset(grad, 0, 3 * n_free * sizeof(double));

    /* --- Bond terms: r = w*(|d| - d0) --- */
    for (int ib = 0; ib < n_bonds; ib++) {
        int si = bond_pairs[2*ib], sj = bond_pairs[2*ib+1];
        const double *ri = coord(x, fixed, n_free, si);
        const double *rj = coord(x, fixed, n_free, sj);
        double d[3];
        sub3(ri, rj, d);
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

    /* --- Angle terms: r = w * 2 * L * sin((theta-theta0)/2) --- */
    for (int ia = 0; ia < n_angles; ia++) {
        int sa = angle_triples[3*ia];
        int sb = angle_triples[3*ia+1];
        int sc = angle_triples[3*ia+2];
        const double *ra = coord(x, fixed, n_free, sa);
        const double *rb = coord(x, fixed, n_free, sb);
        const double *rc = coord(x, fixed, n_free, sc);
        double v1[3], v2[3];
        sub3(ra, rb, v1);
        sub3(rc, rb, v2);
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

        /* dtheta/dra, dtheta/drc */
        double dth_dra[3], dth_drc[3], dth_drb[3];
        for (int k = 0; k < 3; k++) {
            dth_dra[k] = (cos_a * v1h[k] - v2h[k]) / (sin_th * n1);
            dth_drc[k] = (cos_a * v2h[k] - v1h[k]) / (sin_th * n2);
            dth_drb[k] = -(dth_dra[k] + dth_drc[k]);
        }
        /* dr/dra = dr_dtheta * dtheta_dra + dr_dL * 0.5 * v1h */
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

    /* --- Planarity terms: r = w * vol/norm --- */
    for (int ip = 0; ip < n_planar; ip++) {
        int sj = planar_groups[4*ip];
        int sa = planar_groups[4*ip+1];
        int sb = planar_groups[4*ip+2];
        int sc = planar_groups[4*ip+3];
        const double *rj = coord(x, fixed, n_free, sj);
        const double *ra = coord(x, fixed, n_free, sa);
        const double *rb = coord(x, fixed, n_free, sb);
        const double *rc = coord(x, fixed, n_free, sc);
        double v1[3], v2[3], v3[3];
        sub3(ra, rj, v1); sub3(rb, rj, v2); sub3(rc, rj, v3);
        double c23[3], c31[3], c12[3];
        cross3(v2, v3, c23);
        double vol = dot3(v1, c23);
        double n1 = norm3(v1), n2 = norm3(v2), n3 = norm3(v3);
        double nrm = n1 * n2 * n3 + 1e-10;
        double r = w_planar * vol / nrm;
        cost += r * r;

        cross3(v3, v1, c31);
        cross3(v1, v2, c12);
        double inv = 1.0 / nrm;
        double q = vol * inv;

        double dnorm_a[3] = {0,0,0}, dnorm_b[3] = {0,0,0}, dnorm_c[3] = {0,0,0};
        if (n1 > 1e-10) {
            double f = n2 * n3 / n1;
            for (int k = 0; k < 3; k++) dnorm_a[k] = f * v1[k];
        }
        if (n2 > 1e-10) {
            double f = n1 * n3 / n2;
            for (int k = 0; k < 3; k++) dnorm_b[k] = f * v2[k];
        }
        if (n3 > 1e-10) {
            double f = n1 * n2 / n3;
            for (int k = 0; k < 3; k++) dnorm_c[k] = f * v3[k];
        }

        double ga[3], gb[3], gc[3], gj[3];
        double sc2r = 2.0 * r;
        for (int k = 0; k < 3; k++) {
            ga[k] = sc2r * w_planar * (c23[k] * inv - q * dnorm_a[k] * inv);
            gb[k] = sc2r * w_planar * (c31[k] * inv - q * dnorm_b[k] * inv);
            gc[k] = sc2r * w_planar * (c12[k] * inv - q * dnorm_c[k] * inv);
            gj[k] = -(ga[k] + gb[k] + gc[k]);
        }
        scatter(grad, sj, n_free, gj);
        scatter(grad, sa, n_free, ga);
        scatter(grad, sb, n_free, gb);
        scatter(grad, sc, n_free, gc);
    }

    /* --- Chirality terms: r = w * max(0, -sign*vol) --- */
    for (int ic = 0; ic < n_chiral; ic++) {
        int sj = chiral_info[5*ic];
        int sa = chiral_info[5*ic+1];
        int sb = chiral_info[5*ic+2];
        int sc = chiral_info[5*ic+3];
        double sign = (double)chiral_info[5*ic+4];
        const double *rj = coord(x, fixed, n_free, sj);
        const double *ra = coord(x, fixed, n_free, sa);
        const double *rb = coord(x, fixed, n_free, sb);
        const double *rc = coord(x, fixed, n_free, sc);
        double v1[3], v2[3], v3[3], c23[3];
        sub3(ra, rj, v1); sub3(rb, rj, v2); sub3(rc, rj, v3);
        cross3(v2, v3, c23);
        double vol = dot3(v1, c23);
        double raw = -sign * vol;
        if (raw <= 0.0) continue;
        double r = w_chiral * raw;
        cost += r * r;

        double c31[3], c12[3];
        cross3(v3, v1, c31);
        cross3(v1, v2, c12);
        double dr_dvol = -w_chiral * sign;
        double sc2r = 2.0 * r;
        double ga[3], gb[3], gc[3], gj[3];
        for (int k = 0; k < 3; k++) {
            ga[k] = sc2r * dr_dvol * c23[k];
            gb[k] = sc2r * dr_dvol * c31[k];
            gc[k] = sc2r * dr_dvol * c12[k];
            gj[k] = -(ga[k] + gb[k] + gc[k]);
        }
        scatter(grad, sj, n_free, gj);
        scatter(grad, sa, n_free, ga);
        scatter(grad, sb, n_free, gb);
        scatter(grad, sc, n_free, gc);
    }

    /* --- Dihedral terms (shared code for AMSR + EZ) --- */
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
            int s0 = quads[4*id], s1 = quads[4*id+1];
            int s2 = quads[4*id+2], s3 = quads[4*id+3];
            const double *p0 = coord(x, fixed, n_free, s0);
            const double *p1 = coord(x, fixed, n_free, s1);
            const double *p2 = coord(x, fixed, n_free, s2);
            const double *p3 = coord(x, fixed, n_free, s3);

            double actual = measure_torsion(p0, p1, p2, p3);
            double diff = fmod(actual - targets[id] + 540.0, 360.0) - 180.0;
            double r = w * diff;
            cost += r * r;

            /* Blondel-Karplus torsion gradient */
            double b1[3], b2[3], b3[3], n1v[3], n2v[3];
            sub3(p1, p0, b1); sub3(p2, p1, b2); sub3(p3, p2, b3);
            cross3(b1, b2, n1v); cross3(b2, b3, n2v);
            double n1n = norm3(n1v), n2n = norm3(n2v), b2n = norm3(b2);
            if (n1n < 1e-10 || n2n < 1e-10 || b2n < 1e-10) continue;

            /* dt/dp0 = -(b2n/n1n^2) * n1v */
            double f0 = -b2n / (n1n * n1n);
            double f3 = b2n / (n2n * n2n);
            double dt_dp0[3], dt_dp3[3];
            for (int k = 0; k < 3; k++) {
                dt_dp0[k] = f0 * n1v[k];
                dt_dp3[k] = f3 * n2v[k];
            }
            double b1db2 = dot3(b1, b2), b3db2 = dot3(b3, b2);
            double b2sq = b2n * b2n;
            double c1 = b1db2 / b2sq, c2 = b3db2 / b2sq;
            double dt_dp1[3], dt_dp2[3];
            for (int k = 0; k < 3; k++) {
                dt_dp1[k] = -(c1 + 1.0) * dt_dp0[k] + c2 * dt_dp3[k];
                dt_dp2[k] = c1 * dt_dp0[k] - (c2 + 1.0) * dt_dp3[k];
            }
            /* scale: 2*r * w * (180/pi) * dt/dp */
            double sc = 2.0 * w * r * (180.0 / M_PI);
            double g0[3], g1[3], g2[3], g3[3];
            for (int k = 0; k < 3; k++) {
                g0[k] = sc * dt_dp0[k];
                g1[k] = sc * dt_dp1[k];
                g2[k] = sc * dt_dp2[k];
                g3[k] = sc * dt_dp3[k];
            }
            scatter(grad, s0, n_free, g0);
            scatter(grad, s1, n_free, g1);
            scatter(grad, s2, n_free, g2);
            scatter(grad, s3, n_free, g3);
        }
    }

    return cost;
}
