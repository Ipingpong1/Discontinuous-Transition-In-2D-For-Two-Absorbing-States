#!/usr/bin/env python3
"""Pair approximation (PA) and site mean field (MF) for noisy Day & Night.

D&N (B3678/S34678) is self-complementary: a dead cell is born at k in {3,6,7,8} live
neighbours and a live cell dies at k in {0,1,2,5}, i.e. at j = 8 - k in {8,7,6,3}
disagreeing neighbours. So both states flip with the same probability, a function of
j = # disagreeing Moore neighbours only (dnsim.c model_from_rule):

    w(j) = (1-p) alpha [j in {3,6,7,8}] + p j/8 + q [j = qj]          (qj = 2)

PA closure (Bethe lattice, z = 8): given a site's state s, its neighbours are
independent, each disagreeing with prob theta_s = rho / (2 x_s), where x_1 = x (live
density), x_0 = 1 - x and rho = fraction of disagreeing bonds. Variables (m, rho), m = 2x-1.

Continuous time (each site flips at rate w(j); one unit = one dnsim step):
    dx/dt   = (1-x) F_0 - x F_1,                       F_s = sum_j B(j;8,th_s) w(j)
    drho/dt = sum_s [2 n_ss Pc_s - rho Pd_s]
where n_ss = x_s - rho/2 are the concordant-bond fractions and, for a site at one end of
a concordant / discordant bond, the flip probability given its other 7 neighbours is
    Pc_s = sum_k B(k;7,th_s) w(k),   Pd_s = sum_k B(k;7,th_s) w(k+1).
(Equivalent to (2/z) sum_s x_s sum_j B(j;z,th_s) w(j) (z-2j).)

Discrete (synchronous) map, the dnsim update: every site flips at once with prob w(j).
The two ends of a bond flip independently given the bond (other neighbours i.i.d.), so a
bond changes state iff exactly one end flips:
    rho' = rho (1 - Pd_1 - Pd_0 + 2 Pd_1 Pd_0) + sum_s 2 n_ss Pc_s (1 - Pc_s)
    x'   = x - x F_1 + (1-x) F_0

MF: rho slaved to the uncorrelated value rho = 2x(1-x), i.e. th_1 = 1-x, th_0 = x.

Absorbing state (x -> 1, n = 1-x -> 0, theta_0 finite, theta_1 ~ n theta_0), to first order
in n:  dn/dt = n [z w(1) th - F_0(th)],
       dth/dt = (z-2) w(1) th + H(th)/z - th [z w(1) th - F_0(th)],  H = sum B w (z-2j).
Since w(1) = p/8 and F_0(th) >= p th + (rule + noise terms), the growth rate
z w(1) th - F_0(th) < 0 for every p < 1: in PA/MF the absorbing states are always linearly
stable, so an active phase can only appear discontinuously (saddle-node). PA's q* is the
saddle-node (lower spinodal) of the symmetric active fixed point; the simulated q*
(front-velocity zero) should lie at or above it.
"""
import argparse
import sys

import numpy as np
from scipy.optimize import brentq, fsolve
from scipy.special import comb

Z = 8
J8 = np.arange(Z + 1)
J7 = np.arange(Z)
C8 = comb(Z, J8)
C7 = comb(Z - 1, J7)
RULE = np.array([j in (3, 6, 7, 8) for j in J8], float)


def wtable(p, q, alpha, qj=2):
    w = (1 - p) * alpha * RULE + p * J8 / Z
    w[qj] += q
    return np.clip(w, 0, 1)


def b8(th):
    th = np.asarray(th, float)[..., None]
    return C8 * th ** J8 * (1 - th) ** (Z - J8)


def b7(th):
    th = np.asarray(th, float)[..., None]
    return C7 * th ** J7 * (1 - th) ** (Z - 1 - J7)


class Model:
    """w: flip table w[0..8]; disc: synchronous map instead of continuous time."""

    def __init__(self, w, disc=False):
        self.w = np.asarray(w, float)
        if self.w[0] != 0:
            raise ValueError("w(0) must be 0 (absorbing states)")
        self.disc = disc

    def _parts(self, x, rho):
        with np.errstate(divide="ignore", invalid="ignore"):
            th1 = np.clip(np.nan_to_num(rho / (2 * x)), 0.0, 1.0)
            th0 = np.clip(np.nan_to_num(rho / (2 * (1 - x))), 0.0, 1.0)
        w = self.w
        F1, F0 = b8(th1) @ w, b8(th0) @ w
        Pc1, Pc0 = b7(th1) @ w[:-1], b7(th0) @ w[:-1]
        Pd1, Pd0 = b7(th1) @ w[1:], b7(th0) @ w[1:]
        return F1, F0, Pc1, Pc0, Pd1, Pd0

    def rhs(self, m, rho):
        """Continuous: (dm/dt, drho/dt). Discrete: (m'-m, rho'-rho)."""
        x = (1 + m) / 2
        F1, F0, Pc1, Pc0, Pd1, Pd0 = self._parts(x, rho)
        n11, n00 = x - rho / 2, 1 - x - rho / 2
        dx = (1 - x) * F0 - x * F1
        if self.disc:
            drho = (-rho * (Pd1 + Pd0 - 2 * Pd1 * Pd0)
                    + 2 * n11 * Pc1 * (1 - Pc1) + 2 * n00 * Pc0 * (1 - Pc0))
        else:
            drho = 2 * n11 * Pc1 + 2 * n00 * Pc0 - rho * (Pd1 + Pd0)
        return 2 * dx, drho

    def jac(self, m, rho, h=1e-7):
        f = lambda v: np.array(self.rhs(*v))
        J = np.empty((2, 2))
        for i in range(2):
            e = np.zeros(2); e[i] = h
            J[:, i] = (f(np.array([m, rho]) + e) - f(np.array([m, rho]) - e)) / (2 * h)
        return J

    def stable(self, eig):
        eig = np.asarray(eig)
        return bool(np.all(np.abs(1 + eig) < 1)) if self.disc else bool(np.all(eig.real < 0))

    # --- symmetric line m = 0 ---------------------------------------------------------
    def sym_fixed_points(self, ngrid=4000):
        """Fixed points (0, rho*) with rho* in (0, 1]: list of dicts."""
        g = lambda r: self.rhs(0.0, r)[1]
        rs = np.linspace(1e-6, 1.0, ngrid)
        gs = g(rs)
        out = []
        for a, b, ga, gb in zip(rs[:-1], rs[1:], gs[:-1], gs[1:]):
            if ga == 0 or ga * gb < 0:
                r = brentq(g, a, b, xtol=1e-14) if ga != 0 else a
                J = self.jac(0.0, r)
                lm, lr = J[0, 0], J[1, 1]   # decoupled by symmetry at m = 0
                out.append(dict(rho=r, lam_rho=lr, lam_m=lm,
                                stab_rho=self.stable([lr]), stab_m=self.stable([lm])))
        return out

    def asym_fixed_points(self, nm=12, nr=12):
        """Fixed points with m > 0, rho > 0 (ordered active phase), by multistart fsolve."""
        found = []
        for m0 in np.linspace(0.05, 0.95, nm):
            for r0 in np.linspace(0.02, 0.98, nr) * (1 - m0):
                sol, info, ier, _ = fsolve(lambda v: self.rhs(*v), [m0, r0], full_output=True)
                m, r = sol
                if ier != 1 or not (0.01 < m < 0.999 and 1e-6 < r < 1 - m):
                    continue
                if np.max(np.abs(self.rhs(m, r))) > 1e-11:
                    continue
                if any(abs(m - f["m"]) < 1e-6 and abs(r - f["rho"]) < 1e-6 for f in found):
                    continue
                ev = np.linalg.eigvals(self.jac(m, r))
                found.append(dict(m=m, rho=r, eig=ev, stab=self.stable(ev)))
        return found

    # --- slow manifold and front criterion ----------------------------------------------
    def slow_manifold(self, ms, ngrid=1500):
        """For each m, the largest stable root rho_s(m) > 0 of drho = 0 (nan if none)."""
        out = np.full(len(ms), np.nan)
        for i, m in enumerate(ms):
            g = lambda r: self.rhs(m, r)[1]
            rs = np.linspace(1e-9, 1 - abs(m) - 1e-12, ngrid)
            gs = g(rs)
            idx = np.nonzero((gs[:-1] > 0) & (gs[1:] <= 0))[0]   # downward crossings = stable
            if len(idx):
                k = idx[-1]
                out[i] = brentq(g, rs[k], rs[k + 1], xtol=1e-14)
        return out

    def front_integral(self, nm=400):
        """I = int_0^1 dm/dt(m, rho_s(m)) dm along the slow manifold (rho fast). For a
        bistable reaction-diffusion front with constant diffusion, the absorbing phase
        (m = 1) invades the active one (m = 0) iff I > 0."""
        ms = np.linspace(0, 1, nm + 1)[:-1] + 0.5 / nm
        rs = self.slow_manifold(ms)
        if np.any(np.isnan(rs)):
            return np.nan
        return float(np.mean(self.rhs(ms, rs)[0]))

    # --- absorbing state ------------------------------------------------------------
    def absorbing_rate(self, n=1e-9, iters=20000):
        """Growth of a dilute minority (n -> 0) from an isolated-seed shape (theta_0 = 1).
        Continuous: rate lambda (stable < 0). Discrete: per-step factor - 1 (stable < 0)."""
        th = 1.0
        lam = None
        dt = 1.0 if self.disc else 0.05
        for _ in range(iters):
            x, rho = 1 - n, 2 * n * th
            dm, drho = self.rhs(2 * x - 1, rho)
            dn = -dm / 2
            n2, rho2 = n + dt * dn, rho + dt * drho
            th_new = rho2 / (2 * n2)
            lam_new = dn / n
            if lam is not None and abs(lam_new - lam) < 1e-13 and abs(th_new - th) < 1e-12:
                return lam_new, th_new
            lam, th = lam_new, min(max(th_new, 1e-12), 1.0)
        return lam, th


class MeanField:
    """Site MF: theta_1 = 1-x, theta_0 = x. Only m evolves."""

    def __init__(self, w, disc=False):
        self.w, self.disc = np.asarray(w, float), disc

    def dm(self, m):
        x = (1 + m) / 2
        F1, F0 = b8(1 - x) @ self.w, b8(x) @ self.w
        return 2 * ((1 - x) * F0 - x * F1)

    def sym_slope(self, h=1e-6):
        """d(dm)/dm at m = 0: < 0 (cont) / in (-2, 0) (disc) => symmetric state stable."""
        return (self.dm(h) - self.dm(-h)) / (2 * h)

    def sym_stable(self):
        s = self.sym_slope()
        return (-2 < s < 0) if self.disc else s < 0


# --- thresholds -----------------------------------------------------------------------
def pa_status(p, q, alpha, disc):
    M = Model(wtable(p, q, alpha), disc)
    fps = [f for f in M.sym_fixed_points() if f["stab_rho"]]
    full = [f for f in fps if f["stab_m"]]
    return M, fps, full


def _bisect(pred, lo, hi, tol):
    """Smallest q in [lo, hi] with pred(q) True, assuming pred(lo) False, pred(hi) True."""
    while hi - lo > tol * hi:
        mid = np.sqrt(lo * hi) if lo > 0 else hi / 2
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return hi


def q_front(p, alpha, disc=False, qmax=0.2, n=60):
    """Smallest q where the front integral turns negative (active phase invades), by a
    log scan then brentq. None if no + -> - change in (0, qmax]."""
    I = lambda q: Model(wtable(p, q, alpha), disc).front_integral(200)
    qs = np.concatenate([[0.0], np.geomspace(1e-5, qmax, n)])
    prev_q, prev_I = None, None
    for q in qs:
        v = I(q)
        if np.isnan(v):
            break
        if prev_I is not None and prev_I > 0 >= v:
            return brentq(I, prev_q, q, xtol=1e-7)
        if prev_I is None and v <= 0:
            return 0.0
        prev_q, prev_I = q, v
    return None


def thresholds(p, alpha, disc=False, qmax=0.5, tol=1e-4):
    """q_sn: lowest q with a rho-stable symmetric active FP (saddle-node / lower spinodal).
    q_full: lowest q where that FP is also stable against m (no Ising ordering).
    q_mf: lowest q where the MF symmetric state m = 0 is stable.
    q_front: front criterion on the slow manifold (see Model.front_integral).
    rho0: rho* of the symmetric active FP at q = q_front (or q_sn)."""
    qs = np.concatenate([[0.0], np.geomspace(1e-6, qmax, 80)])
    has_rho = lambda q: bool(pa_status(p, q, alpha, disc)[1])
    has_full = lambda q: bool(pa_status(p, q, alpha, disc)[2])
    has_mf = lambda q: MeanField(wtable(p, q, alpha), disc).sym_stable()
    res = {}
    for name, pred in (("q_sn", has_rho), ("q_full", has_full), ("q_mf", has_mf)):
        if pred(0.0):
            res[name] = 0.0
            continue
        vals = [pred(q) for q in qs]
        if not any(vals):
            res[name] = None
        else:
            i = vals.index(True)
            res[name] = _bisect(pred, qs[i - 1], qs[i], tol)
    res["q_front"] = q_front(p, alpha, disc)
    qref = res["q_front"] if res["q_front"] is not None else res["q_sn"]
    f = pa_status(p, qref, alpha, disc)[1] if qref is not None else []
    res["rho0"] = max(x["rho"] for x in f) if f else None
    return res


# --- simulation reference (campaign/README.md) -----------------------------------------
# (p, alpha, q*): Phase 1 E1/E3 fronts (p = 0), E4 fronts and the p = 0.995 strip scan
# (alpha = 0.5); p = 0.1, 0.25 are one-seed local pilots (approximate).
SIM = [
    (0.0, 0.125, 0.0062), (0.0, 0.25, 0.0108), (0.0, 0.5, 0.01837), (0.0, 0.75, 0.038),
    (0.0, 1.0, 0.084),
    (0.1, 0.5, 0.036), (0.25, 0.5, 0.047), (0.5, 0.5, 0.0459), (0.75, 0.5, 0.0276),
    (0.9, 0.5, 0.0119), (0.95, 0.5, 0.0062), (0.98, 0.5, 0.0025), (0.99, 0.5, 0.00127),
    (0.995, 0.5, 0.00062),
]


def fmt(v, spec=".5f"):
    return "—" if v is None else format(v, spec)


def cmd_point(a):
    w = wtable(a.p, a.q, a.alpha)
    print("w(j) =", " ".join(f"{v:.5g}" for v in w))
    for disc in (False, True):
        M = Model(w, disc)
        print(f"\n[{'discrete map' if disc else 'continuous time'}]")
        sym = M.sym_fixed_points()
        if not sym:
            print("  symmetric active FP: none")
        for f in sym:
            print(f"  symmetric FP rho*={f['rho']:.5f}  lam_rho={f['lam_rho']:+.4g} "
                  f"({'stable' if f['stab_rho'] else 'unstable'})  lam_m={f['lam_m']:+.4g} "
                  f"({'stable' if f['stab_m'] else 'UNSTABLE: Ising'})")
        for f in M.asym_fixed_points():
            print(f"  ordered FP m={f['m']:.4f} rho={f['rho']:.5f} eig={np.round(f['eig'], 5)} "
                  f"({'stable' if f['stab'] else 'unstable'})")
        lam, th = M.absorbing_rate()
        print(f"  absorbing: growth {lam:+.4g} at theta0={th:.4f} "
              f"({'stable' if lam < 0 else 'UNSTABLE'})")
        mf = MeanField(w, disc)
        print(f"  MF: d(dm)/dm at m=0 = {mf.sym_slope():+.4g} "
              f"({'m=0 stable' if mf.sym_stable() else 'm=0 unstable'})")


KEYS = ("q_sn", "q_full", "q_mf", "q_front", "rho0")


def cmd_compare(a):
    rows = []
    print(f"{'':22}| continuous time{'':36}| discrete map")
    print(f"{'p':>6} {'alpha':>6} {'q*_sim':>8} |" + 2 * (
        f" {'q_sn':>7} {'q_full':>7} {'q_mf':>7} {'q_front':>8} {'q*/qf':>6} {'rho0':>6} |"))
    for p, al, qs in SIM:
        c, d = thresholds(p, al, False), thresholds(p, al, True)
        rows.append((p, al, qs, c, d))
        line = f"{p:6.3f} {al:6.3f} {qs:8.5f} |"
        for r in (c, d):
            ratio = qs / r["q_front"] if r["q_front"] else None
            line += (f" {fmt(r['q_sn'], '7.4f')} {fmt(r['q_full'], '7.4f')} {fmt(r['q_mf'], '7.4f')}"
                     f" {fmt(r['q_front'], '8.5f')} {fmt(ratio, '6.2f')} {fmt(r['rho0'], '6.3f')} |")
        print(line, flush=True)
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("p,alpha,q_sim," + ",".join(f"{m}_{k}" for m in ("cont", "disc") for k in KEYS) + "\n")
            for p, al, qs, c, d in rows:
                vals = [p, al, qs] + [r.get(k) for r in (c, d) for k in KEYS]
                f.write(",".join("" if v is None else f"{v:.6g}" for v in vals) + "\n")


def cmd_curve(a):
    """q_front along p at fixed alpha, for plotting."""
    ps = np.concatenate([np.linspace(0, 0.9, 19), [0.95, 0.98, 0.99, 0.995]])
    with open(a.csv, "w") as f:
        f.write("p,alpha,cont_q_front,disc_q_front\n")
        for p in ps:
            c, d = q_front(p, a.alpha, False), q_front(p, a.alpha, True)
            f.write(",".join("" if v is None else f"{v:.6g}" for v in (p, a.alpha, c, d)) + "\n")
            print(f"p {p:.3f}: q_front cont {fmt(c)} disc {fmt(d)}", flush=True)


def cmd_selftest(a):
    ok = True
    # voter (p = 1): m conserved, PA rho* = (z-2)/(2(z-1)) (1-m^2), absorbing marginal
    for disc in (False, True):
        M = Model(wtable(1.0, 0.0, 0.5), disc)
        for m in (0.0, 0.3, -0.6):
            dm, _ = M.rhs(m, 0.2)
            ok &= abs(dm) < 1e-14
        if not disc:
            r = brentq(lambda r: M.rhs(0.3, r)[1], 1e-3, 0.9)
            ok &= abs(r - 6 / 14 * (1 - 0.09)) < 1e-10
            lam, th = M.absorbing_rate()
            ok &= abs(lam) < 1e-9 and abs(th - 6 / 7) < 1e-6
    # drho/dt identity: 2 n11 Pc1 - rho Pd1 + ... == (2/z) sum_s x_s sum_j B w (z-2j)
    rng = np.random.default_rng(1)
    for _ in range(20):
        w = np.r_[0, rng.random(8)]
        M = Model(w)
        m = rng.uniform(-0.8, 0.8)
        rho = rng.uniform(0, 1 - abs(m))
        x = (1 + m) / 2
        ref = 2 / Z * sum(xs * b8(rho / (2 * xs)) @ (w * (Z - 2 * J8)) for xs in (x, 1 - x))
        ok &= abs(M.rhs(m, rho)[1] - ref) < 1e-12
        # discrete -> continuous for small rates
        Md, Mc = Model(w * 1e-6, True), Model(w * 1e-6)
        ok &= abs(Md.rhs(m, rho)[1] - Mc.rhs(m, rho)[1]) < 1e-4 * abs(Mc.rhs(m, rho)[1])  # O(w) relative
    # D&N rule is self-complementary in j
    births = {3, 6, 7, 8}
    deaths = {8 - k for k in (0, 1, 2, 5)}
    ok &= births == deaths == {j for j in J8 if RULE[j]}
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("point", help="fixed points and stability at one (p, q, alpha)")
    s.add_argument("--p", type=float, default=0.0)
    s.add_argument("--q", type=float, required=True)
    s.add_argument("--alpha", type=float, default=0.5)
    s.set_defaults(fn=cmd_point)
    s = sub.add_parser("compare", help="PA/MF thresholds vs simulated q*")
    s.add_argument("--csv")
    s.set_defaults(fn=cmd_compare)
    s = sub.add_parser("curve", help="thresholds along p at fixed alpha")
    s.add_argument("--alpha", type=float, default=0.5)
    s.add_argument("--csv", required=True)
    s.set_defaults(fn=cmd_curve)
    s = sub.add_parser("selftest")
    s.set_defaults(fn=cmd_selftest)
    a = ap.parse_args()
    sys.exit(a.fn(a) or 0)


if __name__ == "__main__":
    main()
