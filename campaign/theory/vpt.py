#!/usr/bin/env python3
"""Voter-model perturbation theory for noisy D&N near p = 1 (and for other rules, see below).

Write the dnsim flip table as a perturbation of the (synchronous, Moore) voter model,
    w(j) = j/8 + eps * g(j),   eps = 1 - p,
    g(j) = alpha R(j) - j/8 + kappa [j = 2],   R(j) = [j in {3,6,7,8}],   kappa = q / (1 - p).
To first order in eps the local density u obeys du/dt = (diffusion) + eps f(u) with
    f(u) = E_u[ (1 - xi_0) g(j_0) - xi_0 g(j_0) ]
(origin flips 0 -> 1 at rate g(# live nbrs), 1 -> 0 at rate g(# dead nbrs)), the expectation
taken in the voter state at density u. Through the dual, the 9 sites (origin + Moore
neighbours) carry the colours of the coalescing lineages that reach them (cwalk.c); each
class is independently live with prob u. Hence, for a pattern with classes C_0 (origin),
C_1..C_m holding n_1..n_m neighbours,
    f(u) = sum_patterns P(pattern) sum_colourings u^#live (1-u)^#dead * (+-) g(j_0).
Properties (checked in --selftest): g linear in j gives f = 0 (voter conserves density);
1- and 2-class patterns cancel by symmetry; 3-class patterns give
    u(1-u)(1-2u) [g(n_1) + g(n_2) - g(n_1 + n_2)].
So f = alpha f_R + kappa f_2 and every threshold is kappa* = alpha * (pure number):
  kappa_a      from f'(0) = 0 (absorbing state u = 0 loses stability: GV-like / 'a' = 0),
  kappa_front  from int_0^{1/2} f du = 0 (front between u = 0 and the active state u = 1/2
               is stationary: the first-order transition of the reaction-diffusion limit).
In 2D the lineages all coalesce eventually, so P(pattern) depends on the walk time T ~ 1/eps
(the time the perturbation needs to act); the ratios converge as T -> inf with corrections
in powers of 1/ln T.

Neighbourhood size z (10-30). Everything above holds for a site with z neighbours and any flip
table w(j) = (voter) + g(j), j = 0..z: pattern_poly and drift_basis take z (default 8, Moore),
and `run --z 4` uses cwalk5.c, the dual of the synchronous five-cell (von Neumann + self) voter
model of Molofsky et al., Theor. Popul. Biol. 55, 270 (1999). Their rule is
P(next = 1 | n ones among the 5 cells) = p_n with p_0 = 0, p_5 = 1, p_3 = 1 - p_2, p_4 = 1 - p_1;
as a flip table in the number j of disagreeing neighbours it is w(j) = p_j for both states,
with voter point (p_1, p_2) = (1/5, 2/5), w(j) = j/5. So
    g = (0, d1, d2, -d2, -d1),     d1 = p_1 - 1/5,   d2 = p_2 - 2/5,
and f = d1 G_1 + d2 G_2 with G_i(u) = u(1-u)(1-2u) [a_i + c_i u(1-u)] exactly: 3- and 4-class
patterns give the cubic, the single 5-class pattern (1,1,1,1), weight P_5, gives the mean-field
quintic 5 d1 u(1-u)^4 - ... (c_1 = -5 P_5, c_2 = 10 P_5). With a = d1 a_1 + d2 a_2 and
c = d1 c_1 + d2 c_2: f'(0) = a, f'(1/2) = -(a + c/4)/2, int_0^1/2 f = (6a + c)/192. The
`molofsky` command prints, per T, the rays in the (d1, d2) plane where u = 0 is marginal
(a = 0, on the side c > 0), where the front is stationary (6a + c = 0) and where u = 1/2 is
marginal (a + c/4 = 0); both states are stable between the first and the last. See MOLOFSKY.md.
"""
import argparse
import itertools
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from numpy.polynomial import polynomial as P

HERE = Path(__file__).resolve().parent
Z = 8                                            # default neighbourhood: Moore
J = np.arange(9)
R = np.array([j in (3, 6, 7, 8) for j in J], float)
D2 = (J == 2).astype(float)
H3 = np.array([0.0, 1.0, -3.0, 2.0])             # u(1-u)(1-2u), coeffs low -> high
WALKERS = {8: ("cwalk", "cwalk.c"), 4: ("cwalk5", "cwalk5.c")}
NBRS = {8: [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)],
        4: [(1, 0), (-1, 0), (0, 1), (0, -1)]}
MOVES = {8: NBRS[8], 4: [(0, 0)] + NBRS[4]}      # one dual step: Moore jump / lazy VN step


def read_patterns(path):
    """{T: (trials, [(count, (n_1, ..., n_m)), ...])}"""
    out, cur = {}, None
    for line in open(path):
        a = line.split()
        if a[0] == "T":
            cur = int(a[1]); out[cur] = (int(a[3]), [])
        else:
            out[cur][1].append((int(a[0]), tuple(int(x) for x in a[1].split(","))))
    return out


def check_z(data, z):
    big = max((sum(ns) for _, pats in data.values() for _, ns in pats), default=0)
    if big > z:
        sys.exit(f"pattern file holds {big} neighbours but z = {z}")


def pattern_poly(ns, z=Z):
    """F[j] = polynomial in u (coeffs, low->high) of sum_colourings weight * sign * [j0 = j]
    for one pattern; the origin's class is the extra class."""
    m = len(ns)
    F = np.zeros((z + 1, m + 2))
    for cols in itertools.product((0, 1), repeat=m + 1):
        s0, rest = cols[0], cols[1:]
        j0 = sum(n for n, s in zip(ns, rest) if s != s0)
        nlive = sum(cols)
        w = P.polymul(P.polypow([0, 1], nlive), P.polypow([1, -1], m + 1 - nlive))
        F[j0, :len(w)] += (1 if s0 == 0 else -1) * w
    return F


def drift_basis(pats, trials, z=Z):
    """F[j](u): f(u) = sum_j g(j) F[j](u)."""
    F = np.zeros((z + 1, z + 2))
    for cnt, ns in pats:
        Fp = pattern_poly(ns, z)
        F[:, :Fp.shape[1]] += cnt * Fp
    return F / trials


def f_of(F, g):
    return P.polytrim(g @ F, 1e-300) if np.any(g @ F) else np.zeros(1)


def thresholds(F, alpha=1.0):
    """kappa/alpha at which f'(0) = 0 and at which int_0^1/2 f = 0, plus u = 1/2 stability
    slope there. Returns dict. (D&N rule, z = 8.)"""
    fR, f2 = R @ F, D2 @ F
    d0 = lambda c: P.polyval(0.0, P.polyder(c))
    integ = lambda c: P.polyval(0.5, P.polyint(c))
    res = dict(kappa_a=-d0(fR) / d0(f2), kappa_front=-integ(fR) / integ(f2))
    for key in ("kappa_a", "kappa_front"):
        f = fR + res[key] * f2
        res[key + "_slope_half"] = P.polyval(0.5, P.polyder(f))
    res["fR'(0)"], res["f2'(0)"] = d0(fR), d0(f2)
    return res


def kappa_theta3(pats):
    """kappa/alpha from 3-class patterns only: the finite-T estimate of Cox & Perkins'
    (arXiv:2310.00812, eq. 1.19) Theta_3 = 0, whose K_3 weights are the T -> inf limits of the
    conditional 3-class pattern frequencies."""
    br = sum(c * (R[a] + R[b] - R[a + b]) for c, ns in pats if len(ns) == 2 for a, b in [ns])
    b2 = sum(c * (D2[a] + D2[b] - D2[a + b]) for c, ns in pats if len(ns) == 2 for a, b in [ns])
    return -br / b2 if b2 else float("nan")


def run_cwalk(trials, seed, cks, out, z=Z):
    exe = HERE / WALKERS[z][0]
    src = HERE / WALKERS[z][1]
    if not exe.exists() or exe.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(["cc", "-O3", "-march=native", "-fopenmp", "-o", str(exe), str(src)], check=True)
    with open(out, "w") as f:
        subprocess.run([str(exe), str(trials), str(seed)] + [str(c) for c in cks], stdout=f, check=True)


def exact_t1(z=4):
    """Exact probabilities of the >= 3-class patterns after one dual step, by enumeration of
    all (z+1)^(z+1) move combinations (z = 4 only: 3125 terms)."""
    moves, start = MOVES[z], [(0, 0)] + NBRS[z]
    cnt = Counter()
    for mv in itertools.product(range(len(moves)), repeat=len(start)):
        pos = [(x + moves[m][0], y + moves[m][1]) for (x, y), m in zip(start, mv)]
        others = Counter(p for p in pos[1:] if p != pos[0])
        ns = tuple(sorted(others.values(), reverse=True))
        if len(ns) >= 2:
            cnt[ns] += 1
    tot = len(moves) ** len(start)
    return {ns: c / tot for ns, c in cnt.items()}


# ----------------------------------------------------------------------------- Molofsky rule
def molofsky_g(d1, d2):
    """Flip-table perturbation of the five-cell rule: g(j) = w(j) - j/5, j = 0..4."""
    return np.array([0.0, d1, d2, -d2, -d1])


def molofsky_mf(d1, d2):
    """Mean-field increment x_{t+1} - x_t of Molofsky's map (their eq. 1) at p = voter + g,
    as a polynomial in x: sum_n C(5,n) g(n) x^n (1-x)^(5-n) (the voter part gives x)."""
    from math import comb
    g = molofsky_g(d1, d2)
    f = np.zeros(6)
    for n in range(1, 5):
        w = P.polymul(P.polypow([0, 1], n), P.polypow([1, -1], 5 - n))
        f[:len(w)] += comb(5, n) * g[n] * w
    return f


def molofsky_coeffs(F):
    """f = d1 G_1 + d2 G_2, G_i = h(u) [a_i + c_i u(1-u)], h = u(1-u)(1-2u).
    Returns (a1, c1, a2, c2, resid): resid is the largest violation of that form."""
    out, resid = [], 0.0
    for G in (F[1] - F[4], F[2] - F[3]):
        q, r = P.polydiv(G, H3)
        q = np.pad(q, (0, 3))[:3]
        out += [q[0], q[1]]
        resid = max(resid, float(np.abs(r).max()), abs(q[1] + q[2]))
    return (*out, resid)


def molofsky_rays(a1, c1, a2, c2):
    """Directions (degrees in the (d1, d2) plane, d1 = p_1 - 1/5 along x, d2 = p_2 - 2/5 along y)
    of the three rays: u = 0 marginal (a = 0, taken on the side c > 0 where u = 1/2 is stable),
    front stationary (6a + c = 0) and u = 1/2 marginal (a + c/4 = 0), the last two on the side
    a < 0. Both u = 0 and u = 1/2 are stable between the first and the last."""
    def ray(A, B, side):                      # line A x + B y = 0, direction with side > 0
        x, y = -B, A
        if side(x, y) < 0:
            x, y = B, -A
        return float(np.degrees(np.arctan2(y, x)) % 360)
    a = lambda x, y: a1 * x + a2 * y
    c = lambda x, y: c1 * x + c2 * y
    return (ray(a1, a2, c), ray(6 * a1 + c1, 6 * a2 + c2, lambda x, y: -a(x, y)),
            ray(a1 + c1 / 4, a2 + c2 / 4, lambda x, y: -a(x, y)))


MF_RAYS = molofsky_rays(5.0, -5.0, 0.0, 10.0)                 # T = 0: 90, 111.80, 123.69 deg
DYAD_RAY = float(np.degrees(np.arctan2(3.0, -1.0)))           # 3 p_1 + p_2 = 1: 108.43 deg


def molofsky_rows(data, z=4):
    rows = []
    for T in sorted(data):
        trials, pats = data[T]
        ntot = sum(c for c, _ in pats)
        if ntot == 0:
            continue
        n5 = sum(c for c, ns in pats if len(ns) == 4)
        n3 = sum(c for c, ns in pats if len(ns) == 2)
        a1, c1, a2, c2, resid = molofsky_coeffs(drift_basis(pats, trials, z))
        P5 = n5 / trials
        th = molofsky_rays(a1, c1, a2, c2)
        g1, g2 = molofsky_g(1, 0), molofsky_g(0, 1)          # 3-class brackets only (cf. kappa_theta3)
        b31 = sum(c * (g1[ns[0]] + g1[ns[1]] - g1[ns[0] + ns[1]]) for c, ns in pats if len(ns) == 2)
        b32 = sum(c * (g2[ns[0]] + g2[ns[1]] - g2[ns[0] + ns[1]]) for c, ns in pats if len(ns) == 2)
        rows.append(dict(T=T, p3=ntot / trials, f3=n3 / ntot, f5=n5 / ntot, P5=P5,
                         a1=a1, a2=a2, c1=c1, c2=c2, b1=a1 - 5 * P5, b2=a2, resid=resid,
                         th_a=th[0], th_f=th[1], th_h=th[2], width=th[2] - th[0],
                         th_3=float(np.degrees(np.arctan2(b31, -b32)) % 360) if n3 else float("nan")))
    return rows


def molofsky_boot(data, nb, z=4, seed=1):
    """Multinomial bootstrap of the histogram at each T: sd of the three ray angles."""
    rng = np.random.default_rng(seed)
    out = {}
    for T in sorted(data):
        trials, pats = data[T]
        cnt = np.array([c for c, _ in pats] + [trials - sum(c for c, _ in pats)])
        if cnt[:-1].sum() < 50:
            continue
        Fs = [pattern_poly(ns, z) for _, ns in pats]
        th = []
        for _ in range(nb):
            c = rng.multinomial(trials, cnt / trials)[:-1]
            F = np.zeros((z + 1, z + 2))
            for ci, Fp in zip(c, Fs):
                F[:, :Fp.shape[1]] += ci * Fp
            a1, c1, a2, c2, _ = molofsky_coeffs(F / trials)
            r = molofsky_rays(a1, c1, a2, c2)
            th.append(r + (r[2] - r[0],))
        out[T] = np.std(th, axis=0)
    return out


def cmd_molofsky(a):
    data = read_patterns(a.load)
    check_z(data, 4)
    rows = molofsky_rows(data)
    print("Molofsky five-cell rule: f = d1 G_1 + d2 G_2, G_i = u(1-u)(1-2u)[a_i + c_i u(1-u)]")
    print(f"rays (deg) in the (d1, d2) plane; mean field (T = 0): a = 0 at {MF_RAYS[0]:.2f}, "
          f"front at {MF_RAYS[1]:.2f}, u = 1/2 marginal at {MF_RAYS[2]:.2f}; dyad line "
          f"3 p_1 + p_2 = 1 at {DYAD_RAY:.2f}")
    print(f"{'T':>7} {'P(>=3 cl)':>9} {'P5/P(>=3)':>9} {'P5':>9} {'a1':>8} {'a2':>8} {'b1':>8} "
          f"{'b2':>8} {'c2':>8} {'th_a':>7} {'th_front':>8} {'th_half':>7} {'width':>6} {'th_3cl':>7} {'resid':>8}")
    for r in rows:
        print(f"{r['T']:7d} {r['p3']:9.4g} {r['f5']:9.4f} {r['P5']:9.4g} {r['a1']:8.4f} {r['a2']:8.4f} "
              f"{r['b1']:8.4f} {r['b2']:8.4f} {r['c2']:8.4f} {r['th_a']:7.2f} {r['th_f']:8.2f} "
              f"{r['th_h']:7.2f} {r['width']:6.2f} {r['th_3']:7.2f} {r['resid']:8.1e}")
    if 1 in data:
        trials, pats = data[1]
        ex = exact_t1(4)
        print("\nT = 1 check against exact enumeration (pattern: sampled, exact, pull):")
        for ns, p in sorted(ex.items(), key=lambda kv: -kv[1]):
            n = next((c for c, q in pats if q == ns), 0)
            sd = np.sqrt(p * (1 - p) / trials)
            print(f"  {','.join(map(str, ns)):>8}: {n / trials:.5f} {p:.5f} {(n / trials - p) / sd:+.2f}")
    if a.boot:
        sd = molofsky_boot(data, a.boot)
        print(f"\nbootstrap ({a.boot}): T, sd of th_a, th_front, th_half, width (deg)")
        for T, s in sd.items():
            print(f"{T:7d} {s[0]:.3f} {s[1]:.3f} {s[2]:.3f} {s[3]:.4f}")
    if a.csv:
        sd = sd if a.boot else {}
        with open(a.csv, "w") as f:
            f.write("T,P3,f3,f5,P5,a1,a2,b1,b2,c1,c2,th_a,th_f,th_h,width,th_3,sd_a,sd_f,sd_h,sd_w\n")
            for r in rows:
                s = sd.get(r["T"], (float("nan"),) * 4)
                f.write(",".join(f"{r[k]:.6g}" for k in ("T", "p3", "f3", "f5", "P5", "a1", "a2", "b1",
                                                           "b2", "c1", "c2", "th_a", "th_f", "th_h", "width", "th_3"))
                        + f",{s[0]:.4g},{s[1]:.4g},{s[2]:.4g},{s[3]:.4g}\n")
    return 0


# ------------------------------------------------------------------------------------ commands
def cmd_selftest(a):
    ok = True
    rng = np.random.default_rng(0)
    for _ in range(30):                      # random patterns: voter cancels, symmetry
        m = rng.integers(1, 6)
        ns = tuple(sorted(rng.integers(1, 3, m), reverse=True))
        if sum(ns) > 8:
            continue
        F = pattern_poly(ns)
        ok &= np.allclose(J @ F, 0)
        for g in (R, D2, rng.random(9) * (J > 0)):
            c = g @ F
            us = np.linspace(0, 1, 7)
            ok &= np.allclose(P.polyval(us, c), -P.polyval(1 - us, c))
    for ns in ((1,), (3,), (8,)):            # 2-class patterns cancel
        ok &= np.allclose(R @ pattern_poly(ns), 0)
    for n1, n2 in ((1, 2), (2, 3), (3, 3), (1, 5)):   # 3-class closed form
        g = rng.random(9); g[0] = 0
        c = g @ pattern_poly((n1, n2))
        ref = P.polymul([0, 1, -1], [1, -2]) * (g[n1] + g[n2] - g[n1 + n2])
        ok &= np.allclose(np.pad(c, (0, 5))[:5], np.pad(ref, (0, 5))[:5])
    for ns in ((1, 1, 1), (2, 1, 1), (3, 2, 1), (4, 2, 2)):   # 4-class: same cubic
        g = rng.random(9); g[0] = 0
        c = g @ pattern_poly(ns)
        ref = H3 * (sum(g[n] for n in ns) - g[sum(ns)])
        ok &= np.allclose(np.pad(c, (0, 6))[:6], np.pad(ref, (0, 6))[:6])
    # z = 4 (five-cell rule): the same properties with 4 neighbours
    J4 = np.arange(5)
    pats4 = [(1, 1, 1, 1), (2, 1, 1), (1, 1, 1), (3, 1), (2, 2), (2, 1), (1, 1), (4,), (3,), (2,), (1,)]
    for ns in pats4:
        F = pattern_poly(ns, 4)
        ok &= np.allclose(J4 @ F, 0)
        g = rng.random(5); g[0] = 0
        c = g @ F
        us = np.linspace(0, 1, 7)
        ok &= np.allclose(P.polyval(us, c), -P.polyval(1 - us, c))
        if len(ns) == 1:
            ok &= np.allclose(c, 0)
        if len(ns) in (2, 3):
            ref = H3 * (sum(g[n] for n in ns) - g[sum(ns)])
            ok &= np.allclose(np.pad(c, (0, 6))[:6], np.pad(ref, (0, 6))[:6])
    for d1, d2 in ((0.1, -0.05), (-0.03, 0.2), (0.0, 1.0), (1.0, 0.0)):
        g = molofsky_g(d1, d2)
        ok &= np.allclose(g @ pattern_poly((1, 1, 1, 1), 4), molofsky_mf(d1, d2))   # 5-class = MF
        ok &= np.allclose(molofsky_mf(d1, d2), P.polymul(H3, [5 * d1, 10 * d2 - 5 * d1, 5 * d1 - 10 * d2]))
    a1, c1, a2, c2, resid = molofsky_coeffs(drift_basis([(1, (1, 1, 1, 1))], 1, 4))
    ok &= np.allclose((a1, c1, a2, c2), (5, -5, 0, 10)) and resid < 1e-12       # T = 0 is mean field
    ok &= np.allclose(MF_RAYS, (90.0, 180 - np.degrees(np.arctan(2.5)), 180 - np.degrees(np.arctan(1.5))))
    cnt = rng.integers(1, 1000, len(pats4))                                      # random histogram
    pats = list(zip(cnt, pats4))
    tr = cnt.sum() + 500
    a1, c1, a2, c2, resid = molofsky_coeffs(drift_basis(pats, tr, 4))
    P5 = cnt[0] / tr
    br = {ns: (sum(molofsky_g(1, 0)[n] for n in ns) - molofsky_g(1, 0)[sum(ns)],
               sum(molofsky_g(0, 1)[n] for n in ns) - molofsky_g(0, 1)[sum(ns)]) for ns in pats4 if 2 <= len(ns) <= 3}
    b1 = sum(c * br[ns][0] for c, ns in pats if ns in br) / tr
    b2 = sum(c * br[ns][1] for c, ns in pats if ns in br) / tr
    ok &= resid < 1e-10 and np.allclose((a1, c1, a2, c2), (5 * P5 + b1, -5 * P5, b2, 10 * P5))
    ex = exact_t1(4)                                                             # one dual step
    ok &= abs(sum(ex.values()) - 0.99072) < 1e-12 and abs(ex[(1, 1, 1, 1)] - 1458 / 3125) < 1e-12
    print("selftest", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def cmd_run(a):
    cks = [int(x) for x in a.T.split(",")]
    z = a.z
    if not a.load:
        run_cwalk(a.trials, a.seed, cks, a.out, z)
    data = read_patterns(a.load or a.out)
    check_z(data, z)
    if z != 8:
        print(f"z = {z}: generic pattern statistics; D&N thresholds need z = 8")
        print(f"{'T':>8} {'P(>=3 cl)':>10} {'P(3 cl)/P(>=3)':>14} {'P(z+1 cl)/P(>=3)':>16}")
        for T in sorted(data):
            trials, pats = data[T]
            ntot = sum(c for c, _ in pats)
            if ntot:
                print(f"{T:8d} {ntot / trials:10.4g} {sum(c for c, ns in pats if len(ns) == 2) / ntot:14.4f} "
                      f"{sum(c for c, ns in pats if len(ns) == z) / ntot:16.4f}")
        if z == 4:
            print()
            a.load = a.load or a.out
            a.csv = None
            cmd_molofsky(a)
        return 0
    print(f"{'T':>8} {'P(>=3 cl)':>10} {'P(3 cl)/P(>=3)':>14} {'kappa_a/alpha':>13} "
          f"{'kappa_front/alpha':>17} {'kappa_th3/alpha':>15} {'f_R\'(0)':>10} {'f_2\'(0)':>10}")
    rows = []
    for T in sorted(data):
        trials, pats = data[T]
        n3 = sum(c for c, ns in pats if len(ns) == 2)
        ntot = sum(c for c, _ in pats)
        if ntot == 0:
            continue
        th = thresholds(drift_basis(pats, trials))
        rows.append((T, ntot / trials, th))
        print(f"{T:8d} {ntot / trials:10.4g} {n3 / ntot:14.4f} {th['kappa_a']:13.4f} "
              f"{th['kappa_front']:17.4f} {kappa_theta3(pats):15.4f} {th['fR\'(0)']:10.3g} {th['f2\'(0)']:10.3g}")
    if a.boot:
        boot(data, a.boot)
    return 0


def boot(data, nb):
    """Multinomial bootstrap of the pattern histogram at each T (errors on kappa)."""
    rng = np.random.default_rng(1)
    print(f"\nbootstrap ({nb}): T, kappa_a/alpha +- sd, kappa_front/alpha +- sd")
    for T in sorted(data):
        trials, pats = data[T]
        cnt = np.array([c for c, _ in pats] + [trials - sum(c for c, _ in pats)])
        if cnt[:-1].sum() < 50:
            continue
        Fs = [pattern_poly(ns) for _, ns in pats]
        ka, kf = [], []
        for _ in range(nb):
            c = rng.multinomial(trials, cnt / trials)[:-1]
            F = np.zeros((9, 10))
            for ci, Fp in zip(c, Fs):
                F[:, :Fp.shape[1]] += ci * Fp
            th = thresholds(F / trials)
            ka.append(th["kappa_a"]); kf.append(th["kappa_front"])
        print(f"{T:8d} {np.mean(ka):.4f} +- {np.std(ka):.4f}   {np.mean(kf):.4f} +- {np.std(kf):.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("run")
    s.add_argument("--trials", type=int, default=10**6)
    s.add_argument("--seed", type=int, default=1)
    s.add_argument("--T", default="0,1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384")
    s.add_argument("--z", type=int, default=Z, choices=sorted(WALKERS), help="neighbours: 8 Moore, 4 von Neumann")
    s.add_argument("--out", default=None, help="pattern file (default patterns.txt / patterns5.txt)")
    s.add_argument("--load", help="analyse an existing pattern file instead of running")
    s.add_argument("--boot", type=int, default=0)
    s.set_defaults(fn=cmd_run)
    s = sub.add_parser("molofsky", help="five-cell rule: bistable sector in the (d1, d2) plane vs T")
    s.add_argument("--load", default=str(HERE / "patterns5.txt"))
    s.add_argument("--boot", type=int, default=0)
    s.add_argument("--csv", default=None)
    s.set_defaults(fn=cmd_molofsky)
    s = sub.add_parser("selftest")
    s.set_defaults(fn=cmd_selftest)
    a = ap.parse_args()
    if a.cmd == "run" and a.out is None:
        a.out = str(HERE / ("patterns.txt" if a.z == 8 else "patterns5.txt"))
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
