#!/usr/bin/env python3
"""Front tracking for strip runs (dnsim --init strip:... --profile file).

Reads the time-series header for Lx, Ly and the binary profile records
(int64 t, Ly x uint32 alive per row, Ly x uint32 iface per row) and prints, per output
time, the active width A(t) = sum_y iface_y / (Lx * rho_if), where rho_if is the bulk
interface-cell density of the disordered phase (passed in, or estimated from the
initial random strip), plus the number of rows whose iface density exceeds half of it.
Front velocity per interface ~ dA/dt / 2 (the strip has two fronts on the torus).
"""
import sys
import numpy as np


def read(ts_path, prof_path):
    hdr = {}
    for line in open(ts_path):
        if not line.startswith("#"):
            break
        w = line[1:].split()
        if w and w[0] == "Lx":
            hdr = dict(Lx=int(w[1]), Ly=int(w[3]))
    Lx, Ly = hdr["Lx"], hdr["Ly"]
    raw = np.fromfile(prof_path, dtype=np.uint8)
    rec = 8 + 8 * Ly
    n = raw.size // rec
    raw = raw[: n * rec].reshape(n, rec)
    t = raw[:, :8].copy().view(np.int64).ravel()
    alive = raw[:, 8 : 8 + 4 * Ly].copy().view(np.uint32)
    iface = raw[:, 8 + 4 * Ly :].copy().view(np.uint32)
    return Lx, Ly, t, alive, iface


def main():
    ts_path, prof_path = sys.argv[1], sys.argv[2]
    rho_if = float(sys.argv[3]) if len(sys.argv) > 3 else None
    Lx, Ly, t, alive, iface = read(ts_path, prof_path)
    dens = iface / Lx
    if rho_if is None:
        rho_if = float(np.median(dens[0][dens[0] > 0]))
    A = dens.sum(1) / rho_if
    rows = (dens > 0.5 * rho_if).sum(1)
    print(f"# Lx {Lx} Ly {Ly} rho_if {rho_if:.4f}")
    print("# t A_rows active_rows mean_alive")
    for i in range(len(t)):
        print(f"{t[i]} {A[i]:.2f} {rows[i]} {alive[i].sum() / (Lx * Ly):.5f}")


if __name__ == "__main__":
    main()
