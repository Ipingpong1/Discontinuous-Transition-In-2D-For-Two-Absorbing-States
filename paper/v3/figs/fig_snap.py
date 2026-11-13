"""Fig. snap: configurations at p = 1/2, alpha = 1/2 (local dnsim runs, data/snap/*.ckpt).
Sites with no disagreeing Moore neighbor are white (dead) or gray (alive); sites with at least one
are orange (dead) or red (alive), so the active phase stands out from the two absorbing phases."""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from common import PAGE, ROOT, save

SNAP = os.path.join(ROOT, "data", "snap")


def load(name):
    raw = open(os.path.join(SNAP, name), "rb").read()
    Lx, Ly, t = np.frombuffer(raw[8:32], dtype=np.uint64).astype(int)
    words = np.frombuffer(raw[32:], dtype=np.uint64).reshape(Ly, Lx // 64)
    g = np.unpackbits(words.view(np.uint8).reshape(Ly, Lx // 64, 8), axis=2, bitorder="little")
    return g.reshape(Ly, Lx).astype(np.uint8), t


def classes(g):
    j = np.zeros(g.shape, int)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dx or dy:
                j += np.roll(np.roll(g, dy, 0), dx, 1) != g
    return g + 2 * (j > 0)


cmap = ListedColormap(["#ffffff", "#5a5a5a", "#f6b47c", "#c0262d"])
panels = [("a.ckpt", None, r"(a) $q=0.048195$"),
          ("b7.ckpt", (slice(0, 256), slice(0, 256)), r"(b) $q=0.043605$"),
          ("c.ckpt", (slice(128, 384), slice(0, 256)), r"(c) $q=0.044982$"),
          ("d.ckpt", (slice(128, 384), slice(128, 384)), r"(d) $q=0.05049$")]
fig, axs = plt.subplots(1, 4, figsize=(PAGE, PAGE / 4 + 0.25), gridspec_kw=dict(wspace=0.05))
for ax, (f, crop, title) in zip(axs, panels):
    g, t = load(f)
    c = classes(g)
    if crop:
        c = c[crop]
    ax.imshow(c, cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest")
    ax.set_xticks([]), ax.set_yticks([])
    ax.set_title(title + f", $t={t}$", fontsize=7)
save(fig, "fig_snap")
