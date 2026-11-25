import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

ROOT    = "results"                       # dossier contenant center/ top_left/ top_right/
SOURCES = ["center", "top_left", "top_right"]
KS      = [50, 80, 100, 150]
LEVELS  = [0, 1, 2, 3, 4]
OUT     = "figures"

# ---------- maillage FreeFEM (savemesh) ----------
def read_mesh(path):
    with open(path) as f:
        nv, nt, _ = map(int, f.readline().split()[:3])
        pts = np.array([f.readline().split()[:2] for _ in range(nv)], float)
        tri = np.array([f.readline().split()[:3] for _ in range(nt)], int) - 1
    return pts, tri

# ---------- solution : header "x y re ..." puis x y re ----------
def read_sol(path):
    d = pd.read_csv(path, sep=r"\s+", skiprows=1, header=None).values
    return d[:, 0], d[:, 1], d[:, 2]          # x, y, Re(u)

# ---------- lecture des erreurs ----------
COLS = [
    "level",
    "h",
    "hmax",
    "nt",
    "ndofFull",
    "ndofPhys",
    "errL2",
    "errH1k",
    "relL2",
    "relH1k"
]
def read_err(path):
    df = pd.read_csv(path, sep=r"\s+", comment="#", header=None)
    if not str(df.iloc[0, 0]).replace(".", "", 1).replace("e", "", 1).replace("-", "").isdigit():
        df = df.iloc[1:].astype(float)      # il y avait un en-tête
    df.columns = COLS[:df.shape[1]]
    return df.astype(float).reset_index(drop=True)

def mesh_path(src, k, lev):  return f"{ROOT}/{src}/mesh_for_k{k}lev{lev}.msh"
def sol_path(src, k, lev):   return f"{ROOT}/{src}/sol_appro_k{k}lev{lev}.dat"

for src in SOURCES:
    os.makedirs(f"{OUT}/{src}", exist_ok=True)

    # ===== 1) erreurs vs ndofPhys, une courbe par k =====
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for k in KS:
        e = read_err(f"{ROOT}/{src}/error_k{k}.dat")
        ax[0].loglog(e.ndofPhys, e.relL2,  "o-", label=f"$k={k}$")
        ax[1].loglog(e.ndofPhys, e.relH1k, "s-", label=f"$k={k}$")
    ax[0].set_ylabel(r"Rel Error $L^2$")
    ax[1].set_ylabel(r"eRel Error $H^1_k$")
    for a in ax:
        a.set_xlabel("#dofs "); a.grid(True, which="both", alpha=.3); a.legend()
    plt.tight_layout(); plt.savefig(f"{OUT}/{src}/err.pdf"); plt.close()

    for k in KS:
        # ===== 2) maillages niveau 0 et 4 =====
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.8))
        for a, lev in zip(ax, (0, 4)):
            p, t = read_mesh(mesh_path(src, k, lev))
            a.triplot(p[:, 0], p[:, 1], t, lw=.25, color="k")
            a.set_aspect("equal"); a.set_title(f"$k={k}$, Level {lev} ({len(t)} Triangles)")
        plt.tight_layout(); plt.savefig(f"{OUT}/{src}/mesh_k{k}.pdf"); plt.close()

        # ===== 3) solutions : référence, niveau 0, niveau 4 =====
        files = [("Reference",          f"{ROOT}/{src}/sol_ref_k{k}.dat"),
                 ("Approx, Lev. 0", sol_path(src, k, 0)),
                 ("Approx, Lev 4", sol_path(src, k, 4))]
        data = [(n, *read_sol(f)) for n, f in files]
        vmax = np.abs(data[0][3]).max()
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
        for a, (name, x, y, u) in zip(ax, data):
            tr = mtri.Triangulation(x, y)      # Delaunay sur les dofs (domaine convexe)
            c = a.tricontourf(tr, u, levels=60, cmap="RdBu_r",
                              vmin=-vmax, vmax=vmax, extend="both")
            a.set_aspect("equal"); a.set_title(rf"Re$(u)$ — {name}")
        fig.colorbar(c, ax=ax, shrink=.8)
        plt.savefig(f"{OUT}/{src}/sol_k{k}.pdf", bbox_inches="tight"); plt.close()
print("OK")
