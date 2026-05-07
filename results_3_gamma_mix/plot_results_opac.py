import argparse
import os

import matplotlib.pyplot as plt
import seaborn as sns

from common import load_opacity


def main():
    parser = argparse.ArgumentParser(description="Plot opacity tables.")
    parser.add_argument("--dir", default=".", help="Directory containing opac_k/a/g.txt")
    parser.add_argument("--prefix", default="opac", help="Opacity file prefix")
    args = parser.parse_args()

    cwd = os.getcwd()
    os.chdir(args.dir)
    try:
        opac = load_opacity(args.prefix)
    finally:
        os.chdir(cwd)

    wl = opac["wl"]
    pl = opac["pl"]
    nwl = len(wl)

    fig = plt.figure()
    col = sns.color_palette("husl", nwl)
    for i in range(nwl):
        plt.plot(opac["k_ext"][:, i], pl, c=col[i], label=f"{wl[i]:.2f}")
    plt.yscale("log")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$\kappa_{\rm ext}$")
    plt.ylabel(r"$p$ [bar]")
    plt.legend()

    fig = plt.figure()
    col = sns.color_palette("husl", nwl)
    for i in range(nwl):
        plt.plot(opac["ssa"][:, i], pl, c=col[i], label=f"{wl[i]:.2f}")
    plt.yscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel("single-scattering albedo")
    plt.ylabel(r"$p$ [bar]")
    plt.legend()

    fig = plt.figure()
    col = sns.color_palette("husl", nwl)
    for i in range(nwl):
        plt.plot(opac["g"][:, i], pl, c=col[i], label=f"{wl[i]:.2f}")
    plt.yscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(r"$g$")
    plt.ylabel(r"$p$ [bar]")
    plt.legend()

    plt.show()


if __name__ == "__main__":
    main()
