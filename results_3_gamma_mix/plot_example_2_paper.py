import argparse
import os

import matplotlib.pyplot as plt
import seaborn as sns

from common import load_tracers, compute_bulk


def main():
    parser = argparse.ArgumentParser(description="Paper-style radius/number plot for scalar-q2 gamma output.")
    parser.add_argument("file", nargs="?", default=None)
    parser.add_argument("--dir", default=".", help="Directory containing tracers.txt")
    parser.add_argument("--output", "-o", default=None)
    args = parser.parse_args()

    path = args.file or os.path.join(args.dir, "tracers.txt")
    tr = load_tracers(path)
    bulk = compute_bulk(tr)
    pl = tr["P_bar"]

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twiny()
    col = sns.color_palette("colorblind")

    p_rc = ax1.plot(bulk["r_c"], pl, c=col[0], label=r"$r_{\rm c}$")
    p_nc = ax2.plot(bulk["N_c"], pl, c=col[1], label=r"$N_{\rm c}$", ls="dashed")

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax2.set_xscale("log")
    ax1.invert_yaxis()

    yticks = [100, 10, 1, 0.1, 0.01, 1e-3]
    yticks_lab = ["100", "10", "1", "0.1", "0.01", "10$^{-3}$"]
    ax1.set_yticks(yticks, yticks_lab)
    ax1.set_ylim(300, 3e-3)

    ax1.tick_params(axis="both", which="major", labelsize=14)
    ax2.tick_params(axis="both", which="major", labelsize=14)

    ax1.set_xlabel(r"$r_{\rm c}$ [$\mu$m]", fontsize=16)
    ax2.set_xlabel(r"$N_{\rm c}$ [cm$^{-3}$]", fontsize=16)
    ax1.set_ylabel(r"$p_{\rm gas}$ [bar]", fontsize=16)

    lns = p_rc + p_nc
    labs = [line.get_label() for line in lns]
    ax2.legend(lns, labs, fontsize=10, loc="upper left")

    fig.tight_layout(pad=1.05)
    if args.output:
        fig.savefig(args.output, dpi=200, bbox_inches="tight")
    else:
        plt.show()


if __name__ == "__main__":
    main()
