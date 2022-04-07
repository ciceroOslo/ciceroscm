"""
Two functions to make output plots, saves the plots to a sub-folder in the output-folder
"""
import os

import matplotlib.pyplot as plt
import numpy as np


def fix_plot(axs, ylabel, title, legend=False):
    """
    Add labels to ohc, rib and temp plot
    """
    axs = axs.flatten()
    for i, ylab in enumerate(ylabel):
        axs[i].set_xlabel("Year")
        axs[i].set_ylabel(ylab)
        axs[i].set_title(title[i])
        if legend:
            axs[i].legend()


def plot_output1(pamset, results, nystart, nyend):
    """
    Plot ohc, rib and temp output
    """
    if "output_folder" in pamset:
        outdir = os.path.join(os.getcwd(), pamset["output_folder"])
    else:
        outdir = os.path.join(os.getcwd(), "output")
    plotdir = outdir + "/plots/"
    if not os.path.isdir(plotdir):
        os.system("mkdir " + plotdir)
    indices = np.arange(nystart, nyend + 1)
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(12, 6))
    axs[0].plot(indices, results["OHC700"])
    axs[1].plot(indices, results["OHCTOT"])
    fix_plot(axs, ["OHC [$10^{22}$J]", "OHC [$10^{22}$J]"], ["OHC700", "OHCTOT"])
    fig.suptitle("CICERO SCM simulation, Ocean heat content")
    plt.savefig(plotdir + "ohc.png")
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(12, 6))
    axs[0].plot(indices, results["RIB_glob"], label="RIB_glob")
    axs[1].plot(indices, results["RIB_glob"], label="RIB_glob")
    axs[1].plot(indices, results["RIB_N"], label="RIB_N", linestyle=":")
    axs[1].plot(indices, results["RIB_S"], label="RIB_S", linestyle="--")
    fix_plot(
        axs, ["check", "check"], ["RIB_glob", "RIB_glob"], legend=True,
    )
    fig.suptitle("CICERO SCM simulation, RIB")
    plt.savefig(plotdir + "rib.png")
    fig, axs = plt.subplots(nrows=1, ncols=3, sharex=True, figsize=(14, 6))
    for comp in ["dT_glob", "dT_glob_air", "dT_glob_sea"]:
        axs[0].plot(indices, results[comp], label=comp)
    for comp in [
        "dT_glob",
        "dT_glob_air",
        "dT_glob_sea",
        "dT_NH",
        "dT_NH_air",
        "dT_NH_sea",
        "dT_SH",
        "dT_SH_air",
        "dT_SHsea",
    ]:
        if "NH" in comp:
            axs[1].plot(indices, results[comp], label=comp, linestyle=":")
        elif "SH" in comp:
            axs[1].plot(indices, results[comp], label=comp, linestyle="--")
        else:
            axs[1].plot(indices, results[comp], label=comp)
    for comp in ["dSL(m)", "dSL_thermal(m)", "dSL_ice(m)"]:
        axs[2].plot(indices, results[comp], label=comp)
    fix_plot(
        axs,
        ["Temp. [$^\\circ\\!$C]", "Temp. [$^\\circ\\!$C]", "Sea level rise [m]", ],
        [
            "CICERO SCM simulation, Temperature",
            "CICERO SCM simulation, Temperature",
            "CICERO SCM simulation, SLR \n(!!!!!need to be revised-checked!!!!)",
        ],
        legend=True,
    )
    plt.savefig(plotdir + "temp.png")


def plot_output2(var, df_in, outdir=None):
    """
    Plot concentration, emission and forcing
    """
    if outdir:
        plotdir = outdir + "/plots/"
        if not os.path.isdir(plotdir):
            os.system("mkdir " + plotdir)
    years = df_in["Year"]
    comps = df_in.columns[1:]
    for c, comp in enumerate(comps):
        if c in [0, 16, 32]:
            fig, axs = plt.subplots(nrows=4, ncols=4, sharex=True, figsize=(15, 10))
            if var == "forc":
                title = (
                    "CICERO SCM simulation, Radiative Forcing "
                    + str(int(c / 16 + 1))
                    + " of 3"
                )
                fname = "forc_" + str(int(c / 16 + 1)) + ".png"
                fig.supylabel("RF [Wm$^{-2}$]")
            elif var == "emis":
                title = (
                    "CICERO SCM simulation, Emissions " + str(int(c / 16 + 1)) + " of 3"
                )
                fname = "em_" + str(int(c / 16 + 1)) + ".png"
                fig.supylabel("Emissions [Gg?]")
            else:
                title = (
                    "CICERO SCM simulation, Concentrations "
                    + str(int(c / 16 + 1))
                    + " of 3"
                )
                fname = "conc_" + str(int(c / 16 + 1)) + ".png"
                fig.supylabel("Concentration [?]")
            fig.suptitle(title)
            i = 0
            j = 0
        axs[i, j].plot(years, df_in[comp])
        axs[i, j].set_title(comp)
        if var != "forc":
            axs[i, j].set_ylim(ymin=0)
        if c in [15, 31, len(comps) - 1]:
            if c == len(comps) - 1:
                fig.delaxes(axs[3, 3])
                axs[2, 3].xaxis.set_tick_params(labelbottom=True)
                if var != "forc":
                    fig.delaxes(axs[3, 2])
                    axs[2, 2].xaxis.set_tick_params(labelbottom=True)
            fig.tight_layout(pad=4, h_pad=2, w_pad=1)
            fig.supxlabel("Year")
            fig.savefig(plotdir + fname)
        j += 1
        if j == 4:
            j = 0
            i += 1