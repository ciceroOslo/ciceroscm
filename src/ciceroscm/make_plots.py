"""
Two functions to make output plots, saves the plots to a sub-folder in the output-folder
"""

import os

import matplotlib.pyplot as plt
import numpy as np


def fix_plot(axs, ylabel, title, legend=False):
    """
    Add labels to ohc, rib and temp plot

    Adding Year as xlabel, and a list of ylabels
    and list of titles to each of a list of axes
    If choosen add label
    Used typically  to annotate ohc rib and temp plot

    Parameters
    ----------
    axs : list
       list of plot axes
    ylabel : list
          list of ylabels for each axis
    title : list
          list of ylabels for each axis
    legend : bool
          Whether or not to add a legend to the
          plot axes
    """
    axs = axs.flatten()
    for i, ylab in enumerate(ylabel):
        axs[i].set_xlabel("Year")
        axs[i].set_ylabel(ylab)
        axs[i].set_title(title[i], fontsize=11, fontweight="bold")
        axs[i].set_title(f"{chr(i+97)})", fontsize=11, loc="left")
        if legend:
            axs[i].legend()


def plot_output1(pamset, results, nystart, nyend):
    """
    Plot ohc, rib and temp output

    Parameters
    ----------
    pamset : dict
          pamset possibly containing output_folder path
    results : dict
           Dictionary with results to plot
    nystart : int
           Startyear of dataset
    nyend : int
         Endyear of dataset
    """
    if "output_folder" in pamset:
        outdir = os.path.join(os.getcwd(), pamset["output_folder"])
    else:
        outdir = os.path.join(os.getcwd(), "output")
    if os.path.exists(outdir):
        plotdir = os.path.join(outdir, "plots")
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
    else:
        return
    #    if not os.path.exists("plots"):
    #        os.makedirs("plots")
    indices = np.arange(nystart, nyend + 1)
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, figsize=(12, 6))
    axs[0].plot(indices, results["OHC700"])
    axs[1].plot(indices, results["OHCTOT"])
    fix_plot(axs, ["OHC [$10^{22}$J]", "OHC [$10^{22}$J]"], ["OHC700", "OHCTOT"])
    fig.suptitle("CICERO SCM simulation, Ocean heat content")
    axs[1].yaxis.set_tick_params(labelbottom=True)
    plt.savefig(os.path.join(plotdir, "ohc.png"))
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(12, 6))
    axs[0].plot(indices, results["RIB_glob"], label="RIB_glob")
    axs[1].plot(indices, results["RIB_glob"], label="RIB_glob")
    axs[1].plot(indices, results["RIB_N"], label="RIB_N", linestyle=":")
    axs[1].plot(indices, results["RIB_S"], label="RIB_S", linestyle="--")
    fix_plot(
        axs,
        [r"W $\mathrm{m}^{-2}$", r"W $\mathrm{m}^{-2}$"],
        ["RIB_glob", "RIB_glob"],
        legend=True,
    )
    fig.suptitle("CICERO SCM simulation, RIB")
    plt.savefig(os.path.join(plotdir, "rib.png"))
    fig, axs = plt.subplots(nrows=1, ncols=2, sharex=True, figsize=(14, 6))

    for comp in ["dT_glob", "dT_glob_air", "dT_glob_sea"]:
        colours = {"glob": "tab:green", "air": "tab:orange", "sea": "tab:blue"}
        axs[0].plot(
            indices,
            results[comp],
            label=comp,
            color=colours[comp.rsplit("_", maxsplit=1)[-1]],
        )
    for comp in [
        "dT_NH",
        "dT_NH_air",
        "dT_NH_sea",
        "dT_SH",
        "dT_SH_air",
        "dT_SHsea",
    ]:
        if "NH" in comp:
            colours = {"NH": "tab:green", "air": "tab:orange", "sea": "tab:blue"}
            axs[1].plot(
                indices,
                results[comp],
                label=comp,
                linestyle=":",
                color=colours[comp.rsplit("_", maxsplit=1)[-1]],
            )
        else:
            colours = {"SH": "tab:green", "air": "tab:orange", "SHsea": "tab:blue"}
            axs[1].plot(
                indices,
                results[comp],
                label=comp,
                linestyle="--",
                color=colours[comp.rsplit("_", maxsplit=1)[-1]],
            )
    fix_plot(
        axs,
        ["K", "K"],
        ["Temperature change global", "Temperature change global and hemispheric"],
        legend=True,
    )
    plt.savefig(os.path.join(plotdir, "temp.png"))


def plot_output2(
    var, df_in, outdir, unit=None
):  # pylint: disable=too-many-locals, too-many-branches
    """
    Plot concentration, emission and forcing

    Parameters
    ----------
    var : str
       Denoting whether data is concentrations (conc),
       emissions (emis) or forcing (forc)
    df_in : pd.dataframe
         Dataframe to plot
    outdir : str
          path of where to send plot
    unit : pd.dataframe
        dataframe possibly containing units
    """
    plotdir = os.path.join(outdir, "plots")
    if var == "emis":
        df_in = df_in.drop(labels=unit[unit.values == "X"].index.tolist(), axis=1)
    elif var == "conc":
        df_in = df_in.drop(labels=unit[unit.values == "-"].index.tolist(), axis=1)
    elif var == "forc":
        df_in = df_in.drop(
            labels=["NOx", "CO", "NMVOC", "NH3", "BMB_AEROS_OC", "BMB_AEROS_BC"], axis=1
        )
    years = df_in["Year"]
    comps = df_in.columns[1:]
    for c, comp in enumerate(comps):
        if c in [0, 16, 32]:
            maxrows = 4
            if c == 32:
                maxrows = int(np.ceil((len(comps) - 32) / 4))
            fig, axs = plt.subplots(
                nrows=maxrows, ncols=4, sharex=True, figsize=(15, 10)
            )
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
            else:
                title = (
                    "CICERO SCM simulation, Concentrations "
                    + str(int(c / 16 + 1))
                    + " of 2"
                )
                fname = "conc_" + str(int(c / 16 + 1)) + ".png"
            fig.suptitle(title)
            i = 0
            j = 0
        axs[i, j].plot(years, df_in[comp])
        if comp != "Total_forcing":
            axs[i, j].set_title(comp, fontsize=11, fontweight="bold")
            axs[i, j].set_title(f"{chr(i*4+j+97)})", fontsize=11, loc="left")
        else:
            axs[i, j].set_title("Total anthropogenic", fontsize=11, fontweight="bold")
            axs[i, j].set_title(f"{chr(i*4+j+97)})", fontsize=11, loc="left")
        axs[i, j].xaxis.set_tick_params(labelbottom=True)
        if var != "forc":
            axs[i, j].set_ylim(ymin=0)
            axs[i, j].set_ylabel("[" + unit.loc[comp] + "]")
        if c in [15, 31, len(comps) - 1]:
            if c == len(comps) - 1:
                for j_left in range(j + 1, 4):
                    fig.delaxes(axs[i, j_left])
            fig.tight_layout(pad=4, h_pad=2, w_pad=1)
            fig.supxlabel("Year")
            fig.savefig(os.path.join(plotdir, fname))
        j += 1
        if j == 4:
            j = 0
            i += 1
