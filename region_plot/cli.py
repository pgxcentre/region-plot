"""The main script that do everything."""

# This file is part of region-plot.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


from __future__ import print_function

import os
import sys
import logging
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd

from six.moves import range, zip

from gepyto.utils.genes import ensembl_genes_in_region

from . import utils
from . import __version__
from .error import ProgramError

# importing matplotlib so that it works even if no X11
import matplotlib.pyplot as plt


def main():
    """The main function.

    This is what the pipeline should do:
        1- Find the best hit in an association result
        2- Compute LD with best hit
        3- Produce regional plot

    """
    # Creating the option parser
    desc = "Plots significant regions of GWAS.".format(__version__)
    parser = argparse.ArgumentParser(description=desc)

    # We run the script
    try:
        # Parsing the options
        args = parse_args(parser)
        check_args(args)

        # Adding the logging capability
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=args.log_level,
            handlers=[logging.StreamHandler(),
                      logging.FileHandler(args.log_file, mode="w")]
        )
        logging.info("Logging everything into '{}'".format(args.log_file))

        # We start by reading the association results
        assoc = read_assoc(args.assoc, args)

        # Getting the best hit
        best_hits = get_best_hits(assoc, args)

        # For all hits...
        for best_hit, chrom, start, end in zip(*best_hits):
            # Computing LD with best hit
            ld = compute_ld(args.bfile, assoc, best_hit, chrom, start, end,
                            args)

            # To read the genetic map, we required the chromosomal position of
            # the best hit
            genetic_map = read_genetic_map(chrom, start, end, args.genetic_map,
                                           args)

            # Reading the imputed sites
            imputed_sites = read_imputed_sites(args.imputed_sites)

            # Getting the gene from the region
            gene_list = find_gene_in_region(chrom, start, end, args.build)

            # Plotting the region
            plot_region(best_hit, assoc, ld, genetic_map, imputed_sites, chrom,
                        start, end, gene_list, args)

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    # Catching the ProgramError
    except ProgramError as e:
        parser.error(e.message)


def find_gene_in_region(chrom, start, end, build):
    """Finds the gene in the region."""
    region = "{}:{}-{}".format(chrom, start, end)

    logging.info("Fetching genes in region {}".format(region))
    results = ensembl_genes_in_region(region, bare=True, build=build)

    logging.info("  - {:,d} genes found".format(len(results)))

    # Creating a DataFrame from the list of genes
    genes = []
    for gene in results:
        genes.append((gene.symbol, gene.start, gene.end, gene.strand,
                      gene.biotype))
    genes = pd.DataFrame(genes, columns=["symbol", "start", "end", "strand",
                                         "biotype"])

    # Saving the gene list
    genes.to_csv("genes_in_chr{}_{}_{}.txt".format(chrom, start, end),
                 sep="\t", index=False)

    return genes


def plot_region(best, assoc, ld, genetic_map, imputed_sites, chrom, start, end,
                genes, options):
    """Plots the genomic region."""
    # Getting association data only from required region
    in_region = assoc[options.chr_col] == chrom
    in_region = in_region & (assoc[options.pos_col] >= start)
    in_region = in_region & (assoc[options.pos_col] <= end)
    data = assoc[in_region]

    # Merging association data and LD
    data = pd.merge(data, ld, how="left", left_on=options.snp_col,
                    right_on="SNP_B")

    # Is there NaN values?
    if data.R2.isnull().any():
        logging.debug("NaN values in R2 data")
        data.loc[data.R2.isnull(), "R2"] = 0
        if data.R2.isnull().any():
            raise ProgramError("problem with R2 data")

    # Just checking to be sure...
    if len(data[options.chr_col].unique()) > 1:
        raise ProgramError("more than one chromosome: problem with programmer")

    # Debug
    logging.info("Plotting the region")
    logging.info("  - chr{}:{}-{}".format(data[options.chr_col].unique()[0],
                                          data[options.pos_col].min(),
                                          data[options.pos_col].max()))

    # The colors
    r_colors = ("#0099CC", "#9933CC", "#669900", "#FF8800", "#CC0000")
    i_r_colors = ("#8AD5F0", "#D6ADEB", "#C5E26D", "#FFD980", "#FF9494")

    # the R2 threshold
    r_thresholds = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

    # The figure
    fig = plt.figure(figsize=(11, 4))

    # The axes
    recomb_axe = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    gene_axe = plt.subplot2grid((4, 1), (3, 0), sharex=recomb_axe)
    assoc_axe = recomb_axe.twinx()

    # Setting the recombination axe parameters
    recomb_axe.xaxis.set_ticks_position("none")
    recomb_axe.yaxis.set_ticks_position("right")
    recomb_axe.spines["top"].set_visible(False)
    recomb_axe.spines["left"].set_visible(False)
    recomb_axe.spines["right"].set_position(("outward", 9))
    recomb_axe.spines["bottom"].set_visible(False)
    recomb_axe.axes.get_xaxis().set_visible(False)
    recomb_axe.yaxis.set_label_position("right")

    # Setting the assoc axe parameters
    assoc_axe.xaxis.set_ticks_position("bottom")
    assoc_axe.yaxis.set_ticks_position("left")
    assoc_axe.spines["top"].set_visible(False)
    assoc_axe.spines["bottom"].set_visible(False)
    assoc_axe.spines["right"].set_visible(False)
    assoc_axe.spines["left"].set_position(("outward", 9))
    assoc_axe.yaxis.set_label_position("left")

    # Setting the gene axe parameters
    gene_axe.xaxis.set_ticks_position("bottom")
    gene_axe.yaxis.set_ticks_position("none")
    gene_axe.spines["top"].set_visible(False)
    gene_axe.spines["left"].set_visible(False)
    gene_axe.spines["bottom"].set_position(("outward", 9))
    gene_axe.spines["right"].set_visible(False)
    gene_axe.axes.get_yaxis().set_visible(False)

    # The size of the tick labels
    assoc_axe.tick_params(axis='both', which='major', labelsize=8)
    recomb_axe.tick_params(axis="both", which="major", labelsize=8)
    gene_axe.tick_params(axis="both", which="major", labelsize=8)

    # The title and labels of the figure
    gene_axe.set_xlabel("Position on chr{} (Mb)".format(chrom), fontsize=10,
                        weight="normal")
    assoc_axe.set_ylabel(r"$-\log_{10}(p)$", fontsize=10)
    recomb_axe.set_ylabel("Recombination Rate (cM/Mb)", fontsize=10,
                          weight="normal", rotation=270, va="bottom")

    # Plotting the recombination rate
    recomb_axe.plot(genetic_map[options.genetic_pos_col] / 1e6,
                    genetic_map[options.genetic_rate_col],
                    "-", lw=1, color="black", clip_on=False)

    # Plotting the imputed markers first
    is_imputed = data[options.snp_col].isin(imputed_sites)
    imputed = data[is_imputed]
    logging.info("  - {:,d} imputed markers".format(len(imputed)))
    logging.info("    - {:,d} significant markers".format(
        (imputed[options.p_col] < options.significant).sum(),
    ))
    for i in range(1, len(r_thresholds)):
        # Getting the thresholds
        min_r = r_thresholds[i-1]
        max_r = r_thresholds[i]
        color = i_r_colors[i-1]

        # Plotting the correct threshold
        sub_data = imputed[(min_r < imputed.R2) & (imputed.R2 <= max_r)]
        logging.info("    - {:,d} with r2 <= {}".format(len(sub_data), max_r))

        if len(sub_data) > 0:
            assoc_axe.plot(sub_data[options.pos_col] / 1e6,
                           -np.log10(sub_data[options.p_col]),
                           "D", mec=color, mfc=color, ms=3, clip_on=False,
                           label="_nolegend_")

    # Plotting the genotyped markers
    genotyped = data[~is_imputed]
    logging.info("  - {:,d} genotyped markers".format(len(genotyped)))
    logging.info("    - {:,d} significant markers".format(
        (genotyped[options.p_col] < options.significant).sum()
    ))
    for i in range(1, len(r_thresholds)):
        # Getting the thresholds
        min_r = r_thresholds[i-1]
        max_r = r_thresholds[i]
        color = r_colors[i-1]

        # Plotting the correct threshold
        sub_data = genotyped[(min_r < genotyped.R2) & (genotyped.R2 <= max_r)]
        logging.info("    - {:,d} with r2 <= {}".format(len(sub_data), max_r))

        if len(sub_data) > 0:
            assoc_axe.plot(sub_data[options.pos_col] / 1e6,
                           -np.log10(sub_data[options.p_col]),
                           ".", mec=color, mfc=color, ms=6, clip_on=False,
                           label="_nolegend_")

    # Adding the significant line
    assoc_axe.axhline(-np.log10(options.significant), ls="--",
                      color="#000000", lw=1)

    # Do we need a legend for the assoc axe?
    if len(imputed_sites) > 0:
        assoc_axe.plot([], [], "D", mec="#888888", mfc="#888888",
                       label="Imputed", ms=3)
        assoc_axe.plot([], [], ".", mec="#000000", mfc="#000000",
                       label="Genotyped", ms=6)

    # The R2 legend
    for i, r_color in enumerate(r_colors):
        # The imputed R2 color
        i_r_color = r_colors[i]
        if len(imputed_sites) > 0:
            i_r_color = i_r_colors[i]

        # The fake points
        assoc_axe.plot([], [], "s", mfc=i_r_color, mec=r_color, mew=2,
                       label=r"$r^2 \leq {}$".format(r_thresholds[i+1]))

    # the legend
    assoc_axe.legend(loc="best", fontsize=6, ncol=1, numpoints=1,
                     markerscale=1)

    # Sorting the genes
    genes = genes.sort_values(by=["start", "end"])

    # Getting the figure renderer
    renderer = fig.canvas.get_renderer()

    # The last bbox
    last_t_obj = {}
    last_end = defaultdict(int)

    for i in range(genes.shape[0]):
        gene_start = genes.iloc[i, :].start
        gene_end = genes.iloc[i, :].end
        gene_name = genes.iloc[i, :].symbol

        # Checking the starting position of the gene
        if gene_start < start:
            gene_start = start
        gene_start /= 1e6

        # Checking the ending position of the gene
        if gene_end > end:
            gene_end = end
        gene_end /= 1e6

        # Updating the gene label
        gene_label = None
        if genes.iloc[i, :].strand == 1:
            gene_label = gene_name + r"$\rightarrow$"
        else:
            gene_label = r"$\leftarrow$" + gene_name

        # We find the first j where we can put the line
        j = 0
        while True:
            if last_end[j] < gene_start:
                break
            j -= 1

        # Trying to put the label there
        t = gene_axe.text((gene_start + gene_end) / 2, j - 0.15, gene_label,
                          fontsize=5, ha="center", va="top")

        # Is there a bbox in this location?
        if j in last_t_obj:
            # Getting the bbox
            bb = t.get_window_extent(renderer=renderer)
            last_bb = last_t_obj[j].get_window_extent(renderer=renderer)

            while last_bb.overlaps(bb):
                # BBoxes overlap
                logging.debug("{} overlaps".format(gene_name))
                j -= 1
                t.set_y(j - 0.15)

                # Last j?
                if j not in last_t_obj:
                    break

                # Need to update both bboxes
                bb = t.get_window_extent(renderer=renderer)
                last_bb = last_t_obj[j].get_window_extent(renderer=renderer)

        # Plotting the line
        logging.debug("Putting {} at position {}".format(gene_name, j))
        marker = "-"
        other_param = {}
        if (gene_end - gene_start) < 3e-3:
            # Too small
            marker = "s"
            other_param["ms"] = 1.8
        gene_axe.plot([gene_start, gene_end], [j, j], marker, lw=2,
                      color="#000000", clip_on=False, **other_param)

        # Saving the last position (last end and bbox)
        last_end[j] = gene_end + 3e-3
        last_t_obj[j] = t

    # The limits
    recomb_axe.set_ylim(0, 100)
    recomb_axe.set_xlim(start/1e6, end/1e6)

    # The gene limits should at least be lower than -1
    min_y, max_y = gene_axe.get_ylim()
    if min_y == -1:
        gene_axe.set_ylim(-1.5, max_y)

    # Setting the ticks below the X axis for genes
    gene_axe.get_xaxis().set_tick_params(direction='out')

    # Saving the figure
    o_filename = "chr{}_{}-{}.{}".format(chrom, start, end,
                                         options.plot_format)
    logging.info("  - saving to '{}'".format(o_filename))
    plt.savefig(o_filename, dpi=600, bbox_inches='tight')
    plt.close(fig)


def read_assoc(filename, options):
    """Reads the association file."""
    # Reading the assoc file
    logging.info("Reading assoc file '{}'".format(filename))
    compression = None
    if filename.endswith(".gz"):
        compression = "gzip"
    data = pd.read_csv(filename, delim_whitespace=True,
                       compression=compression)

    # Checking the header
    for column in (options.snp_col, options.chr_col, options.pos_col,
                   options.p_col):
        if column not in data.columns:
            logging.debug(data.columns)
            raise ProgramError("{}: no column named {}".format(filename,
                                                               column))

    # Returning the association data
    logging.info("  - {:,d} markers".format(len(data)))
    data.index = data[options.snp_col]
    return data


def get_best_hits(assoc, args):
    """Gets the IDs of the best hits."""
    # The best hit
    best_hit = assoc[args.p_col].idxmin()

    # Gathering the position
    chrom, pos = assoc.loc[best_hit, [args.chr_col, args.pos_col]]
    start = int(pos - args.region_padding)
    start = start if start > 0 else 0
    end = int(pos + args.region_padding)

    # Do we want the whole region?
    logging.debug("region: chr{}:{}-{}".format(chrom, start, end))
    if args.whole_dataset:
        start = int(assoc[assoc[args.chr_col] == chrom].pos.min())
        end = int(assoc[assoc[args.chr_col] == chrom].pos.max())
        logging.debug("whole region: chr{}:{}-{}".format(chrom, start, end))

    logging.info("Best hit is '{}'".format(best_hit))
    logging.info("  - chr{}:{} (p={})".format(
        chrom,
        pos,
        assoc.loc[best_hit, args.p_col],
    ))

    # Saving in a list
    best_hits = [best_hit]
    chroms = [chrom]
    starts = [start]
    ends = [end]

    if args.whole_dataset:
        return best_hits, chroms, starts, ends

    # Finding the other bests hits
    sub_data = assoc.loc[assoc[args.p_col] < args.plot_p_lower]

    # First exclusions
    exclusion = (
        (sub_data[args.chr_col] == chroms[-1]) &
        (sub_data[args.pos_col] >= starts[-1]) &
        (sub_data[args.pos_col] <= ends[-1])
    )
    sub_data = sub_data.loc[~exclusion, :]

    # Until there are no other points
    while len(sub_data) > 0:
        # Finding the best hit
        best_hit = sub_data[args.p_col].idxmin()

        # Gathering the position
        chrom, pos = sub_data.loc[best_hit, [args.chr_col, args.pos_col]]
        start = int(pos - args.region_padding)
        end = int(pos + args.region_padding)

        # Logging
        logging.info("Secondary hit is '{}'".format(best_hit))
        logging.info("  - chr{}:{} (p={})".format(
            chrom,
            pos,
            sub_data.loc[best_hit, args.p_col],
        ))

        # Saving
        best_hits.append(best_hit)
        chroms.append(chrom)
        starts.append(start)
        ends.append(end)

        # Excluding the previous region
        exclusion = (
            (sub_data[args.chr_col] == chroms[-1]) &
            (sub_data[args.pos_col] >= starts[-1]) &
            (sub_data[args.pos_col] <= ends[-1])
        )
        sub_data = sub_data.loc[~exclusion, :]

    return best_hits, chroms, starts, ends


def compute_ld(prefix, assoc, best_hit, chrom, start, end, args):
    """Compute LD with the best SNP."""
    # Computing the LD for the best hit
    command = [
        "plink",
        "--noweb",
        "--bfile", prefix,
        "--r2",
        "--chr", str(chrom),
        "--from-bp", str(start),
        "--to-bp", str(end),
        "--ld-snp", best_hit,
        "--ld-window-kb", "1000000",
        "--ld-window", "1000000",
        "--ld-window-r2", "0",
        "--out", best_hit,
    ]

    # Launching the command
    if utils.execute_command("PLINK: LD with {}".format(best_hit), command):
        # The task was successful, so we read the LD table
        filename = "{}.ld".format(best_hit)
        if not os.path.isfile(filename):
            raise ProgramError("{}: file wasn't produce by "
                               "plink".format(filename))

        logging.info("Reading LD file '{}'".format(filename))
        ld_data = pd.read_csv(filename, delim_whitespace=True)[["SNP_B", "R2"]]
        logging.info("  - {:,d} LD values".format(len(ld_data)))

        return ld_data


def read_genetic_map(chrom, start, stop, filename, options):
    """Reads the genetic map."""
    compression = None
    if filename.endswith(".gz"):
        compression = "gzip"

    logging.info("Reading genetic map '{}'".format(filename))
    data = pd.read_csv(filename, sep="\t", compression=compression)

    # Checking the column
    for column in (options.genetic_chr_col, options.genetic_pos_col,
                   options.genetic_rate_col):
        if column not in data.columns:
            logging.debug(data.columns)
            raise ProgramError("{}: no column named {}".format(filename,
                                                               column))

    # Sub-setting the data to get a region of X base pair on each side of the
    # hit
    region = data["Chromosome"] == chrom
    region = region & (data["Position(bp)"] >= start)
    region = region & (data["Position(bp)"] <= stop)
    data = data[region]

    logging.info("  - {:,d} data points".format(len(data)))
    return data


def read_imputed_sites(filename):
    """Reads imputed sites (one site per line)."""
    if filename is None:
        logging.info("No imputed sites specified")
        return {}

    logging.info("Reading imputed sites '{}'".format(filename))

    data = None
    with open(filename, "r") as i_file:
        data = {i.rstrip("\r\n") for i in i_file.readlines()}

    logging.info("  - {:,d} imputed sites".format(len(data)))

    return data


def check_args(args):
    """Checks the arguments and options."""
    # Checking the input file
    # The association file
    if not os.path.isfile(args.assoc):
        raise ProgramError("{}: no such file".format(args.assoc))

    # The binary PEDFILE
    for suffix in (".bed", ".bim", ".fam"):
        if not os.path.isfile(args.bfile + suffix):
            raise ProgramError("{}: no such file".format(args.bfile + suffix))

    # Checking the genetic map
    if not os.path.isfile(args.genetic_map):
        raise ProgramError("{}: no such file".format(args.genetic_map))

    # Checking the padding value
    if args.region_padding < 0:
        raise ProgramError("{}: padding should be "
                           "positive".format(args.region_padding))
    if args.region_padding >= 2500000:
        raise ProgramError("{}: padding too large".format(args.region_padding))

    # Checking the imputed sites (if required)
    if args.imputed_sites is not None:
        if not os.path.isfile(args.imputed_sites):
            raise ProgramError("{}: no such file".format(args.imputed_sites))

    # Does the log file ends with '.log'?
    if not args.log_file.endswith(".log"):
        args.log_file += ".log"

    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "--log-level",
        type=str,
        choices=("INFO", "DEBUG"),
        default="INFO",
        help="The logging level [%(default)s]",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        metavar="LOGFILE",
        default="region-plot.log",
        help="The log file [%(default)s]",
    )

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument(
        "--assoc",
        type=str,
        metavar="FILE",
        required=True,
        help="The association file containing the hits",
    )
    group.add_argument(
        "--bfile",
        type=str,
        metavar="PREFIX",
        required=True,
        help="The prefix of the binary PEDFILE to compute LD with best hit",
    )
    group.add_argument(
        "--imputed-sites",
        type=str,
        metavar="FILE",
        help="The file containing the imputed sites (if absent, all points "
             "will have the same darkness)",
    )

    # The association option
    group = parser.add_argument_group("Association Options")
    group.add_argument(
        "--significant",
        type=float,
        metavar="FLOAT",
        default=5e-8,
        help="The significant association threshold [<%(default)e]",
    )
    group.add_argument(
        "--plot-p-lower",
        type=float,
        metavar="FLOAT",
        default=5e-8,
        help="Plot markers with p lower than value [<%(default)e]",
    )
    group.add_argument(
        "--snp-col",
        type=str,
        metavar="COL",
        default="snp",
        help="The name of the SNP column [%(default)s]",
    )
    group.add_argument(
        "--chr-col",
        type=str,
        metavar="COL",
        default="chr",
        help="The name of the chromosome column [%(default)s]",
    )
    group.add_argument(
        "--pos-col",
        type=str,
        metavar="COL",
        default="pos",
        help="The name of the pos column [%(default)s]",
    )
    group.add_argument(
        "--p-col",
        type=str,
        metavar="COL",
        default="p",
        help="The name of the p-value column [%(default)s]",
    )

    # The genetic map option
    group = parser.add_argument_group("Genetic Map Options")
    group.add_argument(
        "--genetic-map",
        type=str,
        metavar="FILE",
        required=True,
        help="The file containing the genetic map",
    )
    group.add_argument(
        "--genetic-chr-col",
        type=str,
        metavar="COL",
        default="chromosome",
        help="The name of chromosome column for the genetic map [%(default)s]",
    )
    group.add_argument(
        "--genetic-pos-col",
        type=str,
        metavar="COL",
        default="position",
        help="The name of the position column for the genetic map "
             "[%(default)s]",
    )
    group.add_argument(
        "--genetic-rate-col",
        type=str,
        metavar="COL",
        default="rate",
        help="The name of the recombination rate column for the genetic map "
             "[%(default)s]",
    )

    # The plot options
    group = parser.add_argument_group("Plot Options")
    group.add_argument(
        "--plot-format",
        type=str,
        choices={"png", "pdf"},
        default="png",
        help="The format of the output file containing the plot (might be "
             "'png' or 'pdf') [%(default)s]",
    )
    group.add_argument(
        "--build",
        type=str,
        choices=("GRCh37", "GRCh38"),
        default="GRCh37",
        help="The build to search the overlapping genes [%(default)s]",
    )
    group.add_argument(
        "--region-padding",
        type=float,
        metavar="FLOAT",
        default=500e3,
        help="The amount of base pairs to pad the region (on each side of the "
             "best hit [%(default).1f]",
    )
    group.add_argument(
        "--whole-dataset",
        action="store_true",
        help="Plot all markers (no padding) (WARNING this might take a lot of "
             "memory)",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
