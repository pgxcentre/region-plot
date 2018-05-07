
# This file is part of region-plot.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import re
import logging
from os import path
from math import floor

import numpy as np
import pandas as pd

from geneparse import parsers, Extractor
from geneparse.utils import compute_ld as geneparse_ld

from .error import ProgramError


logger = logging.getLogger(__name__)


def compute_ld(best, fn, f_format, keep, extract):
    """Compute the LD."""
    # Retrieving the 'keep' list
    samples_to_keep = None
    if keep is not None:
        samples_to_keep = set(keep.read().splitlines())

    # The r2 values
    r2 = pd.Series()

    # Creating the parser
    GenoParser, geno_opts = _get_genotype_parser(fn, f_format)
    logger.debug("parser: {}".format(str(GenoParser)))
    logger.debug("geno_opts: {}".format(str(geno_opts)))
    with GenoParser(**geno_opts) as parser:
        # The mapping info for the sample to keep
        samples = np.array(parser.get_samples(), dtype=str)
        k = _get_sample_select(samples=samples, keep=samples_to_keep)
        logger.debug("{} samples".format(np.sum(k)))

        # Getting the best hit
        best_hit = parser.get_variant_by_name(best)
        if len(best_hit) != 1:
            raise ProgramError("There are {:,d} markers named {}"
                               "".format(len(best_hit), best))
        best_hit = best_hit.pop()
        best_hit.genotypes = best_hit.genotypes[k]

        # Estimating the number of markers to get at a time to not go over 2G
        # of RAM (approximately)
        logger.debug("best hit genotype size is {:,d}"
                     "".format(best_hit.genotypes.nbytes))
        nb_markers = floor(2 / (best_hit.genotypes.nbytes / 1073741824))
        logger.debug("  - Computing {:,d} markers at a time"
                     "".format(nb_markers))

        # Computing the LD per chunks of n markers
        for markers in get_n_markers(Extractor(parser, names=extract), k,
                                     nb_markers):
            r2 = r2.append(geneparse_ld(best_hit, markers, r2=True))
            logger.info("  - {:,d} markers processed".format(r2.shape[0]))

    return r2


def get_n_markers(parser, keep, n):
    """Retrieves n markers in the region using geneparse."""
    # The list of markers in the region
    markers_in_region = []

    # Parsing the genotypes
    for g in parser.iter_genotypes():
        # Changing the name of the variants so that it includes the alleles
        g.variant.name = "{}:{}/{}".format(
            g.variant.name, *sorted([g.reference, g.coded]),
        )

        # Getting the genotypes
        g.genotypes = g.genotypes[keep]
        markers_in_region.append(g)

        if len(markers_in_region) % 1000 == 0:
            logger.debug("Added 1000 markers")

        if len(markers_in_region) >= n:
            logger.info("Yielding")
            yield markers_in_region
            markers_in_region = []

    if len(markers_in_region) > 0:
        yield markers_in_region


def _get_sample_select(samples, keep):
    """Returns a vector of True/False to keep samples."""
    k = np.ones_like(samples, dtype=bool)
    if keep is not None:
        k = np.array([s in keep for s in samples], dtype=bool)
        if np.sum(k) == 0:
            raise ProgramError("No samples matched the keep list")
    return k


def _get_genotype_parser(fn, f_format):
    """Retrieves the correct genotype parsers."""
    # Finding the correct parser
    GenotypesParser = None
    if f_format is not None:
        GenotypesParser = parsers[f_format]
    else:
        # VCF
        if fn.endswith(".vcf") or fn.endswith(".vcf.gz"):
            GenotypesParser = parsers["vcf"]
            f_format = "vcf"

        # BGEN
        elif fn.endswith(".bgen"):
            GenotypesParser = parsers["bgen"]
            f_format = "bgen"

        # IMPUTE2
        elif fn.endswith(".impute2") or fn.endswith(".impute2.gz"):
            GenotypesParser = parsers["impute2"]
            f_format = "impute2"

        # Plink (gave the extension)
        elif fn.endswith(".bed") or fn.endswith(".bim") or fn.endswith(".fam"):
            GenotypesParser = parsers["plink"]
            fn = fn[:-4]
            f_format = "plink"

        # Plink (gave the prefix)
        elif (path.isfile(fn + ".bed") and path.isfile(fn + ".bim") and
              path.isfile(fn + ".fam")):
            GenotypesParser = parsers["plink"]
            f_format = "plink"

        else:
            raise ProgramError("Couldn't guess the genotypes file format.")

    # The parsers options
    parser_options = {}

    if f_format == "impute2":
        # Looking for the sample file
        possible_sample_fn = re.sub(r"\.impute2(\.gz)?$", ".sample", fn)
        if path.isfile(possible_sample_fn):
            parser_options["sample_filename"] = possible_sample_fn
        elif path.isfile(fn + ".sample"):
            parser_options["sample_filename"] = fn + ".sample"
        else:
            raise ProgramError("No sample file found")

    elif f_format == "bgen":
        # Looking for the sample file (might be absent, because it's in the
        # BGEN file itself)
        possible_sample_fn = re.sub(r"\.bgen", ".sample", fn)
        if path.isfile(possible_sample_fn):
            parser_options["sample_filename"] = possible_sample_fn
        elif path.isfile(fn + ".sample"):
            parser_options["sample_filename"] = fn + ".sample"
        parser_options["cpus"] = 2

    if f_format in {"impute2", "bgen", "vcf"}:
        parser_options["filename"] = fn

    elif f_format == "plink":
        parser_options["prefix"] = fn

    return GenotypesParser, parser_options
