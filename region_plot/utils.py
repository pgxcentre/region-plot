"""Utility functions."""

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


class GTFFile(object):
    """Reads a GTF File and allows region based queries.

    Example:

    gtf = GTFFile("Homo_sapiens.GRCh37.87.genes_only.gtf.gz")
    gtf.get_annotations_in_region("16", 56995762, 56995762+5,
                                  ["gene_name", "gene_id"])

    """
    def __init__(self, filename):
        self.df = pd.read_csv(
            filename,
            comment="#",
            names=["seqname", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute"],
            sep="\t",
            dtype={"seqname": str, "start": int, "end": int}
        )

    def get_features_in_region(self, chrom, start, end):
        """Does raw indexing in the underlying GTF file to get annotations."""
        chrom = str(chrom)

        df = self.df

        matching_annotations = df.loc[
            (df.seqname == chrom) &
            (((df.start >= start) & (df.start <= end) |
             ((df.end >= start) & (df.end <= end)))), :
        ]

        return matching_annotations

    def get_annotations_in_region(self, chrom, start, end, preferred_labels):
        """Gets the annotation in the format expected by region-plot."""
        features = self.get_features_in_region(chrom, start, end)

        # Parse the features to return a DF with:
        # start, end, strand, label
        annotations = []
        for _, row in features.iterrows():
            # Parse GTF attributes.
            attributes = {}
            for attr in row.attribute.strip().split(";"):
                if not attr:
                    continue

                splitted = attr.strip().split()
                key = splitted[0]
                value = " ".join(splitted[1:]).strip('"')
                attributes[key] = value

            # Try to find a label in the provided preferred labels.
            label = None
            for try_label in preferred_labels:
                if try_label in attributes:
                    label = attributes[try_label]
                    break

            if label is None:
                label = "Annotation"

            annotations.append((row.start, row.end, row.strand, label))

        return pd.DataFrame(
            annotations,
            columns=["start", "end", "strand", "label"]
        )


def compute_ld(best, fn, f_format, samples_to_keep, extract):
    """Compute the LD."""
    # The r2 values
    r_squared = pd.Series(dtype=float)

    # Creating the parser
    geno_parser, geno_opts = _get_genotype_parser(fn, f_format)
    logger.debug("parser: %s", str(geno_parser))
    logger.debug("geno_opts: %s", str(geno_opts))
    with geno_parser(**geno_opts) as parser:
        # The mapping info for the sample to keep
        samples = np.array(parser.get_samples(), dtype=str)
        keep = _get_sample_select(samples=samples, keep=samples_to_keep)
        logger.debug("%d samples", np.sum(keep))

        # Getting the best hit
        best_hit = parser.get_variant_by_name(best)
        if len(best_hit) != 1:
            raise ProgramError("There are {:,d} markers named {}"
                               "".format(len(best_hit), best))
        best_hit = best_hit.pop()
        best_hit.genotypes = best_hit.genotypes[keep]

        # Estimating the number of markers to get at a time to not go over 2G
        # of RAM (approximately)
        logger.debug("best hit genotype size is {:,d}"
                     "".format(best_hit.genotypes.nbytes))
        nb_markers = floor(2 / (best_hit.genotypes.nbytes / 1073741824))
        logger.debug("  - Computing {:,d} markers at a time"
                     "".format(nb_markers))

        # Computing the LD per chunks of n markers
        for markers in get_n_markers(Extractor(parser, names=extract), keep,
                                     nb_markers):
            r_squared = r_squared.append(
                geneparse_ld(best_hit, markers, r2=True),
            )
            logger.info(
                "  - {:,d} markers processed".format(r_squared.shape[0]),
            )

    return r_squared


def get_n_markers(parser, keep, n):
    """Retrieves n markers in the region using geneparse."""
    # The list of markers in the region
    markers_in_region = []

    # Parsing the genotypes
    for data in parser.iter_genotypes():
        # Changing the name of the variants so that it includes the alleles
        data.variant.name = "{}:{}/{}".format(
            data.variant.name, *sorted([data.reference, data.coded]),
        )

        # Getting the genotypes
        data.genotypes = data.genotypes[keep]
        markers_in_region.append(data)

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
    to_keep = np.ones_like(samples, dtype=bool)
    if keep is not None:
        to_keep = np.array([s in keep for s in samples], dtype=bool)
        if np.sum(to_keep) == 0:
            raise ProgramError("No samples matched the keep list")
    return to_keep


def _get_genotype_parser(fn, f_format):
    """Retrieves the correct genotype parsers."""
    # Finding the correct parser
    genotypes_parser = None

    if f_format is not None:
        genotypes_parser = parsers[f_format]

    else:
        # VCF
        if fn.endswith(".vcf") or fn.endswith(".vcf.gz"):
            genotypes_parser = parsers["vcf"]
            f_format = "vcf"

        # BGEN
        elif fn.endswith(".bgen"):
            genotypes_parser = parsers["bgen"]
            f_format = "bgen"

        # IMPUTE2
        elif fn.endswith(".impute2") or fn.endswith(".impute2.gz"):
            genotypes_parser = parsers["impute2"]
            f_format = "impute2"

        # Plink (gave the extension)
        elif fn.endswith(".bed") or fn.endswith(".bim") or fn.endswith(".fam"):
            genotypes_parser = parsers["plink"]
            fn = fn[:-4]
            f_format = "plink"

        # Plink (gave the prefix)
        elif (path.isfile(fn + ".bed") and path.isfile(fn + ".bim") and
              path.isfile(fn + ".fam")):
            genotypes_parser = parsers["plink"]
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

    return genotypes_parser, parser_options
