# region-plot - A tool to plot significant regions of GWAS

`region-plot` is a bioinformatics pipeline to plot significant regions found by
genome-wide association studies (GWAS). It works on both Python 2 and 3
versions.

The following figure (low resolution) is an example from Tardif *et al.* (2015)
(doi:[10.1161/CIRCGENETICS.114.000663](http://dx.doi.org/10.1161/CIRCGENETICS.114.000663)).
The difference with the original figure is the inclusion of the annotation from
the HAVANA project in this example.

<img src=https://raw.github.com/pgxcentre/region-plot/master/example.png width=728 />

## Dependencies

The tool requires a standard [Python](http://python.org/) installation with the
following packages:

1. [numpy](http://www.numpy.org/) version 1.9.1 or latest
2. [pandas](http://pandas.pydata.org/) version 0.17.0 or latest
3. [six](http://pythonhosted.org/six/) version 1.9.0 or latest
4. [matplotlib](http://matplotlib.org/) version 1.4.3 or latest
5. [gepyto](https://github.com/legaultmarc/gepyto) version 0.9.2 or latest

The tool has been tested on *Linux* only, but should also work on both Mac OSX
and Windows.

## Usage

For Linux users, make sure that the script is executable (using the *chmod*
command).

```console
$ launch-region-plot --help
usage: launch-region-plot [-h] [-v] [--log-level {INFO,DEBUG}] [--log-file LOGFILE] --assoc FILE --genotypes FILE
                          [--imputed-sites FILE] [--annotation-gtf FILE] [--annotation-label LABEL [LABEL ...]]
                          [--genotypes-format FORMAT] [--keep FILE] [--significant FLOAT] [--plot-p-lower FLOAT] [--snp-col COL]
                          [--chr-col COL] [--pos-col COL] [--p-col COL] [--a1-col ALLELE] [--a2-col ALLELE] --genetic-map FILE
                          [--genetic-chr-col COL] [--genetic-pos-col COL] [--genetic-rate-col COL] [--plot-format {png,pdf}]
                          [--build {GRCh37,GRCh38}] [--region-padding FLOAT] [--whole-dataset] [--output-directory DIR]

Plots significant regions of GWAS (0.1.4).

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --log-level {INFO,DEBUG}
                        The logging level. [INFO]
  --log-file LOGFILE    The log file. [region-plot.log]

Input Files:
  --assoc FILE          The association file containing the hits.
  --genotypes FILE      The file containing the genotypes (available format are VCF, IMPUTE2, BGEN or Plink binary files.
  --imputed-sites FILE  The file containing the imputed sites (if absent, all points will have the same darkness).
  --annotation-gtf FILE
                        A GTF file containing annotations.

Annotation Options:
  --annotation-label LABEL [LABEL ...]
                        Labels from the GTF file attributes that will be used as a label in order of preference.

Genotypes Options:
  --genotypes-format FORMAT
                        The genotype file format. If not specified, the tool will try to guess the format and parse the file
                        accordingly.
  --keep FILE           The list of samples to keep for the LD calculation.

Association Options:
  --significant FLOAT   The significant association threshold. [<5.000000e-08]
  --plot-p-lower FLOAT  Plot markers with p lower than value. [<5.000000e-08]
  --snp-col COL         The name of the SNP column. [snp]
  --chr-col COL         The name of the chromosome column. [chr]
  --pos-col COL         The name of the pos column. [pos]
  --p-col COL           The name of the p-value column. [p]
  --a1-col ALLELE       The name of the column containing the first allele. [minor]
  --a2-col ALLELE       The name of the column containing the second allele. [major]

Genetic Map Options:
  --genetic-map FILE    The file containing the genetic map.
  --genetic-chr-col COL
                        The name of chromosome column for the genetic map. [chromosome]
  --genetic-pos-col COL
                        The name of the position column for the genetic map. [position]
  --genetic-rate-col COL
                        The name of the recombination rate column for the genetic map. [rate]

Plot Options:
  --plot-format {png,pdf}
                        The format of the output file containing the plot (might be 'png' or 'pdf'). [png]
  --build {GRCh37,GRCh38}
                        The build to search the overlapping genes. [GRCh37]
  --region-padding FLOAT
                        The amount of base pairs to pad the region (on each side of the best hit. [500000.0]
  --whole-dataset       Plot all markers (no padding) (WARNING this might take a lot of memory).

Output Options:
  --output-directory DIR
                        The output directory. [.]
```
