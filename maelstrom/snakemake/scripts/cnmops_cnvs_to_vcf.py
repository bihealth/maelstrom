"""Convert cn.mops output TSV file to a VCF file."""

import csv
import gzip
import os
import pathlib
import sys
import subprocess
import types
import typing

from logzero import logger
import vcfpy

CONTIGS = (
    {"ID": "1", "length": 249250621},
    {"ID": "2", "length": 243199373},
    {"ID": "3", "length": 198022430},
    {"ID": "4", "length": 191154276},
    {"ID": "5", "length": 180915260},
    {"ID": "6", "length": 171115067},
    {"ID": "7", "length": 159138663},
    {"ID": "8", "length": 146364022},
    {"ID": "9", "length": 141213431},
    {"ID": "10", "length": 135534747},
    {"ID": "11", "length": 135006516},
    {"ID": "12", "length": 133851895},
    {"ID": "13", "length": 115169878},
    {"ID": "14", "length": 107349540},
    {"ID": "15", "length": 102531392},
    {"ID": "16", "length": 90354753},
    {"ID": "17", "length": 81195210},
    {"ID": "18", "length": 78077248},
    {"ID": "19", "length": 59128983},
    {"ID": "20", "length": 63025520},
    {"ID": "21", "length": 48129895},
    {"ID": "22", "length": 51304566},
    {"ID": "X", "length": 155270560},
    {"ID": "Y", "length": 59373566},
)

INFOS = (
    {"ID": "END", "Number": 1, "Type": "Integer", "Description": "End coordinate of this variant"},
    {
        "ID": "SVLEN",
        "Number": 1,
        "Type": "Integer",
        "Description": "Difference in length between REF and ALT alleles",
    },
    {"ID": "SVTYPE", "Number": 1, "Type": "String", "Description": "Type of structural variant"},
    {
        "ID": "SVMETHOD",
        "Number": 1,
        "Type": "String",
        "Description": "Type of approach used to detect SV",
    },
)

ALTS = ({"ID": "CNV", "Description": "Copy number variation"},)

FILTERS = ({"ID": "PASS", "Description": "All filters passed"},)

FORMATS = (
    {"ID": "GT", "Number": 1, "Type": "String", "Description": "All filters passed"},
    {"ID": "CN", "Number": 1, "Type": "Integer", "Description": "Inferred integer copy number"},
    {
        "ID": "cnmops",
        "Number": 1,
        "Type": "Integer",
        "Description": "Whether or not CNV was called by cn.mops",
    },
)


class ChangeDir:
    """Helper context manager for changing the current working directory."""

    def __init__(self, new_path: pathlib.Path):
        self.old_path = None
        self.new_path = new_path

    def __enter__(self):
        logger.info("cd %s", self.new_path)
        self.old_path = pathlib.Path.cwd()
        os.chdir(self.new_path)

    def __exit__(self, _etype, _value, _traceback):
        logger.info("cd %s", self.old_path)
        os.chdir(self.old_path)


def build_grch37_header(samples: typing.Iterable[str]) -> vcfpy.Header:
    """Build ``vcfpy.Header`` for GRCh37 for the given samples."""
    result = vcfpy.Header()
    result.add_line(vcfpy.HeaderLine("fileformat", "VCFv4.2"))
    result.samples = vcfpy.SamplesInfos(list(samples))
    for contig in CONTIGS:
        result.add_contig_line(contig)
    for filter_ in FILTERS:
        result.add_filter_line(filter_)
    for info in INFOS:
        result.add_info_line(info)
    for alt in ALTS:
        result.add_line(vcfpy.AltAlleleHeaderLine.from_mapping(alt))
    for format_ in FORMATS:
        result.add_format_line(format_)
    return result


def main():
    logger.info("performing conversion")
    sample = snakemake.wildcards.sample
    if "WGS" not in sample:
        sample = sample + "-N1-DNA1-WGS1"
    logger.info("processing for sample %s", sample)
    logger.info("constructing header")
    header = build_grch37_header([sample])
    logger.info("opening input file %s", snakemake.input)
    if str(snakemake.input).endswith(".gz"):
        tsvf = gzip.open(str(snakemake.input), "rt")
    else:
        tsvf = open(str(snakemake.input), "rt")
    with tsvf:
        logger.info("opening output file %s", snakemake.output.vcf)
        with vcfpy.Writer.from_path(str(snakemake.output.vcf), header) as writer:
            tsv_header = None
            contig = None
            for row in csv.reader(tsvf, dialect="excel-tab"):
                if tsv_header is None:
                    tsv_header = row
                    continue

                vals = dict(zip(tsv_header, row[1:]))
                if vals["sampleName"] != sample:
                    continue

                if contig != vals["seqnames"]:
                    logger.info("Starting with contig %s", vals["seqnames"])
                    contig = vals["seqnames"]

                # TODO: properly handle sex chromosome!
                cn = int(vals["CN"][2:])
                if vals["CN"] in ("CN0", "CN1"):
                    sv_type = "DEL"
                elif vals["CN"] == "CN2":
                    sv_type = "."
                else:
                    sv_type = "DUP"

                record = vcfpy.Record(
                    CHROM=vals["seqnames"],
                    POS=int(vals["start"]),
                    ID=[],
                    REF="N",
                    ALT=[vcfpy.SymbolicAllele("CNV")],
                    QUAL=None,
                    FILTER=["PASS"],
                    INFO={
                        "END": int(vals["end"]),
                        "SVLEN": abs(int(vals["end"]) - int(vals["start"]) + 1),
                        "SVTYPE": sv_type,
                        "SVMETHOD": "cn.mops",
                    },
                    FORMAT=["GT", "CN", "cnmops"],
                    calls=[
                        vcfpy.Call(
                            sample=sample,
                            data={
                                "GT": "1" if sv_type in ("DEL", "DUP") else "0",
                                "CN": cn,
                                "cnmops": 1 if sv_type in ("DEL", "DUP") else 0,
                            },
                        )
                    ],
                )
                writer.write_record(record)

    logger.info("computing tabix file")
    subprocess.check_call(["tabix", "-f", snakemake.output.vcf])
    logger.info("computing MD5 files")
    path_out = pathlib.Path(str(snakemake.output.vcf))
    with ChangeDir(path_out.parent):
        with open("%s.md5" % path_out.name, "wt") as md5f:
            subprocess.check_call(["md5sum", path_out.name], stdout=md5f)
        with open("%s.tbi.md5" % path_out.name, "wt") as md5f:
            subprocess.check_call(["md5sum", "%s.tbi" % path_out.name], stdout=md5f)

    logger.info("All done. Have a nice day!")


if __name__ == "__main__":
    sys.exit(main())
