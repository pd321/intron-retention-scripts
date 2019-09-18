# ---------------------- get msi  -------------------------------------------
# Get MSI(Mis-splicing index) values for a given bam file to
# introns given in bed format
# ----------------------------------------------------------------------------
import argparse
import logging
import os

import pybedtools


class Intron(object):
    """ Objects of the type intron """

    def __init__(self, id):
        """
        Parameters
        ----------
        bed_file_path: str, filepath
            Path to a .bed file
        """
        self.id = id
        self.a1 = None
        self.a2 = None
        self.exon = None
        self.covfrac = None
        self.msi = None

    def calc_msi(self):
        if self.a1 == 0 and self.a2 == 0 and self.exon == 0:
            self.msi = 0
        else:
            mis_splice_count = self.a1 + self.a2
            self.msi = (mis_splice_count / ((2 * self.exon) + mis_splice_count)) * 100

    def format(self):
        return "{i.id}\t{i.a1}\t{i.a2}\t{i.exon}\t{i.covfrac}\t{i.msi}\n".format(i=self)


def add_unique_id(feature):
    feature.name = "{}|{}|{}|{}".format(feature.chrom, feature.start, feature.stop, feature.strand)
    return feature


def a1(feature):
    feature.stop = feature.start + 3
    feature.start -= 3
    return feature


def a2(feature):
    feature.start = feature.end - 3
    feature.end += 3
    return feature


def exon(feature):
    feature.start -= 3
    feature.end += 3
    return feature


def main(args):
    # Logging
    logging.info("BED: {introns}, BAM: {bam}, GENOMESIZES: {gsizes} OUT: {out}".format(introns=args.introns,
                                                                                       bam=args.bam,
                                                                                       gsizes=args.gsizes,
                                                                                       out=args.out))

    introns_bedtools = pybedtools.BedTool(args.introns)
    bam_bedtools = pybedtools.BedTool(args.bam)

    intron_dict = dict()

    # Add a unique name to each of the introns
    introns_bedtools_named = introns_bedtools.each(add_unique_id).saveas()

    for cov in introns_bedtools_named.each(a1).sort(g=args.gsizes).intersect(bam_bedtools,
                                                                             f=1, split=True, c=True,
                                                                             sorted=True, g=args.gsizes):
        intron_dict.setdefault(cov.name, Intron(cov.name))
        setattr(intron_dict[cov.name], "a1", cov.count)

    for cov in introns_bedtools_named.each(a2).sort(g=args.gsizes).intersect(bam_bedtools,
                                                                             f=1, split=True, c=True,
                                                                             sorted=True, g=args.gsizes):
        setattr(intron_dict[cov.name], "a2", cov.count)

    for cov in introns_bedtools_named.each(exon).sort(g=args.gsizes).intersect(bam_bedtools, f=1, c=True, sorted=True,
                                                                               g=args.gsizes):
        setattr(intron_dict[cov.name], "exon", cov.count)

    for cov_frac in introns_bedtools_named.sort(g=args.gsizes).coverage(bam_bedtools, split=True, sorted=True,
                                                                        g=args.gsizes):
        setattr(intron_dict[cov_frac.name], "covfrac", float(cov_frac.fields[9]))

    with open(args.out, "w") as out_handle:
        out_handle.write("ID\tA1\tA2\tExon\tCovFrac\tMSI\n")
        for intron in intron_dict:
            intron_dict[intron].calc_msi()
            out_handle.write(intron_dict[intron].format())


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)


    def is_valid_file(parser, arg):
        """ Check if file exists """
        if not os.path.isfile(arg):
            parser.error('The file at %s does not exist' % arg)
        else:
            return arg


    epilog = "EXAMPLE: python " + os.path.basename(__file__) + \
             " --bed /path/to/introns.bed  --bam /path/to/reads.bam " \
             "--out /path/to/msi.xls"

    parser = argparse.ArgumentParser(description="Script to calculate MSI for given bam file using an intron bed",
                                     epilog=epilog)

    required_args_group = parser.add_argument_group('required arguments')

    required_args_group.add_argument('-i', '--bed', dest='introns', required=True,
                                     type=lambda x: is_valid_file(parser, x))
    required_args_group.add_argument('-b', '--bam', dest='bam', required=True,
                                     type=lambda x: is_valid_file(parser, x))
    required_args_group.add_argument('-g', '--genomesizes', dest='gsizes', required=True,
                                     type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-o', '--out', dest='out', default="msi.xls")
    args = parser.parse_args()
    # Call up the main
    main(args)
