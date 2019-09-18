# ---------------------- get intron type -------------------------------------
# Classify each intron into U2/U12 type using PWM's from
# [splicerack](http://katahdin.mssm.edu/splice/index.cgi?database=spliceNew)
# Currently hardcoded for human introns
# TODO: expose pwm params
# ----------------------------------------------------------------------------
import argparse
import logging
import os

import pybedtools
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner


class Intron(object):
    """
    Class encapsulating the decision logic explicitly
    """
    def __init__(self, id):
        """
        :param id: Unique id for the intron
        """
        self.id = id
        self.AT_AC_U12_b = None
        self.GT_AG_U12_b = None
        self.AT_AC_U12_d = None
        self.GT_AG_U12_d = None
        self.GT_AG_U2_d = None
        self.GC_AG_U2_d = None
        self.type = None
        self.subtype = None
        self.confidence = None

    def get_type(self):
        """
        Classify as U12/U2 depending on PSM scores.
        """
        # Explicit is better than implicit
        # Clear decision logic which may be reviewed later.
        GT_AG_U12_diff = min(self.GT_AG_U12_d - self.GT_AG_U2_d, self.GT_AG_U12_d - self.GC_AG_U2_d)
        AT_AC_U12_diff = min(self.AT_AC_U12_d - self.GT_AG_U2_d, self.AT_AC_U12_d - self.GC_AG_U2_d)
        GC_AG_U2_diff = min(self.GC_AG_U2_d - self.GT_AG_U12_d, self.GC_AG_U2_d - self.AT_AC_U12_d)
        GT_AG_U2_diff = min(self.GT_AG_U2_d - self.GT_AG_U12_d, self.GT_AG_U2_d - self.AT_AC_U12_d)

        if GT_AG_U12_diff >= 25 and self.GT_AG_U12_d > self.AT_AC_U12_d:
            self.type = "U12"
            self.subtype = "GT_AG_U12"
            self.confidence = "High"
        elif GT_AG_U12_diff >= 10 and self.GT_AG_U12_d > self.AT_AC_U12_d and self.GT_AG_U12_b >= 70:
            self.type = "U12"
            self.subtype = "GT_AG_U12"
            self.confidence = "Mid"
        elif AT_AC_U12_diff >= 25 and self.AT_AC_U12_d > self.GT_AG_U12_d:
            self.type = "U12"
            self.subtype = "AT_AC_U12"
            self.confidence = "High"
        elif AT_AC_U12_diff >= 10 and self.AT_AC_U12_d > self.GT_AG_U12_d and self.AT_AC_U12_b >= 70:
            self.type = "U12"
            self.subtype = "AT_AC_U12"
            self.confidence = "Mid"
        elif GC_AG_U2_diff >= 25 and self.GC_AG_U2_d > self.GT_AG_U2_d:
            self.type = "U2"
            self.subtype = "GC_AG_U2"
            self.confidence = "High"
        elif GT_AG_U2_diff >= 25 and self.GT_AG_U2_d > self.GC_AG_U2_d:
            self.type = "U2"
            self.subtype = "GT_AG_U2"
            self.confidence = "High"
        else:
            self.type = "U2"
            self.subtype = "GT_AG_U2"
            self.confidence = "Low"

    def format(self):
        """
        Format intron for printing
        """
        chrom, start, end, strand = self.id.split("|")
        return "{}\t{}\t{}\t{}\t{i.type}\t{i.subtype}\t{i.confidence}" \
               "\t{i.AT_AC_U12_d}\t{i.GT_AG_U12_d}\t{i.GT_AG_U2_d}" \
               "\t{i.GC_AG_U2_d}\t{i.AT_AC_U12_b}\t{i.GT_AG_U12_b}\n".format(chrom, start, end, strand, i=self)


def add_unique_id(feature):
    """
    Adds a unique id to each feature
    :param feature: pybedtools feature
    :return: updated feature with unique id added
    """
    feature.name = "{}|{}|{}|{}".format(feature.chrom, feature.start, feature.stop, feature.strand)
    return feature


def adj_don(feature):
    """
    Get adjacent donor site for each intron handling strand
    :param feature: pybedtools feature
    :return: updated feature to reflect donor site positions
    """
    if feature.strand == "+":
        feature.stop = feature.start + 10
        feature.start -= 3
    elif feature.strand == "-":
        feature.start = feature.stop - 10
        feature.stop += 3
    return feature


def adj_branch(feature):
    """
    Get adjacent branch site for each intron handling strand
    :param feature: pybedtools feature
    :return: updated feature to reflect branch site positions
    """
    if feature.strand == "+":
        feature.start = feature.stop - 38
        feature.stop = feature.stop - 8
    elif feature.strand == "-":
        feature.stop = feature.start + 38
        feature.start = feature.start + 8
    return feature


def get_fa(intbed, genomefa, type):
    """
    Get fasta sequence for given branch/donor
    :param intbed: intron coords bed file
    :param genomefa: path to ref genome
    :param type: branch/donor
    :return:
    """
    bedtool = pybedtools.BedTool(intbed).each(add_unique_id)
    bedtool_adj = bedtool.each(adj_don) if type == "donor" else bedtool.each(adj_branch)
    return bedtool_adj.sequence(genomefa, name=True)


def get_motif(motif_pwm):
    """
    Extract motifs from a pwm file
    :param motif_pwm: pwm file location
    :return: extracted motifs
    """
    motifs = read_motifs(motif_pwm)
    for motif in motifs:
        motif.pwm_max_score()
        motif.pwm_min_score()
    return motifs


def rescale(org_score, orig_range, new_range):
    """
    rescale the pwm match scores to apply uniform cutoffs
    :param org_score: original pwm score
    :param orig_range: range to adjust within
    :param new_range: range to adjust to
    :return: new score
    """
    delta_orig = orig_range[1] - orig_range[0]
    delta_new = new_range[1] - new_range[0]

    return (delta_new * (org_score - orig_range[0]) / delta_orig) + new_range[0]


def get_motif_scores(fa, motifs):
    s = Scanner()
    s.set_motifs(motifs)
    s.set_threshold(threshold=0.0)
    seqs = Fasta(fa.seqfn)
    for i, result in enumerate(s.scan(seqs, nreport=1)):
        intron_id = seqs.ids[i]
        for m, matches in enumerate(result):
            motif = motifs[m]
            for score, pos, strand in matches:
                if score < 0:
                    score_rescaled = rescale(score, orig_range=[motif.min_score, 0], new_range=[0, 50])
                else:
                    score_rescaled = rescale(score, orig_range=[0, motif.max_score], new_range=[50, 100])
                yield (intron_id, motif.id, score_rescaled)


def main(args):
    logging.info("Received the following args: \n {}".format(args))

    # Vars
    intron_dict = dict()

    logging.info("Getting fa seqeunces for donor and branches")
    don_fa = get_fa(intbed=args.intbed, genomefa=args.genomefa, type="donor")
    branch_fa = get_fa(intbed=args.intbed, genomefa=args.genomefa, type="branch")

    logging.info("Parsing donor and branch pwm files")
    don_motifs = get_motif(args.donorpwm)
    branch_motifs = get_motif(args.branchpwm)

    logging.info("Scoring fa seqeunces for motifs")
    for intron_id, motif_id, score_rescaled in get_motif_scores(don_fa, don_motifs):
        intron_dict.setdefault(intron_id, Intron(id=intron_id))
        setattr(intron_dict[intron_id], "{}_d".format(motif_id), score_rescaled)
    for intron_id, motif_id, score_rescaled in get_motif_scores(branch_fa, branch_motifs):
        setattr(intron_dict[intron_id], "{}_b".format(motif_id), score_rescaled)

    logging.info("Writing output")
    with open(args.out, "w") as out_handle:
        out_handle.write("Chrom\tStart\tEnd\tStrand\tType\tSubType\tConfidence"
                         "\tAT_AC_U12_d\tGT_AG_U12_d\tGT_AG_U2_d\tGC_AG_U2_d\tAT_AC_U12_b\tGT_AG_U12_b\n")
        for intron in intron_dict:
            intron_dict[intron].get_type()
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
             " --bed /path/to/introns.bed  --branch /path/to/branch.pwm --don /path/to/don.pwm" \
             " --genome /path/to/genome.fa --out /path/to/intronType.xls"

    parser = argparse.ArgumentParser(description="Script to classify every intron as U2/U12",
                                     epilog=epilog)

    required_args_group = parser.add_argument_group('required arguments')

    required_args_group.add_argument('-i', '--intbed',
                                     dest='intbed',
                                     required=True,
                                     help="Intron locations in bed format",
                                     type=lambda x: is_valid_file(parser, x))

    required_args_group.add_argument('-b', '--branchpwm',
                                     dest='branchpwm',
                                     required=True,
                                     help="Branch PWM for U12",
                                     type=lambda x: is_valid_file(parser, x))

    required_args_group.add_argument('-d', '--donorpwm',
                                     dest='donorpwm',
                                     required=True,
                                     help="Donor PWM for U12",
                                     type=lambda x: is_valid_file(parser, x))

    required_args_group.add_argument('-g', '--genomefa',
                                     dest='genomefa',
                                     required=True,
                                     help="Genome fasta file",
                                     type=lambda x: is_valid_file(parser, x))

    parser.add_argument('-o', '--out',
                        dest='out',
                        help="Out file for intron types",
                        default="inttype.xls")

    args = parser.parse_args()

    main(args)
