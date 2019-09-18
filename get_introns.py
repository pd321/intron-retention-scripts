# ---------------------- get introns -----------------------------------------
# Script to get a list of clean(starts not overlapping with exons) introns
# ----------------------------------------------------------------------------

import argparse
import logging
import os
import tempfile

import gffutils
import pybedtools


def get_exons(db):
    def gen():
        for exon in db.features_of_type("exon"):
            yield (exon.chrom, exon.start, exon.end)

    return pybedtools.BedTool(gen())


# This will print the introns in the format we require
def get_final_intron_bedtools(final_intron_dict):
    def gen():
        for intron_base_id in final_intron_dict.keys():
            intron_details = tuple(intron_base_id.split("|"))
            # Sorting the intron ids so that we get consistent results every time
            # Adjustment to make is zero or one based
            # zero based
            intron_start = int(intron_details[1]) - 1
            # one based
            # intron_start = intron_details[1]

            yield (intron_details[0], intron_start,
                   intron_details[2],
                   ",".join(sorted(final_intron_dict[intron_base_id])),
                   ".",
                   intron_details[3])

    return pybedtools.BedTool(gen())


def main(args):
    # Logging
    logging.info("GTF: {gtf_db}, WindowSize: {window} Outfile: {out}".format(gtf_db=args.gtf_db, window=args.window,
                                                                             out=args.out))

    # Vars
    window_size = args.window
    actual_window = window_size - 1

    # Read in the GTF db
    db = gffutils.FeatureDB(args.gtf_db)

    # Get the list of introns and add a flank to it
    flanked_introns = list()
    logging.info("Creating a list of introns. This will take a while. Be Patient")
    introns = list(db.create_introns())
    logging.info("Intron List Created")
    for intron in introns:
        # Create a unique id
        intron_id = "{}|{}|{}|{}|{}".format(intron.attributes["transcript_id"][0],
                                            intron.chrom, intron.start, intron.end, intron.strand)
        # Create the left flank
        flanked_introns.append((intron.chrom, intron.start, intron.start + actual_window,
                                intron_id, intron.strand))
        # Create the right flank
        flanked_introns.append((intron.chrom, intron.end - actual_window, intron.end, intron_id, intron.strand))

    intron_bedtool = pybedtools.BedTool(flanked_introns)

    # Get a list of exons
    exons = get_exons(db=db)

    # Find overlapping introns
    logging.info("Finding introns overlapping exons")
    overlapping_introns = intron_bedtool.intersect(exons, wa=True, f=1)
    overlapping_introns_id = set([over_intron.name for over_intron in overlapping_introns])

    # Print the final list of introns
    # Medge introns having the same id for flanks
    final_introns = dict()
    logging.info("Collapsing introns")
    for intron in introns:
        intron_base_id = "{}|{}|{}|{}".format(intron.chrom, intron.start, intron.end, intron.strand)
        intron_crosscheck_id = "{}|{}".format(intron.attributes["transcript_id"][0], intron_base_id)
        # Checking if it was overlapping with any of the exons at the junction
        if intron_crosscheck_id not in overlapping_introns_id:
            final_introns.setdefault(intron_base_id, set())
            final_introns[intron_base_id].add(intron.attributes["gene_id"][0])

    logging.info("Writing output file")
    tmp_dump_file = tempfile.NamedTemporaryFile()
    to_print_bedtools = get_final_intron_bedtools(final_introns)
    # pybedtools has problems directly taking the bedtool for sorting
    # Hence the need to write
    to_print_bedtools.moveto(tmp_dump_file.name).sort().moveto(args.out)
    tmp_dump_file.close()


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)


    def is_valid_file(parser, arg):
        """ Check if file exists """
        if not os.path.isfile(arg):
            parser.error('The file at %s does not exist' % arg)
        else:
            return arg


    epilog = "EXAMPLE: python " + os.path.basename(__file__) + \
             " --gtf /path/to/gtf.db  --out /path/to/introns.bed"

    parser = argparse.ArgumentParser(description="Script to give filtered introns from given gtf db",
                                     epilog=epilog)

    required_args_group = parser.add_argument_group('required arguments')
    required_args_group.add_argument('-d', '--gtfdb', dest='gtf_db', required=True,
                                     type=lambda x: is_valid_file(parser, x))

    parser.add_argument('-w', '--window', dest='window', default=3, help="DEFAULT: 3")
    parser.add_argument('-o', '--out', dest='out', help="DEFAULT: {GTF_BASE}_introns.bed")

    args = parser.parse_args()

    if not args.out:
        args.out = "{base}_introns.bed".format(base=os.path.splitext(args.gtf_db)[0])
    # Call up the main
    main(args)
