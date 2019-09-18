# --------------------------- filter gtf ----------------------------------------
# This script will filter a GTF file to remove any given list of transcript types
# It will print the transcript/gene ids to an output file
# This output file can be then be used in combination with grep to filter the gtf.
# eg. grep -v -f outfile.txt base.gtf
# -------------------------------------------------------------------------------

import argparse
import logging
import os

import gffutils


def main(args):

    logging.info("GTF:{gtf_db}, Transtypes:{transtypes}, outfile:{out}".format(gtf_db=args.gtf_db,
                                                                               transtypes=args.transtypes,
                                                                               out=args.out))

    filtered_ids = set()

    # GTF db to be imported for filtering
    gtf_db = gffutils.FeatureDB(args.gtf_db)

    # Lets go through all the genes in the gtf
    for gene in gtf_db.features_of_type('gene'):
        # For every gene iterate over its transcripts
        for transcript in gtf_db.children(gene, featuretype='transcript'):
            # If it has any of the filtered attributes the process it
            if set(transcript.attributes['transcript_type']).intersection(args.transtypes):
                # If it is a single transcript gene we just remove that gene itself
                if len(list(gtf_db.children(gene, featuretype='transcript'))) == 1:
                    filtered_ids.add(gene.id)
                # Else remove just the transcripts
                else:
                    filtered_ids.add(transcript.id)

    # Write the filtered genes/transcripts to an outfile
    with open(args.out, "w") as out_handle:
        for filt_id in sorted(list(filtered_ids)):
            out_handle.write("{}\n".format(filt_id))


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)


    def is_valid_file(parser, arg):
        """ Check if file exists """
        if not os.path.isfile(arg):
            parser.error('The file at %s does not exist' % arg)
        else:
            return arg


    epilog = "EXAMPLE: python " + os.path.basename(__file__) + \
             " --gtf /path/to/gtf.db  --transtype retained_intron " \
             "--out /path/to/output.txt"

    parser = argparse.ArgumentParser(description="Script to give transcripts/genes of given transcript type",
                                     epilog=epilog)

    required_args_group = parser.add_argument_group('required arguments')
    required_args_group.add_argument('-g', '--gtfdb', dest='gtf_db', required=True,
                                     type=lambda x: is_valid_file(parser, x))
    required_args_group.add_argument('-t', '--transtypes', dest='transtypes', required=True, action="append")
    parser.add_argument('-o', '--out', dest='out', help="DEFAULT: {GTF_BASE}_filt.txt")
    args = parser.parse_args()
    if not args.out:
        args.out = "{base}_filt.txt".format(base=os.path.splitext(args.gtf_db)[0])

    main(args)
