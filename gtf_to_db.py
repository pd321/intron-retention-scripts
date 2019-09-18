# ---------------------- get intron type -------------------------------------
# Convert gtf to an sqlite db which makes filtering and querying easier
# ----------------------------------------------------------------------------

import argparse
import logging
import os

import gffutils


def main(cmds):
    # Logging
    logging.info("GTF:{gtf}, GTFDB:{gtfdb}".format(gtf=args.gtf, gtfdb=args.gtfdb))

    # Now do the actual conversion
    db = gffutils.create_db(args.gtf,
                            args.gtfdb,
                            disable_infer_genes=True,
                            disable_infer_transcripts=True)


if __name__ == '__main__':

    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)


    def is_valid_file(parser, arg):
        """ Check if file exists """
        if not os.path.isfile(arg):
            parser.error('The file at %s does not exist' % arg)
        else:
            return arg


    epilog = "EXAMPLE: python " + os.path.basename(__file__) + \
             " --gtf /path/to/gtf.db  --gtfdb /path/to/gtf.db "

    parser = argparse.ArgumentParser(description="Script to convert gtf to a db for faster processing",
                                     epilog=epilog)
    parser.add_argument('-d', '--gtfdb', dest='gtfdb', help="DEFAULT: {GTF_BASE}.db")
    required_args_group = parser.add_argument_group('required arguments')
    required_args_group.add_argument('-g', '--gtf', dest='gtf', required=True,
                                     type=lambda x: is_valid_file(parser, x))
    args = parser.parse_args()
    if not args.gtfdb:
        args.gtfdb = "{base}.db".format(base=os.path.splitext(args.gtf)[0])

    main(args)
