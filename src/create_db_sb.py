#!/usr/bin/python

import argparse, logging, os, sys

from lib import db

__author__ = 'jmrodriguezc'

def main(args):
    ''' Main function'''

    # check input parameters
    if args.infasta and not args.regex:
        sys.exit("The input database and the regular expression parameters have to be together")
        
    logging.info("create db_creator object")
    w = db.creator(args.species, args.outdir, args.infasta, args.filt)

    logging.info("extract list of identifiers")
    ids = w.extract_identifiers(args.regex)

    # logging.info("download raw file")
    # w.download_raw_dbs(args.filt)

    logging.info("extract categories by "+ args.type)
    output = w.extract_categories(ids, args.type)

    # # 2. REMOVE DUPLICATIONS

    logging.info('print database file')
    w.to_file(output)


if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Create the System Biology database from UniProtKB',
        epilog='''
Examples:
    create_db_sb.py -s pig -o test
    create_db_sb.py -s mouse -f -t gene -o test
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-s',  '--species', required=True, help='First filter based on the species name')
    parser.add_argument('-i',  '--infasta', help='Input database (as FASTA format)')
    parser.add_argument('-r',  '--regex',  help='Regular expression to extract the ID (protein or gene) from the comment of FASTA input file')
    parser.add_argument('-t',  '--type',  default="protein", choices=["protein","gene"], help='Directory where the database will be saved')
    parser.add_argument('-f',  '--filt',  default="no", choices=["no","sw+tr","sw","tr"], help='Directory where the database will be saved')    
    parser.add_argument('-o',  '--outdir', required=True, help='Directory where the database will be saved')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')
