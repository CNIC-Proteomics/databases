#!/usr/bin/python

import argparse, logging, os

from lib import db

__author__ = 'jmrodriguezc'

def main(args):
    ''' Main function'''

    logging.info('create db_creator object')
    w = db.creator(args.species)

    logging.info('download raw database file')
    w.download_raw_db()

    logging.info('parse raw database file')
    w.parse_raw_db()

    logging.info('print database file')
    w.to_file(args.outdir)


if __name__ == "__main__":
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Create the System Biology database from UniProtKB',
        epilog='''
        Example:
            create_db_sb.py -s 'Sus scrofa' -o test
        ''')
    parser.add_argument('-s',  '--species', help='First filter based on the species name')
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

    logging.info('start '+os.path.basename(__file__))
    main(args)
    logging.info('end '+os.path.basename(__file__))