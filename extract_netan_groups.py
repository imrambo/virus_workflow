#!/usr/bin/env python3
"""
Extract node and edge information from NetAn network output.
"""
import re
import os
import argparse
import logging

def connected_components(line):
    """
    Extract connected component members and group number within the network.
    """
    if re.match(r'^- Nodes in CC_\d+:', line):
        try:
            line_list = line.split(':')
            group_number = re.findall('Nodes in CC_(\d+):', line)[0]
        except:
            print('could not extract connected component group number')
        nodes = (group_number, ['{}\t{}\n'.format(group_number, n) for n in line_list[1].strip().split(',')])
        return nodes

    else:
        return None
#=============================================================================
#if __name__ == '__main__':

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, dest='infile', action='store',
help='input NetAn network file')

parser.add_argument('-o', '--output', type=str, dest='output', action='store',
help='output NetAn groups file')

parser.add_argument('-l', '--log', type=str, dest='log', action='store',
help='logfile')


args = parser.parse_args()

#Set up logger
logging_format = '%(levelname)s %(asctime)s - $(message)s'

if args.log:
    #Create the joblog directory if not specified
    if not os.path.exists(os.path.dirname(args.log)):
        os.makedirs(os.path.dirname(args.log))
    else:
        pass
    logging.basicConfig(filename = args.log, level = logging.DEBUG, format = logging_format)
else:
    logfile_default = 'extract_netan_groups.{}.joblog'.format(str(datetime.now().strftime('%Y-%m-%d_%H-%M-%S')))
    logging.basicConfig(filename = logfile_default,
    level = logging.DEBUG, format = logging_format)

logger = logging.getLogger()

try:
    input_handle = open(args.infile, 'r')
    logger.info('opened input handle')
except FileNotFoundError:
    input_error = 'could not open input file {}'.format(args.infile)
    logger.error(input_error)
    print(input_error)

try:
    output_handle = open(args.output, 'w')
    logger.info('opened output handle')
except IOError:
    output_error = 'could not write to output {}'.format(args.output)
    logger.error(output_error)
    print(output_error)

#Loop through the file and extract connected component members
#Write to output
lnum = 0
for line in input_handle:
    line = line.rstrip()
    if line.startswith('- Nodes in CC_'):
        cc = connected_components(line)
        if cc:
            cc_num_msg = 'CC number {}'.format(cc[0])
            cc_mem_msg = '{} members in CC number {}'.format(len(cc[1]), cc[0])
            logger.info(cc_num_msg)
            logger.info(cc_mem_msg)
            for c in cc[1]:
                output_handle.write(c)
        else:
            no_cc_msg = 'Line {}: found "Nodes in CC_", but returned no members'.format(lnum)
            logger.info(no_cc_msg)
    else:
        continue
    lnum += 1

input_handle.close()
logger.info('Closed input file')
output_handle.close()
logger.info('Closed output file')
