"""
Print Bloom filter memory requirements and number of hash functions k for n
items and target maximum false positive rate p.
"""
import argparse
import logging
import os
import sys
from typing import List

from humanfriendly import format_size
from .bloomfilter import BloomFilter, NoMemorySolution

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(os.environ.get('LOG_LEVEL', 'INFO'))
LOGGER.addHandler(logging.StreamHandler())

MSG_PARAMS = 'n=%s; p=%s'
MSG_OUTPUT = 'mem=%s; k=%s'


def parse_cmdargs(args: List[str]) -> argparse.Namespace:
    """
    Returns: Parsed arguments
    """
    parser = argparse.ArgumentParser(
            description=__doc__ + '\n' +
            ('Compare the values of mem and k:\n'
             'memcalc 1000000000 1e-6       # calculate minimum mem required\n'
             'memcalc 1000000000 1e-6 4GiB  # specify 4GiB\n'),
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('n_items',
                        type=int,
                        metavar='N_ITEMS',
                        help='Number of items to store.')
    parser.add_argument('fp_rate',
                        type=float,
                        metavar='FP_RATE',
                        help='Target false positive rate when n items are '
                        'stored.')
    parser.add_argument('mem',
                        nargs='?',
                        metavar='MEM',
                        help='Human-friendly mem size to allow for the Bloom '
                        'filter e.g. "4GiB". If not specified, the memory-'
                        'optimal (minimum) value will be calculated. The '
                        'advantage to specifying more than this is that the '
                        'number of hash functions k required to meet the '
                        'target false-postive rate p is reduced, giving '
                        'better performance. Note that k is very sensitive to '
                        'm around the minimum and as a rule of thumb allowing '
                        'just 1.25x the minimum mem roughly halves the value '
                        'of k. A warning is printed if n, p cannot be '
                        'satisfied with the specified memory.')
    return parser.parse_args(args)

def main() -> None:
    """
    Run as CLI script.
    """
    args = parse_cmdargs(sys.argv[1:])
    if args.mem:
        try:
            m, k = BloomFilter.m_k_mem(args.n_items, args.fp_rate, args.mem)
        except NoMemorySolution as nms:
            LOGGER.warning(nms)
            sys.exit(1)
    else:
        m, k = BloomFilter.m_k_min(args.n_items, args.fp_rate)
    mem = args.mem or format_size(m // 8, keep_width=True,
                                  binary=(m & (m-1) == 0)) 
    LOGGER.info(MSG_PARAMS % (args.n_items, args.fp_rate))
    print(MSG_OUTPUT % (mem, k))


if __name__ == '__main__':
    main()
