"""
Print approximate memory requirement in GB for n items and target maximum
false positive rate p.
"""
import argparse
import sys
from typing import List

from .bloomfilter import BloomFilter


def mem_calc(n: int, p: int) -> float:
    """
    Returns approximate memory requirement in GB for n items and target maximum
    false positive rate p.
    """
    m, _ = BloomFilter.optimal_m_k(n, p)
    return m / 8 / 1024 ** 3


def parse_cmdargs(args: List[str]) -> argparse.Namespace:
    """
    Returns: Parsed arguments
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('n', type=int, help='Number of items to store.')
    parser.add_argument('p', type=float, help='Target false positive rate '
                        'when n items are stored.')
    return parser.parse_args(args)

def main() -> None:
    """
    Run as CLI script.
    """
    args = parse_cmdargs(sys.argv[1:])
    print(f'{mem_calc(args.n, args.p):0.3f}GB')


if __name__ == '__main__':
    main()
