#!/usr/bin/env pypy

from __future__ import print_function

import argparse


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta')
    args = parser.parse_args()

    #
    with open(args.fasta) as fin:
        for id, header, sequence in _parse_fasta(fin):
            if id == 'Y':
                header += ' (XTR masked)'
                sequence = _mask_par_on_chr_y(sequence)

            print('>{}'.format(header))
            print(_wrap(sequence))
            print()


def _parse_fasta(stream):
    id = None
    header = None
    fragments = []

    for line in stream:
        line = line.strip()
        if not line:
            continue

        if line.startswith('>'):
            if header:
                yield id, header, ''.join(fragments)

            id = line[1:].split()[0]
            header = line[1:]
            fragments = []

        else:
            fragments.append(line)

    if header:
        yield id, header, ''.join(fragments)


def _mask_par_on_chr_y(sequence):
    pars_on_chr_y = [
        (10001, 2649520),       # PAR1 (GRCh37)
        (59034050, 59363566),   # PAR2
        (2917959, 6616600)      # XTR
    ]

    original_length = len(sequence)
    for index, (start, end) in enumerate(pars_on_chr_y):
        head = sequence[:start-1]
        par = sequence[start-1:end-1]
        tail = sequence[end-1:]

        if index in (0, 1):
            assert set(par) == set('N')

        sequence = head + ('N' * len(par)) + tail

    assert len(sequence) == original_length
    return sequence


def _wrap(text, width=80):
    return '\n'.join([text[i:i+width] for i in range(0, len(text), width)])


if __name__ == '__main__':
    main()
