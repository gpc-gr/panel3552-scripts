#!/usr/bin/env python3

from __future__ import print_function

import argparse
import json
import os

from ..libs.tqsub import merge_contexts, tqsub


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('template', type=os.path.abspath)
    parser.add_argument('script', type=os.path.abspath)
    parser.add_argument('contexts', metavar='context', nargs='*', type=json.loads)
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--sync', action='store_true')
    parser.add_argument('--ignore-if-already-queued', action='store_true')
    args = parser.parse_args()

    #
    tqsub(
        args.template,
        args.script,
        merge_contexts(*args.contexts),
        dry_run=args.dry_run,
        sync=args.sync,
        ignore_if_already_queued=args.ignore_if_already_queued)
    

if __name__ == '__main__':
    main()

