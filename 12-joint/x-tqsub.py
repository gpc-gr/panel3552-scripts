#!/home/gpc-gr/panel-build37/work/env3/bin/python3

from __future__ import print_function

import argparse
import json
import os
import subprocess

import jinja2


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('template', type=os.path.abspath)
    parser.add_argument('output', type=os.path.abspath)
    parser.add_argument('context', type=json.loads, nargs='+')
    parser.add_argument('--dry-run', action='store_true', default=False)
    parser.add_argument('--sync', action='store_true', default=False)
    args = parser.parse_args()

    #
    merged_context = {
        'script_name': os.path.basename(args.output),
        'script_path': args.output
    }
    for context in args.context:
        merged_context.update(context)

    #
    if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))

    script = _get_template(args.template).render(merged_context)
    with open(args.output, 'w') as fout:
        fout.write(script)

    #
    if args.dry_run:
        return

    command = 'bash' if args.sync else 'qsub'
    subprocess.check_call([command, args.output])


def _get_template(path):
    template_dir, template_name = os.path.split(path)
    
    loader = jinja2.FileSystemLoader(template_dir)
    environment = jinja2.Environment(loader=loader)
    environment.undefined = jinja2.StrictUndefined
    
    return environment.get_template(template_name)


if __name__ == '__main__':
    main()

