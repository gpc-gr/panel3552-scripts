#
#
#

from __future__ import print_function

import json
import os
import re
import subprocess

import jinja2


def tqsub(template_path, script_path, context, dry_run=False, sync=False, ignore_if_already_queued=False):
    #
    if 'script_name' not in context:
        context['script_name'] = os.path.basename(script_path)
    if 'script_path' not in context:
        context['script_path'] = script_path

    #
    script = _get_template(template_path).render(context)

    if ignore_if_already_queued:
        script_name = _get_script_name(script) or os.path.basename(script_path)
        if script_name in _get_queued_script_name():
            return

    #
    if not os.path.exists(os.path.dirname(script_path)):
        os.makedirs(os.path.dirname(script_path))

    with open(script_path, 'w') as fout:
        fout.write(script)

    #
    if dry_run:
        return

    command = 'bash' if sync else 'qsub'
    subprocess.check_call([command, script_path])


def _get_template(path):
    template_dir, template_name = os.path.split(path)
    
    loader = jinja2.FileSystemLoader(template_dir)
    environment = jinja2.Environment(loader=loader)
    environment.undefined = jinja2.StrictUndefined
    
    return environment.get_template(template_name)


def _get_script_name(script):
    for line in script.splitlines():
        if line.startswith('#$'):
            match = re.search('-N ([a-zA-Z0-9_\-\.]+)', line)
            if match:
                return match.group(1)

    return None


def _get_queued_script_name():
    script_names = set()

    output = subprocess.check_output(['qstat', '-xml']).decode('utf-8')
    for line in output.splitlines():
        line = line.strip()
        if line.startswith('<JB_name>') and line.endswith('</JB_name>'):
            script_names.add(line[len('<JB_name>'):-len('</JB_name>')])
   
    return script_names


def merge_contexts(*contexts):
    merged_context = {}
    for context in reversed(contexts):
        merged_context.update(context)

    return merged_context

