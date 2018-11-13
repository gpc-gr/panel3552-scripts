#
#
#

import io
import subprocess


def read_gzip_text(path):
    process = subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, bufsize=-1)
    with io.open(process.stdout.fileno(), closefd=False) as fin:
        for line in fin:
            yield line

    process.wait()
    assert process.returncode == 0

