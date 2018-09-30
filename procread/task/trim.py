#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from ..util.shell_operator import ShellOperator


def trim_adapters(cf, cpus):
    logger = logging.getLogger(__name__)
    output_files = [
        os.path.join(
            cf['paths']['dir']['trim'],
            re.sub(
                r'(\.r[12]\.fastq\.gz)$', r'.trimmed\1',
                os.path.basename(d['{0}_{1}'.format(t, r)])
            )
        )
        for d, t, r in product(cf['paths']['fastq'], ['read1', 'read2'])
    ]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip adapter trimming')
        return
    else:
        os.makedirs(cf['paths']['dir']['trim'], exist_ok=True)

    logger.info('Trim adapter sequences in reads')
    sh = ShellOperator(
        log_txt=os.path.join(
            cf['paths']['dir']['trim'], 'trim_log.txt'
        )
    )
    sh.run('cutadapt --version')

    io = [
        {
            'out': {
                r: os.path.join(
                    cf['paths']['dir']['trim'],
                    re.sub(
                        r'(\.r[12]\.fastq\.gz)$', r'.trimmed\1',
                        os.path.basename(d[r])
                    )
                )
                for r in ['read1', 'read2']
            },
            'in': {
                r: d[r] for r in ['read1', 'read2']
            }
        }
        for d in cf['paths']['fastq']
    ]

    sh.run([
        'cutadapt {0} -o {1} -p {2} {3} {4}'.format(
            cf['cmd_args']['cutadapt'], cpus, f['out']['read1'],
            f['out']['read2'], f['in']['read1'], f['in']['read2']
        ) for f in io
        if not all([os.path.isfile(p) for p in f['out'].values()])
    ])
