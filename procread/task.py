#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from .util import Shell


def do_qc_checks(config, work_dir, cpus):
    logging.info('Start quality control checks for reads')
    qc_dir = os.path.join(work_dir, 'qc')
    os.makedirs(qc_dir, exist_ok=True)
    logging.debug('qc_dir: {}'.format(qc_dir))
    sh = Shell(log_txt=os.path.join(qc_dir, 'fastqc_log.txt'))
    sh.run('fastqc --version')
    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        d = os.path.join(qc_dir, '{0}_{1}'.format(t, r))
        os.makedirs(d, exist_ok=True)
        sh.run(
            arg_str='fastqc {0} --threads {1} --outdir {2} {3}'.format(
                _cmd_arg_str(config, command='fastqc'), cpus, d,
                config['fastq'][t][r]
            )
        )


def trim_adapters(config, work_dir, cpus):
    logging.info('Trim adapter sequences in reads')
    trim_dir = os.path.join(work_dir, 'trim')
    os.makedirs(trim_dir, exist_ok=True)
    logging.debug('trim_dir: {}'.format(trim_dir))
    sh = Shell(log_txt=os.path.join(trim_dir, 'cutadapt_log.txt'))
    sh.run('cutadapt --version')
    for t in ['foreground', 'background']:
        trimmed_fq = {
            k: os.path.join(
                trim_dir,
                re.sub(
                    r'\.fastq\..*$', '.trimmed.{}.fastq.gz'.format(k),
                    os.path.basename(config['fastq'][t][v])
                )
            )
            for k, v in {'r1': 'read1', 'r2': 'read2'}.items()
        }
        sh.run(
            arg_str='cutadapt {0} -a {1} -A {2} -o {3} -p {4} {5} {6}'.format(
                _cmd_arg_str(config, command='cutadapt'),
                config['adapter']['3prime'], config['adapter']['5prime'],
                trimmed_fq['r1'], trimmed_fq['r2'],
                config['fastq'][t]['read1'], config['fastq'][t]['read2']
            )
        )


def _cmd_arg_str(config, command):
    return (
        ' '.join(config['command_args'][command])
        if 'command_args' in config and command in config['command_args']
        else ' '
    )


def map_reads(config, work_dir, cpus):
    pass


def call_variants(config, work_dir, cpus):
    pass
