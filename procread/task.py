#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from .util import Shell, ProcreadRuntimeError


def prepare_paths(config, work_dir, cpus):
    logging.info('Prepare paths to directories and files')
    dir_dict = {
        d: os.path.join(work_dir, d)
        for d in ['input', 'qc', 'trim', 'map', 'call']
    }

    for d in [work_dir, dir_dict['input']]:
        os.makedirs(d, exist_ok=True)
    sh = Shell(log_txt=os.path.join(dir_dict['input'], 'command_log.txt'))
    for c in ['pigz', 'pbzip2']:
        sh.run('{} --version'.format(c))

    fq_dict = {
        '{0}_{1}'.format(t, r):
        os.path.join(
            dir_dict['input'],
            re.sub(
                r'\.(fastq|fq)\.?[^\.]*\.?[^\.]*$',
                '.{}.fastq.gz'.format({'read1': 'r1', 'read2': 'r2'}[r]),
                os.path.basename(config['path']['fastq'][t][r])
            )
        )
        for t, r in product(['foreground', 'background'], ['read1', 'read2'])
    }
    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        fq_gz = fq_dict['{0}_{1}'.format(t, r)]
        if os.path.isfile(fq_gz):
            continue

        fq_src = config['path']['fastq'][t][r]
        fq_src_ext = os.path.splitext(fq_src)[1]
        if fq_src_ext in ['.fastq', '.fq']:
            sh.run('pigz -p {0} -c {1} > {2}'.format(
                cpus, fq_src, fq_gz
            ))
        elif fq_src_ext == '.gz':
            os.symlink(fq_src, fq_gz)
        elif fq_src_ext == '.bz2':
            sh.run('pbzip2 -p# {0} -dc {1} | pigz -p {0} -c - > {2}'.format(
                cpus, fq_src, fq_gz
            ))
        else:
            raise ProcreadRuntimeError(
                'Supported extension for input fastq: '
                '.fastq, .fastq.gz, .fastq.*.gz, .fastq.bz2, .fastq.*.bz2, '
                '.fq, .fq.gz, .fq.*.gz, .fq.bz2, .fq.*.bz2'
            )

    ref_src = config['path']['fasta']['reference']
    ref_fa = os.path.join(
        dir_dict['input'],
        re.sub(
            r'\.(fa|fasta)\.?[^\.]*$', '.ref.fa', os.path.basename(ref_src)
        )
    )
    ref_dict = {'fasta': ref_fa, 'faidx': ref_fa + '.fai'}
    ref_src_ext = os.path.splitext(ref_src)[1]
    if not os.path.isfile(ref_fa):
        if ref_src_ext in ['.fasta', '.fa']:
            os.symlink(ref_src, ref_fa)
        elif ref_src_ext == '.gz':
            sh.run('pigz -p {0} -dc {1} > {2}'.format(
                cpus, ref_src, ref_fa
            ))
        elif ref_src_ext == '.bz2':
            sh.run('pbzip2 -p# {0} -dc {1} > {2}'.format(
                cpus, ref_src, ref_fa
            ))
        else:
            raise ProcreadRuntimeError(
                'Supported extension for input reference fasta: '
                '.fa, .fa.gz, .fa.bz2, .fasta, .fasta.gz, .fasta.bz'
            )

    return {
        'dir': dir_dict,
        'fastq': fq_dict,
        'ref': ref_dict
    }


def make_ref_index(paths):
    if os.path.isfile(paths['ref']['faidx']):
        return

    logging.info('Make reference index files')
    sh = Shell(
        log_txt=os.path.join(paths['dir']['input'], 'command_log.txt'),
        format_log=False
    )
    for c in ['samtools', 'bwa']:
        sh.run(c)
    sh.run('samtools faidx {}'.format(paths['ref']['fasta']))
    sh.run('bwa index -p {0} {1}'.format(
        paths['ref']['faidx'], paths['ref']['fasta']
    ))


def do_qc_checks(config, paths, cpus):
    if os.path.isfile(paths['dir']['qc']):
        return

    logging.info('Do quality control checks for reads')
    os.makedirs(paths['dir']['qc'], exist_ok=True)
    sh = Shell(log_txt=os.path.join(paths['dir']['qc'], 'command_log.txt'))
    sh.run('fastqc --version')

    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        tag = '{0}_{1}'.format(t, r)
        sh.run(
            arg_str='fastqc {0} --threads {1} --outdir {2} {3}'.format(
                _cmd_arg_str(config, command='fastqc'), cpus,
                paths['dir']['qc'], paths['fastq'][tag]
            )
        )


def trim_adapters(config, paths):
    if os.path.isfile(paths['dir']['trim']):
        return

    logging.info('Trim adapter sequences in reads')
    os.makedirs(paths['dir']['trim'], exist_ok=True)
    sh = Shell(log_txt=os.path.join(paths['dir']['trim'], 'command_log.txt'))
    sh.run('cutadapt --version')

    for t in ['foreground', 'background']:
        tag_dict = {
            k: '{0}_{1}'.format(t, r)
            for k, r in {'r1': 'read1', 'r2': 'read2'}.items()
        }
        trimmed_fq = {
            k: os.path.join(
                paths['dir']['trim'],
                re.sub(
                    r'(\.r[12]\.fastq\.gz)$', r'\.trimmed\1',
                    os.path.basename(paths['fastq'][v])
                )
            )
            for k, v in tag_dict.items()
        }
        sh.run(
            arg_str='cutadapt {0} -a {1} -A {2} -o {3} -p {4} {5} {6}'.format(
                _cmd_arg_str(config, command='cutadapt'),
                config['adapter']['3prime'], config['adapter']['5prime'],
                trimmed_fq['r1'], trimmed_fq['r2'],
                paths['fastq'][tag_dict['r1']], paths['fastq'][tag_dict['r1']]
            )
        )


def map_reads(config, paths, cpus):
    if os.path.isfile(paths['dir']['map']):
        return

    logging.info('Map reads to a reference')


def call_variants(config, paths, cpus):
    if os.path.isfile(paths['dir']['call']):
        return

    logging.info('Call SNVs/indels')


def _cmd_arg_str(config, command):
    return (
        ' '.join(config['command_args'][command])
        if 'command_args' in config and command in config['command_args']
        else ' '
    )
