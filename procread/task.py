#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from .util import Shell, ProcreadRuntimeError


def prepare_paths(config, work_dir, cpus):
    logging.info('Prepare paths to directories and files')
    wd = os.path.abspath(work_dir)
    dir_dict = {
        d: os.path.join(wd, d)
        for d in ['input', 'qc', 'trim', 'map', 'call']
    }

    for d in [wd, dir_dict['input']]:
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

        fq_src = os.path.abspath(config['path']['fastq'][t][r])
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

    ref_src = os.path.abspath(config['path']['fasta']['reference'])
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
    if os.path.isdir(paths['dir']['qc']):
        return
    else:
        os.makedirs(paths['dir']['qc'])

    logging.info('Do quality control checks for reads')
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
    if os.path.isdir(paths['dir']['trim']):
        return
    else:
        os.makedirs(paths['dir']['trim'])

    logging.info('Trim adapter sequences in reads')
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


def map_reads(paths, cpus):
    if os.path.isdir(paths['dir']['map']):
        return
    else:
        os.makedirs(paths['dir']['map'])

    logging.info('Map reads to a reference')
    sh = Shell(log_txt=os.path.join(paths['dir']['map'], 'command_log.txt'))

    if os.path.isdir(paths['dir']['trim']):
        fq_dict = {
            k: os.path.join(
                paths['dir']['trim'],
                re.sub(
                    r'(\.r[12]\.fastq\.gz)$', r'\.trimmed\1',
                    os.path.basename(v)
                )
            )
            for k, v in paths['fastq'].items()
        }
    else:
        fq_dict = paths['fastq']

    bam_dict = {
        t: {
            b: os.path.join(paths['dir']['map'], t, '{}.bam'.format(b))
            for b in ['sort', 'fixmate', 'markdup']
        }
        for t in ['foreground', 'background']
    }

    for t in ['foreground', 'background']:
        os.makedirs(os.path.join(paths['dir']['map'], t))
        bd = bam_dict[t]
        sh.run(
            'bwa mem -t {0} {1} {2} {3} '
            '| samtools view -@ {0} -bS - '
            '| samtools sort -@ {0} -o {4} -'.format(
                cpus, paths['ref']['fasta'], fq_dict[t + '_read1'],
                fq_dict[t + '_read2'], bd['sort']
            )
        )
        sh.run(
            'samtools sort -n -@ {0} {1} '
            '| samtools fixmate -m - {2}'.format(
                cpus, bd['sort'], bd['fixmate']
            )
        )
        sh.run(
            'samtools sort -@ {0} {1} '
            '| samtools markdup - {2}'.format(
                cpus, bd['fixmate'], bd['markdup']
            )
        )
        for p in bd.values():
            sh.run('samtools index {}'.format(p))
            sh.run('samtools flagstat {0} > {1}.flagstat.txt'.format(p))


def call_variants(paths):
    if os.path.isdir(paths['dir']['call']):
        return
    else:
        os.makedirs(paths['dir']['call'])

    logging.info('Call SNVs/indels')
    sh = Shell(log_txt=os.path.join(paths['dir']['call'], 'command_log.txt'))

    bam_name = 'markdup.bam'
    raw_bcf = os.path.join(paths['dir']['call'], 'raw.bcf')
    flt_vcf = os.path.join(paths['dir']['call'], 'flt.vcf')
    sh.run(
        'samtools mpileup -uf {0} {1} {2} '
        '| bcftools view -bvcg - > {3}'.format(
            paths['ref']['fasta'],
            os.path.join(paths['dir']['map'], 'foreground', bam_name),
            os.path.join(paths['dir']['map'], 'background', bam_name),
            raw_bcf
        )
    )
    sh.run(
        'bcftools view {0} '
        '| vcfutils.pl varFilter -D100 > {1}'.format(
            raw_bcf, flt_vcf
        )
    )


def _cmd_arg_str(config, command):
    return (
        ' '.join(config['command_args'][command])
        if 'command_args' in config and command in config['command_args']
        else ' '
    )
