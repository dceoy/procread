#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from .util import Shell, ProcreadRuntimeError


def prepare_paths(config, cpus):
    logging.info('Prepare paths to directories and files')
    for d in ['work', 'input']:
        os.makedirs(config['paths']['dir'][d], exist_ok=True)

    sh = Shell(log_txt=os.path.join(
        config['paths']['dir']['input'], 'command_log.txt'
    ))
    for c in ['pigz', 'pbzip2']:
        sh.run('{} --version'.format(c))

    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        fq_gz = config['paths']['fastq']['{0}_{1}'.format(t, r)]
        if not os.path.isfile(fq_gz):
            fq_src = os.path.abspath(config['yml']['path']['fastq'][t][r])
            fq_src_ext = os.path.splitext(fq_src)[1]
            if fq_src_ext in ['.fastq', '.fq']:
                sh.run('pigz -p {0} -c {1} > {2}'.format(
                    cpus, fq_src, fq_gz
                ))
            elif fq_src_ext == '.gz':
                os.symlink(fq_src, fq_gz)
            elif fq_src_ext == '.bz2':
                sh.run(
                    'pbzip2 -p# {0} -dc {1} | pigz -p {0} -c - > {2}'.format(
                        cpus, fq_src, fq_gz
                    )
                )
            else:
                raise ProcreadRuntimeError(
                    'Supported extension for input fastq: '
                    '.fastq, .fastq.gz, .fastq.*.gz, .fastq.bz2, .fastq.*.bz2'
                    ', .fq, .fq.gz, .fq.*.gz, .fq.bz2, .fq.*.bz2'
                )

    ref_src = os.path.abspath(config['yml']['path']['fasta']['reference'])
    ref_fa = config['paths']['ref']['fasta']
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


def make_ref_index(config):
    if os.path.isfile(config['paths']['ref']['faidx']):
        return

    logging.info('Make reference index files')
    sh = Shell(
        log_txt=os.path.join(
            config['paths']['dir']['input'], 'command_log.txt'
        ),
        format_log=False
    )
    sh.run(['samtools --version', 'bwa'], check=False)
    sh.run([
        'samtools faidx {}'.format(config['paths']['ref']['fasta']),
        'bwa index -p {0} {1}'.format(
            config['paths']['ref']['faidx'], config['paths']['ref']['fasta']
        )
    ])


def do_qc_checks(config, cpus):
    if os.path.isdir(config['paths']['dir']['qc']):
        return
    else:
        os.makedirs(config['paths']['dir']['qc'])

    logging.info('Do quality control checks for reads')
    sh = Shell(log_txt=os.path.join(
        config['paths']['dir']['qc'], 'command_log.txt'
    ))
    sh.run('fastqc --version')

    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        sh.run(
            'fastqc {0} --threads {1} --outdir {2} {3}'.format(
                config['cmd_args']['fastqc'], cpus,
                config['paths']['dir']['qc'],
                config['paths']['fastq']['{0}_{1}'.format(t, r)]
            )
        )


def trim_adapters(config, cpus):
    if os.path.isdir(config['paths']['dir']['trim']):
        return
    else:
        os.makedirs(config['paths']['dir']['trim'])

    logging.info('Trim adapter sequences in reads')
    sh = Shell(log_txt=os.path.join(
        config['paths']['dir']['trim'], 'command_log.txt'
    ))
    sh.run('cutadapt --version')

    fq_io = {
        t: {
            'in': {
                k: config['paths']['fastq']['{0}_{1}'.format(t, r)]
                for k, r in {'r1': 'read1', 'r2': 'read2'}.items()
            },
            'out': {
                k: os.path.join(
                    config['paths']['dir']['trim'],
                    re.sub(
                        r'(\.r[12]\.fastq\.gz)$', r'\.trimmed\1',
                        os.path.basename(
                            config['paths']['fastq']['{0}_{1}'.format(t, r)]
                        )
                    )
                )
                for k, r in {'r1': 'read1', 'r2': 'read2'}.items()
            }
        }
        for t in ['foreground', 'background']
    }
    cmds = [
        'cutadapt {0} -a {1} -A {2} ''-o {3} -p {4} {5} {6}'.format(
            config['cmd_args']['cutadapt'], config['yml']['adapter']['3prime'],
            config['yml']['adapter']['5prime'], fq_io[t]['out']['r1'],
            fq_io[t]['out']['r2'], fq_io[t]['in']['r1'], fq_io[t]['in']['r2']
        )
        for t in ['foreground', 'background']
    ]
    if cpus > 1:
        sh.run_parallel(cmds)
    else:
        sh.run(cmds)


def map_reads(config, cpus):
    if os.path.isdir(config['paths']['dir']['map']):
        return
    else:
        os.makedirs(config['paths']['dir']['map'])

    logging.info('Map reads to a reference')
    sh = Shell(log_txt=os.path.join(
        config['paths']['dir']['map'], 'command_log.txt'
    ))

    if os.path.isdir(config['paths']['dir']['trim']):
        fq_dict = {
            k: os.path.join(
                config['paths']['dir']['trim'],
                re.sub(
                    r'(\.r[12]\.fastq\.gz)$', r'\.trimmed\1',
                    os.path.basename(v)
                )
            )
            for k, v in config['paths']['fastq'].items()
        }
    else:
        fq_dict = config['paths']['fastq']

    bam_dict = {
        t: {
            b: os.path.join(
                config['paths']['dir']['map'], t, '{}.bam'.format(b)
            )
            for b in ['sort', 'fixmate', 'markdup']
        }
        for t in ['foreground', 'background']
    }

    for t in ['foreground', 'background']:
        os.makedirs(os.path.join(config['paths']['dir']['map'], t))
        bd = bam_dict[t]
        sh.run([
            'bwa mem -t {0} {1} {2} {3} '
            '| samtools view -@ {0} -bS - '
            '| samtools sort -@ {0} -o {4} -'.format(
                cpus, config['paths']['ref']['faidx'], fq_dict[t + '_read1'],
                fq_dict[t + '_read2'], bd['sort']
            ),
            'samtools sort -n -@ {0} {1} '
            '| samtools fixmate -m - {2}'.format(
                cpus, bd['sort'], bd['fixmate']
            ),
            'samtools sort -@ {0} {1} '
            '| samtools markdup - {2}'.format(
                cpus, bd['fixmate'], bd['markdup']
            )
        ] + [
            'samtools index -@ {0} {1}'.format(cpus, bam_dict[t][b])
            for b in ['sort', 'fixmate', 'markdup']
        ])

    sh.run_parallel([
        'samtools flagstat {0} > {1}.flagstat.txt'.format(bam_dict[t][b])
        for t, b
        in product(
            ['foreground', 'background'], ['sort', 'fixmate', 'markdup']
        )
    ])


def call_variants(config):
    if os.path.isdir(config['paths']['dir']['call']):
        return
    else:
        os.makedirs(config['paths']['dir']['call'])

    logging.info('Call SNVs/indels')
    sh = Shell(log_txt=os.path.join(
        config['paths']['dir']['call'], 'command_log.txt'
    ))

    bam_name = 'markdup.bam'
    raw_bcf = os.path.join(config['paths']['dir']['call'], 'raw.bcf')
    flt_vcf = os.path.join(config['paths']['dir']['call'], 'flt.vcf')
    sh.run([
        'samtools mpileup -uf {0} {1} {2} '
        '| bcftools view -bvcg - > {3}'.format(
            config['paths']['ref']['fasta'],
            os.path.join(
                config['paths']['dir']['map'], 'foreground', bam_name
            ),
            os.path.join(
                config['paths']['dir']['map'], 'background', bam_name
            ),
            raw_bcf
        ),
        'bcftools view {0} '
        '| vcfutils.pl varFilter -D100 > {1}'.format(
            raw_bcf, flt_vcf
        )
    ])
