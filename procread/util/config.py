#!/usr/bin/env python

import logging
import os
import re
import shutil
import yaml


def set_log_config(debug=False):
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(name)-15s %(funcName)s'
        ' - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.DEBUG if debug else logging.WARNING
    )


def write_config_yml(path):
    logger = logging.getLogger(__name__)
    if os.path.exists(path):
        print('The file already exists: {}'.format(path))
    else:
        logger.debug('Write {}'.format(path))
        shutil.copyfile(
            os.path.join(os.path.dirname(__file__), '../static/pread.yml'),
            path
        )
        print('A YAML template was generated: {}'.format(path))


def load_param_config(yml_path, work_dir):
    with open(os.path.expanduser(yml_path)) as f:
        yml_dict = yaml.load(f)
    wd = os.path.abspath(os.path.expanduser(work_dir))
    dir_dict = {
        d: os.path.join(wd, d)
        for d in ['input', 'qc', 'trim', 'map', 'call']
    }
    return {
        'yml': yml_dict,
        'paths': {
            'dir': dict([('work', wd)] + list(dir_dict.items())),
            'fastq': [
                dict(
                    [('id', d['id'])] +
                    [
                        (
                            r,
                            os.path.join(
                                dir_dict['input'],
                                re.sub(
                                    r'\.(fastq|fq)\.?[^\.]*\.?[^\.]*$',
                                    '.{}.fastq.gz'.format(
                                        {'read1': 'r1', 'read2': 'r2'}[r]
                                    ),
                                    os.path.basename(d[r])
                                )
                            )
                        )
                        for r in ['read1', 'read2']
                        if r in d
                    ]
                )
                for d in yml_dict['path']['fastq']
            ],
            'ref': {
                f: os.path.join(
                    dir_dict['input'],
                    re.sub(
                        r'\.(fa|fasta)\.?[^\.]*$', e,
                        os.path.basename(
                            yml_dict['path']['fasta']['reference']
                        )
                    )
                )
                for f, e
                in {'fasta': '.ref.fa', 'faidx': '.ref.fa.fai'}.items()
            }
        },
        'cmd_args': {
            c: (
                ' '.join(yml_dict['command_args'][c])
                if 'command_args' in yml_dict and c in yml_dict['command_args']
                else ' '
            )
            for c in ['fastqc', 'cutadapt']
        }
    }
