#!/usr/bin/env python

from itertools import product
import logging
import os
import re
import shutil
import subprocess
import yaml


class Shell:
    def __init__(self, log_txt=None, quiet=False, format_log=False,
                 executable='/bin/bash'):
        self.executable = executable
        self.log_txt = log_txt
        self.quiet = quiet
        if log_txt:
            if format_log:
                logging.debug('Write a log file: {}'.format(log_txt))
                with open(log_txt, 'w') as f:
                    f.write('# shell: {0}{1}'.format(executable, os.linesep))
            self.post_proc = (
                ' >> {} 2>&1'.format(log_txt) if quiet
                else ' 2>&1 | tee -a {}'.format(
                    log_txt
                ) + ' && exit ${PIPESTATUS[0]}'
            )
        else:
            self.post_proc = ' > /dev/null 2>&1' if quiet else ''

    def run(self, args, check=True, cwd=None, prompt=None):
        pp = prompt or '[{}] $'.format(os.getcwd())
        for a in (args if isinstance(args, list) else [args]):
            cmd = a + self.post_proc
            logging.debug('shell:{0}{1} {2}'.format(os.linesep, pp, cmd))
            if self.log_txt:
                with open(self.log_txt, 'a') as f:
                    f.write('{0}{1} {2}{0}'.format(os.linesep, pp, cmd))
            subprocess.run(
                cmd, executable=self.executable, stdout=None, stderr=None,
                shell=True, check=check, cwd=cwd
            )

    def run_parallel(self, args, check=True, cwd=None, prompt=None):
        pp = prompt or '[{}] $'.format(os.getcwd())
        if self.log_txt:
            tmp_log_txts = [
                self.log_txt + '.{}'.format(i) for i, a in enumerate(args)
            ]
            cmds = [
                (
                    a + ' >> {} 2>&1'.format(l)
                    if self.quiet else
                    a + ' 2>&1 | tee -a {}'.format(l) +
                    ' && exit ${PIPESTATUS[0]}'
                )
                for a, l in zip(args, tmp_log_txts)
            ]
            for c, l in zip(cmds, tmp_log_txts):
                with open(l, 'w') as f:
                    f.write('{0}{1} {2}{0}'.format(os.linesep, pp, c))
        else:
            cmds = [
                a + (' > /dev/null 2>&1' if self.quiet else '')
                for a in args
            ]
        logging.debug('shell:{0}{1}'.format(
            os.linesep, [(pp + ' ' + c + os.linesep) for c in cmds]
        ))
        procs = [
            subprocess.Popen(
                c, executable=self.executable, stdout=None, stderr=None,
                shell=True, cwd=cwd
            )
            for i, c in enumerate(cmds)
        ]
        try:
            for p, l in zip(procs, tmp_log_txts):
                e = p.communicate()[1]
                with open(self.log_txt, 'a') as lt:
                    with open(l, 'r') as tlt:
                        lt.write(tlt.read())
                os.remove(l)
                if p.returncode != 0:
                    logging.error(e)
                    raise subprocess.CalledProcessError(
                        'Command \'{0}\' returned non-zero exit status '
                        '{1}.'.format(p.args, p.returncode)
                    )
        except ProcreadRuntimeError as err:
            for p, l in zip(procs, tmp_log_txts):
                if not p.returncode:
                    p.kill()
                if os.path.isfile(l):
                    os.remove(l)
            raise err


class ProcreadRuntimeError(Exception):
    pass


def set_log_config(debug=False):
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG if debug else logging.WARNING)


def check_files_exist(paths):
    return all([os.path.isfile(p) for p in paths])


def write_config_yml(path):
    if os.path.exists(path):
        print('The file already exists: {}'.format(path))
    else:
        logging.debug('Write {}'.format(path))
        shutil.copyfile(
            os.path.join(os.path.dirname(__file__), 'pread.yml'), path
        )
        print('A YAML template was generated: {}'.format(path))


def generate_param_config(yml_path, work_dir):
    with open(yml_path) as f:
        yml_dict = yaml.load(f)
    wd = os.path.abspath(work_dir)
    dir_dict = {
        d: os.path.join(wd, d)
        for d in ['input', 'qc', 'trim', 'map', 'call']
    }

    return {
        'yml': yml_dict,
        'paths': {
            'dir': dict([('work', wd)] + list(dir_dict.items())),
            'fastq': [
                {
                    '{0}_{1}'.format(t, r):
                    os.path.join(
                        dir_dict['input'],
                        re.sub(
                            r'\.(fastq|fq)\.?[^\.]*\.?[^\.]*$',
                            '.{}.fastq.gz'.format(
                                {'read1': 'r1', 'read2': 'r2'}[r]
                            ),
                            os.path.basename(d[t][r])
                        )
                    )
                    for t, r
                    in product(
                        ['foreground', 'background'], ['read1', 'read2']
                    )
                }
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
