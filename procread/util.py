#!/usr/bin/env python

import logging
import os
import shutil
import subprocess
import yaml


class Shell:
    def __init__(self, log_txt=None, quiet=False, executable='/bin/bash'):
        self.executable = executable
        self.log_txt = log_txt
        self.quiet = quiet
        if log_txt:
            logging.debug('Write a log file: {}'.format(log_txt))
            with open(log_txt, 'w') as f:
                f.write('# shell: {0}{1}'.format(executable, os.linesep))
            self.post_proc = (
                ' >> {} 2>&1'.format(log_txt) if quiet
                else ' 2>&1 | tee -a {}'.format(log_txt)
            )
        else:
            self.post_proc = ' > /dev/null 2>&1' if quiet else ''

    def run(self, arg_str, prompt=None):
        cmd = arg_str + self.post_proc
        pr = prompt or '[{}] $'.format(os.getcwd())
        logging.debug('shell:{0}{1} {2}'.format(os.linesep, pr, cmd))
        if self.log_txt:
            with open(self.log_txt, 'a') as f:
                f.write('{0}{1} {2}{0}'.format(os.linesep, pr, cmd))
        return subprocess.run(
            cmd, executable=self.executable, stdout=None, stderr=None,
            shell=True, check=True
        )


class ProcreadError(Exception):
    pass


def set_log_config(debug=False):
    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG if debug else logging.WARNING)


def read_yaml(path):
    with open(path) as f:
        d = yaml.load(f)
    return d


def dump_yaml(dict, flow=False):
    return yaml.dump(dict, default_flow_style=flow)


def write_config_yml(path):
    if os.path.exists(path):
        print('The file already exists: {}'.format(path))
    else:
        logging.debug('Write {}'.format(path))
        shutil.copyfile(
            os.path.join(os.path.dirname(__file__), 'pread.yml'), path
        )
        print('A YAML template was generated: {}'.format(path))
