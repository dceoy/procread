#!/usr/bin/env python

import logging
import os
import shutil
import subprocess
import yaml


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


def sh_c(args, stdout=None, stderr=None, executable='/bin/bash'):
    logging.debug('sh_c:{0}{1}'.format(os.linesep, args))
    return subprocess.run(
        args, stdout=stdout, stderr=stderr, shell=True, check=True,
        executable=executable
    )
