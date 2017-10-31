#!/usr/bin/env python
"""
Run read-to-variant pipelines for DNA-seq analyses

Usage:
    procread init [--debug] [--file=<yaml>]
    procread run [--debug] [--file=<yaml>]
    procread -h|--help
    procread -v|--version

Options:
    -h, --help      Print help and exit
    -v, --version   Print version and exit
    --debug         Execute a command with debug messages
    --file=<yaml>   Set a path to a YAML for configurations [$PROCREAD_YML]

Commands:
    init            Generate a YAML template for configuration
    run             Run a variant calling pipeline
"""

import logging
import os
from docopt import docopt
from . import __version__
from .util import set_log_config, set_config_yml, write_config_yml, read_yaml


def main():
    args = docopt(__doc__, version='procread {}'.format(__version__))
    set_log_config(debug=args['--debug'])
    logging.debug('args:{0}{1}'.format(os.linesep, args))
    config_yml = set_config_yml(path=args['--file'])

    if args['init']:
        logging.debug('Initiation')
        write_config_yml(path=config_yml)
    elif args['run']:
        logging.debug('config_yml: {}'.format(config_yml))
        config = read_yaml(path=config_yml)
        logging.debug('config:{0}{1}'.format(os.linesep, config))
