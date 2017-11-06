procread
========

Processing tool for read-to-variant workflows in DNA-seq analyses

[![wercker status](https://app.wercker.com/status/7fb703c37b40bee7d216d7d86f8eb7f2/m/master "wercker status")](https://app.wercker.com/project/byKey/7fb703c37b40bee7d216d7d86f8eb7f2)

Dependencies:
  - FastQC
  - Cutadapt
  - SAMtools
  - BWA
  - BCFtools

Docker image
------------

```sh
# build using docker
$ docker image build -t procread .

# build using docker-compose
$ docker-compose build
```

Usage
-----

```sh
$ procread --help
Run read-to-variant pipelines for DNA-seq analyses

Usage:
    procread init [--debug] [--file=<yaml>] [--work=<dir>]
    procread qc [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread trim [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread map [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread call [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread run [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread -h|--help
    procread -v|--version

Options:
    -h, --help      Print help and exit
    -v, --version   Print version and exit
    --debug         Execute a command with debug messages
    --file=<yaml>   Set a path to a configuration YAML [default: ./pread.yml]
    --cpus=<int>    Limit CPU cores for use
    --work=<dir>    Set a working directory [default: .]

Commands:
    init            Generate a YAML template for configuration
    qc              Do quality control checks
    trim            Trim adapter sequences in reads
    map             Map reads to a reference
    call            Call SNVs/indels
    run             Run a variant calling pipeline
```
