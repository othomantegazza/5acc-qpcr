#!/usr/bin/env python3


singularity_container = ('shub://TomHarrop/'
                         'singularity-containers:r_3.5.2'
                         '@f4e1feb1a9776e3e255c8fafeb033219')

rule target:
    input:
        'fig/Figure_S4.pdf',
        'fig/Figure_S10.pdf'

rule qpcr_confirms_sampling:
    input:
        chip3 = "data/chip3-normalized.Rds",
        chip4 = "data/chip4-normalized.Rds"
    output:
        f1 = 'fig/Figure_S4.pdf'
    log:
        'logs/fig_S4_qpcr_confirms_sampling.log'
    singularity:
        singularity_container
    script:
        '03-qpcr-confirms-sampling.R'


rule ap2:
    input:
        chip3 = "data/chip3-normalized.Rds",
        chip4 = "data/chip4-normalized.Rds"
    output:
        f1 = "fig/Figure_S10.pdf"
    log:
        'logs/fig_S8_ap2.log'
    singularity:
        singularity_container
    script:
        '02-ap2.R'

rule normalise:
    input:
        chip3 = "data/table-SX-fluidigm-chip3.csv",
        chip4 = "data/table-SX-fluidigm-chip4.csv"
    output:
        chip3 = "data/chip3-normalized.Rds",
        chip4 = "data/chip4-normalized.Rds"
    log:
        'logs/normalise.log'
    singularity:
        singularity_container
    script:
        '01-normalize.R'

