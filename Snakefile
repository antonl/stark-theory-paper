import yaml

with open('config.yaml') as f:
    cfg = yaml.load(f)

TMPDIR=cfg['tmpdir']
SCRATCHDIR=cfg['scratchdir']
OUTPUTDIR=cfg['outputdir']

rule make_spectral_densities:
    input:
        "spectral-densities/oscillator-modes.csv"
    output:
        "spectral-densities/Ji-spd.txt",
        "spectral-densities/Jii-spd.txt",
        "spectral-densities/Jlowfreq-spd.txt",
        "spectral-densities/Jnjp-spd.txt",
        "spectral-densities/sims-spectral-densities.png"
    shell:
        "cd spectral-densities;"
        "python make-spectral-densities.py"

rule copy_spd:
    input:
        "spectral-densities/{spdname}.txt"
    output:
        "{simdir}/{spdname}.txt"
    shell:
        "cp {input} {output}"

rule prep_ddess_sim:
    input:
        "{simdir}/Ji-spd.txt",
        "{simdir}/Jii-spd.txt",
        "{simdir}/Jlowfreq-spd.txt",
        "{simdir}/Jnjp-spd.txt",
    output:
        "{simdir}/template-cfg.yaml",
        "{simdir}/metacfg.yaml",
        touch("{simdir}/.made-2dess")
    shell:
        "cd {wildcards.simdir};"
        "PYQCFP_TMPDIR={TMPDIR}/{wildcards.simdir} "
        "PYQCFP_SCRATCHDIR={SCRATCHDIR}/{wildcards.simdir} "
        "PYQCFP_OUTPUTDIR={OUTPUTDIR}/{wildcards.simdir} "
        "python make-cfg.py"

rule run_2dess_sim:
    input:
        "{simdir}/template-cfg.yaml",
        "{simdir}/metacfg.yaml",
        "{simdir}/.made-2dess"
    output:
        "{simdir}/pump-probe.h5"
    shell:
        "cd {wildcards.simdir}; "
        "python simulation-meta.py; "
        "rm .made-2dess"

rule prep_linear_sim:
    input:
        "{simdir}/template-cfg.yaml",
        "{simdir}/metacfg.yaml",
        "{simdir}/.made-2dess"
    output:
        touch("{simdir}/.made-linear")
    shell:
        "cd {wildcards.simdir}; "
        "python make-cfg-linear.py; "
        "rm .made-2dess"

rule run_linear_sim:
    input:
        "{simdir}/template-cfg.yaml",
        "{simdir}/metacfg.yaml",
        "{simdir}/.made-linear"
    output:
        "{simdir}/absorption.h5",
    shell:
        "cd {wildcards.simdir}; "
        "python simulation-meta.py; "
        "rm .made-linear"
