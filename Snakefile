import yaml
import pathlib

with open('config.yaml') as f:
    cfg = yaml.load(f)

TMPDIR=cfg['tmpdir']
SCRATCHDIR=cfg['scratchdir']
OUTPUTDIR=cfg['outputdir']
THREADS=cfg['threads']

SIMSCRIPT_PATH=str(pathlib.Path('scripts/simulation-meta.py').absolute())
MKLIN_PATH=str(pathlib.Path('scripts/make-cfg-linear.py').absolute())
PLOTSCRIPT_PATH=str(pathlib.Path('scripts/make-plots.py').absolute())

rule make_spectral_densities:
    input:
        "rawdata/oscillator-modes.csv"
    output:
        "spectral-densities/Ji-spd.txt",
        "spectral-densities/Jii-spd.txt",
        "spectral-densities/Jlowfreq-spd.txt",
        "spectral-densities/Jnjp-spd.txt",
        "spectral-densities/sims-spectral-densities.png"
    shell:
        "cd spectral-densities;"
        "python ../scripts/make-spectral-densities.py"

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
    output:
        "{simdir}/pump-probe.h5",
    threads: THREADS
    shell:
        "cd {wildcards.simdir}; "
        "python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml template-cfg.yaml; "

rule prep_linear_sim:
    input:
        "{simdir}/template-cfg.yaml",
        "{simdir}/metacfg.yaml",
    output:
        "{simdir}/template-cfg-linear.yaml"
    shell:
        "cd {wildcards.simdir}; "
        "python {MKLIN_PATH}; "

rule prepare_complex_dephasing:
    input:
        "simulation/compare-dephasing/{simdir}/Ji-spd.txt",
        "simulation/compare-dephasing/{simdir}/Jii-spd.txt",
        "simulation/compare-dephasing/{simdir}/Jlowfreq-spd.txt",
        "simulation/compare-dephasing/{simdir}/Jnjp-spd.txt",
    output:
        "simulation/compare-dephasing/{simdir}/template-cfg-complex.yaml",
        "simulation/compare-dephasing/{simdir}/template-cfg-real.yaml"
        "simulation/compare-dephasing/{simdir}/metacfg.yaml",
    run:
        path = "simulation/compare-dephasing/{simdir}"
        shell("cd {path}", path=wildcards.path)
        shell("python make-cfg.py")

rule run_compare_complex_dephasing:
    input:
        "simulation/compare-dephasing/{simdir}/template-cfg-complex.yaml",
        "simulation/compare-dephasing/{simdir}/template-cfg-real.yaml"
        "simulation/compare-dephasing/{simdir}/metacfg.yaml",
    output:
        "simulation/compare-dephasing/{simdir}/pump-probe-complex.h5"
        "simulation/compare-dephasing/{simdir}/pump-probe-real.h5"
    threads: 2
    run:
        path = "simulation/compare-dephasing/{simdir}"
        shell("cd {path}", path=wildcards.path)
        shell("python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml "
              "template-cfg-complex.yaml;"
              "mv pump-probe.h5 pump-probe-complex.h5;")
        shell("python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml "
              "template-cfg-real.yaml;"
              "mv pump-probe.h5 pump-probe-real.h5;")

rule run_linear_sim:
    input:
        "{simdir}/template-cfg-linear.yaml",
        "{simdir}/metacfg.yaml",
    output:
        "{simdir}/absorption.h5",
    threads: THREADS
    shell:
        "cd {wildcards.simdir}; "
        "python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml template-cfg-linear.yaml; "

rule plot_sim_results:
    input:
        "{simdir}/absorption.h5",
        "{simdir}/pump-probe.h5"
    output:
        expand("{{simdir}}/figures/{name}",
            name=['eigen-energies.info',
                  '2d-reference.png',
                  '2d-fieldon.png',
                  '2d-fieldoff.png',
                  '2d-stark.png',
                  'linear-reference.png',
                  'linear-fieldoff.png',
                  'linear-fieldon.png',
                  'linear-stark.png',
                  'linear-projections.png',
                  'linear-stark-projections.png'])
    threads: 6
    shell:
        "cd {wildcards.simdir}; "
        "python {PLOTSCRIPT_PATH} . -c {threads} --limits 14.25 15.25 --fudge-factor 6.6; "
