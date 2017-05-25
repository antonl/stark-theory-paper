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
DEPHASINGPLOTSCRIPT_PATH=str(pathlib.Path('scripts/compare-dephasing-flag.py').absolute())

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
        "python {MKLIN_PATH} template-cfg.yaml; "

rule prepare_complex_dephasing:
    input:
        "simulations/compare-dephasing/{simdir}/Ji-spd.txt",
        "simulations/compare-dephasing/{simdir}/Jii-spd.txt",
        "simulations/compare-dephasing/{simdir}/Jlowfreq-spd.txt",
        "simulations/compare-dephasing/{simdir}/Jnjp-spd.txt",
    wildcard_constraints:
        simdir="[\d\w\-+]+"
    output:
        "simulations/compare-dephasing/{simdir}/template-cfg-complex.yaml",
        "simulations/compare-dephasing/{simdir}/template-cfg-real.yaml",
        "simulations/compare-dephasing/{simdir}/template-cfg-complex-linear.yaml",
        "simulations/compare-dephasing/{simdir}/template-cfg-real-linear.yaml",
        "simulations/compare-dephasing/{simdir}/metacfg.yaml",
    shell:
        "cd simulations/compare-dephasing/{wildcards.simdir};"
        "PYQCFP_TMPDIR={TMPDIR}/simulations/compare-dephasing/{wildcards.simdir} "
        "PYQCFP_SCRATCHDIR={SCRATCHDIR}/simulations/compare-dephasing/{wildcards.simdir} "
        "PYQCFP_OUTPUTDIR={OUTPUTDIR}/simulations/compare-dephasing/{wildcards.simdir} "
        "python make-cfg.py;"
        "python {MKLIN_PATH} template-cfg-real.yaml;"
        "python {MKLIN_PATH} template-cfg-complex.yaml;"

rule run_compare_complex_dephasing:
    input:
        "simulations/compare-dephasing/{simdir}/template-cfg-complex.yaml",
        "simulations/compare-dephasing/{simdir}/template-cfg-real.yaml",
        "simulations/compare-dephasing/{simdir}/metacfg.yaml",
    wildcard_constraints:
        simdir="[\d\w\-+]+"
    output:
        "simulations/compare-dephasing/{simdir}/pump-probe-complex.h5",
        "simulations/compare-dephasing/{simdir}/pump-probe-real.h5",
        "simulations/compare-dephasing/{simdir}/absorption-complex.h5",
        "simulations/compare-dephasing/{simdir}/absorption-real.h5",
    threads: 2
    shell:
        "cd simulations/compare-dephasing/{wildcards.simdir};"
        "python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml template-cfg-complex.yaml;"
        "mv pump-probe.h5 pump-probe-complex.h5;"
        "python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml template-cfg-complex-linear.yaml;"
        "mv absorption.h5 absorption-complex.h5;"
        "python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml template-cfg-real.yaml;"
        "mv pump-probe.h5 pump-probe-real.h5;"
        "python {SIMSCRIPT_PATH} -c {threads} metacfg.yaml template-cfg-real-linear.yaml;"
        "mv absorption.h5 absorption-real.h5;"

rule plot_compare_complex_dephasing:
    wildcard_constraints:
        simdir="[\d\w\-+]+"
    input:
        "simulations/compare-dephasing/{simdir}/pump-probe-complex.h5",
        "simulations/compare-dephasing/{simdir}/pump-probe-real.h5",
        "simulations/compare-dephasing/{simdir}/absorption-complex.h5",
        "simulations/compare-dephasing/{simdir}/absorption-real.h5"
    output:
        "simulations/compare-dephasing/{simdir}/figures/2d-real.png",
        "simulations/compare-dephasing/{simdir}/figures/2d-complex.png",
    threads: 2
    shell:
        "cd simulations/compare-dephasing/{wildcards.simdir};"
        "python {DEPHASINGPLOTSCRIPT_PATH} . -c {threads} --limits 14.25 15.25 --fudge-factor 0.0; "

ruleorder: 
    prepare_complex_dephasing > prep_ddess_sim

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
        "python {PLOTSCRIPT_PATH} . -c {threads} --limits 14.25 15.25 --fudge-factor 0.0; "

rule clean_sim:
    input:
        "{simdir}/Ji-spd.txt"
    output:
        temp("{simdir}/.clean")
    shell:
        "rm -f {wildcards.simdir}/*.{{wrk,txt,yaml,h5}};"
        "rm -rf {wildcards.simdir}/figures;"
        "touch {wildcards.simdir}/.clean"

rule clean_all:
    input:
        "simulations/quick/Ji-monomer-mu/.clean",
        "simulations/quick/Ji-monomer-alpha/.clean",
        "simulations/quick/Ji-dimer-mu/.clean",
        "simulations/quick/Ji-dimer-mu-uncoupled/.clean",
        "simulations/compare-dephasing/Ji-monomer-mu/.clean",


