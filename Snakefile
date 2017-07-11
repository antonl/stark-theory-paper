import yaml
import pathlib

with open('config.yaml') as f:
    cfg = yaml.load(f)

TMPDIR=cfg['tmpdir']
SCRATCHDIR=cfg['scratchdir']
OUTPUTDIR=cfg['outputdir']
THREADS=cfg['threads']
VD_COUNT=cfg['voltage-dependence-count']
VD_RANGE_MIN, VD_RANGE_MAX=cfg['voltage-dependence-range']

SIMSCRIPT_PATH=str(pathlib.Path('scripts/simulation-meta.py').absolute())
MKLIN_PATH=str(pathlib.Path('scripts/make-cfg-linear.py').absolute())
PLOTSCRIPT_PATH=str(pathlib.Path('scripts/make-plots.py').absolute())
DEPHASINGPLOTSCRIPT_PATH=str(pathlib.Path('scripts/compare-dephasing-flag.py').absolute())
MKANALYTICPLOTSCRIPT_PATH=str(pathlib.Path('scripts/make-analytic-plot.py').absolute())
MKVOLT_PATH=str(pathlib.Path('scripts/make-cfg-voltage-dependence.py').absolute())
MKDDESSPLOT_PATH=str(pathlib.Path('scripts/make-ddess-plot.py').absolute())
MKLINEARPLOT_PATH=str(pathlib.Path('scripts/make-linear-plot.py').absolute())
TTTSIM_PATH=str(pathlib.Path('scripts/simulation-ttt.py').absolute())

PLOT_FILES_DDESS = ['eigen-energies.info',
                    '2d-reference.png',
                    '2d-fieldon.png',
                    '2d-fieldoff.png',
                    '2d-stark.png',]

PLOT_FILES_LINEAR = ['linear-reference.png',
                     'linear-fieldoff.png',
                     'linear-fieldon.png',
                     'linear-stark.png',]
                     #'linear-projections.png',
                     #'linear-stark-projections.png']

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
        "{simdir}/V-spd.txt",
        "{simdir}/Jso-spd.txt",
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
        "{simdir}/{cfgname}.yaml",
        "{simdir}/metacfg.yaml",
    output:
        "{simdir}/{cfgname}-linear.yaml"
    shell:
        "cd {wildcards.simdir}; "
        "python {MKLIN_PATH} {wildcards.cfgname}.yaml; "

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
        "simulations/compare-dephasing/{simdir}/figures/linear-real.png",
        "simulations/compare-dephasing/{simdir}/figures/linear-complex.png",
    threads: 2
    shell:
        "cd simulations/compare-dephasing/{wildcards.simdir};"
        "python {DEPHASINGPLOTSCRIPT_PATH} . -c {threads} --limits 14.25 15.25 --fudge-factor 6.6; "

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

rule prepare_voltage_dependence_ddess:
    input:
        "simulations/voltage-dependence/{simdir}/template-cfg.yaml",
        "simulations/voltage-dependence/{simdir}/metacfg.yaml",
    wildcard_constraints:
        simdir="[\d\w\-+]+"
    output:
        "simulations/voltage-dependence/{simdir}/voltagecfg.yaml",
        expand("simulations/voltage-dependence/{{simdir}}/template-cfg-{field_id:03d}.yaml",
            field_id=range(VD_COUNT))
    shell:
        "cd simulations/voltage-dependence/{wildcards.simdir}; "
        "python {MKVOLT_PATH} --range {VD_RANGE_MIN} {VD_RANGE_MAX} --count {VD_COUNT} template-cfg.yaml; "

rule run_voltage_dependence_ddess:
    input:
        "simulations/voltage-dependence/{simdir}/metacfg.yaml",
        "simulations/voltage-dependence/{simdir}/template-cfg-{field_id}.yaml",
    output:
        "simulations/voltage-dependence/{simdir}/pump-probe-{field_id}.h5",
    threads: THREADS
    shell:
        "python {SIMSCRIPT_PATH} -c {threads} {input[0]} {input[1]};"
        "mv pump-probe.h5 {output[0]};"

rule run_voltage_dependence_linear:
    input:
        "simulations/voltage-dependence/{simdir}/metacfg.yaml",
        "simulations/voltage-dependence/{simdir}/template-cfg-{field_id}-linear.yaml",
    output:
        "simulations/voltage-dependence/{simdir}/absorption-{field_id}.h5",
    threads: THREADS
    shell:
        "python {SIMSCRIPT_PATH} -c {threads} {input[0]} {input[1]};"
        "mv absorption.h5 {output[0]};"

rule plot_voltage_dependence_ddess:
    input:
        "simulations/voltage-dependence/{simdir}/pump-probe-{field_id}.h5",
    output:
        expand("simulations/voltage-dependence/{{simdir}}/figures/pump-probe-{{field_id}}/{file}",
            file=PLOT_FILES_DDESS)
    threads: 6
    shell:
        #"cd simulations/voltage-dependence/{wildcards.simdir};"
        "python {MKDDESSPLOT_PATH} {input} -c {threads} --limits 14.25 15.25 --fudge-factor 0.0; "

rule plot_voltage_dependence_linear:
    input:
        "simulations/voltage-dependence/{simdir}/absorption-{field_id}.h5",
    output:
        expand("simulations/voltage-dependence/{{simdir}}/figures/absorption-{{field_id}}/{file}",
            file=PLOT_FILES_LINEAR)
    threads: 6
    shell:
        #"cd simulations/voltage-dependence/{wildcards.simdir};"
        "python {MKLINEARPLOT_PATH} {input} -c {threads} --limits 14.25 15.25 --fudge-factor 0.0; "

rule collect_voltage_dependence_ddess_dir:
    input:
        expand("simulations/voltage-dependence/{{simdir}}/pump-probe-{field_id:03d}.h5",
            field_id=range(VD_COUNT))
    output:
        temp("simulations/voltage-dependence/{simdir}/.ddess-voltage-dep-done")
    shell:
        "touch simulations/voltage-dependence/{wildcards.simdir}/.ddess-voltage-dep-done"

rule collect_voltage_dependence_linear_dir:
    input:
        expand("simulations/voltage-dependence/{{simdir}}/absorption-{field_id:03d}.h5",
            field_id=range(VD_COUNT))
    output:
        temp("simulations/voltage-dependence/{simdir}/.linear-voltage-dep-done")
    shell:
        "touch simulations/voltage-dependence/{wildcards.simdir}/.linear-voltage-dep-done"

rule plot_voltage_dependence_ddess_dir:
    input:
        expand("simulations/voltage-dependence/{{simdir}}/figures/pump-probe-{field_id:03d}/{file}",
            field_id=range(VD_COUNT), file=PLOT_FILES_DDESS)
    output:
        temp("simulations/voltage-dependence/{simdir}/.plot-ddess-voltage-dep-done")
    shell:
        "touch simulations/voltage-dependence/{wildcards.simdir}/.plot-ddess-voltage-dep-done"

rule plot_voltage_dependence_linear_dir:
    input:
        expand("simulations/voltage-dependence/{{simdir}}/figures/absorption-{field_id:03d}/{file}",
            field_id=range(VD_COUNT), file=PLOT_FILES_LINEAR)
    output:
        temp("simulations/voltage-dependence/{simdir}/.plot-linear-voltage-dep-done")
    shell:
        "touch simulations/voltage-dependence/{wildcards.simdir}/.plot-linear-voltage-dep-done"

rule make_voltage_dependence_ddess_video:
    input:
        expand("simulations/voltage-dependence/{{simdir}}/figures/pump-probe-{field_id:03d}/{{plotfile}}.png",
            field_id=range(VD_COUNT))
    output:
        "simulations/voltage-dependence/{simdir}/figures/{plotfile}.mp4"
    shell:
        "cd simulations/voltage-dependence/{wildcards.simdir}/figures; "
        "ffmpeg -framerate 1 -i pump-probe-%03d/{wildcards.plotfile}.png -c:v libx264 -r 30 "
        "-pix_fmt yuv420p {wildcards.plotfile}.mp4"

rule make_voltage_dependence_linear_video:
    input:
        expand("simulations/voltage-dependence/{{simdir}}/figures/absorption-{field_id:03d}/{{plotfile}}.png",
            field_id=range(VD_COUNT))
    output:
        "simulations/voltage-dependence/{simdir}/figures/{plotfile}.mp4"
    shell:
        "cd simulations/voltage-dependence/{wildcards.simdir}/figures; "
        "ffmpeg -framerate 1 -i absorption-%03d/{wildcards.plotfile}.png -c:v libx264 -r 30 "
        "-pix_fmt yuv420p {wildcards.plotfile}.mp4"

rule run_ttt_sim:
    input:
        "simulations/ttt-sim/{simdir}/template-cfg.yaml",
    output:
        expand("simulations/ttt-sim/{{simdir}}/{files}.{exts}",
            files=['rephasing', 'nonrephasing'], exts=['inp', 'outp', 'text']),
        expand("simulations/ttt-sim/{{simdir}}/{files}-{wn}.png",
            files=['absorptive', 'rephasing', 'nonrephasing', 'rephasing-windowed',
            'nonrephasing-windowed'], wn=['rect', 'gauss'])
    run:
        files = expand("simulations/ttt-sim/{simdir}/{files}.png",
                files=['absorptive', 'rephasing', 'nonrephasing', 'rephasing-windowed',
                'nonrephasing-windowed'], simdir=wildcards.simdir)

        shell("cd simulations/ttt-sim/{wildcards.simdir}; "
        "python {TTTSIM_PATH} template-cfg.yaml --limits 14.25 15.25; ")
        #"python {TTTSIM_PATH} template-cfg.yaml;")

        for p in [pathlib.Path(n) for n in files]:
            s = p.with_suffix('')
            s = s.with_name(s.name + '-rect').with_suffix('.png')
            shell("mv {inp} {outp}", inp=str(p), outp=str(s))

        shell("cd simulations/ttt-sim/{wildcards.simdir}; "
        "python {TTTSIM_PATH} template-cfg.yaml --window --limits 14.25 15.25;")
        #"python {TTTSIM_PATH} template-cfg.yaml --window;")

        for p in [pathlib.Path(n) for n in files]:
            s = p.with_suffix('')
            s = s.with_name(s.name + '-gauss').with_suffix('.png')
            shell("mv {inp} {outp}", inp=str(p), outp=str(s))

rule run_all_ttt:
    input:
        expand("simulations/ttt-sim/{simdir}/rephasing.outp",
            simdir=[
            'Ji-monomer-mu', 'Ji-dimer-mu', 'Ji-dimer-ct-mu',
            'Jii-monomer-mu', 'Jii-dimer-mu', 'Jii-dimer-ct-mu'])

rule plot_grid_size:
    input:
        "simulations/scan-grid-size/{simdir}/pump-probe.h5"
    output:
        "simulations/scan-grid-size/{simdir}/figures/2d-reference.png"
    shell:
        "cd simulations/scan-grid-size/{wildcards.simdir};"
        "python {MKANALYTICPLOTSCRIPT_PATH} pump-probe.h5 --limits 14.25 15.25"

rule plot_all_grid_size:
    input:
        expand("simulations/scan-grid-size/{simdir}-{size}/figures/2d-reference.png",
            simdir=['Ji-monomer-mu', 'Ji-dimer-mu', 'Ji-dimer-ct-mu'],
            size=['small', 'medium', 'large'])

rule plot_all_complex_dephasing:
    input:
        expand("simulations/compare-dephasing/{simdir}/figures/2d-real.png",
            simdir=['Ji-monomer', 'Ji-dimer'])

rule plot_all_quick:
    input:
        "simulations/quick/Ji-monomer-mu/figures/2d-reference.png",
        "simulations/quick/Ji-monomer-alpha/figures/2d-reference.png",
        "simulations/quick/Ji-dimer-mu/figures/2d-reference.png",
        "simulations/quick/Ji-dimer-alpha/figures/2d-reference.png",
        "simulations/quick/Ji-dimer-mu-uncoupled/figures/2d-reference.png",
        "simulations/quick/Ji-dimer-ct-mu/figures/2d-reference.png",
        "simulations/quick/Ji-dimer-ct-alpha/figures/2d-reference.png",

rule plot_all_large_mesh:
    input:
        "simulations/large-mesh/Ji-monomer-mu/figures/2d-reference.png",
        "simulations/large-mesh/Ji-monomer-alpha/figures/2d-reference.png",
        "simulations/large-mesh/Ji-dimer-mu/figures/2d-reference.png",
        "simulations/large-mesh/Ji-dimer-alpha/figures/2d-reference.png",
        "simulations/large-mesh/Ji-dimer-mu-uncoupled/figures/2d-reference.png",
        "simulations/large-mesh/Ji-dimer-ct-mu/figures/2d-reference.png",
        "simulations/large-mesh/Ji-dimer-ct-alpha/figures/2d-reference.png",

rule plot_all_voltage_dependence:
    input:
        "simulations/voltage-dependence/Ji-monomer-mu/figures/2d-stark.mp4",
        "simulations/voltage-dependence/Ji-monomer-mu/figures/2d-fieldon.mp4",
        "simulations/voltage-dependence/Ji-monomer-mu/figures/linear-stark.mp4",
        "simulations/voltage-dependence/Ji-monomer-mu/figures/linear-fieldon.mp4",
        "simulations/voltage-dependence/Ji-dimer-mu/figures/2d-stark.mp4",
        "simulations/voltage-dependence/Ji-dimer-mu/figures/2d-fieldon.mp4",
        "simulations/voltage-dependence/Ji-dimer-mu/figures/linear-stark.mp4",
        "simulations/voltage-dependence/Ji-dimer-mu/figures/linear-fieldon.mp4",
        "simulations/voltage-dependence/Ji-dimer-ct-mu/figures/2d-stark.mp4",
        "simulations/voltage-dependence/Ji-dimer-ct-mu/figures/2d-fieldon.mp4",
        "simulations/voltage-dependence/Ji-dimer-ct-mu/figures/linear-stark.mp4",
        "simulations/voltage-dependence/Ji-dimer-ct-mu/figures/linear-fieldon.mp4",

rule clean_sim:
    input:
        "{simdir}/Ji-spd.txt"
    output:
        temp("{simdir}/.clean")
    shell:
        "rm -f {wildcards.simdir}/*.{{wrk,txt,yaml,h5,inp,outp,png,text}};"
        "rm -rf {wildcards.simdir}/figures;"
        "touch {wildcards.simdir}/.clean"

rule clean_voltage_dependence:
    input:
        "simulations/voltage-dependence/Ji-monomer-mu/.clean",
        "simulations/voltage-dependence/Ji-dimer-mu/.clean",
        "simulations/voltage-dependence/Ji-dimer-ct-mu/.clean",

rule clean_ttt:
    input:
        expand("simulations/ttt-sim/{simdir}/.clean",
            simdir=['Ji-monomer-mu', 'Ji-dimer-mu', 'Ji-dimer-ct-mu',
                    'Jii-monomer-mu', 'Jii-dimer-mu', 'Jii-dimer-ct-mu'])

rule clean_scan_grid:
    input:
        expand("simulations/scan-grid-size/{simdir}-{size}/.clean",
            simdir=['Ji-monomer-mu', 'Ji-dimer-mu', 'Ji-dimer-ct-mu'],
            size=['small', 'medium', 'large'])

rule clean_quick:
    input:
        expand("simulations/quick/{simdir}/.clean",
            simdir=[
                'Ji-monomer-mu', 
                'Ji-monomer-alpha',
                'Ji-dimer-mu', 
                'Ji-dimer-mu-uncoupled', 
                'Ji-dimer-alpha', 
                'Ji-dimer-ct-mu',
                'Ji-dimer-ct-alpha',
                'Ji-dimer-ct-mu']),

rule clean_compare_dephasing:
    input:
        expand("simulations/compare-dephasing/{simdir}/.clean",
            simdir=['Ji-monomer', 'Ji-dimer',])

rule clean_large_mesh:
    input:
        expand("simulations/large-mesh/{simdir}/.clean",
            simdir=[
                'Ji-monomer-mu', 
                'Ji-monomer-alpha',
                'Ji-dimer-mu', 
                'Ji-dimer-mu-uncoupled', 
                'Ji-dimer-alpha', 
                'Ji-dimer-ct-mu',
                'Ji-dimer-ct-alpha',
                'Ji-dimer-ct-mu']),

rule clean_all:
    input:
        # quick simulations
        rules.clean_quick.input,
        # compare no complex dephasing flag
        rules.clean_compare_dephasing.input,
        # large mesh simulations
        rules.clean_large_mesh.input,
        # time-domain sims
        rules.clean_ttt.input,
        # grid scan
        rules.clean_scan_grid.input,
        # voltage dependence
        rules.clean_voltage_dependence.input,

rule run_all:
    input:
        rules.plot_all_quick.input,
        rules.plot_all_large_mesh.input,
        rules.plot_all_voltage_dependence.input,
        rules.plot_all_complex_dephasing.input,
        rules.run_all_ttt.input,

ruleorder:
    prepare_complex_dephasing > prep_ddess_sim

ruleorder:
    plot_grid_size > plot_sim_results
