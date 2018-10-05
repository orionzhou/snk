#!python

#!/usr/bin/env python3
import os
import sys

from snakemake.utils import read_job_properties

jobscript = sys.argv[1]
c = read_job_properties(jobscript)
cc = c['cluster']
cp = c['params']

os.system("qsub -q {q} 
    -l nodes={nodes}:ppn={ppn},walltime={runtime}:00:00,mem={mem}gb \
    -N {N} -M {M} -m {m} -r {r} -o {o} -e {e} {script}".format(
        q = cc['q'],
        nodes = cc['nodes'],
        ppn = c['threads'],
        runtime = cc['runtime'],
        mem = cc['mem'],
        N = cc['N'],
        M = cc['M'],
        m = cc['m'],
        r = cc['r'],
        o = cc['o'],
        e = cc['e'],
        script = jobscript))
