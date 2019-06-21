#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess

from snakemake.utils import read_job_properties, makedirs

parser=argparse.ArgumentParser(add_help=False)
parser.add_argument("--depend", help="Space separated list of ids for jobs this job should depend on.")
parser.add_argument("-a", help="Declare the time when the job becomes eligible for execution.")
parser.add_argument("-A", help="Define the account string.")
parser.add_argument("-b", help="PBS Server timeout.")
parser.add_argument("-c", help="Checkpoint options.")
parser.add_argument("-C", help="Directive prefix in script file.")
parser.add_argument("-d", help="Working directory to be used (default: ~). PBS_O_INITDIR")
parser.add_argument("-D", help="Root directory to be used. PBS_O_ROOTDIR")
parser.add_argument("-e", help="standard error path.")
parser.add_argument("-f", help="Fault tolerant.",action="store_true")
parser.add_argument("-h", help="Apply user hold at submission time",action="store_true")
parser.add_argument("-j", help="Merge standard error and standard out. (oe or eo)")
parser.add_argument("-l", help="Resource list.")
parser.add_argument("-m", help="Mail options.")
parser.add_argument("-M", help="Mail users.")
parser.add_argument("-N", help="Name for the job.")
parser.add_argument("-o", help="standard output path.")
parser.add_argument("-p", help="Set job priority.")
parser.add_argument("-P", help="Proxy user for job.")
parser.add_argument("-q", help="Set destination queue.")
parser.add_argument("-t", help="Array request.")
parser.add_argument("-u", help="Set user name for job.")
parser.add_argument("-v", help="Environment variables to export to the job.")
parser.add_argument("-V", help="Export all environment variables.",action="store_true")
parser.add_argument("-w", help="Set working directory. PBS_O_WORKDIR")
parser.add_argument("-W", help="Additional attributes.")
parser.add_argument("--help", help="Display help message.",action="store_true")

parser.add_argument("positional",action="append",nargs="?")
args = parser.parse_args()

if args.help :
    parser.print_help()
    sys.exit(0)

jobscript = sys.argv[-1]

job_properties = read_job_properties(jobscript)
if job_properties.get("params", {}).get("j") == None:
    print("no job parameters params[j]", job_properties)
    sys.exit(2)

atime, acc_string, pbs_time, chkpt, pref, \
dd, rd, se, ft, hold, \
j, resource, mail, mailuser, jname, \
so, priority, proxy, q, ar, \
user, ev, eall, wd, add, \
depend, resourceparams, extras = [''] * 28

nodes, ppn, mem, walltime, jobe, jobo = [''] * 6

cc = job_properties['cluster']
if 'q' in cc: q = " -q " + cc['q']
if 'm' in cc: mail = " -m " + cc['m']
if 'M' in cc: mailuser = " -M " + cc['M']
if 'N' in cc: jname = " -N " + cc['N']
if 'o' in cc:
    so = " -o " + cc['o']
    jobo = cc['o']
if 'e' in cc:
    se = " -e " + cc['e']
    jobe = cc['e']
#if 'r' in cc: r = " -r " + cc['r']
if 'nodes' in cc: nodes = "nodes=" + str(cc["nodes"])
if 'ppn' in cc: ppn = "ppn=" + str(cc["ppn"])
if 'mem' in cc: mem = "mem=" + str(cc["mem"]) + 'gb'
if 'runtime' in cc: walltime = "walltime=" + str(cc["runtime"]*3600)

if args.depend:
    for m in set(args.depend.split(" ")):
        depend = depend + ":" + m
if depend:
    depend = " -W \"depend=afterok" + depend + "\""

if args.positional:
    for m in args.positional:
        extras = extras + " " + m

if args.a: atime = " -a " + args.a
if args.A: acc_string = " -A " + args.A
if args.b: pbs_time = " -b " + args.b
if args.c: chkpt = " -c " + args.c
if args.C: pref = " -C " + args.C
if args.d: dd = " -d " + args.d
if args.D: rd = " -D " + args.D
if args.e: se = " -e " + args.e
if args.f: ft = " -f"
if args.h: hold = " -h"
if args.j: j = " -j " + args.j
if args.l: resource = " -l " + args.l
if args.m: mail = " -m " + args.m
if args.M: mailuser = " -M " + args.M
if args.N: jname = " -N " + args.N
if args.o: so = " -o " + args.o
if args.p: priority = " -p " + args.p
if args.P: proxy = " -P " + args.P
if args.q: q = " -q " + str(qmap[args.q])
if args.t: ar = " -t " + args.ar
if args.u: user = " -u " + args.u
if args.v: ev = " -v " + args.v
if args.V: eall = " -V"
if args.w: wd = " -w " + args.w
if args.W: add= " -W \"" + args.W + "\""


params = job_properties["params"]
if 'N' in params: jname = " -N " + str(params['N'])
if 'o' in params:
    so = " -o " + str(params['o'])
    jobo = params['o']
if 'e' in params:
    se = " -e " + str(params['e'])
    jobe = params['e']

retry = 0
if "resources" in job_properties:
    retry = job_properties["resources"]['attempt'] - 1

params = params['j']
if 'q' in params: q = " -q %s" % params['q']
if 'm' in params: mail = " -m %s" % str(params['m'])
if 'M' in params: mailuser = " -M %s" % str(params['M'])
if "nodes" in params: nodes="nodes=%d" % params["nodes"]
if 'ppn' in params: ppn = "ppn=%d" % (params["ppn"] + retry * params['appn'])
if ppn and not nodes : nodes="nodes=1"
if "mem" in params: mem="mem=%dgb" % (params["mem"] + retry * params['amem'])
if "runtime" in params: walltime="walltime=%d:00:00" % (params["runtime"] + retry * params['aruntime'])

print('  '.join((jname, jobo, jobe)), file=sys.stderr)
print("  ".join((q, ppn, mem, walltime)), file=sys.stderr)
#sys.exit(2)

for jdir in set([os.path.dirname(p) for p in [jobe, jobo]]):
    if not os.path.isdir(jdir):
        makedirs(jdir)

if nodes or ppn or mem or walltime: resourceparams = " -l \""
if nodes: resourceparams = resourceparams + nodes
if nodes and ppn: resourceparams = resourceparams + ":" + ppn
if nodes and mem: resourceparams = resourceparams + ","
if mem: resourceparams = resourceparams + mem
if walltime and (nodes or mem): resourceparams = resourceparams + ","
if walltime: resourceparams = resourceparams + walltime
if nodes or mem or walltime: resourceparams = resourceparams + "\""

cmd = "qsub {a}{A}{b}{c}{C}{d}{D}{e}{f}{h}{j}{l}{m}{M}{N}{o}{p}{P}{q}{t}{u}{v}{V}{w}{W}{rp}{dep}{ex}".format(\
    a=atime,A=acc_string,b=pbs_time,c=chkpt,C=pref,d=dd,D=rd,e=se,f=ft,h=hold,j=j,l=resource,m=mail,M=mailuser,\
    N=jname,o=so,p=priority,P=proxy,q=q,t=ar,u=user,v=ev,V=eall,w=wd,W=add,rp=resourceparams,dep=depend,ex=extras)

try:
    res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
except subprocess.CalledProcessError as e:
    raise e

res = res.stdout.decode()
print(res.strip())
