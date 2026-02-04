import sys, json
from robio import analyze
from plots import gc_plot
from gel import virtual_gel

fasta=sys.argv[1]
job=sys.argv[2]

res=analyze(fasta,job)

gc_plot(res["dna"],f"python-runner/results/{job}_gc.png")

cuts=[]
for k,v in res["enzymes"].items():
    cuts+=v

virtual_gel(cuts,res["length"],f"python-runner/results/{job}_gel.png")

del res["dna"]

with open(f"python-runner/results/{job}.json","w") as f:
    json.dump(res,f,indent=2)
