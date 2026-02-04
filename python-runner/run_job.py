import sys,base64,json,os
from robio_engine import analyze_fasta
from plots import gc_plot
from gel import virtual_gel

fasta=base64.b64decode(sys.argv[1]).decode()
job=sys.argv[2]

result=analyze_fasta(fasta)

os.makedirs("python-runner/results",exist_ok=True)
os.makedirs("python-runner/images",exist_ok=True)

dna=fasta.splitlines()[1].replace("U","T")

gc_plot(dna,f"python-runner/images/{job}_gc.png")

cuts=[]
for s in result["restriction_sites"]:
    cuts+=s["positions"]

virtual_gel(len(dna),cuts,f"python-runner/images/{job}_gel.png")

result["gc_plot"]=f"{job}_gc.png"
result["gel_plot"]=f"{job}_gel.png"

with open(f"python-runner/results/{job}.json","w") as f:
    json.dump(result,f,indent=2)
