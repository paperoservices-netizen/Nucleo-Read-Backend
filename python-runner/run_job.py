import sys, json, os
from robio import analyze

fasta=sys.argv[1]
job_id=sys.argv[2]

os.makedirs("python-runner/results",exist_ok=True)

result=analyze(fasta,job_id)

with open(f"python-runner/results/{job_id}.json","w") as f:
    json.dump(result,f,indent=2)

print("Nucleoread job complete")
