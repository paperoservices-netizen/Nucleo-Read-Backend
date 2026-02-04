import sys, os, json, base64
from robio import analyze

fasta = base64.b64decode(sys.argv[1]).decode("utf-8")
job_id = sys.argv[2]

os.makedirs("python-runner/results", exist_ok=True)

result = analyze(fasta, job_id)

with open(f"python-runner/results/{job_id}.json", "w") as f:
    json.dump(result, f, indent=2)

print("Nucleoread job complete")
