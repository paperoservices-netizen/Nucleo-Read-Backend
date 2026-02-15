# python-runner/run_job.py
import sys, json, os
from robio_engine import analyze_fasta


fasta_path = sys.argv[1]
job = sys.argv[2]

with open(fasta_path) as f:
    fasta = f.read()

print("🔬 Processing job", job)
result = analyze_fasta(fasta, job)

os.makedirs("python-runner/results", exist_ok=True)
with open(f"python-runner/results/{job}.json", "w") as f:
    json.dump(result, f, indent=2)
