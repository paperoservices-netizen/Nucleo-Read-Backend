import sys
import json
import os
from robio_engine import analyze_fasta

if len(sys.argv) < 3:
    print("Usage: python run_job.py <fasta_file> <job_id>")
    sys.exit(1)

fasta_file = sys.argv[1]
job_id = sys.argv[2]

with open(fasta_file) as f:
    fasta = f.read()

print(f"🔬 Processing job {job_id}")
print(f"FASTA length: {len(fasta)}")

result = analyze_fasta(fasta)

os.makedirs("python-runner/results", exist_ok=True)
out_file = f"python-runner/results/{job_id}.json"

with open(out_file, "w") as f:
    json.dump({"job_id": job_id, **result}, f, indent=2)

print(f"✅ Result written to {out_file}")
