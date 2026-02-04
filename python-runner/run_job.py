import sys, os, json
from robio import analyze

job_id = sys.argv[1]

input_file = f"input/{job_id}.fasta"

if not os.path.exists(input_file):
    print("Input FASTA not found:", input_file)
    exit(1)

with open(input_file, "r") as f:
    fasta = f.read()

os.makedirs("python-runner/results", exist_ok=True)

result = analyze(fasta, job_id)

with open(f"python-runner/results/{job_id}.json", "w") as f:
    json.dump(result, f, indent=2)

print("Nucleoread job complete")
