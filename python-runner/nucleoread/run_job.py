import sys
import json
import os
from .analysis import run_analysis

def main():
    if len(sys.argv) < 3:
        print("Usage: python -m python-runner.nucleoread.run_job <fasta_file> <job_id>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    job_id = sys.argv[2]

    with open(fasta_file, "r") as f:
        fasta = f.read()

    try:
        result = run_analysis(fasta)
    except Exception as e:
        result = {
            "status": "failed",
            "error": str(e)
        }

    os.makedirs("results", exist_ok=True)
    out_path = f"results/{job_id}.json"

    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"✅ Result written to {out_path}")

if __name__ == "__main__":
    main()
