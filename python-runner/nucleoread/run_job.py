import sys
import json
from .analysis import run_full_analysis

def main():
    if len(sys.argv) < 2:
        print("Usage: python -m python-runner.nucleoread.run_job <job_id>")
        sys.exit(1)

    job_id = sys.argv[1]

    # Read FASTA from STDIN
    fasta = sys.stdin.read()

    if not fasta.strip():
        raise ValueError("FASTA input is empty (stdin)")

    result = run_full_analysis(fasta, job_id)

    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
