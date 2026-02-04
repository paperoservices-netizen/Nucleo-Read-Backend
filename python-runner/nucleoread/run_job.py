import sys, json, os
from .analysis import run_full_analysis

def main():
    fasta = sys.argv[1]
    job_id = sys.argv[2]

    os.makedirs("results", exist_ok=True)

    result = run_full_analysis(fasta, job_id)

    with open(f"results/{job_id}.json", "w") as f:
        json.dump(result, f, indent=2)

    print(f"✅ Result written to results/{job_id}.json")

if __name__ == "__main__":
    main()
