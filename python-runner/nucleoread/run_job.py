import sys
import os
import json
from .analysis import run_full_analysis


def main():
    if len(sys.argv) < 3:
        print("Usage: run_job <fasta> <job_id>")
        sys.exit(1)

    fasta = sys.argv[1]
    job_id = sys.argv[2]

    os.makedirs("results", exist_ok=True)

    try:
        result = run_full_analysis(fasta, job_id)
        out_file = f"results/{job_id}.json"
        with open(out_file, "w") as f:
            json.dump(result, f, indent=2)

        print(f"✅ Result written to {out_file}")

    except Exception as e:
        err = {
            "status": "failed",
            "error": str(e)
        }
        out_file = f"results/{job_id}.json"
        with open(out_file, "w") as f:
            json.dump(err, f, indent=2)
        print("❌ Analysis failed:", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
