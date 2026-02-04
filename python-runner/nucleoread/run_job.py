import sys, json, os
from .analysis import run_full_analysis

def main():
    if len(sys.argv) < 3:
        print("Usage: python -m python_runner.nucleoread.run_job <fasta_file> <job_id>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    job_id = sys.argv[2]

    with open(fasta_file, "r") as f:
        fasta_text = f.read()

    result = run_full_analysis(fasta_text)

    os.makedirs("results", exist_ok=True)
    output_path = f"results/{job_id}.json"
    with open(output_path, "w") as f:
        json.dump(result, f, indent=4)

    print(f"✅ Result saved to {output_path}")

if __name__ == "__main__":
    main()
