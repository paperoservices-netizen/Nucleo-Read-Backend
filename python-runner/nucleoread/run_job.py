import sys
import json
import os

from .analysis import run_analysis
from .plots import save_gc_plot

if len(sys.argv) < 3:
    print("Usage: python -m python-runner.nucleoread.run_job <fasta> <job_id>")
    sys.exit(1)

fasta = sys.argv[1]
job_id = sys.argv[2]

result = run_analysis(fasta)

RESULT_DIR = "python-runner/results"
os.makedirs(RESULT_DIR, exist_ok=True)

gc_plot_path = f"{RESULT_DIR}/gc_{job_id}.png"
plot = save_gc_plot(fasta, gc_plot_path)

if plot:
    result["plots"] = {"gc_plot": f"gc_{job_id}.png"}

out_file = f"{RESULT_DIR}/{job_id}.json"
with open(out_file, "w") as f:
    json.dump(result, f, indent=2)

print(f"✅ NucleoRead result written to {out_file}")
