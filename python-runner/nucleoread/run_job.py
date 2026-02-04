import sys, json, os
from .analysis import run_analysis
from .plots import save_gc_plot

fasta = sys.argv[1]
job_id = sys.argv[2]

result = run_analysis(fasta)

os.makedirs("python-runner/results", exist_ok=True)

json_path = f"python-runner/results/{job_id}.json"
gc_path = f"python-runner/results/gc_{job_id}.png"

save_gc_plot(fasta, gc_path)

result["plots"] = {
    "gc_plot": f"gc_{job_id}.png"
}

with open(json_path, "w") as f:
    json.dump(result, f, indent=2)

print(f"✅ NucleoRead completed: {json_path}")
