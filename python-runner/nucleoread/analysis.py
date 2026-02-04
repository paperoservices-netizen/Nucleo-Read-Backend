import re
from .config import *

def parse_fasta(text):
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("FASTA header missing")
    return "".join(lines[1:]).upper()

def detect_type(seq):
    return "RNA" if "U" in seq else "DNA"

def find_orfs(seq):
    seq = seq.replace("U","T")
    orfs = []
    for frame in range(3):
        for i in range(frame, len(seq)-2, 3):
            if seq[i:i+3] == START_CODON:
                for j in range(i+3, len(seq)-2, 3):
                    if seq[j:j+3] in STOP_CODONS:
                        orfs.append({
                            "frame": frame+1,
                            "start": i,
                            "end": j+3,
                            "length": j+3-i
                        })
                        break
    return orfs

def translate(seq):
    seq = seq.replace("U","T")
    return "".join(GENETIC_CODE.get(seq[i:i+3], "X")
                   for i in range(0, len(seq)-2, 3))

def design_primers(seq, length=20):
    fwd = seq[:length]
    rev = seq[-length:].translate(str.maketrans("ATGC","TACG"))[::-1]
    return {"forward": fwd, "reverse": rev}

def find_restriction_sites(seq):
    sites = []
    for name, pat in RESTRICTION_ENZYMES.items():
        hits = [m.start()+1 for m in re.finditer(f'(?={pat})', seq)]
        if hits:
            sites.append({"enzyme": name, "positions": hits})
    return sites

def run_analysis(fasta):
    seq = parse_fasta(fasta)
    dna = seq.replace("U","T")
    orfs = find_orfs(dna)
    best = max(orfs, key=lambda x: x["length"]) if orfs else None

    protein = translate(dna[best["start"]:best["end"]]) if best else ""

    return {
        "sequence": {
            "type": detect_type(seq),
            "length": len(seq),
            "gc_percent": round((dna.count("G")+dna.count("C"))/len(dna)*100,2)
        },
        "orf": best,
        "protein": protein,
        "primers": design_primers(dna),
        "restriction_sites": find_restriction_sites(dna)
    }
