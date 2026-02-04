import re
from .constants import *
from .plots import gc_plot, virtual_gel

import re

def parse_fasta(fasta: str) -> str:
    lines = fasta.strip().splitlines()

    if not lines or not lines[0].startswith(">"):
        raise ValueError("Invalid FASTA: missing header")

    # Join all non-header lines, strip spaces
    seq = "".join(lines[1:]).replace(" ", "").upper()

    if not seq:
        raise ValueError("Invalid FASTA: empty sequence")

    if not re.fullmatch(r"[ACGTN]+", seq):
        raise ValueError("Invalid FASTA characters")

    return seq



def find_orfs(seq):
    seq = seq.replace("U", "T")
    orfs = []

    for frame in range(3):
        i = frame
        while i <= len(seq) - 3:
            if seq[i:i+3] == START_CODON:
                for j in range(i+3, len(seq)-2, 3):
                    if seq[j:j+3] in STOP_CODONS:
                        orfs.append({
                            "frame": frame + 1,
                            "start": i,
                            "end": j + 3,
                            "length": j + 3 - i
                        })
                        break
                i += 3
            else:
                i += 3
    return orfs


def translate(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3], "X")
                   for i in range(0, len(seq)-2, 3))


def design_primers(seq, length=20):
    fwd = seq[:length]
    rev = seq[-length:].translate(str.maketrans("ATGC","TACG"))[::-1]

    def tm(p):
        gc = p.count("G") + p.count("C")
        return round(64.9 + 41*(gc-16.4)/len(p), 1)

    return {
        "forward": fwd,
        "reverse": rev,
        "tm_forward": tm(fwd),
        "tm_reverse": tm(rev)
    }


def restriction_sites(seq):
    hits = []
    for name, site in RESTRICTION_ENZYMES.items():
        pos = [m.start()+1 for m in re.finditer(f"(?={site})", seq)]
        if pos:
            hits.append({"enzyme": name, "site": site, "positions": pos})
    return hits


def run_full_analysis(fasta, job_id):
    seq = parse_fasta(fasta)
    dna = seq.replace("U", "T")

    gc = dna.count("G") + dna.count("C")
    at = dna.count("A") + dna.count("T")
    amb = dna.count("N")

    comp = dna.translate(str.maketrans("ATGC","TACG"))
    rev_comp = comp[::-1]

    orfs = find_orfs(dna)
    best = max(orfs, key=lambda x: x["length"]) if orfs else None

    coding = protein = None
    if best:
        coding = dna[best["start"]:best["end"]]
        protein = translate(coding)

    primers = design_primers(coding) if coding else None
    enzymes = restriction_sites(dna)
    all_cuts = [p for e in enzymes for p in e["positions"]]

    gc_img = gc_plot(dna, f"results/{job_id}_gc.png")
    gel_img = virtual_gel(len(dna), all_cuts, f"results/{job_id}_gel.png")

    return {
        "core_analysis": {
            "sequence_type": "RNA" if "U" in seq else "DNA",
            "length": len(seq),
            "gc_percent": round(gc/len(seq)*100, 2),
            "at_percent": round(at/len(seq)*100, 2),
            "ambiguous_bases": amb
        },
        "transformations": {
            "dna_to_rna": dna.replace("T","U"),
            "complement": comp,
            "reverse_complement": rev_comp
        },
        "orf_analysis": {
            "total_orfs": len(orfs),
            "longest_orf": best,
            "coding_sequence": coding,
            "protein": protein
        },
        "primers": primers,
        "restriction_analysis": {
            "total_cuts": len(all_cuts),
            "enzymes": enzymes
        },
        "plots": {
            "gc_plot": f"{job_id}_gc.png",
            "virtual_gel": f"{job_id}_gel.png"
        }
    }
