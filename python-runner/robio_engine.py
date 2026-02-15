import re
import matplotlib.pyplot as plt
import numpy as np
import os

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

RESTRICTION_ENZYMES = {
    "EcoRI":"GAATTC","BamHI":"GGATCC","HindIII":"AAGCTT","NotI":"GCGGCCGC",
    "XhoI":"CTCGAG","KpnI":"GGTACC","AluI":"AGCT","HaeIII":"GGCC","TaqI":"TCGA"
}

def parse_fasta(text):
    lines = text.strip().splitlines()
    # allow plain sequence or FASTA with header
    if lines and lines[0].startswith(">"):
        seq = "".join(lines[1:]).upper()
    else:
        seq = "".join(lines).upper()
    if not seq or not re.fullmatch("[ATUGCN]+", seq):
        raise ValueError("Invalid FASTA")
    return seq.replace("U", "T")

def translate(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3], "X")
                   for i in range(0, len(seq)-2, 3))

# 3‑frame ORFs with metadata
def find_orfs_with_meta(seq):
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

def plot_gc(seq, out):
    window = 50
    if len(seq) <= window:
        return
    vals = [
        (seq[i:i+window].count("G") + seq[i:i+window].count("C")) / window * 100
        for i in range(len(seq) - window)
    ]
    plt.figure(figsize=(8, 3))
    plt.plot(vals)
    plt.title("GC Content (Sliding Window)")
    plt.xlabel("Position")
    plt.ylabel("GC %")
    plt.tight_layout()
    plt.savefig(out, dpi=150)
    plt.close()

def draw_gel(length, cuts, out, circular=True):
    cuts = sorted(set(cuts))

    if not cuts:
        fragments = [length]
    elif circular:
        fragments = [cuts[i+1] - cuts[i] for i in range(len(cuts)-1)]
        fragments.append((length - cuts[-1]) + cuts[0])
    else:
        full = [0] + cuts + [length]
        fragments = [full[i+1] - full[i] for i in range(len(full)-1)]

    fragments = sorted([f for f in fragments if f > 0], reverse=True)

    ladder = [100,200,300,400,500,600,700,800,900,
              1000,1500,2000,3000,5000,10000]

    plt.figure(figsize=(5, 8))
    ax = plt.gca()
    ax.set_facecolor("#121212")

    # ladder lane background
    plt.axvline(x=0.3, color="gray", lw=60, alpha=0.1)
    # sample lane background
    plt.axvline(x=0.7, color="gray", lw=60, alpha=0.1)

    # ladder bands
    for l in ladder:
        if l <= length * 1.5:
            plt.hlines(y=l, xmin=0.2, xmax=0.4, color="white", lw=1)
            plt.text(0.02, l, f"{l}", color="white", fontsize=8, va="center")

    # sample fragments
    last_y = None
    for frag in fragments:
        plt.hlines(y=frag, xmin=0.6, xmax=0.8, color="cyan", lw=3)
        if last_y is None:
            y = frag
        else:
            # avoid overlapping labels on log scale
            y = frag if abs(np.log10(last_y) - np.log10(frag)) > 0.05 else last_y * 0.85
        plt.text(0.85, y, f"{frag}bp", color="deepskyblue",
                 fontsize=9, fontweight="bold", va="center")
        last_y = y

    plt.yscale("log")
    plt.title("Virtual Gel Simulation", fontsize=14)
    plt.text(0.3, 13000, "Ladder", color="white", ha="center", fontweight="bold")
    plt.text(0.7, 13000, "Sample", color="white", ha="center", fontweight="bold")
    plt.ylabel("Size (bp)")
    plt.ylim(10, 15000)
    plt.xticks([])
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()

def design_primers(seq, length=20):
    seq = seq.replace("U", "T")
    if len(seq) < length * 2:
        return {"fwd": "N/A", "rev": "N/A", "fwd_tm": 0, "rev_tm": 0}

    fwd = seq[:length]
    rev = seq[-length:].translate(str.maketrans("ATGC", "TACG"))[::-1]

    def calc_tm(p):
        gc = p.count("G") + p.count("C")
        return round(64.9 + 41 * (gc - 16.4) / len(p), 1)

    return {
        "fwd": fwd,
        "rev": rev,
        "fwd_tm": calc_tm(fwd),
        "rev_tm": calc_tm(rev)
    }

def find_restriction_sites(seq):
    seq = seq.replace("U", "T")
    sites = []
    for name, site in RESTRICTION_ENZYMES.items():
        positions = [m.start() + 1 for m in re.finditer(f"(?={site})", seq)]
        if positions:
            sites.append({
                "name": name,
                "site": site,
                "positions": positions
            })
    return sites

def analyze_fasta(fasta, job):
    os.makedirs("python-runner/images", exist_ok=True)
    os.makedirs("python-runner/results", exist_ok=True)

    dna = parse_fasta(fasta)
    length = len(dna)

    gc_count = dna.count("G") + dna.count("C")
    at_count = dna.count("A") + dna.count("T")
    amb = dna.count("N")

    gc = gc_count / length * 100
    at = at_count / length * 100

    rna = dna.replace("T", "U")
    comp = dna.translate(str.maketrans("ATGC", "TACG"))
    rev = comp[::-1]

    # 6‑frame ORFs
    all_orfs = find_orfs_with_meta(dna)
    for o in find_orfs_with_meta(rev):
        o["frame"] = -o["frame"]
        all_orfs.append(o)
    best = max(all_orfs, key=lambda x: x["length"]) if all_orfs else None

    coding_seq = ""
    protein = ""
    best_orf_info = None
    primers = {"fwd": "N/A", "rev": "N/A", "fwd_tm": 0, "rev_tm": 0}

    if best:
        target = rev if best["frame"] < 0 else dna
        coding_seq = target[best["start"]:best["end"]]
        protein = translate(coding_seq)
        primers = design_primers(coding_seq)
        best_orf_info = {
            "frame": best["frame"],
            "start": best["start"],
            "end": best["end"],
            "length": best["length"],
            "strand": "Sense" if best["frame"] > 0 else "Antisense"
        }

    # restriction mapping
    sites = find_restriction_sites(dna)
    all_cuts = [p for s in sites for p in s["positions"]]

    # plots
    plot_gc(dna, f"python-runner/images/{job}_gc.png")
    draw_gel(len(dna), all_cuts, f"python-runner/images/{job}_gel.png", circular=True)

    # ---- nice console output (optional but matches your examples) ----
    print(f"🔬 CORE ANALYSIS")
    print("=" * 60)
    print(f"📏 Length                : {length} bases")
    print(f"📊 GC content (%)        : {round(gc, 2)}")
    print(f"📊 AT/AU content (%)     : {round(at, 2)}")
    print(f"✅ {amb} ambiguous bases detected")

    if best_orf_info:
        print(f"\n🧬 LONGEST ORF DETAILS (6-FRAME SCAN)")
        print("=" * 60)
        print(f"Frame     : {best_orf_info['frame']} ({best_orf_info['strand']})")
        print(f"Start-End : {best_orf_info['start']} – {best_orf_info['end']}")
        print(f"Length    : {best_orf_info['length']} bp")

    if sites:
        print(f"\n✂️ RESTRICTION MAPPING")
        print("=" * 60)
        for s in sites:
            print(f" 🔹 {s['name']} found at: {s['positions']}")

    # -----------------------------------------------------------------

    return {
        "length": length,
        "gc": round(gc, 2),
        "at": round(at, 2),
        "ambiguous": amb,
        "rna": rna,
        "complement": comp,
        "reverse_complement": rev,
        "best_orf": best_orf_info,
        "protein": protein,
        "primers": primers,
        "restriction_sites": sites
    }

