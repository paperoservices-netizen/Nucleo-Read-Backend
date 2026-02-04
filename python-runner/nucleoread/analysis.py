import re
import matplotlib.pyplot as plt
import numpy as np
import base64
from io import BytesIO

# ================= CONSTANTS =================
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
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC", "XhoI": "CTCGAG", "KpnI": "GGTACC",
    "AluI": "AGCT", "HaeIII": "GGCC", "TaqI": "TCGA"
}

# ================= CORE FUNCTIONS =================
def parse_fasta(text):
    lines = text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("FASTA header missing (must start with >)")
    seq = "".join(lines[1:]).upper().replace(" ", "")
    if not re.fullmatch("[ATUGCN]+", seq):
        raise ValueError("Invalid characters detected (Allowed: A,T,U,G,C,N)")
    return seq

def detect_type(seq):
    return "RNA" if "U" in seq else "DNA"

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
                            "frame": frame+1,
                            "start": i,
                            "end": j+3,
                            "length": j+3-i
                        })
                        break
                i += 3
            else:
                i += 3
    return orfs

def translate_dna(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3], "X")
                   for i in range(0, len(seq)-2, 3))

def design_primers(seq, length=20):
    seq = seq.replace("U", "T")
    if len(seq) < length * 2:
        return {"forward":"N/A","reverse":"N/A","tm_forward":0,"tm_reverse":0}

    fwd = seq[:length]
    rev = seq[-length:].translate(str.maketrans("ATGC","TACG"))[::-1]

    def calc_tm(p):
        gc = p.count("G") + p.count("C")
        return round(64.9 + 41 * (gc - 16.4) / len(p), 1)

    return {"forward":fwd, "reverse":rev,
            "tm_forward":calc_tm(fwd), "tm_reverse":calc_tm(rev)}

def find_restriction_sites(seq):
    seq = seq.replace("U", "T")
    sites = []
    for name, site in RESTRICTION_ENZYMES.items():
        positions = [m.start()+1 for m in re.finditer(f'(?={site})', seq)]
        if positions:
            sites.append({"enzyme":name,"site":site,"positions":positions})
    return sites

# ================= PLOTTING =================
def gc_plot_base64(seq, window=50):
    seq = seq.replace("U", "T")
    if len(seq) < window: return None
    vals = [(seq[i:i+window].count("G")+seq[i:i+window].count("C"))/window*100
            for i in range(len(seq)-window)]
    plt.figure(figsize=(10,4))
    plt.plot(vals, color="#2c3e50")
    plt.title("GC Content (Sliding Window)")
    plt.xlabel("Position")
    plt.ylabel("GC %")
    buf = BytesIO()
    plt.savefig(buf, dpi=150, bbox_inches="tight")
    plt.close()
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def virtual_gel_base64(seq_length, cut_positions):
    cuts = sorted(set(cut_positions))
    fragments = [seq_length] if not cuts else [cuts[i+1]-cuts[i] for i in range(len(cuts)-1)] + [(seq_length - cuts[-1]) + cuts[0]]
    fragments = sorted([f for f in fragments if f>0], reverse=True)
    plt.figure(figsize=(5,8))
    for i, frag in enumerate(fragments):
        plt.hlines(y=frag, xmin=0.6, xmax=0.8, color="cyan", lw=3)
    plt.yscale("log")
    buf = BytesIO()
    plt.savefig(buf, dpi=150, bbox_inches="tight")
    plt.close()
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

# ================= MAIN FUNCTION =================
def run_full_analysis(fasta_text):
    seq = parse_fasta(fasta_text)
    seq_type = detect_type(seq)
    dna_seq = seq.replace("U","T")
    complement = dna_seq.translate(str.maketrans("ATGC","TACG"))
    reverse_complement = complement[::-1]

    gc_count = dna_seq.count("G") + dna_seq.count("C")
    at_count = dna_seq.count("A") + dna_seq.count("T")
    ambiguous_count = dna_seq.count("N")

    # ORF
    orfs = find_orfs(dna_seq)
    best_orf = max(orfs, key=lambda x:x["length"]) if orfs else None
    coding_sequence = translate_protein = primers = None
    if best_orf:
        coding_sequence = dna_seq[best_orf["start"]:best_orf["end"]]
        translate_protein = translate_dna(coding_sequence)
        primers = design_primers(coding_sequence)

    # Restriction sites
    restriction_sites = find_restriction_sites(dna_seq)
    all_cuts = [p for s in restriction_sites for p in s["positions"]]

    # Plots
    gc_plot_b64 = gc_plot_base64(dna_seq)
    gel_plot_b64 = virtual_gel_base64(len(dna_seq), all_cuts)

    return {
        "core_analysis": {
            "sequence_type": seq_type,
            "length": len(seq),
            "gc_percent": round(gc_count/len(seq)*100,2),
            "at_percent": round(at_count/len(seq)*100,2),
            "ambiguous_bases": ambiguous_count
        },
        "transformations": {
            "dna_to_rna": dna_seq.replace("T","U"),
            "complement": complement,
            "reverse_complement": reverse_complement
        },
        "orf_analysis": {
            "total_orfs": len(orfs),
            "longest_orf": best_orf,
            "coding_sequence": coding_sequence,
            "protein": translate_protein
        },
        "primers": primers,
        "restriction_analysis": {
            "total_cuts": len(all_cuts),
            "enzymes": restriction_sites
        },
        "plots": {
            "gc_plot": gc_plot_b64,
            "virtual_gel": gel_plot_b64
        }
    }
