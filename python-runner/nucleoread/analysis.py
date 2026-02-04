import re

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

# ---------------- FASTA PARSER ----------------
def parse_fasta(text: str) -> str:
    lines = [l.strip() for l in text.strip().splitlines() if l.strip()]

    if not lines or not lines[0].startswith(">"):
        raise ValueError("Invalid FASTA: missing header")

    seq = "".join(lines[1:]).upper()

    if not seq:
        raise ValueError("Invalid FASTA: empty sequence")

    if not re.fullmatch("[ATUGCN]+", seq):
        raise ValueError("Invalid FASTA: illegal characters")

    return seq

# ---------------- ORF FINDER ----------------
def find_orfs(seq: str):
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

# ---------------- TRANSLATION ----------------
def translate(seq: str) -> str:
    seq = seq.replace("U", "T")
    return "".join(
        GENETIC_CODE.get(seq[i:i+3], "X")
        for i in range(0, len(seq) - 2, 3)
    )

# ---------------- MAIN ANALYSIS ----------------
def run_analysis(fasta_text: str) -> dict:
    seq = parse_fasta(fasta_text)
    dna = seq.replace("U", "T")

    length = len(dna)
    gc_count = dna.count("G") + dna.count("C")
    at_count = dna.count("A") + dna.count("T")

    gc_percent = round((gc_count / length) * 100, 2) if length > 0 else 0.0

    # ORFs (sense + antisense)
    orfs = find_orfs(dna)

    rev_comp = dna.translate(str.maketrans("ATGC", "TACG"))[::-1]
    for o in find_orfs(rev_comp):
        o["frame"] = -o["frame"]
        orfs.append(o)

    best_orf = max(orfs, key=lambda x: x["length"]) if orfs else None

    protein = None
    if best_orf:
        target = rev_comp if best_orf["frame"] < 0 else dna
        coding = target[best_orf["start"]:best_orf["end"]]
        protein = translate(coding)

    return {
        "status": "success",
        "sequence_type": "RNA" if "U" in seq else "DNA",
        "length": length,
        "gc_percent": gc_percent,
        "at_percent": round((at_count / length) * 100, 2) if length > 0 else 0.0,
        "orfs_found": len(orfs),
        "longest_orf": best_orf,
        "protein": protein
    }
