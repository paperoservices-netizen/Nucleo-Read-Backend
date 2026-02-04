import os
import json
import math
import matplotlib.pyplot as plt

DNA_COMPLEMENT = str.maketrans("ATGCatgc", "TACGtacg")


def parse_fasta(fasta: str):
    lines = [l.strip() for l in fasta.strip().splitlines() if l.strip()]
    if not lines or not lines[0].startswith(">"):
        raise ValueError("Invalid FASTA: missing header")

    seq = "".join(lines[1:]).upper()
    if not seq:
        raise ValueError("Invalid FASTA: empty sequence")

    return lines[0][1:], seq


def core_stats(seq):
    length = len(seq)
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    ambiguous = length - (gc + at)

    return {
        "sequence_type": "DNA",
        "length": length,
        "gc_content": round(gc / length * 100, 2),
        "at_content": round(at / length * 100, 2),
        "ambiguous_bases": ambiguous,
    }


def transformations(seq):
    rna = seq.replace("T", "U")
    complement = seq.translate(DNA_COMPLEMENT)
    rev_complement = complement[::-1]

    return {
        "dna_to_rna": rna,
        "complement": complement,
        "reverse_complement": rev_complement,
    }


def find_longest_orf(seq):
    stops = {"TAA", "TAG", "TGA"}
    longest = {"length": 0}

    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            codon = seq[i : i + 3]
            if codon == "ATG":
                j = i
                while j < len(seq) - 2:
                    stop = seq[j : j + 3]
                    if stop in stops:
                        length = j + 3 - i
                        if length > longest.get("length", 0):
                            longest = {
                                "frame": frame + 1,
                                "start": i + 1,
                                "end": j + 3,
                                "length": length,
                                "sequence": seq[i : j + 3],
                            }
                        break
                    j += 3
            i += 3
    return longest


CODON_TABLE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}


def translate(seq):
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        protein += CODON_TABLE.get(seq[i : i + 3], "X")
    return protein


def primer_design(seq):
    fwd = seq[:20]
    rev = seq[-20:].translate(DNA_COMPLEMENT)[::-1]

    def tm(s):
        return round(2 * (s.count("A") + s.count("T")) + 4 * (s.count("G") + s.count("C")), 1)

    return {
        "forward_primer": {"sequence": fwd, "tm": tm(fwd)},
        "reverse_primer": {"sequence": rev, "tm": tm(rev)},
    }


def gc_plot(seq, job_id):
    window = 20
    values = []
    for i in range(len(seq) - window + 1):
        frag = seq[i : i + window]
        values.append((frag.count("G") + frag.count("C")) / window * 100)

    plt.figure()
    plt.plot(values)
    plt.xlabel("Window position")
    plt.ylabel("GC %")
    plt.title("GC Content Sliding Window")
    path = f"results/{job_id}_gc.png"
    plt.savefig(path)
    plt.close()
    return path


def run_full_analysis(fasta, job_id):
    header, seq = parse_fasta(fasta)

    core = core_stats(seq)
    trans = transformations(seq)
    orf = find_longest_orf(seq)
    protein = translate(orf["sequence"]) if orf else ""
    primers = primer_design(seq)
    gc_img = gc_plot(seq, job_id)

    return {
        "header": header,
        "core_analysis": core,
        "transformations": trans,
        "longest_orf": orf,
        "protein": protein,
        "primers": primers,
        "plots": {
            "gc_content": gc_img
        }
    }
