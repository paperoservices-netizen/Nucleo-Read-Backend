import re
import matplotlib.pyplot as plt
import numpy as np
import os

START_CODON="ATG"
STOP_CODONS={"TAA","TAG","TGA"}

GENETIC_CODE={
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

RESTRICTION_ENZYMES={
 "EcoRI":"GAATTC","BamHI":"GGATCC","HindIII":"AAGCTT","NotI":"GCGGCCGC",
 "XhoI":"CTCGAG","KpnI":"GGTACC","AluI":"AGCT","HaeIII":"GGCC","TaqI":"TCGA"
}

def parse_fasta(text):
    lines=text.strip().splitlines()
    seq="".join(l for l in lines if not l.startswith(">")).upper()
    if not seq or not re.fullmatch("[ATUGCN]+",seq):
        raise ValueError("Invalid FASTA")
    return seq.replace("U","T")

def find_orfs(seq):
    orfs=[]
    for frame in range(3):
        for i in range(frame,len(seq)-2,3):
            if seq[i:i+3]==START_CODON:
                for j in range(i+3,len(seq)-2,3):
                    if seq[j:j+3] in STOP_CODONS:
                        orfs.append((i,j+3))
                        break
    return orfs

def translate(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3],"X") for i in range(0,len(seq)-2,3))

def plot_gc(seq,out):
    window=50
    vals=[(seq[i:i+window].count("G")+seq[i:i+window].count("C"))/window*100 for i in range(len(seq)-window)]
    plt.figure(figsize=(8,3))
    plt.plot(vals)
    plt.savefig(out,dpi=150)
    plt.close()

def draw_gel(length, cuts, out, circular=True):
    cuts = sorted(cuts)
    
    # calculate fragments
    if not cuts:
        fragments = [length]
    elif circular:
        fragments = [cuts[0]] + [cuts[i+1]-cuts[i] for i in range(len(cuts)-1)] + [length-cuts[-1]]
    else:
        full = [0] + cuts + [length]
        fragments = [full[i+1]-full[i] for i in range(len(full)-1)]
    
    fragments = sorted([f for f in fragments if f>0], reverse=True)

    # ladder sizes
    ladder = [100,200,300,400,500,600,700,800,900,
              1000,1500,2000,3000,5000,10000]

    plt.figure(figsize=(5,8))
    ax = plt.gca()
    ax.set_facecolor("#121212")

    # draw ladder
    for l in ladder:
        if l <= length*1.5:
            plt.hlines(y=l, xmin=0.2, xmax=0.4, color="white", lw=1)
            plt.text(0.02, l, f"{l}", color="white", fontsize=8, va="center")

    # draw sample fragments
    for frag in fragments:
        plt.hlines(y=frag, xmin=0.6, xmax=0.8, color="cyan", lw=3)
        plt.text(0.82, frag, f"{frag}bp", color="deepskyblue",
                 fontsize=9, fontweight="bold", va="center")

    plt.yscale("log")
    plt.title("Virtual Gel Simulation")
    plt.ylabel("Size (bp)")
    plt.xticks([])
    plt.grid(False)
    plt.savefig(out,dpi=150,bbox_inches="tight")
    plt.close()

def design_primers(seq, length=20):
    seq = seq.replace("U", "T")
    if len(seq) < length * 2:
        return {"fwd":"N/A","rev":"N/A","fwd_tm":0,"rev_tm":0}

    fwd = seq[:length]
    rev = seq[-length:].translate(str.maketrans("ATGC","TACG"))[::-1]

    def calc_tm(p):
        gc = p.count("G") + p.count("C")
        return round(64.9 + 41 * (gc - 16.4) / len(p), 1)

    return {
        "fwd": fwd,
        "rev": rev,
        "fwd_tm": calc_tm(fwd),
        "rev_tm": calc_tm(rev)
    }


def analyze_fasta(fasta, job):
    os.makedirs("python-runner/images", exist_ok=True)
    os.makedirs("python-runner/results", exist_ok=True)

    dna = parse_fasta(fasta)
    gc = (dna.count("G")+dna.count("C"))/len(dna)*100
    rna = dna.replace("T","U")
    comp = dna.translate(str.maketrans("ATGC","TACG"))
    rev = comp[::-1]

    orfs = find_orfs(dna)
    start, end = max(orfs, key=lambda x: x[1]-x[0])
    coding_seq = dna[start:end]
    protein = translate(coding_seq)

    # restriction cuts
    cuts = []
    for enz, site in RESTRICTION_ENZYMES.items():
        for m in re.finditer(site, dna):
            cuts.append(m.start())

    # plots
    plot_gc(dna,f"python-runner/images/{job}_gc.png")
    draw_gel(len(dna), cuts, f"python-runner/images/{job}_gel.png", circular=True)

    # primers
    primers = design_primers(coding_seq)

    return {
        "length": len(dna),
        "gc": round(gc,2),
        "rna": rna,
        "complement": comp,
        "reverse_complement": rev,
        "orf": [start,end],
        "protein": protein,
        "primers": primers
    }

