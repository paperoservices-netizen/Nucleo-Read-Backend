import re, json
import matplotlib.pyplot as plt
import numpy as np

START_CODON = "ATG"
STOP_CODONS = {"TAA","TAG","TGA"}

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

def parse_fasta(text):
    lines=text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("FASTA header missing")
    seq="".join(lines[1:]).upper().replace(" ","")
    return seq

def detect_type(seq):
    return "RNA" if "U" in seq else "DNA"

def complement(seq):
    return seq.translate(str.maketrans("ATGC","TACG"))

def reverse_complement(seq):
    return complement(seq)[::-1]

def find_orfs(seq):
    seq=seq.replace("U","T")
    best=None
    for frame in range(3):
        i=frame
        while i<len(seq)-2:
            if seq[i:i+3]=="ATG":
                for j in range(i+3,len(seq)-2,3):
                    if seq[j:j+3] in STOP_CODONS:
                        length=j+3-i
                        if not best or length>best["length"]:
                            best={"frame":frame+1,"start":i,"end":j+3,"length":length}
                        break
            i+=3
    return best

def translate(seq):
    seq=seq.replace("U","T")
    p=""
    for i in range(0,len(seq)-2,3):
        p+=GENETIC_CODE.get(seq[i:i+3],"X")
    return p

def primers(seq,n=20):
    fwd=seq[:n]
    rev=complement(seq[-n:])[::-1]
    def tm(p):
        gc=p.count("G")+p.count("C")
        return round(64.9+41*(gc-16.4)/len(p),1)
    return {"fwd":fwd,"rev":rev,"fwd_tm":tm(fwd),"rev_tm":tm(rev)}

def restriction_sites(seq):
    out=[]
    for name,site in RESTRICTION_ENZYMES.items():
        pos=[m.start()+1 for m in re.finditer(f"(?={site})",seq)]
        if pos:
            out.append({"name":name,"positions":pos})
    return out

def analyze_fasta(fasta):
    seq=parse_fasta(fasta)
    dna=seq.replace("U","T")

    gc=(dna.count("G")+dna.count("C"))/len(dna)*100

    orf=find_orfs(dna)
    coding=dna[orf["start"]:orf["end"]]
    protein=translate(coding)
    prim=primers(coding)
    sites=restriction_sites(dna)

    return {
        "length":len(dna),
        "gc":round(gc,2),
        "rna":dna.replace("T","U"),
        "complement":complement(dna),
        "reverse_complement":reverse_complement(dna),
        "orf":[orf["start"],orf["end"]],
        "protein":protein,
        "primers":prim,
        "restriction_sites":sites
    }
