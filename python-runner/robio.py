import re, matplotlib.pyplot as plt

START="ATG"
STOP={"TAA","TAG","TGA"}

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

def parse_fasta(txt):
    lines=txt.strip().splitlines()
    if not lines[0].startswith(">"): raise ValueError("Invalid FASTA")
    return "".join(lines[1:]).replace(" ","").upper()

def find_orfs(seq):
    out=[]
    for f in range(3):
        for i in range(f,len(seq)-2,3):
            if seq[i:i+3]==START:
                for j in range(i+3,len(seq)-2,3):
                    if seq[j:j+3] in STOP:
                        out.append({"frame":f+1,"start":i,"end":j+3,"length":j+3-i})
                        break
    return out

def translate(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3],"X") for i in range(0,len(seq)-2,3))

def gc_plot(seq,path):
    w=50
    vals=[(seq[i:i+w].count("G")+seq[i:i+w].count("C"))/w*100 for i in range(len(seq)-w)]
    plt.plot(vals)
    plt.savefig(path)
    plt.close()

def analyze(fasta,job_id):
    dna=parse_fasta(fasta).replace("U","T")
    orfs=find_orfs(dna)
    best=max(orfs,key=lambda x:x["length"]) if orfs else None

    result={
        "length":len(dna),
        "gc":round((dna.count("G")+dna.count("C"))/len(dna)*100,2),
        "orfs":orfs,
        "best_orf":best
    }

    if best:
        coding=dna[best["start"]:best["end"]]
        result["protein"]=translate(coding)
        gc_plot(dna,f"python-runner/results/{job_id}_gc.png")

    return result
