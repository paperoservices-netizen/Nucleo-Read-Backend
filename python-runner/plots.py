import matplotlib.pyplot as plt

def gc_plot(dna,out):
    w=20
    vals=[(dna[i:i+w].count("G")+dna[i:i+w].count("C"))/w*100 for i in range(len(dna)-w)]
    plt.plot(vals)
    plt.title("GC Content")
    plt.ylabel("%GC")
    plt.xlabel("Position")
    plt.savefig(out)
    plt.close()
