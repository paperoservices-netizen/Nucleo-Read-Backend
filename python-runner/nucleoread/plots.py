import matplotlib.pyplot as plt
import numpy as np

def gc_plot(seq, out_path):
    window = 50
    if len(seq) < window:
        return None

    vals = [
        (seq[i:i+window].count("G") + seq[i:i+window].count("C")) / window * 100
        for i in range(len(seq) - window)
    ]

    plt.figure(figsize=(10,4))
    plt.plot(vals)
    plt.title("GC Content (Sliding Window)")
    plt.xlabel("Position")
    plt.ylabel("GC %")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    return out_path


def virtual_gel(sequence_length, cuts, out_path):
    if not cuts:
        fragments = [sequence_length]
    else:
        cuts = sorted(cuts)
        fragments = [cuts[i+1]-cuts[i] for i in range(len(cuts)-1)]
        fragments.append(sequence_length - cuts[-1] + cuts[0])

    fragments = sorted(fragments, reverse=True)

    ladder = [100,200,300,500,700,1000,1500,2000,3000,5000,10000]

    plt.figure(figsize=(5,8))
    ax = plt.gca()
    ax.set_facecolor("#121212")

    for l in ladder:
        if l <= sequence_length * 1.5:
            plt.hlines(l, 0.1, 0.3, color="white")
            plt.text(0.02, l, str(l), color="white", fontsize=8)

    for f in fragments:
        plt.hlines(f, 0.6, 0.9, color="cyan", lw=3)
        plt.text(0.92, f, f"{f}bp", color="cyan", fontsize=9)

    plt.yscale("log")
    plt.ylabel("Base pairs")
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    return out_path
