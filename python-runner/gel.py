import matplotlib.pyplot as plt
import numpy as np

def virtual_gel(length,cuts,out):
    if not cuts:
        frags=[length]
    else:
        cuts=sorted(cuts)
        frags=[cuts[i+1]-cuts[i] for i in range(len(cuts)-1)]
        frags.append(length-cuts[-1]+cuts[0])

    plt.figure(figsize=(5,8))
    for f in frags:
        plt.hlines(f,0.5,1)
    plt.yscale("log")
    plt.savefig(out,dpi=150)
    plt.close()
