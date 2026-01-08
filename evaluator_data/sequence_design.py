###Pulling sequences and measurements from BigWig files
#May 13th
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pyfaidx
import pysam
import random
import numpy as np
random.seed(10)
from Bio import SeqIO
from Bio.Seq import Seq


enhancer_data = pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/final_APIs_toPOST/DeepBICCN2_EnhancerLabel_Evaluator/evaluator_data/biccn_enhancers.csv")
print(enhancer_data)
fasta_path_human = "/arc/project/st-cdeboer-1/iluthra/hg38.fa"
genome_human = pysam.FastaFile(fasta_path_human)

fasta_path_mouse = '/arc/project/st-cdeboer-1/Genomes/mm10.fa'
genome_mm10 = pysam.FastaFile(fasta_path_mouse)

enhancer_data["enhancer_sequence"] = None


for i in range(0, enhancer_data.shape[0]):
    chrom = enhancer_data.iloc[i]["chrom"]
    start = enhancer_data.iloc[i]["start"]
    end   = enhancer_data.iloc[i]["end"]

    if enhancer_data["species"].iloc[i] == "mouse":

        seq = genome_mm10.fetch(chrom, start,  end)
        seq = seq.upper()
        print(len(seq))
        enhancer_data.loc[i, "enhancer_sequence"] = seq
    if enhancer_data["species"].iloc[i] == "human":
        seq = genome_human.fetch(chrom, start,  end)
        seq = seq.upper()
        print(len(seq))
        enhancer_data.loc[i, "enhancer_sequence"] = seq


enhancer_data.to_csv("/scratch/st-cdeboer-1/iluthra/game_apis/final_APIs_toPOST/DeepBICCN2_EnhancerLabel_Evaluator/evaluator_data/biccn_enhancers_withSequence.csv", index = False, sep = '\t')
