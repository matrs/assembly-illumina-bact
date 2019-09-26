import pandas as pd

configfile: "config.yaml"
df = pd.read_csv("samples.tsv", dtype=str, sep="\t").set_index("name", drop=False)
sample_name = df.index.to_list()[0]

def get_fastqs(wildcard):
    """Get raw FASTQ files from sample sheet. Returns a dictionary
    with keys `r1` and `r2`"""

    return df.loc[wildcard, ["r1", "r2"]].dropna().to_dict()

def get_reads_stats(file):
    with open(file, 'r') as fh: 
        stats = dict()
        for line in fh: 
            if line.startswith("min_len"):
                max_l, avg_l = line.strip().split(';')[1:3]
                # print(max_l, avg_l)
                stats["max_len"] = int(int(max_l.split(':')[1]) + 0.5)
                stats["avg_len"] = int(float(avg_l.split(':')[1])+ 0.5)
                stats["min_len"] = int(stats["avg_len"] / 2 + 0.5)
            if line.startswith("ALL"):
                stats["total_bp"] = int(line.strip().split('\t')[1])*2
    return stats

def read_params(file):
    df = pd.read_csv(file, sep='\t')
    records = df.to_dict('records')[0] #Returns a dict (original is a list of dict/s)
    # print('records',records)
    return records


def estimating_kmers(file):
    stats = read_params(file)
    # print(stats)
    max_k = min(config["max_k"], int(config["kmer_read_frac"]*stats["avg_len"]))
    min_k = 21 if stats["avg_len"] < 75 else config["min_k"]
    ks = max(5, int((max_k - min_k)/(config["kn"]-1)))
    ks = ks + 1 if ks % 2 == 1 else ks
    k_vals = [k + ks for k in range(min_k, max_k+1, ks) if (k+ks) <= max_k]
    k_vals.insert(0, min_k)
    # print("k_vals", k_vals)
    return ",".join( repr(e) for e in k_vals ) 

    
