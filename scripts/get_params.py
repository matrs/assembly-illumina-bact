import pandas as pd

def get_reads_stats(file):
    with open(file, 'r') as fh: 
        stats = dict()
        for line in fh: 
            if line.startswith("min_len"):
                max_l, avg_l = line.strip().split(';')[1:3]
                stats["max_len"] = int(int(max_l.split(':')[1]) + 0.5)
                stats["avg_len"] = int(float(avg_l.split(':')[1])+ 0.5)
                stats["min_len"] = int(stats["avg_len"] / 2 + 0.5)
            if line.startswith("ALL"):
                stats["total_bp"] = int(line.strip().split('\t')[1])*2
    return stats

def genome_size(file):
    with open(file, "r") as fh:
        g_size = dict()
        for line in fh:
            if "size" in line:
                g_size["gsize"] = int(float(line.strip().split(" ")[3]))
    #print(g_size)
    return g_size

# gs  =genome_size(file)
# stats= get_reads_stats(file)
all_vals = {**genome_size(snakemake.input.gsize), **get_reads_stats(snakemake.input.stats)}
df = pd.DataFrame(all_vals, index=[0]) #Index has to be a collection of some sort

df.to_csv(snakemake.output[0], sep='\t')



# df.to_dict('records')[0] 