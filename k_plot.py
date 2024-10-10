import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('quantiles_70-500_0.1.csv', index_col=0)

k_values = df.index

for gene in df.columns:
    
    plt.subplots_adjust()
    plt.plot(k_values, df[gene], marker='o', linestyle='-', label=gene)
    plt.xlabel('K values')
    plt.ylabel('Gamma values')
    plt.title(f'Plot for {gene}')
    #plt.legend()
    #plt.show()
    plt.title("K influence gamma values in gene: {}\n".format(gene))
    plt.savefig(f'k_plot/{gene}.png')

    plt.clf()
