import logging

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.INFO)

data = pd.read_csv('intermediate_results.csv', header = None)
data.columns = ['seed', 'doublet_rate', 'cell_count', 'score', 'sample_identity_is_pathological', 'v_score']
sns.set_theme(style="whitegrid")

data = data.loc[data['cell_count'].isin([1000, 2000, 3000, 4000])]

# Customize the appearance of the plots
plt.rcParams.update({
    'axes.titlesize': 26,
    'axes.labelsize': 26,
    'xtick.labelsize': 20,
    'ytick.labelsize': 26,
    'legend.fontsize': 16,
    'figure.titlesize': 18
})
fig, axes = plt.subplots(1, 5, figsize=(30, 6), sharey=True)

doublet_rates = [0, 0.02, 0.04, 0.06, 0.08]

for ax, rate in zip(axes, doublet_rates):
    sns.boxplot(x='cell_count', y='score', data=data.loc[(~data['sample_identity_is_pathological']) & (data['doublet_rate'] == rate)], ax=ax, color='#C97B84',flierprops=dict(marker='o', color='red', markersize=5, markerfacecolor='black'))
    ax.set_title(f'Doublet Rate: {rate}')
    ax.set_xlabel('Cell Count')

axes[0].set_ylabel('Sample assignment score')
plt.tight_layout()

logging.info("Saving figure4a.png")
plt.savefig("figure4a.png", dpi = 300)

plt.close()


plt.figure(figsize=(10, 8))
plt.scatter(x='v_score', y='score', data = data, color = '#C97B84')

# Customize the plot
plt.grid(True)
plt.xlabel('V-Score')
plt.ylabel('Sample assignment score')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_color('#444444')
plt.gca().spines['bottom'].set_color('#444444')
plt.gca().xaxis.label.set_color('#444444')
plt.gca().yaxis.label.set_color('#444444')
plt.gca().tick_params(axis='x', colors='#444444')
plt.gca().tick_params(axis='y', colors='#444444')

logging.info("Saving figure4b.png")
plt.savefig('figure4b.png', dpi = 300)
plt.close()


plt.figure(figsize=(10, 8))
sns.boxplot(x='sample_identity_is_pathological', y='score', data=data, color='#C97B84', flierprops=dict(marker='o', color='red', markersize=5, markerfacecolor='black'))

# Customize the plot
plt.grid(True)
plt.xlabel('Incorrect sample demultiplexing')
plt.ylabel('Sample assignment score')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_color('#444444')
plt.gca().spines['bottom'].set_color('#444444')
plt.gca().xaxis.label.set_color('#444444')
plt.gca().yaxis.label.set_color('#444444')
plt.gca().tick_params(axis='x', colors='#444444')
plt.gca().tick_params(axis='y', colors='#444444')

logging.info("Saving figure4c.png")
plt.savefig('figure4c.png', dpi = 300)
plt.close()
# Show the plot

