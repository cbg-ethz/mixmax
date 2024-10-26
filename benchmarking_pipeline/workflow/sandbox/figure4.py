import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

data = pd.read_csv('intermediate_results.csv', header = None)
data.columns = ['seed', 'doublet_rate', 'cell_count', 'score', 'sample_identity_is_pathological', 'v_score']
sns.set_theme(style="whitegrid")

# Customize the appearance of the plots
plt.rcParams.update({
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.titlesize': 18
})
fig, axes = plt.subplots(1, 6, figsize=(30, 6), sharey=True)

doublet_rates = [0, 0.02, 0.04, 0.06, 0.08]

for ax, rate in zip(axes, doublet_rates):
    sns.boxplot(x='cell_count', y='score', data=data.loc[(~data['sample_identity_is_pathological']) & (data['doublet_rate'] == rate)], ax=ax, color='#C97B84',flierprops=dict(marker='o', color='red', markersize=5, markerfacecolor='black'))
    ax.set_title(f'Doublet Rate: {rate}')
    ax.set_xlabel('Cell Count')
axes[0].set_ylabel('Score')
plt.tight_layout()
plt.savefig("figure3.png")
