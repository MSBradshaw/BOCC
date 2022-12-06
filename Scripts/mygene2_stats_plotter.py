import pandas as pd
import matplotlib.pyplot as plt
import statistics


def clear_ax(ax, top=False, bottom=False, left=False, right=False):
    ax.spines['top'].set_visible(top)
    ax.spines['bottom'].set_visible(bottom)
    ax.spines['left'].set_visible(left)
    ax.spines['right'].set_visible(right)
    # ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(True)
    ax.yaxis.set_tick_params(width=0.0, labelsize=8)
    ax.xaxis.set_tick_params(width=0.0, labelsize=8)


df = pd.read_csv('my_gene_2_variantes_by_family_tables.csv')

# distribution of HPOs per family
df.groupby('family_id')
hpos = list(df.groupby('family_id')['phenotype_hpo_id'].nunique())

fig, ax = plt.subplots()
ax.hist(hpos, bins=20)
ax.set_xlabel('# HPOs per patient', size=10)
ax.set_ylabel('Count', size=10)
clear_ax(ax)
plt.tight_layout()
plt.savefig('Figures/mg2_hpo_dist.png')
plt.show()

# distribution of Genes per family
genes = list(df.groupby('family_id')['gene'].nunique())

fig, ax = plt.subplots()
ax.hist(genes, bins=30)
ax.set_xlabel('# Variants per patient', size=10)
ax.set_ylabel('Count', size=10)
clear_ax(ax)
plt.tight_layout()
plt.savefig('Figures/mg2_variant_dist.png')
plt.show()

# scatter of genes vs hpos, points are families
fig, ax = plt.subplots()
ax.scatter(genes,hpos,alpha=.3)
ax.set_xlabel('# Variants per patient', size=10)
ax.set_ylabel('# HPOs per patient', size=10)
clear_ax(ax)
plt.tight_layout()
plt.savefig('Figures/mg2_variant_hpo_scatter.png')
plt.show()

# dist of number of pairs per patient
pairs = [genes[i] * hpos[i] for i in range(len(genes))]

fig, ax = plt.subplots()
ax.hist(pairs, bins=100)
ax.set_xlabel('# Pairs per patient', size=10)
ax.set_ylabel('Count', size=10)
clear_ax(ax)
plt.tight_layout()
plt.yscale('log')
plt.savefig('Figures/mg2_pairs_dist.png')
plt.show()

# pairs without the outlier
filtered_pairs = [x for x in pairs if x < 1000]
fig, ax = plt.subplots()
ax.hist(filtered_pairs, bins=100)
ax.set_xlabel('# Pairs per patient', size=10)
ax.set_ylabel('Count', size=10)
clear_ax(ax)
plt.tight_layout()
plt.yscale('log')
plt.savefig('Figures/mg2_filtered_pairs_dist.png')
plt.show()


# pairs non-resolved patients
non_resolved_pairs = [genes[i] * hpos[i] for i in range(len(genes)) if genes[i] > 1]

fig, ax = plt.subplots()
ax.hist(non_resolved_pairs, bins=100)
ax.set_xlabel('# Pairs per non-resolved patient', size=10)
ax.set_ylabel('Count', size=10)
clear_ax(ax)
plt.tight_layout()
# plt.yscale('log')
plt.savefig('Figures/mg2_nonresolved_pairs_dist.png')
plt.show()

# pairs non-resolved patients and no outlier
filtered_non_resolved_pairs = [x for x in non_resolved_pairs if x < 1000]

fig, ax = plt.subplots()
ax.hist(filtered_non_resolved_pairs, bins=100)
ax.set_xlabel('# Pairs per non-resolved patient', size=10)
ax.set_ylabel('Count', size=10)
clear_ax(ax)
plt.tight_layout()
# plt.yscale('log')
plt.savefig('Figures/mg2_filtered_nonresolved_pairs_dist.png')
plt.show()

print('Mode:', str(statistics.mode(filtered_non_resolved_pairs)))
print('Median:', str(statistics.median(filtered_non_resolved_pairs)))
print('Mean:', str(statistics.mean(filtered_non_resolved_pairs)))