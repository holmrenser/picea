# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Sequences
#
# picea represents a single biological sequence as a {class}`~picea.Sequence`, and groups of sequences as a
# {class}`~picea.SequenceCollection` (unaligned) or a {class}`~picea.MultipleSequenceAlignment` (aligned).

# %%
import numpy as np
from matplotlib import pyplot as plt

from picea import MultipleSequenceAlignment, Sequence, SequenceCollection

# %% [markdown]
# ## Single sequences
#
# A sequence has a header and a sequence string. When no alphabet is given, picea picks the best matching one (DNA
# or amino acid).

# %%
dna = Sequence("my_gene", "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGATAA")
dna

# %%
dna.alphabet.name, len(dna)

# %%
Sequence("peptide", "MKVLAAGIVGLL").alphabet.name

# %% [markdown]
# Transformations return new sequence objects, so they can be chained.

# %%
dna.reverse_complement.sequence

# %%
dna.amino_acids.sequence

# %%
dna[:12].lowercase.sequence

# %%
print(dna.to_fasta(linewidth=20))

# %% [markdown]
# ## Sequence collections
#
# Collections are read from and written to fasta or json. This file contains protein sequences of
# hydroxycinnamoyl transferase (HCT) homologs from several plant species.

# %%
seqs = SequenceCollection.from_fasta(filename="data/HCT.fasta")
len(seqs), seqs.headers[:3]

# %% [markdown]
# Index by header to get a {class}`~picea.Sequence`, iterate to get all of them, or use `iloc` to select a subset by
# position.

# %%
seqs["AT5G48930.1"].sequence[:60]

# %%
lengths = [len(seq) for seq in seqs]
min(lengths), max(lengths)

# %%
subset = seqs.iloc[:3]
print(subset.to_fasta(linewidth=60))

# %% [markdown]
# Headers and sequences can be modified in place.

# %%
subset.rename_inplace(lambda header: header.split(".")[0])
subset.headers

# %%
subset.add(Sequence("my_protein", "MASKGEELFTG"))
subset.headers

# %%
print(subset.to_json(indent=2)[:200])

# %% [markdown]
# ## Multiple sequence alignments
#
# A {class}`~picea.MultipleSequenceAlignment` stores sequences of equal length (gaps are `-`).

# %%
msa = MultipleSequenceAlignment.from_fasta(filename="data/multiple_sequence_alignment.fasta")
msa.n_seqs, msa.n_chars

# %% [markdown]
# Alignment columns are easy to analyse with numpy, for example the fraction of gaps per alignment position.

# %%
columns = np.array([list(sequence) for sequence in msa.sequences])
gap_fraction = (columns == "-").mean(axis=0)

fig, ax = plt.subplots(figsize=(10, 2))
ax.plot(gap_fraction)
ax.set_xlabel("Alignment position")
ax.set_ylabel("Gap fraction");

# %% [markdown]
# Unaligned collections can be aligned with an external aligner that reads fasta from stdin, such as
# [MAFFT](https://mafft.cbrc.jp/alignment/software/) (not run here):
#
# ```python
# msa = seqs.align(method="mafft")
# ```
