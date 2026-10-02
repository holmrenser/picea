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
# # Sequence annotation
#
# Gene models are read from GFF3 (or GTF) into a {class}`~picea.SequenceAnnotation`: a directed acyclic graph of
# {class}`~picea.SequenceInterval` objects (genes, mRNAs, exons, CDSs, ...) linked by their `Parent` attributes.

# %%
from picea import SequenceAnnotation, SequenceInterval

# %% [markdown]
# ## Single intervals
#
# One GFF3 line corresponds to one interval. The first eight columns become attributes, and so does every key in
# the ninth column. Attribute keys are lowercased (except `ID`) and attribute values are always lists.

# %%
interval = SequenceInterval.from_gff_line("ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN")
interval

# %%
interval.seqid, interval.start, interval.end, interval.strand, interval.name

# %%
interval.to_gff_line()

# %% [markdown]
# ## Gene models
#
# This file contains a single gene from the *Medicago truncatula* genome annotation.

# %%
annotation = SequenceAnnotation.from_gff(filename="data/genemodel.gff3")
len(annotation)

# %% [markdown]
# Intervals are accessed by ID, and can be grouped or filtered with any function of an interval.

# %%
by_type = annotation.groupby(lambda interval: interval.interval_type)
{interval_type: len(intervals) for interval_type, intervals in by_type.items()}

# %%
long_exons = annotation.filter(lambda interval: interval.interval_type == "exon" and interval.end - interval.start > 500)
[(exon.ID, exon.end - exon.start) for exon in long_exons]

# %% [markdown]
# `children` and `parents` are **transitive**: they contain all descendants and all ancestors of an interval. They
# return a new {class}`~picea.SequenceAnnotation`, so `groupby` and `filter` work on them as well.

# %%
gene = annotation["gene:MtrunA17Chr1g0184451"]
for child in gene.children:
    print(f"{child.interval_type:16}{child.start:>10}{child.end:>10}  {child.ID}")

# %%
cds = by_type["CDS"].elements[0]
[(parent.interval_type, parent.ID) for parent in cds.parents]

# %%
gene.children.groupby(lambda interval: interval.interval_type)["exon"].elements

# %% [markdown]
# ## Output
#
# Intervals and annotations can be written to GFF3, GTF, and json.

# %%
gene.to_gff_line()

# %%
print(gene.to_json(indent=2))

# %% [markdown]
# GTF files have no interval IDs or `Parent` attributes: intervals are linked by their `gene_id` and `transcript_id`
# attributes instead. Reading GTF recreates the gene model.

# %%
gtf = annotation.to_gtf()
print("\n".join(gtf.split("\n")[:3]))

# %%
from_gtf = SequenceAnnotation.from_gtf(string=gtf)
[(interval.interval_type, interval.ID) for interval in from_gtf["gene:MtrunA17Chr1g0184451"].children][:4]
