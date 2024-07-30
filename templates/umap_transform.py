#!/usr/bin/env python

import umap
import pandas as pd

df = pd.read_csv("$kmer_freqs", delimiter="\t")

#UMAP
motifs = [x for x in df.columns.values if x not in ["read", "length"]]
X = df.loc[:,motifs]
X_embedded = umap.UMAP(n_neighbors=15, min_dist=0.1, verbose=2).fit_transform(X)

df_umap = pd.DataFrame(X_embedded, columns=["D1", "D2"])
umap_out = pd.concat([df["read"], df["length"], df_umap], axis=1)

umap_out.to_csv("umap.output.tsv", sep="\t", index=False)
