# KEGG-modules-MMSEQS2-search

From KEGG database, all module id's downloaded. In module page all orthologs of that module and all genes retrieved using KEGG API. 
Modules are filtered to contain only prokaryotic genes from all orthologs belonging to that module.
Fasta sequences belonging to each gene is retrieved.

Previous database containin predicted proteins from metagenomic samples (with predicted taxonomies) are filtered into 2 categories: Eukaryotic and Unknown. Rest of the proteins are not included in the further search.

KEGG genes from prokaryotic modules searched against predicted Eukaryotic and predicted Unknown protein datasets using MMSEQS2 program.

Hit results from MMSEQS2 further investigated.
Output file converted into data frame, hits that belong to same contig, same module and same ortholog removed. (Only the hit wiht lowest E-value kept in the data frame).
"1 ortholog contains many gene sequences from different species which results in the same alignment therefore removing the duplicates was necessary"

Further grouping applied to the data frame:
Hits in 1 contig from 1 module (to see how many of the proteins in the same KEGG pathway (belonging the same module) gave a hit with same contig.
Hits in 1 sample from 1 module (to see how many of the proteins in the same KEGG pathway (belonging the same module) gave a hit with same sample.
+How many of these hits were from the different proteins in that contig.
