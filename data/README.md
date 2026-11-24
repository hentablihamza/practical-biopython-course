# Datasets

All files in this folder are small, real, public records, committed so that every notebook runs offline.
Nothing was invented or edited by hand: each file is the unmodified response of the source service, or
(where stated) a direct extraction from such a response by [`scripts/refresh_cache.py`](../scripts/refresh_cache.py).
Retrieval dates, queries, and tool versions are recorded in [`provenance.json`](provenance.json).
Rebuild everything with `NCBI_EMAIL=you@example.org python scripts/refresh_cache.py`.

## `sequences/` — UniProt (notebooks 01–03, 06, 09, 12)

Downloaded from the [UniProt REST API](https://rest.uniprot.org) (UniProtKB/Swiss-Prot, reviewed entries).

- `P00533_EGFR_human.fasta` — human EGFR, the course's running example (UniProt FASTA, unmodified).
- `P00533_uniprot_features.json` — the EGFR entry's signal peptide, topology, transmembrane helix, kinase domain,
  binding-site, active-site, region, and glycosylation features (a subset of the entry's JSON `features` list).
- `kinase_panel_full.fasta` — full-length sequences of 15 kinases: the human ErbB family (EGFR P00533, ERBB2 P04626,
  ERBB3 P21860, ERBB4 Q15303), EGFR orthologues (mouse Q01279, fruit fly P04412, nematode LET-23 P24348), and drug
  targets from other kinase families (ABL1 P00519, SRC P12931, BRAF P15056, VEGFR2/KDR P35968, BTK Q06187,
  CDK2 P24941, MET P08581, ALK Q9UM73).
- `kinase_domains.fasta` — the "Protein kinase" domain of each panel protein, cut at the boundaries UniProt annotates.
- `kinase_panel.csv` — accession, entry name, gene, protein name, organism, length, kinase-domain start/end,
  and UniProt entry version.

**Cite:** The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025.
*Nucleic Acids Res.* **2025**, *53*, D609–D617. https://doi.org/10.1093/nar/gkae1010
UniProt data are freely available under the CC BY 4.0 licence (https://www.uniprot.org/help/license).

## `entrez_cache/` — NCBI Entrez (notebooks 02–05)

Raw responses from the NCBI [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) via `Bio.Entrez`.
Notebook 05 queries NCBI live and only falls back to these files when NCBI is unreachable or `COURSE_OFFLINE=1`.

- `NM_005228.gb`, `NM_005228.fasta` — RefSeq human EGFR mRNA, transcript variant 1 (`NM_005228.5`), efetch GenBank and FASTA.
- `esearch_nucleotide_EGFR.xml`, `esummary_NM_005228.xml` — the search `EGFR[Gene Name] AND Homo sapiens[Organism] AND refseq_select[filter]` and its summary.
- `esearch_pubmed_osimertinib.xml`, `esummary_pubmed_osimertinib.xml` — PubMed search `osimertinib AND EGFR AND resistance AND review[pt]` (top 10 by relevance) and their summaries.

**Source:** NCBI RefSeq and PubMed, National Center for Biotechnology Information, U.S. National Library of Medicine.
See NCBI's policies at https://www.ncbi.nlm.nih.gov/home/about/policies/. **Cite:** Sayers, E. W. *et al.*
Database resources of the National Center for Biotechnology Information in 2025. *Nucleic Acids Res.* **2025**, *53*, D20–D29.
https://doi.org/10.1093/nar/gkae979

## `blast_cache/` — BLAST (notebooks 08, 12)

- `EGFR_kinase_domain_query.fasta` — the query: EGFR residues 712–979 (UniProt kinase domain).
- `EGFR_kinase_blastp_swissprot.xml` — BLAST XML (outfmt 5) for `blastp` against Swiss-Prot, E-value ≤ 1e-10,
  50 best hits, BLOSUM62.

The remote NCBI web BLAST queue (`NCBIWWW.qblast`) did not return this job within 40 minutes when the cache was built,
so the cached result was produced with the stand-alone **BLAST+ 2.17.0** `blastp` against NCBI's own preformatted
`swissprot` database (https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz; database date in `provenance.json`).
The XML format is identical to a web BLAST result. Re-run with `python scripts/refresh_cache.py blast` (remote) or
`blast-local` (BLAST+).

**Cite:** Camacho, C. *et al.* BLAST+: architecture and applications. *BMC Bioinformatics* **2009**, *10*, 421.
https://doi.org/10.1186/1471-2105-10-421

## `alignments/` — multiple sequence alignment (notebooks 07, 09, 11, 12)

- `kinase_domains_famsa.aln` (Clustal format) and `kinase_domains_famsa.fasta` — alignment of the 15 sequences in
  `sequences/kinase_domains.fasta`, computed with FAMSA (via `pyfamsa`, UPGMA guide tree). Notebook 07 recomputes it
  and checks that the result is identical.

**Cite:** Deorowicz, S., Debudaj-Grabysz, A. & Gudyś, A. FAMSA: Fast and accurate multiple sequence alignment of huge
protein families. *Sci. Rep.* **2016**, *6*, 33964. https://doi.org/10.1038/srep33964

## `structures/` — Protein Data Bank (notebooks 10, 12)

- `1M17.cif`, `1M17.pdb` — PDB entry **1M17**, EGFR tyrosine kinase domain in complex with erlotinib (ligand AQ4),
  X-ray, 2.6 Å, downloaded unmodified from RCSB (https://files.rcsb.org/download/). Notebooks 10 and 12 download the
  entry live and fall back to these files offline.

**Cite the structure:** Stamos, J., Sliwkowski, M. X. & Eigenbrot, C. Structure of the epidermal growth factor receptor
kinase domain alone and in complex with a 4-anilinoquinazoline inhibitor. *J. Biol. Chem.* **2002**, *277*, 46265–46272.
https://doi.org/10.1074/jbc.M207135200
**Cite the PDB:** Berman, H. M. *et al.* The Protein Data Bank. *Nucleic Acids Res.* **2000**, *28*, 235–242.
https://doi.org/10.1093/nar/28.1.235. PDB data are freely available under CC0 1.0 (https://www.rcsb.org/pages/usage-policy).
