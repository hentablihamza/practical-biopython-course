"""Rebuild the cached data in data/ from the original public sources.

The notebooks query NCBI, UniProt, and RCSB live when they can and fall back to these
cached files when they cannot (or when COURSE_OFFLINE=1). Run this script to refresh them:

    NCBI_EMAIL=you@university.edu python scripts/refresh_cache.py            # everything
    NCBI_EMAIL=you@university.edu python scripts/refresh_cache.py blast      # one step

Steps: uniprot, entrez, blast, blast-local, msa, pdb (default: all except blast-local).
Remote BLAST (`blast`) can wait in NCBI's queue for a long time; `blast-local` runs the same
search with BLAST+ and NCBI's preformatted swissprot database. The msa step needs `pip install pyfamsa`.
"""
import datetime as dt
import json
import os
import sys
import time
from pathlib import Path

import requests
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW

DATA = Path(__file__).resolve().parents[1] / "data"
TODAY = dt.date.today().isoformat()
Entrez.email = os.environ.get("NCBI_EMAIL", "your.name@example.org")
Entrez.tool = "practical-biopython-course"

# UniProt accessions: the human ErbB family, EGFR orthologues, and other drug-target kinases.
PANEL = ["P00533", "P04626", "P21860", "Q15303", "Q01279", "P04412", "P24348",
         "P00519", "P12931", "P15056", "P35968", "Q06187", "P24941", "P08581", "Q9UM73"]


def meta_update(key, value):
    f = DATA / "provenance.json"
    meta = json.loads(f.read_text()) if f.exists() else {}
    meta[key] = value
    f.write_text(json.dumps(meta, indent=2) + "\n")


def uniprot():
    seqdir = DATA / "sequences"
    seqdir.mkdir(exist_ok=True)
    full, doms, rows = [], [], []
    for acc in PANEL:
        r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=30)
        r.raise_for_status()
        d = r.json()
        gene = d["genes"][0]["geneName"]["value"]
        org = d["organism"]["scientificName"]
        seq = d["sequence"]["value"]
        kd = [f for f in d["features"] if f["type"] == "Domain" and f.get("description") == "Protein kinase"][0]
        s, e = kd["location"]["start"]["value"], kd["location"]["end"]["value"]
        name = d["proteinDescription"]["recommendedName"]["fullName"]["value"]
        full.append(f">{acc}|{d['uniProtkbId']} {name} OS={org} GN={gene}\n{seq}\n")
        doms.append(f">{d['uniProtkbId']} {acc} kinase domain {s}-{e} OS={org} GN={gene}\n{seq[s-1:e]}\n")
        rows.append(f"{acc},{d['uniProtkbId']},{gene},\"{name}\",{org},{len(seq)},{s},{e},"
                    f"{d['entryAudit']['entryVersion']}")
        if acc == "P00533":
            keep = [f for f in d["features"] if f["type"] in ("Signal", "Topological domain", "Transmembrane", "Domain",
                                                         "Region", "Binding site", "Active site", "Glycosylation")]
            (seqdir / "P00533_uniprot_features.json").write_text(json.dumps(
                {"accession": acc, "entryAudit": d["entryAudit"], "features": keep}, indent=1) + "\n")
        time.sleep(0.3)
    (seqdir / "kinase_panel_full.fasta").write_text("".join(full))
    (seqdir / "kinase_domains.fasta").write_text("".join(doms))
    (seqdir / "kinase_panel.csv").write_text(
        "accession,entry_name,gene,protein_name,organism,length,kinase_start,kinase_end,uniprot_entry_version\n"
        + "\n".join(rows) + "\n")
    r = requests.get("https://rest.uniprot.org/uniprotkb/P00533.fasta", timeout=30)
    r.raise_for_status()
    (seqdir / "P00533_EGFR_human.fasta").write_text(r.text)
    meta_update("uniprot", {"retrieved": TODAY, "source": "https://rest.uniprot.org/uniprotkb/", "accessions": PANEL})


def entrez():
    c = DATA / "entrez_cache"
    c.mkdir(exist_ok=True)
    term = "EGFR[Gene Name] AND Homo sapiens[Organism] AND refseq_select[filter]"
    with Entrez.esearch(db="nucleotide", term=term, retmax=20) as h:
        (c / "esearch_nucleotide_EGFR.xml").write_bytes(h.read())
    time.sleep(0.5)
    with Entrez.efetch(db="nucleotide", id="NM_005228", rettype="gb", retmode="text") as h:
        (c / "NM_005228.gb").write_text(h.read())
    time.sleep(0.5)
    with Entrez.efetch(db="nucleotide", id="NM_005228", rettype="fasta", retmode="text") as h:
        (c / "NM_005228.fasta").write_text(h.read())
    time.sleep(0.5)
    with Entrez.esummary(db="nucleotide", id="NM_005228") as h:
        (c / "esummary_NM_005228.xml").write_bytes(h.read())
    time.sleep(0.5)
    pterm = "osimertinib AND EGFR AND resistance AND review[pt]"
    with Entrez.esearch(db="pubmed", term=pterm, retmax=10, sort="relevance") as h:
        raw = h.read()
    (c / "esearch_pubmed_osimertinib.xml").write_bytes(raw)
    time.sleep(0.5)
    ids = Entrez.read(Entrez.esearch(db="pubmed", term=pterm, retmax=10, sort="relevance"))["IdList"]
    time.sleep(0.5)
    with Entrez.esummary(db="pubmed", id=",".join(ids)) as h:
        (c / "esummary_pubmed_osimertinib.xml").write_bytes(h.read())
    meta_update("entrez", {"retrieved": TODAY, "nucleotide_term": term, "pubmed_term": pterm,
                           "records": ["NM_005228"]})


def _blast_query():
    c = DATA / "blast_cache"
    c.mkdir(exist_ok=True)
    egfr = SeqIO.read(DATA / "sequences" / "P00533_EGFR_human.fasta", "fasta")
    query = egfr.seq[711:979]  # UniProt kinase domain 712-979
    (c / "EGFR_kinase_domain_query.fasta").write_text(f">EGFR_HUMAN_kinase_712-979\n{query}\n")
    return c, query


def blast():
    """Remote BLAST at NCBI (can wait in the queue for a long time)."""
    c, query = _blast_query()
    t0 = time.time()
    h = NCBIWWW.qblast("blastp", "swissprot", str(query), hitlist_size=50, expect=1e-10)
    (c / "EGFR_kinase_blastp_swissprot.xml").write_text(h.read())
    meta_update("blast", {"retrieved": TODAY, "how": "NCBIWWW.qblast (remote NCBI BLAST)", "program": "blastp",
                          "database": "swissprot", "query": "P00533 residues 712-979", "hitlist_size": 50,
                          "expect": 1e-10, "runtime_s": round(time.time() - t0)})


def blast_local():
    """Same search with a local BLAST+ install and NCBI's preformatted swissprot database.

    Needs `blastp` on PATH and BLASTDB pointing to a folder with swissprot.* from
    https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
    """
    import subprocess
    c, _ = _blast_query()
    out = c / "EGFR_kinase_blastp_swissprot.xml"
    subprocess.run(["blastp", "-query", str(c / "EGFR_kinase_domain_query.fasta"), "-db", "swissprot",
                    "-evalue", "1e-10", "-max_target_seqs", "50", "-matrix", "BLOSUM62", "-outfmt", "5",
                    "-out", str(out)], check=True)
    version = subprocess.run(["blastp", "-version"], capture_output=True, text=True).stdout.splitlines()[0]
    info = subprocess.run(["blastdbcmd", "-db", "swissprot", "-info"], capture_output=True, text=True).stdout
    db_date = next((l.split("\t")[0].replace("Date: ", "") for l in info.splitlines() if l.startswith("Date:")), "?")
    meta_update("blast", {"retrieved": TODAY, "how": f"local {version} against NCBI's preformatted swissprot database",
                          "database": "swissprot (https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz)",
                          "database_date": db_date, "program": "blastp", "query": "P00533 residues 712-979",
                          "max_target_seqs": 50, "expect": 1e-10, "matrix": "BLOSUM62"})


def msa():
    """Align the kinase domains locally with FAMSA (pip install pyfamsa)."""
    import pyfamsa
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import AlignIO

    recs = list(SeqIO.parse(DATA / "sequences" / "kinase_domains.fasta", "fasta"))
    seqs = [pyfamsa.Sequence(r.id.encode(), str(r.seq).encode()) for r in recs]
    aln = pyfamsa.Aligner(guide_tree="upgma").align(seqs)
    order = {r.id: i for i, r in enumerate(recs)}
    rows = sorted(aln, key=lambda s: order[s.id.decode()])
    msa_obj = MultipleSeqAlignment([SeqRecord(Seq(s.sequence.decode()), id=s.id.decode(), description="")
                                    for s in rows])
    (DATA / "alignments").mkdir(exist_ok=True)
    AlignIO.write(msa_obj, DATA / "alignments" / "kinase_domains_famsa.aln", "clustal")
    AlignIO.write(msa_obj, DATA / "alignments" / "kinase_domains_famsa.fasta", "fasta")
    meta_update("msa", {"retrieved": TODAY, "tool": f"FAMSA via pyfamsa {pyfamsa.__version__}, UPGMA guide tree",
                        "input": "data/sequences/kinase_domains.fasta",
                        "outputs": ["kinase_domains_famsa.aln", "kinase_domains_famsa.fasta"]})


def pdb():
    s = DATA / "structures"
    s.mkdir(exist_ok=True)
    for ext in ("pdb", "cif"):
        r = requests.get(f"https://files.rcsb.org/download/1M17.{ext}", timeout=60)
        r.raise_for_status()
        (s / f"1M17.{ext}").write_text(r.text)
    meta_update("pdb", {"retrieved": TODAY, "source": "https://files.rcsb.org/download/", "ids": ["1M17"]})


STEPS = {"uniprot": uniprot, "entrez": entrez, "blast": blast, "blast-local": blast_local, "msa": msa, "pdb": pdb}

if __name__ == "__main__":
    for name in sys.argv[1:] or [s for s in STEPS if s != "blast-local"]:
        print("refreshing", name, "...", flush=True)
        STEPS[name]()
    print("done")
