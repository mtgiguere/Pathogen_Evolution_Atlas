from __future__ import annotations

import argparse
import os
from pathlib import Path

from Bio import Entrez


def build_query(min_len: int, max_len: int) -> str:
    # SARS-CoV-2 + "complete genome" + length filter
    return (
        '"Severe acute respiratory syndrome coronavirus 2"[Organism] '
        'AND "complete genome"[Title] '
        f'AND {min_len}:{max_len}[SLEN]'
    )


def fetch_accessions(query: str, email: str, retmax: int) -> list[str]:
    Entrez.email = email
    with Entrez.esearch(db="nuccore", term=query, retmax=retmax) as handle:
        res = Entrez.read(handle)
    # IdList here are GI/UIDs; we ask Entrez to translate to accessions via efetch later,
    # but in practice nuccore IdList works fine with efetch rettype=acc.
    ids = res.get("IdList", [])
    if not ids:
        return []

    with Entrez.efetch(db="nuccore", id=",".join(ids), rettype="acc", retmode="text") as handle:
        text = handle.read()

    accessions = [line.strip() for line in text.splitlines() if line.strip()]
    # De-dup while preserving order
    seen = set()
    out: list[str] = []
    for a in accessions:
        if a not in seen:
            seen.add(a)
            out.append(a)
    return out


def write_accessions(accessions: list[str], out_path: str | Path) -> None:
    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(accessions) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description="Build a SARS-CoV-2 complete-genome accession list from NCBI.")
    ap.add_argument("--out", default="data/accessions/covid_accessions.txt", help="Output .txt path")
    ap.add_argument("--n", type=int, default=300, help="How many accessions to write")
    ap.add_argument("--min-len", type=int, default=29000, help="Minimum sequence length (bp)")
    ap.add_argument("--max-len", type=int, default=31000, help="Maximum sequence length (bp)")
    args = ap.parse_args()

    email = os.getenv("NCBI_EMAIL")
    if not email:
        raise SystemExit("NCBI_EMAIL environment variable not set (required by NCBI).")

    query = build_query(args.min_len, args.max_len)
    accessions = fetch_accessions(query=query, email=email, retmax=args.n)

    write_accessions(accessions, args.out)
    print(f"Wrote {len(accessions)} accessions -> {args.out}")
    print(f"Query: {query}")


if __name__ == "__main__":
    main()
