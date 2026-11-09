# Reference data

The databases are distributed separately (size and redistribution licensing).
Create the layout below and point `server/configs/config2.yaml` (or the GUI) at
it.

```
data/
├── virus_taxids_high_pathogenic.txt              # target taxids, one per line
├── hvn/                                          # high pathogenic virus, nucleotide
│   ├── hvn.fa                                    # BLASTN database source
│   └── hvn.fa.{nhr,nin,nsq}                      # makeblastdb output
├── hvp/                                          # high pathogenic virus, protein
│   └── hvp.dmnd                                  # DIAMOND pre-screen database
└── databases/
    └── nr.dmnd                                   # NCBI nr (comprehensive pass)
host_idx/
└── host_idx.{1..4,rev.1..4}.bt2                   # bowtie2 host index
```

## 1. Target taxids

`virus_taxids_high_pathogenic.txt` lists the taxids hvFinder reports. One taxid per
line; blank lines and lines starting with `#` are ignored. Genus-level taxids are
allowed when a whole genus is in scope. Resolve names to taxids with the NCBI
Taxonomy web search or the E-utilities API.

## 2. `hvn` — high pathogenic virus nucleotide database (BLASTN)

`hvn` (hv-**n**ucleotide) is a curated nucleotide database of highly pathogenic
viruses. `hvn.fa` holds the reference sequences; build the BLAST database next
to it:

```bash
makeblastdb -in data/hvn/hvn.fa -dbtype nucl -out data/hvn/hvn.fa
```

The FASTA header order/format matters for reporting — keep one accession per
record with a stable identifier.

## 3. `hvp` — high pathogenic virus protein database (`hvp.dmnd`)

`hvp` (hv-**p**rotein) is the DIAMOND database used by the pre-screen. Build it with
**taxonomy attached**, so that `staxids` can be filtered:

```bash
diamond makedb --in hvp.faa --db data/hvp/hvp.dmnd --threads 16 \
               --taxonmap prot.accession2taxid \
               --taxonnodes nodes.dmp --taxonnames names.dmp
```

`nodes.dmp` / `names.dmp` / `prot.accession2taxid` come from the NCBI Taxonomy and
accession-mapping distributions. Record the taxdump release date alongside the
database — it determines how taxids are resolved.

## 4. Comprehensive database (`nr.dmnd`) — NCBI nr

```bash
diamond makedb --in nr.faa --db nr.dmnd --threads 16 --taxonmap ... --taxonnodes ... --taxonnames ...
```

`nr.dmnd` is large (hundreds of GB); the detection stage spends most of its wall
clock loading it, which is why the pre-screen exists.

## 5. Host index

Any bowtie2 index of the host genome works; pass its **prefix** as
`bowtie2.params.index_path`.

```bash
bowtie2-build host.fa host_idx/host_idx
```
