# intron-retention-scripts

*under active development*

Component scripts of an intron retention analysis toolkit. These scripts will be integrated as addons for intron retention analysis a snakemake based rnaseq worklfow. 

## Usage

1. `gtf_to_db.py`

converts a given gtf file to an sqlite database using gffutils. 

```bash
python gtf_to_db.py --gtf /path/to/input.gtf  --gtfdb /path/to/output_gtf.db
```

2. `filter_gtf.py`

optionally filters a given gtf file to remove transcripts of a particular type e.g. `retained_intron`

```bash
python filter_gtf.py --gtf /path/to/gtf.db  --transtype retained_intron -out /path/to/output.txt
```

3. `get_introns.py`

get a list of clean(starts not overlapping with exons) introns from a gtf db

```bash
python get_introns.py --gtf /path/to/gtf.db  --out /path/to/introns.bed
```

4. `get_intron_type.py`

classify each intron into U2/U12 type using PWM's from [splicerack](http://katahdin.mssm.edu/splice/index.cgi?database=spliceNew)

```bash
python get_intron_type.py --bed /path/to/introns.bed  --branch /path/to/branch.pwm --don /path/to/don.pwm --genome /path/to/genome.fa --out /path/to/intronType.xls
```

5. `get_msi.py`

calculate mis splicing index(MSI) for a given sample bam file for introns supplied in a bed format

```bash
python --bed /path/to/introns.bed  --bam /path/to/reads.bam --out /path/to/msi.xls
```

