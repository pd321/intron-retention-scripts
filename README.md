# intron-retention-scripts

**Not under active development**

Component scripts of a snakemake based pipeline for intron retention analysis.

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

get a list of clean(starts not overlapping with exons) introns from a gtf db. The `data` directory stores intron output bed files for mouse (Gencode M21) and human (Gencode 30).

```bash
python get_introns.py --gtf /path/to/gtf.db  --out /path/to/introns.bed
```

4. `get_intron_type.py`

classify each intron into U2/U12 type using PWM's from [splicerack](http://katahdin.mssm.edu/splice/index.cgi?database=spliceNew). The `data` directory stores gzipped intron type output xls files for mouse (Gencode M21) and human (Gencode 30).

```bash
python get_intron_type.py --bed /path/to/introns.bed  --branch /path/to/branch.pwm --don /path/to/don.pwm --genome /path/to/genome.fa --out /path/to/intronType.xls
```

5. `get_msi.py`

calculate mis splicing index(MSI) for a given sample bam file for introns supplied in a bed format

```bash
python get_msi.py --bed /path/to/introns.bed  --bam /path/to/reads.bam --out /path/to/msi.xls
```

6. `deltaMSI.R`

calculate deltaMSI values for when given a treatment and control MSI output file from `get_msi.py`

```bash
Rscript --vanilla deltaMSI.R --trt /path/to/trt_msi.xls --cnt /path/to/cnt_msi.xls --trtname trt --cntname cnt --out /path/to/deltaMSI.xls
```

7. `addAnnotation.R`

Add gene and intron type annotation to the raw dmsi file produced above. The data directory contains gene info files for human (hg38, hg19) and mouse (mm10) genomes.

```bash
Rscript --vanilla addAnnotation.R --gene /path/to/introns.bed --type /path/to/intronType.xls --geneinfo /path/to/geneinfo.xls --deltamsi /path/to/raw_dmsi.xls --out /path/to/processed_dmsi.xls
```

8. `plotdeltaMSI.R`

plot the deltaMSI values for a given comparison as a dot/density plot post removing outliers.

```bash
Rscript --vanilla plotdeltaMSI.R --deltaMSI /path/to/deltaMSI.xls --outdot /path/to/dotplot.pdf --outdensity /path/to/densityplot.pdf
```