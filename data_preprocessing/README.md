# Data preprocessing

The following steps were used to preprocess raw ONT data:  

Concatenate raw files
```
zcat ${dir}/fastq_pass/*.fastq.gz | gzip -c > ${prefix}.raw.fastq.gz
```
Filter concatenated files by read length ([NanoFilt](https://github.com/wdecoster/nanofilt)), visualize reads ([NanoPlot](https://github.com/wdecoster/NanoPlot))
```
gunzip -c ${prefix}.raw.fastq.gz | NanoFilt -l 40000 | gzip > ${prefix}.40kbp.raw.fastq.gz
NanoPlot --huge --outdir ${prefix}.40kbp_nanoplot --prefix ${prefix}.40kbp_ --fastq ${prefix}.40kbp.raw.fastq.gz
```
HERRO correction (requires GPU) [Dorado](https://github.com/nanoporetech/dorado/?tab=readme-ov-file#read-error-correction)
```
gunzip -c ${prefix}.40kbp.raw.fastq.gz > ${prefix}.40kbp.raw.fastq
dorado correct ${prefix}.40kbp.raw.fastq > ${prefix}.40kbp.corrected.fastq
gzip ${prefix}.40kbp.corrected.fastq
```
