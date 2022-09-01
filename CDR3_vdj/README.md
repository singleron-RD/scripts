Convert `{sample}_cell_confident.tsv` to `filtered_contig_annotations.csv` which can be used as input of scRepertoire.
https://ncborcherding.github.io/vignettes/vignette.html


```
perl filter.vdjfile.1.pl \
 -inputlist {sample}_cell_confident.tsv \
 -namelist {sample} \
 -VDJtype TCR \
 -Datatype sgr \
 -outdir ./ \
 -chain_pair  NO \
 -prefix {sample}
```
