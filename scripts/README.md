## Assembly Concordance script

### Example Usage
```
$ python assembly_concordance.py -v HG002.PAW71238.haplotagged.strkit.vcf -a hg002_mat_asm5_hg38.sorted.bam hg002_pat_asm5_hg38.sorted.bam -t strkit -bed test-isolated-vc-catalog.strkit.bed
```

### Example Output

```
#catalog_locus  motif   vcf_locus       coord_change    length_difference       tool_genotype   assembly_genotype       tool    status
chr1-103260-103280      TAAAA   chr1-103259-103280      CHANGED NA      NA      NA      strkit  MISSING
chr1-434350-434825      GACAAACACGTGGATACATGGAGGGGAACAACACACACCAGGGCCTCTCAGGGGGACAGGGGGTAGGAGACCATCAG   chr1-434349-434889      CHANGED NA      NA      NA      strkit  MISSING
chr1-862076-862088      GATG    chr1-862068-862088      CHANGED 0,0     20,20   20,20   strkit  MATCH
chr1-905012-905114      GCGGGGGAGGCTGTTGGGGACGTTCGTG    chr1-905011-905120      CHANGED 0,0     109,109 109,109 strkit  MATCH
chr1-912714-912732      GCCCCTCCT       chr1-912713-912732      CHANGED 0,0     19,19   19,19   strkit  MATCH
chr1-930081-930152      GCCCCACCTTCCTCTCCTCCT   chr1-930080-930152      CHANGED 0,0     72,72   72,72   strkit  MATCH
chr1-969958-970156      CAGGTTCTGTGCCTGTGGGGGGCTGTTCTTCACGTATGTGTTGTGTGGCTGTGCA chr1-969957-970178      CHANGED 0,0     221,221 221,221 strkit  MATCH
chr1-1091706-1091782    GGGGCTGTGGTGGAGGGTGGGGCCAAATGGAAGTGGGC  chr1-1091705-1091821    CHANGED 0,0     116,116 116,116 strkit  MATCH
chr1-1141113-1141137    CCCCATCCCCAC    chr1-1141112-1141162    CHANGED 0,0     50,50   50,50   strkit  MATCH
chr1-1158353-1158375    GGGCTGG chr1-1158352-1158375    CHANGED 0,0     23,23   23,23   strkit  MATCH
chr1-1168009-1168125    GGGCGGAGGGCCGAGCGGGGCCAGCAGACGGGTGA     chr1-1168008-1168147    CHANGED 0,0     104,139 104,139 strkit  MATCH
chr1-1221641-1221664    AATTAAA chr1-1221640-1221664    CHANGED 0,0     24,24   24,24   strkit  MATCH
#Total regions in VCF: 12038
#Regions missed (no assembly genotype): 175
#Regions called (with assembly genotype): 11863
#Matching genotypes: 11239
#One-off genotypes: 386
#Concordance: 94.74 %
```