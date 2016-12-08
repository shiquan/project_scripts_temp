
HLA typing from Sanger sequencing ABI file
------------------------------------------



Recently I got an AB1 file which is the Sanger sequence result of the HLA-A gene (partly). I tried to use this file to find out the exactly HLA-A haplotypes, but I didn't find any open sources or free software to help me to do typing. So I decided to do everything from scratch. Here is the full practice. And of course you could copy and modify these steps freely for any purpose.



**Design**

First of all, we should convert the Sanger sequence signals in ab1 file into a plain sequence, because it is much easier to analysis a sequence than using an ab1 file.
The difference among MHC haplotypes is used to do typing and consist of wide ranges of genetic variantions. My protocol is to find these variantions and compare these differences with the HLA-A gene databases, to find the best hit(s).


**Convert AB1 to Fasta**

Most converter with GUI or batch mode could convert the AB1 file into fasta or fastq format, however I didn't find any program could also keep the variantions as well as they only export the best choice of each position.  I decided to use IUPAC notations instead of only use (ATCG)s to contruct the DNA sequence. I wrote a program named [ABI_convertor](https://github.com/shiquan/small_projects_collections/blob/master/projects/abi/ab1_convert.c) could handle this well. Here is the command:



> ABI_convertor HLA_seq.ab1 −o seq.txt 



This is the orignial AB1 file.

![](https://github.com/shiquan/small_projects_collections/blob/master/projects/hla_typing/sanger_ab1_demo.png)

And this is the converted sequence. Remeber to trim the first and last few bases because for sanger sequencing the tails of the sequence usually come along with low quality and bias. And try to google IUPAC for more details if you don't fullly understand what's the difference between [WYKSR] and [ATCG]s.

>\>seq
>
>GGGACGAGGAGACASGGAAWGTGAAGGCCCACTCACAGAYTGACCGAGWG
>RACCTGSGGAYCSTGCKCGGCTACTACAACCAGAGCGAGGCCGGTGAGTG
>ACCCCRGCCCGGGGCGCAGGTCACGACCTCTCATCCCCCACGGACGGGCC
>RGGTCRCCCACAGTCTCCGGGTCCGAGATCCACCCCGAAGCCGCGGGACC
>CCGAGACCCTTGCCCCGGGAGAGGCCCAGGCGCCTTWACCCGGTTTCATT
>TTCAGTTTAGGCCAAAAATCCCCCCGGGTTGGTCGGGGCCGGACGGGGCT
>CGGGGGACTGGGCTGACCGYGGGGTCGGGGCCAGGTTCTCACACCMTCCA
>GATGATGTATGGCTGCGACGTGGGGTCGGACGGGCGCTTCCTCCGCGGGT
>ACCASCAGKACGCCTACGACGGCAAGGATTACATCG



**align the sequence against database**

Before to align the sequence against the HLA databses, we should download the databases from [ebi](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/). And please notice that ebi database only give you haplotypes of genes in  MHC region. But we all know that human being is diplotype. So we need reconstruct the diplotype database base on the haplotype sequences. Here I worte a program named [database_construct](https://github.com/shiquan/small_projects_collections/blob/master/projects/hla_typing/database_construct.c).

Download the hapotype databases, here we only download HLA-A gene, for other genes try to change the utl.

> wget −c ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/A_gen.txt	

Before to reconstruct the databases, try to use mafft to realign all the subtypes with each other, to find the most common regions.

> mafft A_gen.fasta | database_construct > database.fa		

Now align the Sanger	sequence against our database. Remember there are IUPAC code in the sequence, so we might not use the current stat NGS aligner like bwa and bowtie2 etc. Here we use blat in protein alignment mode in conversion, but in long term or big projects we should consider to develop a program to support IUPAC. Here is the command.

> blat −out=blast −t=prot −q=prot −stepSize=5 −repMatch=2253 −minMatch=11 −minScore=90 −minIdentity=0.99 database.fa seq.txt /dev/stdout > result.txt	



***Results***

The result file show several possible genotypes. Now we just pick the best hits. Because this fragement is short we may not find the exactly genotypes but the result should be improved if we have more fragement and sequences.



>Reference:  Kent, WJ. (2002) BLAT - The BLAST-like alignment tool
>
>...
>
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA08602_A*24:215_3116_bp           1085   0.0
>HLA:HLA12600_A*31:01:02:02_3000_bp,HLA:HLA06368_A*23:38N_3020_bp     1085   0.0
>HLA:HLA02651_A*31:14N_3090_bp,HLA:HLA14804_A*24:03:01:02_3072_bp     1085   0.0
>HLA:HLA12600_A*31:01:02:02_3000_bp,HLA:HLA08602_A*24:215_3116_bp     1085   0.0
>HLA:HLA12600_A*31:01:02:02_3000_bp,HLA:HLA11397_A*24:280_3072_bp     1085   0.0
>HLA:HLA12600_A*31:01:02:02_3000_bp,HLA:HLA14804_A*24:03:01:02_3072_bp  1085   0.0
>HLA:HLA04466_A*31:01:04_2918_bp,HLA:HLA05800_A*24:152_3176_bp        1085   0.0
>HLA:HLA04466_A*31:01:04_2918_bp,HLA:HLA06368_A*23:38N_3020_bp        1085   0.0
>HLA:HLA04466_A*31:01:04_2918_bp,HLA:HLA08602_A*24:215_3116_bp        1085   0.0
>HLA:HLA04466_A*31:01:04_2918_bp,HLA:HLA11397_A*24:280_3072_bp        1085   0.0
>HLA:HLA04466_A*31:01:04_2918_bp,HLA:HLA14804_A*24:03:01:02_3072_bp   1085   0.0
>HLA:HLA02651_A*31:14N_3090_bp,HLA:HLA11397_A*24:280_3072_bp          1085   0.0
>HLA:HLA02651_A*31:14N_3090_bp,HLA:HLA08602_A*24:215_3116_bp          1085   0.0
>HLA:HLA12433_A*31:01:24_2918_bp,HLA:HLA14804_A*24:03:01:02_3072_bp   1085   0.0
>HLA:HLA02651_A*31:14N_3090_bp,HLA:HLA06368_A*23:38N_3020_bp          1085   0.0
>HLA:HLA12433_A*31:01:24_2918_bp,HLA:HLA11397_A*24:280_3072_bp        1085   0.0
>HLA:HLA12433_A*31:01:24_2918_bp,HLA:HLA08602_A*24:215_3116_bp        1085   0.0
>HLA:HLA12433_A*31:01:24_2918_bp,HLA:HLA06368_A*23:38N_3020_bp        1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA03183_A*23:19N_3105_bp           1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA05800_A*24:152_3176_bp           1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA06112_A*24:02:01:03_3075_bp      1085   0.0
>HLA:HLA12433_A*31:01:24_2918_bp,HLA:HLA05800_A*24:152_3176_bp        1085   0.0
>HLA:HLA12600_A*31:01:02:02_3000_bp,HLA:HLA05800_A*24:152_3176_bp     1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA14895_A*23:01:19_2943_bp         1085   0.0
>HLA:HLA02651_A*31:14N_3090_bp,HLA:HLA05800_A*24:152_3176_bp          1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA14804_A*24:03:01:02_3072_bp      1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA14803_A*24:02:01:08_2902_bp      1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA14802_A*24:02:01:07_2902_bp      1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA14800_A*24:02:01:06_2902_bp      1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA13786_A*24:02:01:05_2979_bp      1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA12634_A*24:293_2902_bp           1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA06368_A*23:38N_3020_bp           1085   0.0
>HLA:HLA05799_A*31:46_3075_bp,HLA:HLA11397_A*24:280_3072_bp           1085   0.0



And let's just pick the first several genotypes and find the mismatches, and we could see this mismatches come from the bias of ABI_convertor, but I think use imporved algorithm and more sequences will fix this problem.

![](https://github.com/shiquan/small_projects_collections/blob/master/projects/hla_typing/sanger_align_demo.png)



Recently, I tried to genotype the next generation sequencing data, and I think I could just write a program to convert the vcf file to fasta sequences and apply this pipeline to do it directly. I will update this protocol soon, please let me know if you have any comments or questions.