
HLA typing from Sanger sequencing ABI file
------------------------------------------



Recently I got a AB1 file which is the sequence result of partly HLA gene. And tried to use this file to genotype this sample, but I didn't find any open sources or free software to help me to do that. So I decide to do this from scratch. And of course you could copy and modify these steps freely for any purpose.



**Design**

First of all we should convert the sequence signals in ab1 file into a sequence file, so we could use current-stat-of-art programs to align and call variants and genotype. And second, because there are  a lot of variants usually saps in this region, so we much keep this information and so we could use it for typing.  Last we need comparing out sequence with HLA database, to find the best hit(s). Here are the details.



**Convert AB1 to Fasta**

Most converter with GUI or batch mode could convert the AB1 file into fasta or fastq format, however I didn't find any program could also keep the variantions as well as they only export the best choice of each position.  I decided to use IUPAC notations instead of only use (ATCG)s to contruct the DNA sequence. I wrote a program named [ABI_convertor](https://github.com/shiquan/small_projects_collections/blob/master/projects/abi/ab1_convert.c) could handle this well. Here is the command:



> ABI_convertor HLA_seq.ab1 −o seq.txt 



![](https://github.com/shiquan/small_projects_collections/blob/master/projects/hla_typing/sanger_ab1_demo.png)



**align the sequence against database**

Before to align the sequence against the HLA databses, we should download the databases from [ebi](ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/). And please notice that ebi database only give you haplotypes of genes in  MHC region. But we all know that human being is diplotype. So we need reconstruct the diplotype database base on the haplotype sequences. Here I worte a program named [database_construct](https://github.com/shiquan/small_projects_collections/blob/master/projects/hla_typing/database_construct.c).

Download the hapotype databases, here we only download HLA-A gene, for other genes try to change the utl.

> wget −c ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/A_gen.txt	

Before to reconstruct the databases, try to use mafft to realign all the subtypes with each other, to find the most common regions.

> mafft A_gen.fasta | database_construct > database.fa		

Now align the Sanger	sequence against our database. Remember there are IUPAC code in the sequence, so we might not use the current stat NGS aligner like bwa and bowtie2 etc. Here we use blat in protein alignment mode in conversion, but in long term or big projects we should consider to develop a program to support IUPAC. Here is the command.

> blat −out=blast −t=prot −q=prot −stepSize=5 −repMatch=2253 −minMatch=11 −minScore=90 −minIdentity=0.99 database.fa seq.txt /dev/stdout > result.txt	




​				
​			
​		
​	







