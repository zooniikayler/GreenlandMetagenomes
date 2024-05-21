### Commands for Greenland Metagenomes Thesis


### Versions
# trimmomatic v0.32
# bwa 0.7.17
# samtools 1.13-30
# megahit 1.2.9
# metabat 2.15
# maxbin 2.2.7
# dastools 1.1.4
# checkm 1.2.1
# prodigal 2.6.3
# gtdbtk 2.1.0
# kegg decoder 1.3

###Step 1 Raw
##RC Run time <4 hours - hmem 24GB

m=sub-glacial-Leverett
n=sub-glacial-Russell-East
x=sub-glacial-Russell-West
v=sub-glacial-Kangaussarssup
a=sub-glacial

for i in {$m,$n,$x,$v}; do trimmomatic PE -threads 8 -phred33 -trimlog "$i"-trimlog "$i"_rawdata/"$i"-R1-raw.fastq "$i"_rawdata/"$i"-R2-raw.fastq -baseout "$i"-trimmed.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 MINLEN:40; 
mkdir "$i"_trimmed; mv "$i"-trimmed* "$i"_trimmed; done

cat "$m"_trimmed/*R1P "$n"_trimmed/*R1P "$x"_trimmed/*R1P "$v"_trimmed/*R1P > "$a"-trimmed-R1P
cat "$m"_trimmed/*R2P "$n"_trimmed/*R2P "$x"_trimmed/*R2P "$v"_trimmed/*R2P > "$a"-trimmed-R2P

###Step 2 Assembly - Run time 2 days - cpu 100GB - hmem

megahit -1 "$a"-trimmed-R1P -2 "$a"-trimmed-R2P -t 24 -m 200000000 --presets meta-large -o "$a"-megahit-assembly

###Step 2.2 - Mapping
mkdir mapping/
cp "$a"-megahit_assembly/final.contigs.fa  mapping/"$a"-final_contigs.fa
##Create index for contigs FASTA
bwa index mapping/"$a"-final_contigs.fa
##Align reads against contigs and sort BAM 
bwa mem -t 16 mapping/"$a"-final_contigs.fa "$a"-trimmed-R1P "$a"-trimmed-R2P > mapping/"$a"-aln.sam | samtools sort -o mapping/"$a"-aln.bam  
for i in {$m,$n,$x,$v}; do bwa mem -t 16 mapping/"$a"-final_contigs.fa "$i"_trimmed/*R1P "$i"_trimmed/*R2P > mapping/"$i"-aln.sam | samtools sort -o mapping/"$i"-aln.bam
##index BAM file
for i in {$m,$n,$x,$v,$a}; do samtools index mapping/"$i"-aln.bam

###Step 3 Contigs to MAGs

##Binning
## 3.1 Metabat 
runMetaBat.sh "$a"-megahit-assembly/"$a"-final_contigs.fa mapping/"$a"-aln.bam -o "$a"-metabat
mkdir "$a"-metabat_bins
mv "$a"-metabat* "$a"-metabat_bins/

## 3.2 Maxbin
run_MaxBin.pl -contig mapping/"$a"-final_contigs.fa -out "$a"-maxbin -reads "$a"-trimmed-R1P -reads2 "$a"-trimmed-R2P -thread 12
mkdir "$a"-max_bins
mv "$a"-maxbin* "$a"-max_bins/

#maxbin creates extra files in the output folder that must be removed before running das_input.py 
#run das_input.py to create dastool input tsv files from bins 

## 3.3 DasTool 
DAS_Tool -i "$a"-dastools/das_input_maxbin.tsv,"$a"-dastools/das_input_metabat.tsv -c mapping/"$a"-final_contigs.fa -o "$a"-dastool -l maxbin,metabat2 --threads=12 --write_bins --write_unbinned --write_bin_evals

## 3.4 Checkm Quality 
#runs tree, tree_qa, lineage_set, analyze, qa 
mkdir "$a"-checkm
checkm lineage_wf "$a"-dastool/dastool_"$a"_DASTool_bins/ "$a"-checkm/ -t 24 -x fa --tab_table
mkdir "$a"-quality_bins
#run bin_quality.py; moves bins with >50% completeness, <10% contanimation to "$a"-quality_bins

### Step 4 Sequence to function 

#directory to write prodigal outputs
mkdir $a-prodigal
mkdir $a-prodigal/proteins
mkdir $a-prodigal/dna
mkdir $a-prodigal/gff

#run prodigal on all bins in the dastools output directory, save aa, and gff output

for i in "$a"-quality_bins/*; do prodigal -i $i -f gff "$a"-prodigal/gff/$i ; done
for i in "$a"-quality_bins/*; do prodigal -i $i -a $a-prodigal/proteins/$i -o $a-prodigal/dna/$i ; done

#move unbinned.fa out of consensus bin directory & run make_GK_input_fasta.py 

### Step 5 Structure and Function 
### GK web tool output ###

KEGG-decoder -i "$a"-ko.txt -o "$a"-functions.tsv 
#runs identify, align, classify   
gtdbtk classify_wf -x fa --cpus 8 --genome_dir "$a"-quality_bins/ --out_dir "$a"-gtdbtk




