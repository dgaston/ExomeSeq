#!/usr/bin/perl
use strict;
use warnings;

my $sample_id = shift;
my $outfile = "run_alignment_" . $sample_id . ".sh";
my $paired1 = "";
my $paired2 = "";

my $GATK_PATH = "/sb/project/hmr-143-aa/dgaston/GATK";
my $reference = "$GATK_PATH/ucsc.hg19.fasta";
my $dbsnp = "$GATK_PATH/dbsnp_137.hg19.vcf";
my $indel_ref1 = "$GATK_PATH/1000G_phase1.indels.hg19.vcf";
my $indel_ref2 = "$GATK_PATH/Mills_and_1000G_gold_standard.indels.hg19.vcf";

my $sam = $sample_id . ".sam";
my $bam = $sample_id . ".bam";
my $added = $sample_id . ".rg.bam";
my $sorted = $sample_id . ".sorted.bam";
my $duplicates_removed = $sample_id . ".dup_removed.sorted.bam";
my $realigned = $sample_id . ".realigned.sorted.bam";
my $recal = $sample_id . ".recalibrated.sorted.bam";
my $reduced = $sample_id . ".reduced.bam";

my $metrics_file = $sample_id . "duplicate.metrics.txt";
my $realignment_targets = $sample_id . ".targets.intervals";
my $recal_config = $sample_id . ".recal.grp";
my $recal_pdf = $sample_id . ".recal.pdf";

open(OUT, ">$outfile") or die "Could not open $outfile. Exiting...\n\n";

#Header
print OUT "#!/bin/bash\n";
print OUT "#PBS -l nodes=1:ppn=12\n";
print OUT "#PBS -l walltime=96:00:00\n";
print OUT "#PBS -o $sample_id.outfile\n";
print OUT "#PBS -e $sample_id.errorfile\n";
print OUT "#PBS -V\n";
print OUT "#PBS -N " . $sample_id . "_pipeline\n";

print OUT "cd /sb/project/hmr-143-aa/dgaston/Exomes/\n";

#Body
print OUT "bowtie2 -p 12 --sensitive -x /sb/project/hmr-143-aa/dgaston/Bowtie2Indexes/genome --phred33 --local --rdg 3,2 --rfg 3,2 -1 $paired1 -2 $paired2 -S $sam\n";
print OUT "java -Xmx4g -jar ~/bin/SamFormatConverter.jar INPUT=$sam OUTPUT=$bam CREATE_MD5_FILE=true VALIDATION_STRINGENCY=LENIENT\n";
print OUT "java -Xmx4g -jar ~/bin/AddOrReplaceReadGroups.jar INPUT=$bam OUTPUT=$added CREATE_MD5_FILE=true RGLB=$sample_id RGPL=illumina RGPU=$sample_id RGSM=$sample_id SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT\n";
print OUT "java -Xmx4g -jar ~/bin/ReorderSam.jar CREATE_INDEX=true CREATE_MD5_FILE=true INPUT=$added OUTPUT=$sorted REFERENCE=$reference VALIDATION_STRINGENCY=LENIENT\n";
print OUT "java -Xmx4g -jar ~/bin/MarkDuplicates.jar CREATE_INDEX=true INPUT=$sorted OUTPUT=$duplicates_removed METRICS_FILE=$metrics_file REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT\n";
print OUT "java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 12 -R $reference --known $indel_ref1 --known $indel_ref2 -I $duplicates_removed -o $realignment_targets\n";
print OUT "java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner -I $duplicates_removed -o $realigned -targetIntervals $realignment_targets -R $reference\n";
print OUT "java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T BaseRecalibrator -I $realigned -o $recal_config -R $reference --knownSites $dbsnp --plot_pdf_file $recal_pdf\n";
print OUT "java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T PrintReads -I $realigned -o $recal -R $reference -BQSR $recal_config\n";
#print OUT "java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T ReduceReads -I $recal -o $reduced -R $reference\n";

close(OUT);