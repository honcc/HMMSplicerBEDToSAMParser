#!/usr/bin/perl/ -w
$|++;
use strict;
use File::Path;
use Time::HiRes qw( time );

######################################################################################################################################################
#
#	Description
#		This is a perl script to parse the non-colpased BED file generated from HMMSplicer into SAM alignment. It will first extract the successfully 
#	aligned (unspliced) sequences from the initial SAM file generated from HMMSplicer (have to turn on the -S option in bowtie command of the HMMSPLicer python script), 
#	and it read the BED file for spliced reads, and takes the sequence from the SAM file (unaligned ones). User have to specify paired-end or single-end. 
#	In case of paired-end, the read name has to be ended with 1 and 2. The aim of this script is to generate a SAM file that is comparable to that of Tophat 
#	so that it can be loaded to Cufflinks.
#
#	Input	
#
#		--samFile=					the samFile found in the tmp folder of HMMsplicer output, used to extract the sequences, normally named "bowtie.txt". You have to have to turn on the -S option in bowtie command of the HMMSPLicer python script;
#		--multipleCutoff=			the minimum score for junctions with more than 1 supported read, default = 400;
#		--singleCutoff=				the minimum score for junctions with only 1 supported read=, default = 600;
#		--uniqueIndivBEDFile= 		the uncollapsed BED file contains the BED info for all unfiltered reads, normally it is "junctionAgain.bed";
#		--duplicatedIndivBEDFile=	the BED file that contains the junctions with duplicated flanking sequence (since uniqueIndivBEDFile doesn't contain duplicated junctions), normally "junction_withDups.bed";
#		--topHatCrossValidBED=		the bed file output from topHat to filter the duplicate junctions, as topHat takes into account "exon island" and pair end read information; "no" to turn off;
#		--refFasta=					the reference genome for extracting the flanking sequences;
#		--SAMHeader=				to preserve the SAM header. "yes" or "no", if "no", the header will be print as a separate file;
#		--justCollapse=				"yes" or "no"; if "yes", the program will quit after collapsing the BEDs and will not generate the spliced samfile;
#		--nonCanOutSam=				"yes" or "no"; if "yes", the non-cannoical junctions will be output in the spliced samfile; default=no
#		--outNonSplicedAligned		"yes" or "no"; if "yes",  to output aligned but unspliced reads into the spliced sam; use "no" to to prevent the unclear bowtie aligning parameters defined in HMMSplicer to mess thing up; default = no
#		--outDir= 					the directory for output;
#
#	Output
#		1. unfilter raw junction info ---> $outDir/HMMSplicer/junctionInfo/HMMSplicer.raw.JnctnInfo.txt
#		2. filtered all junction info ---> $outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.JnctnInfo.txt
#		3. filtered unique canonical junction info ---> $outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.JnctnInfo.txt
#		4. filtered unique non-canonical junction info ---> $outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.noncan.unique.JnctnInfo.txt
#		5. samFile of all mapped read, including both spliced and unscpliced ---> $outDir/HMMSplicer/finalSAM/HMMSplicer.spliced.sam
#		6. unaligned read in fastq ---> $outDir/HMMSplicer/finalSAM/HMMSplicer.unaligned.fastq
#
#	Usage example
#		perl HMMSplicerBEDToSAMParser_v0.5.pl --samFile=/Users/chung/Desktop/NGS/test/results/tmp/bowtie.txt --multipleCutoff=400 --singleCutoff=600 --uniqueIndivBEDFile=/Users/chung/Desktop/NGS/test/results/tmp/junctionAgain.bed --duplicatedIndivBEDFile=/Users/chung/Desktop/NGS/test/results/tmp/junction_withDups.bed --topHatCrossValidBED=no --refFasta=/Volumes/CCHON1TBA/softwareForNGS/myPerlScripts/pipeLines/allTheWayHomeRunner/resources/genome/EHI_v13.fa --SAMHeader=no --outDir=./
#
#	Assumption 
#
#	Versions
#
#		v0.1
#		-debut;
#
#		v0.2
#		-the bed files contains the duplicated junction will be collaped;
#
#		v0.3
#		-the single pair end option is discarded. Now it only support single end. (Focus on single end, HMMSplicer's advantage is sensitivity, my plan is to use Tophat's advanatge of pairend to control the speficity of duplicated junctions in HMMSplicer)
#		-use the topHat bed output as a filter for duplicated junctions
#
#		v0.4
#		-the flanking sequence will be extracted, and the junction will be analyzed as cannonical and non-cannonical;
#		-outDir option added
#
#		v0.5
#		-"no" tophat bed file option added.
#		-BED files will be generated along side the junction info.txt;
#
#		V0.6
#		-added option to not produced the sam but just Collapse Junct (--justCollapse=)
#
#		v0.7
#		-added option to choose to output non-cannonical junctions in sam or not; default=no;
#		-added option to output aligned but unspliced reads into the spliced sam or not; The purpose is to prevent the unclear bowtie aligning parameters defined in HMMSplicer to mess thing up
#
######################################################################################################################################################

#==========================================================Main body starts==========================================================================#


#1----------Read the parameters----------#
use vars qw ($samFile $multipleCutoff $singleCutoff $uniqueIndivBEDFile  $duplicatedIndivBEDFile $topHatCrossValidBED $refFasta $SAMHeader $justCollapse $nonCanOutSam $outNonSplicedAligned $outDir);
my ($samFile, $multipleCutoff, $singleCutoff, $uniqueIndivBEDFile, $duplicatedIndivBEDFile, $topHatCrossValidBED, $refFasta, $SAMHeader, $justCollapse, $nonCanOutSam, $outNonSplicedAligned, $outDir) = readParameters();
printCMDLogOrFinishMessage("CMDLog");

#----------Read the reference genome----------#
my $cntgSeqHsh_ref = readMultiFastaFile($refFasta);

#----------Read the tophat BED-----------#
my $tophatJnctnScoreHsh_ref = readTophatBEDFile();

#2----------read the two final BED files----------#
my (%hitBEDLineHsh);
open (RAWJNC, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.raw.JnctnInfo.txt");
open (FLTRJNC, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.JnctnInfo.txt");
open (FLTRCNUNQJNC, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.JnctnInfo.txt");
open (FLTRCNALLJNC, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.JnctnInfo.txt");
open (FLTRNONCNUNQJNC, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.noncan.unique.JnctnInfo.txt");

open (RAWBED, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.raw.bed");
print RAWBED "track name=rawHMMSplicer description=rawHMMSplicer useScore=1\n";
open (FLTRBED, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.bed");
print FLTRBED "track name=allFilteredHMMSplicer description=allFilteredHMMSplicer useScore=1\n";
open (FLTRCNUNQBED, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.bed");
print FLTRCNUNQBED "track name=uniqueFilteredCanonicalHMMSplicer description=uniqueFilteredCanonicalHMMSplicer useScore=1\n";
open (FLTRCNALLBED, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.can.unique.dup.bed");
print FLTRCNALLBED "track name=uniqueAndDupFilteredCanonicalHMMSplicer description=uniqueAndDupFilteredCanonicalHMMSplicer useScore=1\n";
open (FLTRNONCNUNQBED, ">$outDir/HMMSplicer/junctionInfo/HMMSplicer.filter.noncan.unique.bed");
print FLTRNONCNUNQBED "track name=uniqueFilteredNoncanonicalHMMSplicer description=uniqueFilteredNoncanonicalHMMSplicer useScore=1\n";

my $hitBEDLineHsh_ref = \%hitBEDLineHsh;
$hitBEDLineHsh_ref = \%hitBEDLineHsh;
my %jnctnStrFilterHsh; #---empty hash, as unique jnctn doesnt need filter
my $jnctnStrFilterHsh_ref = \%jnctnStrFilterHsh;

$hitBEDLineHsh_ref = collapseBEDFile($hitBEDLineHsh_ref, $uniqueIndivBEDFile, "unique", $jnctnStrFilterHsh_ref, $cntgSeqHsh_ref);
#---tophatJnctnScoreHsh_ref will be an empty hash if no tophat BED was specified
$hitBEDLineHsh_ref = collapseBEDFile($hitBEDLineHsh_ref, $duplicatedIndivBEDFile, "duplicate", $tophatJnctnScoreHsh_ref, $cntgSeqHsh_ref);

close RAWJNC;
close FLTRJNC;
close FLTRCNALLJNC;
close FLTRCNUNQJNC;
close FLTRNONCNUNQJNC;

close RAWBED;
close FLTRBED;
close FLTRCNALLBED;
close FLTRCNUNQBED;
close FLTRNONCNUNQBED;

if ($justCollapse ne "yes") {
	my $SAMFlagTableHsh_ref = defineSAMFlagTable();

	my $hitBEDLineNum = 0;
	open (SPLICEDSAM, ">$outDir/HMMSplicer/finalSAM/HMMSplicer.spliced.sam");	
	open (UNALIGNEDFQ, ">$outDir/HMMSplicer/finalSAM/HMMSplicer.unaligned.fastq");
	open (SAMHEADER, ">$outDir/HMMSplicer/finalSAM/SAMHeader.txt") if ($SAMHeader eq "no");	
	readBowtieSAMAndPrintSplicedSAM($hitBEDLineHsh_ref, $SAMFlagTableHsh_ref);
	close SPLICEDSAM;
	close UNALIGNEDFQ;
	close SAMHEADER if ($SAMHeader eq "no");
}

printCMDLogOrFinishMessage("finishMessage");
exit;
#========================================================= Main body ends ===========================================================================#
################################################################## readParameters ####################################################################
sub readParameters {
	
	$singleCutoff = 600;
	$multipleCutoff = 400;
	$topHatCrossValidBED = "no";
	$SAMHeader = "no";
	$justCollapse = "no";
	$nonCanOutSam = "no";
	$outNonSplicedAligned = "no";

	foreach (@ARGV) {
		if ($_ =~ m/--samFile=/) {$samFile = substr ($_, index ($_, "=")+1);} 
		elsif ($_ =~ m/--multipleCutoff=/) {$multipleCutoff = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--singleCutoff=/) {$singleCutoff = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--uniqueIndivBEDFile=/) {$uniqueIndivBEDFile = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--duplicatedIndivBEDFile=/) {$duplicatedIndivBEDFile = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--topHatCrossValidBED=/) {$topHatCrossValidBED = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--refFasta=/) {$refFasta = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--justCollapse=/) {$justCollapse = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--nonCanOutSam=/) {$nonCanOutSam = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--outNonSplicedAligned=/) {$outNonSplicedAligned = substr ($_, index ($_, "=")+1);}
		elsif ($_ =~ m/--outDir=/) {$outDir = substr ($_, index ($_, "=")+1);}
	}
	
	system ("mkdir -p -m 777 $outDir/HMMSplicer/finalSAM/"); #---to make sure the output
	system ("mkdir -p -m 777 $outDir/HMMSplicer/junctionInfo/"); #---to make sure the output
	
	return ($samFile, $multipleCutoff, $singleCutoff, $uniqueIndivBEDFile, $duplicatedIndivBEDFile, $topHatCrossValidBED, $refFasta, $SAMHeader, $justCollapse, $nonCanOutSam, $outNonSplicedAligned, $outDir);
}
################################################################## readParameters ####################################################################
sub readTophatBEDFile {

	#	track name=junctions description="TopHat junctions"
	#	DS571731	2015	2246	JUNC00000001	20	+	2015	2246	255,0,0	2	92,94	0,137
	#	DS571736	1442	1672	JUNC00000002	10	+	1442	1672	255,0,0	2	94,77	0,153
	#	DS571148	19849	20136	JUNC00000003	1	-	19849	20136	255,0,0	2	92,8	0,279
	#	DS571148	19857	20226	JUNC00000004	22	-	19857	20226	255,0,0	2	86,96	0,273
	#	DS571148	25913	26149	JUNC00000005	23	+	25913	26149	255,0,0	2	92,95	0,141
	#	DS571148	31882	32108	JUNC00000006	6	-	31882	32108	255,0,0	2	81,87	0,139
	#	DS571148	31872	32047	JUNC00000007	2	-	31872	32047	255,0,0	2	91,9	0,166
	#	DS571148	31874	32262	JUNC00000008	10	-	31874	32262	255,0,0	2	89,85	0,303
	#	DS571148	32030	32265	JUNC00000009	9	-	32030	32265	255,0,0	2	92,88	0,147
	#	DS571148	40786	41020	JUNC00000010	3	+	40786	41020	255,0,0	2	96,92	0,142
	#	DS571148	41887	42099	JUNC00000011	13	+	41887	42099	255,0,0	2	75,87	0,125

	my %tophatJnctnScoreHsh;
	
	if ($topHatCrossValidBED ne "no") {
		open (INFILE, "$topHatCrossValidBED") || die "$topHatCrossValidBED can be read\n";
		print "Reading all junctions reads from $topHatCrossValidBED\n";

		while (my $theLine = <INFILE>) {
			next if ($theLine =~ m/^track name/);
			my @theLineSplt = split /\t/, $theLine;
			my $cntg = $theLineSplt[0];
			my $cntgStart = $theLineSplt[1];
			my $cntgEnd = $theLineSplt[2];
			my $juncID = $theLineSplt[3];
			my $score =$theLineSplt[4];
			my $strd = $theLineSplt[5];
			my $blkSize = $theLineSplt[10];
			my @blkSizeAry = split /,/, $blkSize;
			my $blkSize1 = $blkSizeAry[0];
			my $blkSize2 = $blkSizeAry[1];
			my $intronStart = $cntgStart + $blkSize1 + 1;
			my $intronEnd = $cntgEnd - $blkSize2;
			my $jnctnStr = $cntg.",".$intronStart.",".$intronEnd;
			$tophatJnctnScoreHsh{$jnctnStr} = $score;
		}
		close INFILE;
		my $topHatJuncNum = keys %tophatJnctnScoreHsh;
		print "Totally $topHatJuncNum tophat junctions have been stored.\n";
	} else {
		print  "NO tophat bed file will be read.\n";
	}

	return \%tophatJnctnScoreHsh;
}
###########################################§######### collapse and filter the BED files contains the individual spliced reads ########################
sub collapseBEDFile {

#track name='HMMSplicer Junctions' description='HMMSplicer Junctions' useScore=1
#DS571286	30173	30382	EHI_026410_0_r34_f238_DS571286_30077/30176/cnt/100M/+_30215/30366/msp/junc_2030/53M52N47M/-_/2|junc=24	1081.30673742	-	30173	30382	0	2	94,63,	0,146,
#DS571286	31440	31666	EHI_026430_0_r8_f254_DS571286_31687/31786/cnt/100M/-_31485/31632/msp/junc_2032/44M48N56M/+_/2|junc=11	1091.21669189	+	31440	31666	0	2	88,90,	0,136,
#DS571286	30743	30989	EHI_026420_0_r57_f285_DS571286_30790/30949/msp/junc_2031/48M60N52M/+_31035/31134/cnt/100M/-_/1|junc=24	1165.62798932	-	30743	30989	0	2	94,92,	0,154,
#DS571286	31762	31988	EHI_026430_0_r13_f261_DS571286_31794/31944/msp/junc_2034/57M51N43M/-_31633/31732/cnt/100M/+_/1|junc=8	1053.33607763	-	31762	31988	0	2	88,87,	0,139,
#DS571286	31768	32181	EHI_026430_1_r11_f258_DS571286_31793/32181/msp/junc_2035/58M289N42M/-_31635/31734/cnt/100M/+_/1|junc=4	913.556319795	-	31768	32181	0	2	82,42,	0,371,
#DS571286	32537	32950	EHI_026440_0_r33_f238_DS571286_32439/32538/cnt/100M/+_32577/32910/msp/junc_2037/49M234N51M/-_/2|junc=20	1124.95804055	+	32537	32950	0	2	88,91,	0,322,
	
	my %hitBEDLineHsh = %{$_[0]};
	my $theBEDFile = $_[1];
	my $BEDFileTag = $_[2];
	my %jnctnStrFilterHsh = %{$_[3]};
	my %cntgSeqHsh = %{$_[4]};

	my %allJnctnHsh;

	my $junctNum = 0;
	my $IDHead = "unq";
	$IDHead = "dup" if ($BEDFileTag eq "duplicate");
	my $rawJunct = 0;
	
	print "Reading all spliced reads from $theBEDFile\n";
	open (INFILE, "$theBEDFile") || die "Can't read $theBEDFile\n";
	while (my $theLine = <INFILE>) {
		next if ($theLine =~ m/^track name/);
		$rawJunct++;
		print "$rawJunct raw junctions stored.                       \r";
		my @theLineSplt = split /\t/, $theLine;
		my $cntg = $theLineSplt[0];
		my $cntgStart = $theLineSplt[1];
		my $cntgEnd = $theLineSplt[2];
		my $seqID = $theLineSplt[3];
		my $score = sprintf ("%.1f", $theLineSplt[4]);
		my $strd = $theLineSplt[5];
		my $blkSize = $theLineSplt[10];
		my @blkSizeAry = split /,/, $blkSize;
		my $blkSize1 = $blkSizeAry[0];
		my $blkSize2 = $blkSizeAry[1];
		my $intronStart = $cntgStart + $blkSize1;
		my $intronEnd = $cntgEnd - $blkSize2 - 1;		
		my $jnctnStr = $cntg.",".$intronStart.",".$intronEnd;
		push @{$allJnctnHsh{$jnctnStr}}, $seqID.",:,".$score.",:,".$blkSize1.",:,".$blkSize2.",:,".$strd.",:,".$cntgStart.",:,".$cntgEnd;
	}
	close INFILE;

	print "\n";
	
	my $discardJnctnNum = my $belowCutoffJnctnNum = my $validJnctnNum = 0;

	my $jnctnNum = keys %allJnctnHsh;
	print "$jnctnNum junctions have been stored. Start collapsing and filtering.\n";
	
	#---collapse these junctions
	my $procJunctNum = 0;
	foreach my $jnctnStr (sort {$a cmp $b} keys %allJnctnHsh) {
		$procJunctNum++;
		print "$procJunctNum / $jnctnNum junctions collapsed. junct : $jnctnStr                    \r";
		my $supportSeqNum = @{$allJnctnHsh{$jnctnStr}};
		my $scoreCutOff = $multipleCutoff; #---multiple
		$scoreCutOff = $singleCutoff if ($supportSeqNum == 1); #---single 
		
		#---sort the array by scores, to make sure the pariwise comparison starts from the highest scoring read
		my $highestScore = -99999;
		my $indexWithHighestScore;
		my %tmpScoreAryIndexForSorting;
		for (my $i=0; $i < @{$allJnctnHsh{$jnctnStr}}; $i++) {
			my @theJnctnInfoStrSplt = split /,:,/, ${$allJnctnHsh{$jnctnStr}}[$i]; 
			my $curntScore = $theJnctnInfoStrSplt[1];
			$tmpScoreAryIndexForSorting{$i} = $curntScore;
			if ($curntScore >= $highestScore) {
				$indexWithHighestScore = $i;
				$highestScore = $curntScore;
			}
		}
		
		#---get the info of the read with lowest score as the base 
		my $baseJnctnInfoStr = ${$allJnctnHsh{$jnctnStr}}[$indexWithHighestScore];
		my @baseJnctnInfoStrSplt = split /,:,/, $baseJnctnInfoStr; 
		my $extendScore = $baseJnctnInfoStrSplt[1];
		my $extendBlk1Size = $baseJnctnInfoStrSplt[2];
		my $extendBlk2Size = $baseJnctnInfoStrSplt[3];
		my %dirCountHsh;
		
		$dirCountHsh{"+"} = 0;
		$dirCountHsh{"-"} = 0;
		
		#	copied from the original PLoS One 
		#
		#	All reads creating the same intron are collapsed into a single junction with the score for these reads increased in relation to number of 
		#   additionally covered bases. Distinct reads covering the same junction add significantly to a its potential to be real, but two identical reads 
		#   may be from the same source, such as PCR amplification artifacts. To follow the previous example, a 35/10 split (35 bp on the first exon, 10 bp on 
		#   the second exon) combined with another 35/10 split would not increase the score, but the 35/10 split plus a 10/35 split would yield a substantial 
		#   boost to the score because the covered bases would now be now 35/35. To be exact, imagine the 10/35 junction read has a score of 800 and the 35/10 
		#   junction read has a score of 600. The higher score read is considered first, then the second read is collapsed onto it. In this case, the new junction adds 25 
		#   bases out of a total of, now, 70 bases covered, so a value of (25 / 70) * 600 is added to the original score of 800, yielding a collapsed score of 1214.2.		
		
		#---the protocol below reproduces exactly the collapsed output of HMMSplicer
		
		foreach my $aryIndex (sort {$tmpScoreAryIndexForSorting{$b} <=> $tmpScoreAryIndexForSorting{$a}} keys %tmpScoreAryIndexForSorting) {#---sort the arry index based on score
			my $otherJnctnInfoStr = ${$allJnctnHsh{$jnctnStr}}[$aryIndex];		
			my @otherJnctnInfoStrSplt = split /,:,/, $otherJnctnInfoStr; 
			my $curntScore = $otherJnctnInfoStrSplt[1];
			my $curntBlk1Size = $otherJnctnInfoStrSplt[2];
			my $curntBlk2Size = $otherJnctnInfoStrSplt[3];
			my $strd = $otherJnctnInfoStrSplt[4];
			$dirCountHsh{$strd}++;

			my $smallerScore = $curntScore; #---smaller score will be the base score to add to extendscore
			if ($curntScore >= $extendScore) {#---$curntScore will be the extendScore
				$smallerScore = $extendScore;
				$extendScore = $curntScore;
			}
			
			#---assumming all reads are of the same length, if would only be either of the following case
			my $extendLength = 0;
			if ($curntBlk1Size > $extendBlk1Size) {#---blk 1 is extended
				$extendLength = $curntBlk1Size - $extendBlk1Size;				
				$extendBlk1Size = $curntBlk1Size;
			} elsif ($curntBlk2Size > $extendBlk2Size) {#---blk 2 is extended
				$extendLength = $curntBlk2Size - $extendBlk2Size;				
				$extendBlk2Size = $curntBlk2Size;
			} 			
			
			my $extendBothBlkSize = $extendBlk1Size + $extendBlk2Size;
			$extendScore += $smallerScore*($extendLength/$extendBothBlkSize); #---wont add anything if there's only 1 seq or no extension as $extendLength = 0;
		}
		$extendScore = sprintf ("%.1f", $extendScore);

		#---get the flanking sequence of the junction
		my @jnctnStrSplt = split /,/, $jnctnStr;
		my $cntg = $jnctnStrSplt[0];
		my $intronStart = $jnctnStrSplt[1];
		my $intronEnd = $jnctnStrSplt[2];
		my $cntgSeq = $cntgSeqHsh{$cntg};
		my $flankLength = 30;
		my $leftFlankSeq = substr $cntgSeq, $intronStart-$flankLength, ($flankLength*2)+2;
		my $rightFlankSeq = substr $cntgSeq, $intronEnd-$flankLength-1, ($flankLength*2)+2;
		my $left2ntJnctnSeq = substr $cntgSeq, $intronStart, 2;
		my $right2ntJnctnSeq = substr $cntgSeq, $intronEnd-1, 2;
		my $jnctnSeq = $left2ntJnctnSeq.$right2ntJnctnSeq;
		
		my $BEDBlkLeftSeq = substr $cntgSeq, $intronStart-$extendBlk1Size, $extendBlk1Size;
		my $BEDBlkRightSeq = substr $cntgSeq, $intronEnd+1, $extendBlk2Size;
		
		my ($canonicality, $canJnctnDir);
		
		if ($jnctnSeq =~ m/GTAG/i) {#---plus strand cannoical
			
			$canonicality = "cannoical";
			$canJnctnDir = "+";
			
		} elsif ($jnctnSeq =~ m/CTAC/i) {#---minius strand cannoical

			$jnctnSeq = "GTAG";
			$canonicality = "cannoical";
			$canJnctnDir = "-";
			($leftFlankSeq, $rightFlankSeq) = ($rightFlankSeq, $leftFlankSeq);
			$leftFlankSeq = reverse $leftFlankSeq;
			$rightFlankSeq = reverse $rightFlankSeq;
			$leftFlankSeq =~ tr/ACGTacgt/TGCAtgca/;
			$rightFlankSeq =~ tr/ACGTacgt/TGCAtgca/;

			($BEDBlkLeftSeq, $BEDBlkRightSeq) = ($BEDBlkRightSeq, $BEDBlkLeftSeq);
			$BEDBlkLeftSeq = reverse $BEDBlkLeftSeq;
			$BEDBlkRightSeq = reverse $BEDBlkRightSeq;
			$BEDBlkLeftSeq =~ tr/ACGTacgt/TGCAtgca/;
			$BEDBlkRightSeq =~ tr/ACGTacgt/TGCAtgca/;

		} else {#---non-cannoical

			$canonicality = "non-cannoical";
			$canJnctnDir = "*";
		}

		#---calculate the info for bed file
		my $bedStart = $intronStart - $extendBlk1Size;
		my $bedEnd = $intronEnd + $extendBlk2Size + 1;
		my $blk1Start = 0;
		my $blk2Start = $intronEnd - $intronStart + $extendBlk1Size + 1;
		
		$junctNum++;
		my $junctID = $IDHead.$junctNum;
		
		#---defined the bedLine info and the junction txt info
		my $BEDLine = $cntg."\t".$bedStart."\t".$bedEnd."\t"."read=$supportSeqNum|score=$extendScore|$junctID"."\t".$extendScore."\t".$canJnctnDir."\t".$bedStart."\t".$bedEnd."\t"."0"."\t"."2"."\t"."$extendBlk1Size,$extendBlk2Size,"."\t"."$blk1Start,$blk2Start,";
		my $juncInfoLine = $canonicality."\t".$BEDFileTag."\t".$jnctnStr."\t".$extendScore."\t".$supportSeqNum."\t".$dirCountHsh{"+"}."\t".$dirCountHsh{"-"}."\t".$leftFlankSeq."\t".$rightFlankSeq."\t".$jnctnSeq."\t".$canJnctnDir."\t".$junctID."\t".$BEDBlkLeftSeq."\t".$BEDBlkRightSeq;
		
		print RAWJNC $juncInfoLine."\n";
		print RAWBED $BEDLine."\n";

		#---store if passed the cut-off, pass the spliced seq info into a hash
		if ($extendScore >= $scoreCutOff) {
			
			#---if not tophat filter or have filter and present
			if (((keys %jnctnStrFilterHsh) <= 1) or (exists $jnctnStrFilterHsh{$jnctnStr})) {
				
				$validJnctnNum++;
				
				print FLTRJNC $juncInfoLine."\n";
				print FLTRCNALLJNC $juncInfoLine."\n" if ($canonicality eq "cannoical");
				print FLTRCNUNQJNC $juncInfoLine."\n" if (($canonicality eq "cannoical") and ($BEDFileTag eq "unique"));
				print FLTRNONCNUNQJNC $juncInfoLine."\n" if (($canonicality eq "non-cannoical") and ($BEDFileTag eq "unique"));

				print FLTRBED $BEDLine."\n";
				print FLTRCNALLBED $BEDLine."\n" if ($canonicality eq "cannoical");
				print FLTRCNUNQBED $BEDLine."\n" if (($canonicality eq "cannoical") and ($BEDFileTag eq "unique"));
				print FLTRNONCNUNQBED $BEDLine."\n" if (($canonicality eq "non-cannoical") and ($BEDFileTag eq "unique"));
				
				if (($nonCanOutSam ne "no") or ($canonicality eq "cannoical")) {#---store info only for the cannoical splicing and skip non-cannoical when nonCanOutSam is "no"
				
					foreach my $splicedSeqInfo (@{$allJnctnHsh{$jnctnStr}}) {
						#	original variable ->>> push @{$allJnctnHsh{$jnctnStr}}, $seqID.",:,".$score.",:,".$blkSize1.",:,".$blkSize2.",:,".$dir.",:,".$cntgStart.",:,".$cntgEnd;	
	
						#---get the info	
						my @splicedSeqInfoSplt = split /,:,/, $splicedSeqInfo;
						my $seqID = $splicedSeqInfoSplt[0];
						my $blkSize1 = $splicedSeqInfoSplt[2];
						my $blkSize2 = $splicedSeqInfoSplt[3];
						my $strd = $splicedSeqInfoSplt[4];
						my $cntgStart = $splicedSeqInfoSplt[5];
						my $cntgEnd = $splicedSeqInfoSplt[6];	
		
						#---convert the info to SAM needed info and store in a Hash
						my $samStart = $cntgStart +1;
						my $samEnd = $cntgEnd;
						my $cigarHead = $blkSize1;
						my $cigarIntron = $intronEnd - $intronStart + 1;
						my $cigarTail = $blkSize2;
						my $cigarStr = $cigarHead."M".$cigarIntron."N".$cigarTail."M";
		
						my $seqIDCount = (keys %{$hitBEDLineHsh{$seqID}}) + 1; #---multiple hit;
						${$hitBEDLineHsh{$seqID}}{$seqIDCount} = $cntg.",:,".$samStart.",:,".$strd.",:,".$cigarStr; #---could be both mates spliced in a pair
					}
				}#---end of foreach my $splicedSeqInfo (@{$allJnctnHsh{$jnctnStr}})
			
			} else {
			
				$discardJnctnNum++;
			}
		
		} else {#---end of if ($extendScore >= $scoreCutOff)
			
			$belowCutoffJnctnNum++;
			
		}
		
		delete $allJnctnHsh{$jnctnStr}; #---delete the hsh to release memory
		
	}#---End of foreach my $jnctnStr (sort {$a cmp $b} keys %allJnctnHsh)
	
	print "\n";
	
	my $filteredSplicedSeqNum = keys %hitBEDLineHsh;
	print "Totally $filteredSplicedSeqNum spliced sequences have been stored.\n";
	print "validJnctnNum = $validJnctnNum.\n";
	print "discardJnctnNum = $discardJnctnNum.\n";
	print "belowCutoffJnctnNum  = $belowCutoffJnctnNum.\n";
	
	return \%hitBEDLineHsh;
}
####################################### read the UNspliced Bowtie SAM and convert it to spliced SAM using the BED informations #######################
sub readBowtieSAMAndPrintSplicedSAM {

	my %hitBEDLineHsh = %{$_[0]};
	my %SAMFlagTableHsh = %{$_[1]};


	my %seqTypeHash;

	my $unsplicedAlginedReadNum = my $splicedAlignedReadNum = my $totalReadNum = my $allAlignedReadNum = my $unalignedReadNum = 0;

	#---read the samfile
	print "Scanning the samfile for the BED features.\n";
	open (SAMFILE, "$samFile")  || die "Can't read $samFile\n";
	while (my $theLine = <SAMFILE>) {#--check and print the header
		chomp $theLine;
		
		if ($theLine =~ m/^@/) {
		 print SAMHEADER $theLine."\n" if ($SAMHeader eq "no"); 
		 print SPLICEDSAM $theLine."\n" if ($SAMHeader eq "yes"); 
		 next;
		}
		
		$totalReadNum++;
		
		my @theLineSplt = split /\t/, $theLine;
		my $seqID = $theLineSplt[0];
		
		#---store all SAM bits in a Hsh
		my $SAMFlag = $theLineSplt[1];
		my $SAMBitStr = $SAMFlagTableHsh{$SAMFlag};
		my @SAMBitAry = split /\+/, $SAMBitStr;
		my %SAMBitHsh;
		foreach my $SAMBit (@SAMBitAry) {$SAMBitHsh{$SAMBit}++;}
		
		#if (not exists $seqTypeHash{$seqID}) { 
		#	${$seqTypeHash{$seqID}}{"aligned"} = 0;
		#	${$seqTypeHash{$seqID}}{"spliced"} = 0;
		#	${$seqTypeHash{$seqID}}{"all"} = 0;
		#}

		#---unaligned reads as the bit 4 is present
		if ((exists $SAMBitHsh{4}) and (not(exists $hitBEDLineHsh{$seqID}))) {#---unaligned and NOT spliced
			
				print UNALIGNEDFQ "@".$theLineSplt[0]."\n";
				print UNALIGNEDFQ $theLineSplt[9]."\n";
				print UNALIGNEDFQ "+".$theLineSplt[0]."\n";
				print UNALIGNEDFQ $theLineSplt[10]."\n";
				$unalignedReadNum++;

		} else {#---aligned, maybe spliced or maynot be splcied
		
			$allAlignedReadNum++;
			
			#---check if unspliced and aligned
			if (not exists $SAMBitHsh{4}) {#---aligned, maybe unspliced (unique) and maybe spliced (if it is a duplicated read it is possible)
				
				if ($outNonSplicedAligned ne "no") {
					print SPLICEDSAM $theLine."\n";#---go straight to SPLICEDSAM
					$unsplicedAlginedReadNum++;
				} else {
					print UNALIGNEDFQ "@".$theLineSplt[0]."\n";
					print UNALIGNEDFQ $theLineSplt[9]."\n";
					print UNALIGNEDFQ "+".$theLineSplt[0]."\n";
					print UNALIGNEDFQ $theLineSplt[10]."\n";
					$unalignedReadNum++;
				}
			}

			#---check if the read is spliced
			if (exists $hitBEDLineHsh{$seqID}) {#---this SAM Line contains the hit seqID, may be the spliced seq itself or the mate

				foreach my $seqIDCount (keys %{$hitBEDLineHsh{$seqID}}) {
	
					$splicedAlignedReadNum++;
					
					#${$seqTypeHash{$seqID}}{"spliced"}++;
					#${$seqTypeHash{$seqID}}{"all"}++;
					my $BEDInfo = ${$hitBEDLineHsh{$seqID}}{$seqIDCount};
					my @BEDInfoSplt = split /,:,/, $BEDInfo;
					my $BEDCntg = $BEDInfoSplt[0];
					my $BEDSamStart = $BEDInfoSplt[1];
					my $BEDStrd = $BEDInfoSplt[2];
					my $BEDCigarStr = $BEDInfoSplt[3];				
					
				 	my ($flag, $seq, $qual);
					if ($BEDStrd eq "+") {
						$flag = 0;
						$seq = $theLineSplt[9];
						$qual = $theLineSplt[10];
							
					} elsif ($BEDStrd eq "-") {
						$flag = 16;
						$seq = reverse $theLineSplt[9]; #---reverse
						$seq =~ tr/ACGTacgt/TGCAtgca/; #---complement
						$qual = reverse $theLineSplt[10];
					}
						
					my @SAMAry;
					$SAMAry[0] = $seqID;
					$SAMAry[1] = $flag;
					$SAMAry[2] = $BEDCntg;
					$SAMAry[3] = $BEDSamStart;
					$SAMAry[4] = "255";
					$SAMAry[5] = $BEDCigarStr;
					$SAMAry[6] = "*";
					$SAMAry[7] = 0;
					$SAMAry[8] = 0;
					$SAMAry[9] = $seq;
					$SAMAry[10] = $qual;
					$SAMAry[11] = "XA:i:0\tXS:A:$BEDStrd"; #---end user tag by bowite, meaning aligned;
					my $outSAMLine = join "\t", @SAMAry;
					
					print SPLICEDSAM $outSAMLine."\n";
				
				}#--- end of foreach my $seqIDCount (keys %{$hitBEDLineHsh{$seqID}}) {

			}#---end of if (exists $hitBEDLineHsh{$seqID}) {

		}#---end of if ((exists $SAMBitHsh{4}) and (not(exists $hitBEDLineHsh{$seqID}))) 

	}#--end of while (my $theLine = <SAMFILE>) {
	close SAMFILE;
	
	open (LOGFILE, ">$outDir/HMMSplicer/junctionInfo/HMMSplicerBEDToSAMParser.log.txt");
	print LOGFILE "total number of read input in HMMSplicer = $totalReadNum\n";
	print LOGFILE "number of unspliced but aligned read = $unsplicedAlginedReadNum\n";
	print LOGFILE "number of spliced and aligned read = $splicedAlignedReadNum\n";
	print LOGFILE "number of unaligned read = $unalignedReadNum\n";
	close LOGFILE;	

	print "\ntotal number of read input in HMMSplicer = $totalReadNum\n";
	print "number of unspliced but aligned read = $unsplicedAlginedReadNum\n";
	print "number of spliced and aligned read = $splicedAlignedReadNum\n";
	print "number of unaligned read = $unalignedReadNum\n";

}
######################################################################################################################################################
sub defineSAMFlagTable {
	
#
#copied from http://bioinformatics.;.edu/chuanglab/wiki/index.php/SAM_pairwise_flag_translation_table
#
#0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair 
#0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment)  
#0x0004 4 the query sequence itself is unmapped 
#0x0008 8 the mate is unmapped  
#0x0010 16 strand of the query (0 for forward; 16 for reverse strand) 
#0x0020 32 strand of the mate  (0 for forward; 32 for reverse strand) 
#0x0040 64 the read is the ﬁrst read in a pair  
#0x0080 128 the read is the second read in a pair 
#0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records) 
	
	my %SAMFlagTableHsh = (

		0 => "0",
		1 => "1",
		2 => "2",
		3 => "1+2",
		4 => "0+4",
		5 => "1+4",
		6 => "0+2+4",
		7 => "1+2+4",
		8 => "0+8",
		9 => "1+8",
		10 => "0+2+8",
		11 => "1+2+8",
		12 => "0+4+8",
		13 => "1+4+8",
		14 => "0+2+4+8",
		15 => "1+2+4+8",
		16 => "0+16",
		17 => "1+16",
		18 => "0+2+16",
		19 => "1+2+16",
		20 => "0+4+16",
		21 => "1+4+16",
		22 => "0+2+4+16",
		23 => "1+2+4+16",
		24 => "0+8+16",
		25 => "1+8+16",
		26 => "0+2+8+16",
		27 => "1+2+8+16",
		28 => "0+4+8+16",
		29 => "1+4+8+16",
		30 => "0+2+4+8+16",
		31 => "1+2+4+8+16",
		32 => "0+32",
		33 => "1+32",
		34 => "0+2+32",
		35 => "1+2+32",
		36 => "0+4+32",
		37 => "1+4+32",
		38 => "0+2+4+32",
		39 => "1+2+4+32",
		40 => "0+8+32",
		41 => "1+8+32",
		42 => "0+2+8+32",
		43 => "1+2+8+32",
		44 => "0+4+8+32",
		45 => "1+4+8+32",
		46 => "0+2+4+8+32",
		47 => "1+2+4+8+32",
		48 => "0+16+32",
		49 => "1+16+32",
		50 => "0+2+16+32",
		51 => "1+2+16+32",
		52 => "0+4+16+32",
		53 => "1+4+16+32",
		54 => "0+2+4+16+32",
		55 => "1+2+4+16+32",
		56 => "0+8+16+32",
		57 => "1+8+16+32",
		58 => "0+2+8+16+32",
		59 => "1+2+8+16+32",
		60 => "0+4+8+16+32",
		61 => "1+4+8+16+32",
		62 => "0+2+4+8+16+32",
		63 => "1+2+4+8+16+32",
		64 => "0+64",
		65 => "1+64",
		66 => "0+2+64",
		67 => "1+2+64",
		68 => "0+4+64",
		69 => "1+4+64",
		70 => "0+2+4+64",
		71 => "1+2+4+64",
		72 => "0+8+64",
		73 => "1+8+64",
		74 => "0+2+8+64",
		75 => "1+2+8+64",
		76 => "0+4+8+64",
		77 => "1+4+8+64",
		78 => "0+2+4+8+64",
		79 => "1+2+4+8+64",
		80 => "0+16+64",
		81 => "1+16+64",
		82 => "0+2+16+64",
		83 => "1+2+16+64",
		84 => "0+4+16+64",
		85 => "1+4+16+64",
		86 => "0+2+4+16+64",
		87 => "1+2+4+16+64",
		88 => "0+8+16+64",
		89 => "1+8+16+64",
		90 => "0+2+8+16+64",
		91 => "1+2+8+16+64",
		92 => "0+4+8+16+64",
		93 => "1+4+8+16+64",
		94 => "0+2+4+8+16+64",
		95 => "1+2+4+8+16+64",
		96 => "0+32+64",
		97 => "1+32+64",
		98 => "0+2+32+64",
		99 => "1+2+32+64",
		100 => "0+4+32+64",
		101 => "1+4+32+64",
		102 => "0+2+4+32+64",
		103 => "1+2+4+32+64",
		104 => "0+8+32+64",
		105 => "1+8+32+64",
		106 => "0+2+8+32+64",
		107 => "1+2+8+32+64",
		108 => "0+4+8+32+64",
		109 => "1+4+8+32+64",
		110 => "0+2+4+8+32+64",
		111 => "1+2+4+8+32+64",
		112 => "0+16+32+64",
		113 => "1+16+32+64",
		114 => "0+2+16+32+64",
		115 => "1+2+16+32+64",
		116 => "0+4+16+32+64",
		117 => "1+4+16+32+64",
		118 => "0+2+4+16+32+64",
		119 => "1+2+4+16+32+64",
		120 => "0+8+16+32+64",
		121 => "1+8+16+32+64",
		122 => "0+2+8+16+32+64",
		123 => "1+2+8+16+32+64",
		124 => "0+4+8+16+32+64",
		125 => "1+4+8+16+32+64",
		126 => "0+2+4+8+16+32+64",
		127 => "1+2+4+8+16+32+64",
		128 => "0+128",
		129 => "1+128",
		130 => "0+2+128",
		131 => "1+2+128",
		132 => "0+4+128",
		133 => "1+4+128",
		134 => "0+2+4+128",
		135 => "1+2+4+128",
		136 => "0+8+128",
		137 => "1+8+128",
		138 => "0+2+8+128",
		139 => "1+2+8+128",
		140 => "0+4+8+128",
		141 => "1+4+8+128",
		142 => "0+2+4+8+128",
		143 => "1+2+4+8+128",
		144 => "0+16+128",
		145 => "1+16+128",
		146 => "0+2+16+128",
		147 => "1+2+16+128",
		148 => "0+4+16+128",
		149 => "1+4+16+128",
		150 => "0+2+4+16+128",
		151 => "1+2+4+16+128",
		152 => "0+8+16+128",
		153 => "1+8+16+128",
		154 => "0+2+8+16+128",
		155 => "1+2+8+16+128",
		156 => "0+4+8+16+128",
		157 => "1+4+8+16+128",
		158 => "0+2+4+8+16+128",
		159 => "1+2+4+8+16+128",
		160 => "0+32+128",
		161 => "1+32+128",
		162 => "0+2+32+128",
		163 => "1+2+32+128",
		164 => "0+4+32+128",
		165 => "1+4+32+128",
		166 => "0+2+4+32+128",
		167 => "1+2+4+32+128",
		168 => "0+8+32+128",
		169 => "1+8+32+128",
		170 => "0+2+8+32+128",
		171 => "1+2+8+32+128",
		172 => "0+4+8+32+128",
		173 => "1+4+8+32+128",
		174 => "0+2+4+8+32+128",
		175 => "1+2+4+8+32+128",
		176 => "0+16+32+128",
		177 => "1+16+32+128",
		178 => "0+2+16+32+128",
		179 => "1+2+16+32+128",
		180 => "0+4+16+32+128",
		181 => "1+4+16+32+128",
		182 => "0+2+4+16+32+128",
		183 => "1+2+4+16+32+128",
		184 => "0+8+16+32+128",
		185 => "1+8+16+32+128",
		186 => "0+2+8+16+32+128",
		187 => "1+2+8+16+32+128",
		188 => "0+4+8+16+32+128",
		189 => "1+4+8+16+32+128",
		190 => "0+2+4+8+16+32+128",
		191 => "1+2+4+8+16+32+128",
		192 => "0+64+128",
		193 => "1+64+128",
		194 => "0+2+64+128",
		195 => "1+2+64+128",
		196 => "0+4+64+128",
		197 => "1+4+64+128",
		198 => "0+2+4+64+128",
		199 => "1+2+4+64+128",
		200 => "0+8+64+128",
		201 => "1+8+64+128",
		202 => "0+2+8+64+128",
		203 => "1+2+8+64+128",
		204 => "0+4+8+64+128",
		205 => "1+4+8+64+128",
		206 => "0+2+4+8+64+128",
		207 => "1+2+4+8+64+128",
		208 => "0+16+64+128",
		209 => "1+16+64+128",
		210 => "0+2+16+64+128",
		211 => "1+2+16+64+128",
		212 => "0+4+16+64+128",
		213 => "1+4+16+64+128",
		214 => "0+2+4+16+64+128",
		215 => "1+2+4+16+64+128",
		216 => "0+8+16+64+128",
		217 => "1+8+16+64+128",
		218 => "0+2+8+16+64+128",
		219 => "1+2+8+16+64+128",
		220 => "0+4+8+16+64+128",
		221 => "1+4+8+16+64+128",
		222 => "0+2+4+8+16+64+128",
		223 => "1+2+4+8+16+64+128",
		224 => "0+32+64+128",
		225 => "1+32+64+128",
		226 => "0+2+32+64+128",
		227 => "1+2+32+64+128",
		228 => "0+4+32+64+128",
		229 => "1+4+32+64+128",
		230 => "0+2+4+32+64+128",
		231 => "1+2+4+32+64+128",
		232 => "0+8+32+64+128",
		233 => "1+8+32+64+128",
		234 => "0+2+8+32+64+128",
		235 => "1+2+8+32+64+128",
		236 => "0+4+8+32+64+128",
		237 => "1+4+8+32+64+128",
		238 => "0+2+4+8+32+64+128",
		239 => "1+2+4+8+32+64+128",
		240 => "0+16+32+64+128",
		241 => "1+16+32+64+128",
		242 => "0+2+16+32+64+128",
		243 => "1+2+16+32+64+128",
		244 => "0+4+16+32+64+128",
		245 => "1+4+16+32+64+128",
		246 => "0+2+4+16+32+64+128",
		247 => "1+2+4+16+32+64+128",
		248 => "0+8+16+32+64+128",
		249 => "1+8+16+32+64+128",
		250 => "0+2+8+16+32+64+128",
		251 => "1+2+8+16+32+64+128",
		252 => "0+4+8+16+32+64+128",
		253 => "1+4+8+16+32+64+128",
		254 => "0+2+4+8+16+32+64+128",
		255 => "1+2+4+8+16+32+64+128",
	);
	 
	 return \%SAMFlagTableHsh;
}
#################################################Read and return the sequence of the reference #######################################################
sub readMultiFastaFile {

	my $fastaFile = $_[0];
	my ($seq, $seqName, %fastaHsh);
	my $i = 0;
	open (INFILE, $fastaFile) || die "Can't read $fastaFile.\n";
	print "Reading $fastaFile into a hash.\n";	
	chomp (my $theCurntLine = <INFILE>); #get the first line
	while (my $theNextLine = <INFILE>) {
		chomp $theNextLine;
		
		#---Only two types of line in current line, the header or seq
		if ($theCurntLine =~ m/^>/) {#-- header line
			my @theLineSplt = split (/\|/, $theCurntLine);
			$seqName = $theLineSplt[0]; #---get the second tag
			$seqName =~ s/ //g; #---remove space
			$seqName =~ s/>//g; #---remove space
		} else {#--seq line
			$seq = $seq.$theCurntLine;
		}
		
		#---check if next line has a > or that's the end of file
		if ($theNextLine =~ m/^>/) {
			$fastaHsh{$seqName} = $seq;
			$seq = "";
		} elsif (eof(INFILE)) {#---this is the last line
			$seq = $seq.$theNextLine;
			$fastaHsh{$seqName} = $seq;
		}
		
		#---next line becomes current line
		$theCurntLine = $theNextLine;
	}
	close INFILE;
	
	return (\%fastaHsh);
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];

	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptPathXext = $0;
		$scriptPathXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptPathXext.cmd.log.txt"); #---append the CMD log file
		open (OUTDIRCMDLOG, ">$outDir/run.cmd.log.txt"); #---make the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;
		print CMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print OUTDIRCMDLOG "[".$runTime."]\t"."perl $0 ".(join " ", @ARGV)."\n";
		close OUTDIRCMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}

}
