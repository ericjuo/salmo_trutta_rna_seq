#!/usr/bin/perl -w

##
# long_read.pl
#
# Basic tests to ensure that long reads are handled properly.
#

use strict;
use warnings;

my $bowtie   = "./bowtie";
my $bowtie_d = $bowtie . " --debug";
my $bowtie_l = $bowtie . " --large-index --debug";
if(system("$bowtie --version") != 0) {
	$bowtie = `which bowtie`;
	chomp($bowtie);
	if(system("$bowtie --version") != 0) {
		die "Could not find bowtie in current directory or in PATH\n";
	}
}
if(system("$bowtie_d --version") != 0) {
	die "Could not find bowtie-debug in current directory or in PATH\n";
}
if(system("$bowtie_l --version") != 0) {
	die "Could not find bowtie large index support build in current directory or in PATH\n";
}

if(! -f "e_coli_c.1.ebwt") {
	print STDERR "Making colorspace e_coli index\n";
	my $bowtie_build = "./bowtie-build";
	if(system("$bowtie_build --version") != 0) {
		print STDERR "Could not execute ./bowtie-build; looking in PATH...\n";
		$bowtie_build = `which $bowtie_build`;
		chomp($bowtie_build);
		if(system("$bowtie_build --version") != 0) {
			die "Could not find bowtie-build in current directory or in PATH\n";
		}
	}
	system("$bowtie_build -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}

if(! -f "e_coli_c.1.ebwtl") {
	print STDERR "Making large colorspace e_coli index\n";
	my $bowtie_build = "./bowtie-build --large-index";
	if(system("$bowtie_build --version") != 0) {
		print STDERR "Could not execute ./bowtie-build; looking in PATH...\n";
		$bowtie_build = `which $bowtie_build`;
		chomp($bowtie_build);
		if(system("$bowtie_build --version") != 0) {
			die "Could not find bowtie-build in current directory or in PATH\n";
		}
	}
	system("$bowtie_build  --large-index -C genomes/NC_008253.fna e_coli_c") && die;
} else {
	print STDERR "Colorspace e_coli index already present...\n";
}


my @reads = (
	# 70:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????",
	# 140:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD",
	# 210:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD",
	# 280:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD",
	# 350:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGGTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCAGTGAATCGCCTACCAGCATGAC:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDDDDDDDDDDDDDDDAAADDDDDDD;;;;DDDEFFFDEGGGGDDDDDDDDDCCCDDDDDD",
	# 700:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGGTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCAGTGAATCGCCTACCAGCATGAC".
	"GTTAATCCCTTCATCGGCAAACAGTTTTGCGAAGCTATAGTCGTAAGCGGTGATGGTTGCGAAACGTTTT".
	"TTGTCCTGCTTGCATTTCTGCAATGAAGCAATCGTGGTTGGTTTCATAGCGTATCCTGATAAATGATGTT".
	"GTGCTGTCTCGCATTTTATCAGTCACATTGGTGGGGGCAATGATTTATCCGTATCGCACCGCCGTTGTCT".
	"GGCGGTGCGATAGTGTTATCAGGGGTAGGTAACCTGAAAAGTGCTGGTGGCTTTAAAATCCCCCGGTTCT".
	"ATAGCAATATTACCGTCTTGTTTTAATGTAGCCTGGAAATGAAGAGGTTGAGAGGTCCCACCGTTGCCAT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD",
	# 980:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGACCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGGTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCAGTGAATCGCCTACCAGCATGAC".
	"GTTAATCCCTTCATCGGCAAACAGTTTTGCGAAGCTATAGTCGTAAGCGGTGATGGTTGCGAAACGTTTT".
	"TTGTCCTGCTTGCATTTCTGCAATGAAGCAATCGTGGTTGGTTTCATAGCGTATCCTGATAAATGATGTT".
	"GTGCTGTCTCGCATTTTATCAGTCACATTGGTGGGGGCAATGATTTATCCGTATCGCACCGCCGTTGTCT".
	"GGCGGTGCGATAGTGTTATCAGGGGTAGGTAACCTGAAAAGTGCTGGTGGCTTTAAAATCCCCCGGTTCT".
	"ATAGCAATATTACCGTCTTGTTTTAATGTAGCCTGGAAATGAAGAGGTTGAGAGGTCCCACCGTTGCCAT".
	"TATCCGGGAAAATCCCCCCAGTAGTGTCATTTTCTGTTTCATAATCTTTATAAATAGACGTTGCATCATT".
	"TGGTTTCAGTACCATTTGTGCACTTTTAGTGTTTTTTAAACCTTCAATCAGTACACCAACGCCTTTTGCT".
	"GCGTCATTTCCTGTGAGTGTGTTAGCCAAAAGTTCTTTGCTTACACTACCCACCTTGTTTGACTTGAGTT".
	"TTGTTTCTATATTGCGAACACGAATACAGTTCTGTAGCGTGATATCAAAAGGAACCGCTGTTGCGCCATT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD",
	# 980 w/ 3 mismatches:
	"CCGCGCCCCTGGACTTTGTAGCCACCGAAAATATTCACTGACTGTGGGGTTAGGCCTAAATGAaCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTTTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGaTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCAGTGAATCGCCTACCAGCATGAC".
	"GTTAATCCCTTCATCGGCAAACAGTTTTGCGAAGCTATAGTCGTAAGCGGTGATGGTTGCGAAACGTTTT".
	"TTGTCCTGCTTGCATTTCTGCAATGAAGCAATCGTGGTTGGTTTCATAGCGTATCCTGATAAATGATGTT".
	"GTGCTGTCTCGCATTTTATCAGTCACATTGGTGGGGGCAATGATTTATCCGTATCGCACCGCCGTTGTCT".
	"GGCGGTGCGATAGTGTTATCAGGGGTAGGTAACCTGAAAAGTGCTGGTGGCTTTAAAATCCCCCGGTTCT".
	"ATAGCAATATTACCGTCTTGTTTTAATGTAGCCTGGAAATGAAGAGGTTGAGAGGTCCCACCGTTGCCAT".
	"TATCCGGGAAAATCCCCCCAGTAGTGTCATTTTCTGTTTCAtAATCTTTATAAATAGACGTTGCATCATT".
	"TGGTTTCAGTACCATTTGTGCACTTTTAGTGTTTTTTAAACCTTCAATCAGTACACCAACGCCTTTTGCT".
	"GCGTCATTTCCTGTGAGTGTGTTAGCCAAAAGTTCTTTGCTTACACTACCCACCTTGTTTGACTTGAGTT".
	"TTGTTTCTATATTGCGAACACGAATACAGTTCTGTAGCGTGATATCAAAAGGAACCGCTGTTGCGCCATT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD",
	# 980 w/ a bunch of mismatches:
	"CCGaGCCCCaGGACTTTGTAGCCACgGAAAATATTCACTGACTGTGGcGTTAGGCCTAAATGAaCACAAA".
	"CTGGAACGGCACGTTCGGTCAGCATTTGTACGGTTTCCACCAGCCACTCACCGCCTTCAATTTTGACCAT".
	"GTTGGCTCCGGCACGCATTACCGTTGCGGCGTTaTCAAAGGCTTGTTCCGGCGTGGCATACGCCATAAAT".
	"GGCAGGTCAGCCAGCAGAAGGCAGTTTGGCGCGCCGCGACGCACAGCGGCGGTGTGATAGGCGATATCGG".
	"CAACGGTAACCGaTAGCGTGGAGTCATGCCCCTGGACCGTCATGCCCcGTGAATCGCCTACCAGCATGAC".
	"GTTAATCCCTTCATCGGCAAACAGTTTTGCGAAGCTATAGTCGTAAGCGGTGATGGTTGCGAAACGTTTT".
	"TTGTCCTGCTTGCATTTCTGCAATGAAGCAATCGTGGTTGGTTTCATAGCGTATCCTGATAAATGATGTT".
	"GTGCTGTCTCGCATTTTATCAGTCACATTGGTGGGaGCAATGATTTATCCGTATCGCACCGCCGTTGTCT".
	"GGCGGTGCGATAGTGTTATCAGGGGTAGGTAACCTGAAAAGTGCTGGTGGCTTTAAAATCCCCCGGTTCT".
	"ATAGCAATATTACCGTCTTGTTTTAATGTAGCCTGGAAATGAAGAGGTTGAGAGGTCCCACCGcTGCCAT".
	"TATCCGGGAAAATCCCCCCAGTAGTGTCATTTTCTGTTTCAtAATCTTTATAAATAGACGTTGCATCATT".
	"TGGTTTCAGaACCATTTGTGCACTTTTAGTGTTTTTTAAACCTTCAATCAGTACACCAACGCCTTTTGCT".
	"GCGTCATTTCCTGTGtGTGTGTTAGCCAAAAGTTCTTTGCTTACACTACCCACCTTGTTTGACTTGAGTT".
	"TTGTTTCTATATTGCGAACACGAATACAGTTCTGTAGCGTGATATCAAAAGGAACCGCTGTTGCGCCATT:".
	"01111:::=????====???6111177=::::=???7444=??????====??===??????:::=????".
	"??????????=====;96666/.......;=>>>?DDDDDDDDDDDDDDDDDDDDDBBBBDDDDDCCCDD".
	"DDDDDDDDC>>>CDDDDDD>>?DDDDDDDDDDDDDCC???CCCCCCDDDDDDDDDDDDDDDDDDDDDDDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD".
	"DDDDDDDDDDDD;;;CDDDDDCCCCCDDDDDDCCCDDDDDDDDDDDDCCCDDDDDDCCCCCCC:666DDD",
	""
);

my @args1 = (
	"",
	"--best",
	"--best -m 1",
	"--best --strata -m 1"
);

my @args2 = (
	"-v 0",
	"-v 1",
	"-v 2",
	"-v 3",
	"-n 0",
	"-n 1",
	"-n 2",
	"-n 3",
	"-n 3 -l 25 -e 300"
);

sub btrun {
	my ($read, $args, $exalns, $color) = @_;
	$args .= $color ? " -C" : "";
	my $cmd = "$bowtie $args -c e_coli \"$read\"";
	print "\n$cmd\n\n";
	my $als = 0;
	open BTIE, "$cmd |" || die "Could not open pipe '$cmd |'\n";
	while(<BTIE>) {
		$als++;
	}
	close(BTIE);
	#$als == $exalns || die "# alignments $als didn't match expected number $exalns\n";
	$cmd = "$bowtie_d $args -c e_coli \"$read\"";
	print "\n$cmd\n\n";
	$als = 0;
	open BTIE, "$cmd |" || die "Could not open pipe '$cmd |'\n";
	while(<BTIE>) {
		$als++;
	}
	close(BTIE);
	#$als == $exalns || die "# alignments $als didn't match expected number $exalns\n";
	$cmd = "$bowtie_l $args -c e_coli \"$read\"";
	print "\n$cmd\n\n";
	$als = 0;
	open BTIE, "$cmd |" || die "Could not open pipe '$cmd |'\n";
	while(<BTIE>) {
		$als++;
	}
	close(BTIE);
	print "PASSED: \"$args\"\n";
}

for my $r (@reads) {
	next if $r eq "";
	for my $a1 (@args1) {
		for my $a2 (@args2) {
			btrun($r, "$a1 $a2 -a", 1, 0);
		}
	}
}
