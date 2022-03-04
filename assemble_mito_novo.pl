#!/usr/bin/perl -w
#This script is written by Wenjun Chen(LMSE, Shanghai Ocean University) to assemble the mitochondrial genome of vertebrates for LMSE gene-capture sequencing results.
#

use strict;
use warnings;
use Getopt::Long;
use data;

my $query;                     #to store the filename of pair-end sequences;
my $blastref;                  #to store the filename of BLAST reference;
my $help;                      #show usage;

my $opt = GetOptions( 'query:s', \$query,
                      'reference:s', \$blastref,
                      'help!', \$help);
					  



if (!($opt && $query) || $help) {                                         #check for the required inputs;
   print STDERR "\nExample usage:\n";
   print STDERR "\nperl $0 -query taxon.txt -reference ref.fas\n\n";
   print STDERR "Options:\n";
   print STDERR "        -query = name of the sequencing data";
   print STDERR "        -reference = a fasta file with sequences of 11 sites for BLAST.\n";
   print STDERR "\n\nPlease check readme.txt for further information.\n";
   
   exit;
}


my $queryfile = $query;                                          #get the input file names;
my $blastreffile = $blastref;
my @sitenames;
my $BLASTREF;                                                             #get the site names to assemble the blast results respectively;

open ($BLASTREF,"<$blastreffile") or die "Can't open the blast reference file!";

while (my $line = readline ($BLASTREF)){
	
	chomp $line;
	if ($line =~ /^>(\S+)/){
		
		my $sitename = substr ($line,1);
		push (@sitenames,$sitename);
		
	}
	
}

close $BLASTREF or die "Can't close the blast reference file!";


my @taxon;
my $TAXONIN;

open ($TAXONIN,"<$queryfile") or die "Can't open the query file!";        #split the taxon names through the taxon file;

while (my $taxonname = readline($TAXONIN)){
	
	chomp $taxonname;
	push (@taxon,$taxonname);
	
}

close $TAXONIN or die "Can't close the query file";


foreach my $taxon (@taxon){                                               #we assemble each taxon respectively, but the functions can do all the taxons alone;
	
	my @taxon_rmdup = rmdup($taxon);                                      #first we need to remove the duplications and make database for blast and locating;
    my @taxon_blast = blastn($blastreffile,@taxon_rmdup);                 #then blast the 11sites to the pair-end reads;
    my $taxon_blast = $taxon_blast[0];                                    #get the filename;
	my @taxon_blastout = blastout($taxon_blast,@taxon_rmdup);             #output the sequences of blast results;
	cleandb($taxon);                                                      #delete the database file which we won't use any longer;
	
	my @trinity_result;
	
	foreach my $sitename (@sitenames){                                    #first we assemble each blast result by Trinity and collect the results;
		
		my $R1 = $sitename . "_R1.fq";
		my $R2 = $sitename . "_R2.fq";
		my $trinity_output = "Trinity" . $sitename;
		
		`Trinity --seqType fq --CPU 12 --max_memory 10G --left $R1 --right $R2 --output $trinity_output`;
		`rm -rf $R1 $R2`;
		
		push (@trinity_result,$trinity_output);
		
	}
	
	my $merged_sequences = extract(@trinity_result);                      #merge all of the assemble sequences;
	my @database_name = makedb($merged_sequences);                        #make a new database to select the best sequences for each site;
	@taxon_blast = blastn($blastreffile,@database_name);                  #blast all the sites to the assemble results to select the best match sequences;
	my @ref = blastsel($merged_sequences,@taxon_blast);
	my $ref = $ref[0];
	cleandb($merged_sequences);                                           #after get the reference for mapping, we delete the files which we won't use any longer;
	$merged_sequences = $merged_sequences . ".fas";
	
	my @clean_files = (@taxon_blast,$merged_sequences);
	
	foreach my $clean (@clean_files){
		
		`rm -rf $clean`;
		
	}
	
	my $novoconfig = mknovoconfig($taxon);                               #make config files for seeds of each taxon
	my $novoassemble = novoassemble($taxon,$ref,$novoconfig);            #use novoplasty to extend
}


exit 0;




#this function is to remove the duplication sequences;
sub rmdup{
	
#input taxon name like "denio_rerio" and return the taxon names which has done rmdup；
	
	my @taxon = @_;
	my @taxon_rmdup;
	
	foreach my $taxon (@taxon) {
		
		#get the input and output filenames for each taxon;
		
		my $infile1 = $taxon . "_R1.fq";
		my $infile2 = $taxon . "_R2.fq";
		my $newfile = $taxon . ".fas";
		my $newfile1 = $taxon . "_rmdup_R1.fq";
		my $newfile2 = $taxon . "_rmdup_R2.fq";
		my $index1 = $taxon . "_rmdup_R1.index";
		my $index2 = $taxon . "_rmdup_R2.index";
		
		my ($INFILE1, $INFILE2, $NEWFILE, $NEWFILE1, $NEWFILE2);
		
		#open the input and output files;

        open $INFILE1, "<$infile1" or die ("Cannot open $infile1 for reading ($!)");
        open $INFILE2, "<$infile2" or die ("Cannot open $infile2 for reading ($!)");
        open $NEWFILE, ">$newfile" or die ("Cannot open $newfile for writing ($!)");
        open $NEWFILE1, ">$newfile1" or die ("Cannot open $newfile for writing ($!)");
        open $NEWFILE2, ">$newfile2" or die ("Cannot open $newfile for writing ($!)");
		
		my ($id1,$id2);                                                             #to store the sequence names of pair-end sequences;
		
		my %uniseq;                                                                 #hash for storing unique sequence, hash key: 20bp first read + 20bp second read;
		my ($key1, $key2, $key, $seq1, $seq2, $quality1, $quality2, $quality);      #to store the sequences and qualities to compare;
		
		while (my $line1 = readline($INFILE1)){
			
			chomp $line1;
			if ($line1 =~ /^@\S+\s+\S+/){                                           #if the line is starting with @;
				if ($id1){                                                          #if there is a sequence name;
					$key = "$key1" . "$key2";                                       #merge the starting and ending sequences of the pair-end reads;
					$quality = "$quality1" . "$quality2";                           #merge the quality of the pair-end reads;
					if ($uniseq{$key}){                                             #if the merged sequence exsits,which means these two sequences are duplicated;
						my $qold = &data::averq($uniseq{$key}->{quality});          #we use the data.pm module to calculate the average quality of each sequence;
						my $qnew = &data::averq($quality);
						if ($qold < $qnew) {                                        #we keep the one which has a higher quality;
							$uniseq{$key} = { id => $id1,
							                  quality => $quality,
											  seq1 => $seq1,
											  seq2 => $seq2
						                    };
						}
					}else{                                                          #if there is not a duplication,and the two sequences from two files are pair-end, keep them;
						if ($id1 eq $id2) {
							$uniseq{$key} = { id => $id1,
							                  quality => $quality,
											  seq1 => $seq1,
											  seq2 => $seq2
						                    };
						}else{                                                      #if the sequences are not pair-end, die and print the error;
							die "The $taxon two reads files are not consistene!\n";
						}
					}
				}
				($id1) = $line1 =~ /^@(\S+)/;                                       #get the next name of R1 sequence before the first space;
				$seq1 = readline ($INFILE1);                                        #get the sequence;
				chomp $seq1;
				$key1 = substr($seq1,0,20);                                         #cut the first 20bp to further compare;
				readline ($INFILE1);                                                #skip the '+';
				$quality1 = readline ($INFILE1);                                    #get the quality of the sequence;
				chomp $quality1;
			}
			
			my $line2 = readline ($INFILE2);                                        #get the same information of R2;
			chomp $line2;
			if ($line2 =~ /^@\S+\s+\S+/){
				($id2) = $line2 =~ /^@(\S+)/;
				$seq2 = readline($INFILE2);
				chomp $seq2;
				$key2 = substr($seq2,0,20);
				readline ($INFILE2);
				$quality2 = readline ($INFILE2);
				chomp $quality2;
			}
		}
		
		
		#now write to the output files;
		foreach my $key (sort (keys (%uniseq))){
			my $id = $uniseq{$key}->{id};
			my $seq1 = $uniseq{$key}->{seq1};
			my $seq2 = $uniseq{$key}->{seq2};
			
			my $q = $uniseq{$key}->{quality};
			my $q1 = substr($q, 0, length($seq1));
			my $q2 = substr($q, -length($seq2));
			
			print $NEWFILE ">$id\n$seq1\n$seq2\n";
			print $NEWFILE1 "\@$id 1\n$seq1\n+\n$q1\n";
			print $NEWFILE2 "\@$id 2\n$seq2\n+\n$q2\n";
			
		}
		
		
        close ($INFILE1) or die "Can't close the $taxon new file!!!";
        close ($INFILE2) or die "Can't close the $taxon new file!!!";
        close ($NEWFILE) or die "Can't close the $taxon new file!!!";
        close ($NEWFILE1) or die "Can't close the $taxon new file!!!";
        close ($NEWFILE2) or die "Can't close the $taxon new file!!!";
		
		&data::indexfasta($taxon);                                                  #make an index file for the new combined fasta file;

        `makeblastdb -dbtype nucl -in $newfile -out $taxon -max_file_sz 2GB`;       #make a blast database for the new fasta file;
	
	    &data::indexfastq($newfile1, $index1);                                      #index the read1 fastq;
    
        &data::indexfastq($newfile2, $index2);                                      #index the read2 fastq;
		
		push (@taxon_rmdup,$taxon);                                                 #push the taxon name whose duplication been removed;
	
	}
	
	
	return @taxon_rmdup;                                                            #return the taxon names;
	
}


#this function is to call BLASTN;
sub blastn{
	
#input the reference file name like "human.fas" and database name like "danio_rerio", and return the file name of blast result like "human.fas.danio_rerio.blast.txt";
	
	my $queryfile = shift @_;
	my @database = @_;
	my @taxon_blasted;
	foreach my $database (@database) {
		
		my $blastout = "$queryfile.$database.blast.txt";
		`blastn -query "$queryfile" -task blastn -db "$database" -out "$blastout" -word_size 7 -gapopen 5 -gapextend 2 -penalty -1 -reward 1 -evalue 0.0001 -outfmt 6`;
		push (@taxon_blasted,$blastout);
		
	}
	
	return @taxon_blasted;	

}


#this function is to output the sequences of the pair-end reads after the first blastn;
sub blastout{
	
#input the reference filename like "human.fas.danio_rerio.blast.txt" and the taxon name like "denio_rerio", and return the taxon name which has done blastout like "danio_rerio";
#will remove the useless index files;
	
	my $ref = shift @_;                                                             #in this script, we use only one reference;
	my @taxon_blasted = @_;
	my @taxon_return;
	
    foreach my $taxon_blastout (@taxon_blasted){                                    #we output the results for each single taxon;
	
	    my ($SEQ_FILE1, $INDEX_FILE1);                                              #get the file names of pair-end sequences and their indexes;
	    my $seqfile1 = $taxon_blastout . "_rmdup_R1.fq";
	    my $indexfile1 = $taxon_blastout . "_rmdup_R1.index";
	    open ($SEQ_FILE1,"<$seqfile1") or die "Can't open $seqfile1 for blastout!($!)";
	    open ($INDEX_FILE1,"<$indexfile1") or die "Can't open $indexfile1 for blastout!($!)";
	
	    my %index1;                                                                 #hash to store the sequence id and its position;
    	my $id1;
	    while (my $line = readline ($INDEX_FILE1)){                                 #readline and hash it;
	    	chomp $line;
		    my ($id1, $position1) = split /\t/, $line;
		    $index1{$id1} = $position1;
	    }
	
	    my ($SEQ_FILE2, $INDEX_FILE2);                                              #do the same thing to indexfile2;
	    my $seqfile2 = $taxon_blastout . "_rmdup_R2.fq";
	    my $indexfile2 = $taxon_blastout . "_rmdup_R2.index";
	    open ($SEQ_FILE2,"<$seqfile2") or die "Can't open $seqfile2 for blastout!($!)";
	    open ($INDEX_FILE2,"<$indexfile2") or die "Can't open $indexfile2 for blastout!($!)";
	
	    my %index2;
	    my $id2;
	    while (my $line = readline ($INDEX_FILE2)){
	    	chomp $line;
		    my ($id2, $position2) = split /\t/, $line;
		    $index2{$id2} = $position2;
	    }
	
	
	    my $BLASTOUT;
	    open ($BLASTOUT,"<$ref") or die "Can't open the blast_result file of $taxon_blastout!";
	    my ($geneidlag, $hitidlag, $OUTFILE1, $OUTFILE2);
	
	    my $line = readline ($BLASTOUT);                                            #get the first line of blast result;
	    my ($geneid, $hitid) = $line =~ /^(\S+)\s+(\S+)/;                           #get the site name and the hit sequence name;
	    my $outfile1 = "$geneid" . "_R1." . "fq";	                                #then open the output file to print pair-end sequences;
		my $outfile2 = "$geneid" . "_R2." . "fq";
	    open ($OUTFILE1,">$outfile1") or die "Can't open the blast_output file1 of $taxon_blastout!";
	    open ($OUTFILE2,">$outfile2") or die "Can't open the blast_output file2 of $taxon_blastout!";
		
	    $geneidlag = $geneid;                                                       #store the last site id sequence id;
	    $hitidlag = $hitid;
	
	    while (my $line = readline ($BLASTOUT)){
		
		    ($geneid, $hitid) = $line =~ /^(\S+)\s+(\S+)/;                          #read next line;
		
		    if ($hitid ne $hitidlag){                                               #if the sequence id are different;
			
			    seek $SEQ_FILE1, $index1{$hitidlag},0;                              #seek the sequence in the sequences file by using the index;
			    my $first = readline ($SEQ_FILE1);                                  #merge the information of this sequence and output it;
			    my $second = readline ($SEQ_FILE1);
			    my $third = readline ($SEQ_FILE1);
			    my $fourth = readline ($SEQ_FILE1);
			    my $seq = "$first" . "$second" . "$third" . "$fourth";
			    print $OUTFILE1 "$seq";
				
			    seek $SEQ_FILE2, $index2{$hitidlag},0;                              #and do the same thing in R2 file;
			    $first = readline ($SEQ_FILE2);
			    $second = readline ($SEQ_FILE2);
			    $third = readline ($SEQ_FILE2);
			    $fourth = readline ($SEQ_FILE2);
			    $seq = "$first" . "$second" . "$third" . "$fourth";
			    print $OUTFILE2 "$seq";
		    	
		    }
		
		    if ($geneid ne $geneidlag){                                             #if the site has changed;
			
		    	close ($OUTFILE1) or die "Can't close the blast_output file1 of $taxon_blastout!";
				close ($OUTFILE2) or die "Can't close the blast_output file2 of $taxon_blastout!";
		    	$outfile1 = $geneid . "_R1." . "fq";
			    $outfile2 = $geneid . "_R2." . "fq";
			    open ($OUTFILE1,">$outfile1") or die "Can't open the blast_output file1 of $taxon_blastout!";
				open ($OUTFILE2,">$outfile2") or die "Can't open the blast_output file2 of $taxon_blastout!";
			    $geneidlag = $geneid;                 
				#go on outputting;
		    }
		    
		    $hitidlag = $hitid;
		
	    }
	
        close ($SEQ_FILE1) or die "Can't close the R1 file of $taxon_blastout while outputting!";
        close ($INDEX_FILE1) or die "Can't close the index1 file of $taxon_blastout!";
        close ($SEQ_FILE2) or die "Can't close the R2 file of $taxon_blastout while outputting!";
        close ($INDEX_FILE2) or die "Can't close the index2 file of $taxon_blastout!";
		
		my $dbfile = $taxon_blastout . ".fas";
		`rm -rf $dbfile $indexfile1 $indexfile2 $ref`;                                      #delete the index file which we won't use any longer;
	    push (@taxon_return,$taxon_blastout);
		
    }
	
	return @taxon_return;
	
}


#this function is to make blastn database file of pair-end sequences and return databasefile name;
sub makedb{
	
#input what you chose to name the database like "danio_rerio" and return the database name like "danio_rerio";
	
	my @taxon = @_;
	my @database_name;
	
	foreach my $taxon(@taxon){
		
		my $filename = $taxon . ".fas";
		
		`makeblastdb -dbtype nucl -in $filename -out $taxon	-max_file_sz 2GB`;
		
		push (@database_name,$taxon);
		
	}
	
	return @database_name;
	
}


#this function is to put the assemble results of Trinity together and rename each sequence;
sub extract{
	
#input the filename like "danio_rerio.fas" and return the filename "merged";
#will remove the input directions;
	
	my $output = "merged.fas";
	my $return = "merged";
	
	my $OUT;
	open ($OUT,">$output") or die "Can't open the $output!";
	my $count = 1;
	
	foreach my $dir (@_){
		
		my $file = $dir . "/" . "Trinity.fasta";
		my $file_exist = -e $file;
		if ($file_exist){
	        my $TRINITY_RESULT;
	        open ($TRINITY_RESULT,"<$file") or die "Can't open the $file!";
		    my $firstline = readline($TRINITY_RESULT);
		    print $OUT ">sequence$count\n";
		    $count ++;
		
		    while (my $line = readline($TRINITY_RESULT)){
			
		    	chomp $line;
			
		    	if ($line =~ /^>.+/){
				
			    	print $OUT "\n>sequence$count\n";
			    	$count ++;
				
			    }else{
				
			    	print $OUT "$line";
				
			    }
			
		    }
		
		    close $TRINITY_RESULT or die "Can't close $file!";
		
		    print $OUT "\n";
		
		    `rm -rf $dir`;
		
	    }else{
			
			`rm -rf $dir`;
			next;
			
		}
		
	}
	
	close $OUT or die "Can't close the $output!";
	
	return $return;
	
}


#this function is to output the best match of each sites in the blast_output file for the first time;
sub blastsel{
	
#input the database name like "merged" and the filename of blast results like "denio_rerio.fas.merged.blast.txt", and return the filename like "mergedref.fas";
#will remove the "merged.fas" and blast result file;
	
	my $sequence_name = shift @_;                                                          #the first file name is the merged sequence file;
	my $sequence_file = $sequence_name . ".fas";
	my @taxon_blasted = @_;                                                                #then each taxon file;
	my @blast_outfile;                                                                     #the output file of each taxon;                                                                     
	
	foreach my $blast_output(@taxon_blasted){
		
		my %blast_select;                                                                  #to store the sequence's name of each site;
		my $IN;
		open ($IN,"<$blast_output") or die "Can't open the blast file $blast_output!";
		my $firstline = readline ($IN);                                                    #almost the same steps like blastout, but we only need the first sequence name of each site;
		my ($site, $id) = $firstline =~ /^(\S+)\s+(\S+)/;
		$blast_select{$id} = $site;
		
		while (my $line = readline ($IN)){
			
			my ($site_next, $id_next) = $line =~ /^(\S+)\s+(\S+)/;
			if ($site_next ne $site){
				
				if (not exists $blast_select{$id_next}){
					
					$site = $site_next;
					$blast_select{$id_next} = $site;
					
				}else{
					
					$site = $site_next;
					
				}
				
			}
			
		}
		
		close $IN or die "Can't close the blast file $blast_output!";                      #now all the names of sequences are selected, we need to output the sequences;
		
		my $MERGED;
		
		my $blast_outfile = $sequence_name . "ref" . ".fas";                               #make an output file;
		my $OUT;
		open ($OUT,">$blast_outfile") or die "Can't open $blast_outfile!";
		
		foreach my $key (sort keys %blast_select){
			
			open ($MERGED,"<$sequence_file") or die "Can't open $sequence_file!";          #open the sequence file;
			
			while (my $line = readline($MERGED)){
				
				my ($id) = $line =~ /^>(\S+)/;
				if ($id eq $key){
					
					my $seq = readline($MERGED);
					print $OUT "$line$seq";
					
				}else{
					
					readline($MERGED);
					
				}
				
			}
			
			close $MERGED or die "Can't close $sequence_file!";
			
		}
		
		
		close $OUT or die "Can't open $blast_outfile!";
		
		push (@blast_outfile,$blast_outfile);                                              #store the filenames;
		`rm -rf $sequence_file $blast_output`;
		
	}
	
	
	return @blast_outfile;                                                                 #return the filenames;
	
	
}


#this function is to remove the database files;
sub cleandb{
	
#input the database names and return the names of databases which have been removed;
	
	my @db_delete;
	my @db_cleaned;
	foreach my $database(@_){
		
		my $dbindex = $database . ".index";
		my $dbnhr = $database . ".nhr";
		my $dbnin = $database . ".nin";
		my $dbnsq = $database . ".nsq";
		push (@db_delete,$dbindex,$dbnhr,$dbnin,$dbnsq);
		push (@db_cleaned,$database);
	}
	
	foreach my $filename(@db_delete){
		
		`rm -rf $filename`;
		
	}
	
	return @db_cleaned;
	
}

#this function is to make config file;
sub mknovoconfig{
	
	my $filename = shift @_;
	my $novoR1 = $filename . "_rmdup_R1.fq";
	my $novoR2 = $filename . "_rmdup_R2.fq";
	my $configfile = "config_" . $filename . ".txt";
	my $outdir = "$filename" . "_result";
	`mkdir $outdir`;
	open (CONFIG,">$configfile") or die "can not open $configfile!";
	print CONFIG "Project:\n-----------------------\n";
	print CONFIG "Project name          = batch:batch_file_$filename.txt\n";
	print CONFIG "Type                  = mito\n";
	print CONFIG "Genome Range          = 14000-20000\n";
	print CONFIG "K-mer                 = 31\n";
	print CONFIG "Max memory            = 8\n";
	print CONFIG "Extended log          = 0\n";
	print CONFIG "Save assembled reads  = no\n";
	print CONFIG "Seed Input            = batch\n";
	print CONFIG "Extend seed directly  = no\n";
	print CONFIG "Reference sequence    = \n";
	print CONFIG "Variance detection    = no\n";
	print CONFIG "Chloroplast sequence  = \n";
	print CONFIG "Dataset 1:\n-----------------------\n";
	print CONFIG "Read Length           = 145\n";
	print CONFIG "Insert size           = 300\n";
	print CONFIG "Platform              = illumina\n";
	print CONFIG "Single/Paired         = PE\n";
	print CONFIG "Combined reads        = \n";
	print CONFIG "Forward reads         = $novoR1\n";
	print CONFIG "Reverse reads         = $novoR2\n";
	print CONFIG "Store Hash            =\n";
	print CONFIG "\n";
	print CONFIG "Heteroplasmy:\n-----------------------\n";
	print CONFIG "MAF                   = \n";
	print CONFIG "HP exclude list       = \n";
	print CONFIG "PCR-free              = yes\n";
	print CONFIG "\n";
	print CONFIG "Optional:\n-----------------------\n";
	print CONFIG "Insert size auto      = yes\n";
	print CONFIG "Use Quality Scores    = no\n";
	print CONFIG "Output path           = $outdir/\n";
	
	close CONFIG or die "can not close CONFIGFILE";
	return $configfile;
	
}

#this function is to extend mitochondrial sequence by novoplasty, please check the path of the software before using;
sub novoassemble{
	
	my $taxon = shift @_;
	my $seqfl = shift @_;
	my $cfg = shift @_;
	my $batchfile = "batch_file_" . $taxon . ".txt";
	my @seeds;
	
	open (SEQ,"<$seqfl") or die "can not open $seqfl!";
	
	while (my $line = readline(SEQ)){
		
		if ($line =~ /^>(\S+)/){
			
			chomp $line;
			$line = substr($line,1);
			my $seed = $line . ".fas";
			my $sq = readline(SEQ);
			chomp $sq;
			open (SEED,">$seed") or die "can not open $seed!";
			
			print SEED ">$line\n$sq\n";
			push(@seeds,$seed);
			
			close SEED;
		}
		
	}

	close SEQ;
	
	open (BAT,">$batchfile") or die "can not open $batchfile!";
	
	my $project_number = 1;
	foreach my $seed (@seeds){
		
		print BAT "Project$project_number\n$seed\n";
		$project_number ++;
		
	}
	
	close BAT;
	
	`perl /mnt/disk6/Lab_Users/cwj/software/NOVOPlasty-master/NOVOPlasty4.3.pl -c $cfg`;
    `rm -rf $seqfl $batchfile $cfg`;
	foreach my $sd(@seeds){
		
		`rm -rf $sd`;
		
	}
	
}
