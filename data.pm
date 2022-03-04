package data;

###################################################################################################
#
# Moduel name: data.pm
#
# This module include sub routine for data treatment:
#     indexfastq (indexing fastq data file),
#     indexfasta (indexing fasta data file),
#     averq (calculating the average of q-score),
#     retrievefasta (retrieve fasta sequence),
#     retrievefastq (retrieve fastq sequence).
#     
#
# Written by 
#                 Chenhong Li
#                 Shanghai Ocean University, China
#                 Created on Oct 2013
#                 Last modified on
#
###################################################################################################

use strict;
use warnings;


###################################################################################################
#
# sub routine: indexfastq
#
# This sub is make index file of fastq files.
#
# Input data include infile name and the index outfile name.
#
# Output file is the index file.
#
###################################################################################################

sub indexfastq {

    my $infile = shift(@_); #get the name of the infile
    my $index = shift(@_); #get the output index file

    my ($INFILE, $INDEX);

    open $INFILE, "<$infile" or die ("Cannot open $infile for reading ($!)");
    open $INDEX, ">$index" or die ("Cannot open $index for writing ($!)");

    my ($id);
    my $file_pointer = 0;
    my $file_pointerlag = 0;
    
    while (my $line = readline ($INFILE)) {
        chomp $line;
        if ($line =~ /^@\S+\s+\S+/) { # if we find @
            if ($id) { # if id exists we first need to write old information
        	    print $INDEX "$id\t$file_pointer\n";
        	    $file_pointer = $file_pointerlag;
            }
            ($id) = $line =~ /^@(\S+)/; #extract the new id
        }
        $file_pointerlag += length("$line\n"); #add value to file_pointer1
    }
    # We need to store the last sequence as well
    print $INDEX "$id\t$file_pointer\n";
        
    
    close ($INFILE) or die "Can't close the $infile file!!!";
    close ($INDEX) or die "Can't close the $index file!!!";
    
    return(1);

}

###################################################################################################
#
# sub routine: indexfasta
#
# This sub is make index file of fasta files.
#
# Imput data include infile name and the index outfile name.
#
# Output file is the index file.
#
###################################################################################################

sub indexfasta {

    my $infile = shift(@_); #get the name of the infile

    my ($OLD_FILE, $INDEX_FILE);

    my ($old, $index);
    
    $old = "$infile.fas";
    $index = "$infile.index";

    open $OLD_FILE, "<$old" or die ("Cannot open $old for reading ($!)");
    open $INDEX_FILE, ">$index" or die ("Cannot open $index for writing ($!)");

    my ($id, $seq);
    my $file_pointer = 0;
    while (my $line = readline ($OLD_FILE)) {
        chomp $line;
        if ($line =~ /^>/) { # if we find >
            if ($id) { # if id exists we first need to write old information
        		print $INDEX_FILE "$id\t$file_pointer\n";
                $file_pointer += length(">$id\n$seq");
            } 
        	# .. And then extract the new id and set sequence to empty
        	($id) = $line =~ /^>(\S+)/;
        	$seq = "";
        } else {
        		$seq = $seq . $line. "\n";
        }
    }
    # We need to store the last sequence as well
    print $INDEX_FILE "$id\t$file_pointer\n";

    close ($OLD_FILE) or die "Can't close the old file!!!";
    close ($INDEX_FILE) or die "Can't close the new file!!!";
    
    
    return(1);

}

###################################################################################################
#
# sub routine: averq
#
# This sub is used to calculate the average of q-score.
#
# Imput data include string of ascII q-score.
#
# Return is the average value of the q-score.
#
###################################################################################################

sub averq {

    my $string = shift(@_); #get the $string from main program
    
    my @array = split(//, $string);
    my $sum = 0;
    my $average;
    my $numofchar = scalar(@array);
    
    foreach my $char (@array) {
        $sum = $sum + (ord($char) - 33);
    }
    
    $average = $sum/$numofchar;
    
    return ($average);

}


###################################################################################################
#
# sub routine: retrieve fasta sequence
#
# This sub is used to retrieve the sequences from a fasta file.
#
# Imput data include the name of the fasta file and associated index file.
#
# Return is the sequence in fasta formate.
#
###################################################################################################

sub retrievefasta {

    my $seqid = shift(@_); #get the id of the sequence from the main program
    my $seqfile =  shift(@_); #get the name of the fasta sequence file from the main program
    my $indexfile =  shift(@_); #get the name of the index file from the main program
    

    my ($SEQ_FILE, $INDEX_FILE);
    open ($SEQ_FILE, "<$seqfile") or die "Cannot open $seqfile for reading ($!)";
    open ($INDEX_FILE, "<$indexfile") or die "Cannot open $indexfile for reading ($!)";


    # First read the index file
    my %index;
    my $id;
    while (my $line = readline ($INDEX_FILE)) {
    	chomp $line;
    	my ($id, $position) = split /\t/, $line;# the id in index file is longer than in embl file
    	$index{$id} = $position;
    }


    # and now we have superfast access to our sequences:
    seek $SEQ_FILE, $index{$seqid}, 0;  # let’s jump to our seq
    my $first = readline($SEQ_FILE);  # first line is header, skip
    my $second = readline($SEQ_FILE);	# second line is the whole sequence
    my $seq = "$first" . "$second";

    close ($SEQ_FILE) or die "Can't close the new file!!!";
    close ($INDEX_FILE) or die "Can't close the index file!!!";
    
    return ($seq);

}



###################################################################################################
#
# sub routine: retrieve fastq sequence
#
# This sub is used to retrieve the sequences from a fastq file.
#
# Imput data include the name of the fastq file and associated index file.
#
# Return is the sequence in fastq formate.
#
###################################################################################################

sub retrievefastq {

    my $seqid = shift(@_); #get the id of the sequence from the main program
    my $seqfile =  shift(@_); #get the name of the fasta sequence file from the main program
    my $indexfile =  shift(@_); #get the name of the index file from the main program
    

    my ($SEQ_FILE, $INDEX_FILE);
    open ($SEQ_FILE, "<$seqfile") or die "Cannot open $seqfile for reading ($!)";
    open ($INDEX_FILE, "<$indexfile") or die "Cannot open $indexfile for reading ($!)";


    # First read the index file
    my %index;
    my $id;
    while (my $line = readline ($INDEX_FILE)) {
    	chomp $line;
    	my ($id, $position) = split /\t/, $line;# the id in index file is longer than in embl file
    	$index{$id} = $position;
    }


    # and now we have superfast access to our sequences:
    seek $SEQ_FILE, $index{$seqid}, 0;  # let’s jump to our seq
    my $first = readline($SEQ_FILE);  # first line is header, skip
    my $second = readline($SEQ_FILE);	# second line is the whole sequence
    my $third = readline($SEQ_FILE);	# third line is "+"
    my $fourth = readline($SEQ_FILE);	# fourth line is the quality score
    my $seq = "$first" . "$second" . "$third" . "fourth";

    close ($SEQ_FILE) or die "Can't close the new file!!!";
    close ($INDEX_FILE) or die "Can't close the index file!!!";
    
    return ($seq);

}


1