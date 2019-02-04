#!/usr/local/bin/perl

##################################################################
=head1 # SEQUENCE RANDOMIZER #####################################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program accepts an input file filled with nucleotide
or protein sequences (it can actually accept any string of characters).
It then randomizes the characters within the sequences and outputs the
new randomized sequences to an output file.

The script allows the user to generate one or multiple randomized output
sequence for each input sequence. 

=head2 INPUT:

After the executable, one should type the path and filename of the
input file as an argument. The input file is a tab-delineated file,
where the first column is a unique ID for the gene/protein and the
second column is the sequence. 

=head2 REQUIREMENTS:

None

=head2 DEFAULT ARGUMENT VALUES

=cut

use warnings;
use strict;

my ($output_file_name, $output_seq_per_input_seq);

$output_file_name = 'OUTPUT_FILE_NAME.tab';
$output_seq_per_input_seq = 50;

=head2 ARGUMENT DESCRIPTIONS

$output_file_name = Path and filename of the output file containing all the
newly-calculated randomized sequences. This file is a tab-delineated file,
where the first column contains the original unique ID to which the listed
randomized sequence matches and the second column contains the newly
randomized sequence. 

$output_seq_per_input_seq = This the number of of randomized output
sequences that are generated for each input sequence. Must be an integer
value.

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($input_file_line,
    $original_sequence,
    $sequence_name,    
    $randomized_sequence,
    $sequence_number_count,
    );
  
open(RANDOM_SEQ_OUTPUT, ">", $output_file_name) || die("Cannot open the $output_file_name file: $!\n");

# We will now randomize each sequence.

while (defined( $input_file_line = <> ))
    {
    chomp $input_file_line;
    
    $input_file_line =~ m/^(\w+)\t([a-z]+)$/gi;
    
    $sequence_name = $1;
    $original_sequence = $2;
    
    for ($sequence_number_count = 1; $sequence_number_count <= $output_seq_per_input_seq; $sequence_number_count++)
        {
        $randomized_sequence = &randomizer($original_sequence);    
        
        # Now that the sequence is randomized, we can now print it
        # to the output file.
        
        print RANDOM_SEQ_OUTPUT "$sequence_name\t$randomized_sequence\n";
        };
    
    print "$sequence_name\n";    
    };

close RANDOM_SEQ_OUTPUT;

exit;

####################################################################################
# SUBROUTINES ######################################################################
####################################################################################

# RANDOMIZER

# The following subroutine will randomize the input sequence. That is,
# it will randomly order the characters that are present in the
# input string. The randomizing engine is the perl "rand" command. 

sub randomizer
    {
    
    # Input Variables
    
    my  ($original_sequence) = @_;
    
    # Subroutine Variables
    
    my  (@original_seq_array,
        $remaining_aa_number,
        $index,
        @randomized_seq_array, 
        $randomized_sequence,
        );
    
    @original_seq_array = split('', $original_sequence);
    
    # We can now begin randomizing the sequence.
    
    until ( scalar(@original_seq_array) == 0 )
        {    
        $remaining_aa_number = scalar(@original_seq_array);
        
        # The value that is output from rand() ranges from 0 to
        # LESS THAN the number input into rand() (in this case
        # $remaining_aa_number). Therefore, although we are using
        # the "scalar" command instead of "$#", we are still
        # giving rand() the correct range.
        
        $index = int( rand($remaining_aa_number) );
        
        # We will now place the amino acid at $index into the
        # new randomized array
        
        push( @randomized_seq_array, $original_seq_array[$index] );
        
        splice(@original_seq_array, $index, 1);
        };
    
    undef @original_seq_array;
    
    $randomized_sequence = join('', @randomized_seq_array);
    
    undef @randomized_seq_array;

    return $randomized_sequence;
    
    };

########################################################################################  
