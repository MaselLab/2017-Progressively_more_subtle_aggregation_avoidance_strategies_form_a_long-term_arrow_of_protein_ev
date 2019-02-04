#!/usr/local/bin/perl

#################################################################################
=head1 # HYDROPHOBIC DISPERSION CONTROL POLYPEPTIDE GENERATOR ###################
#################################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will calculate a scrambled polypeptide sequence that has
a similar hydrophobicity profile to that of a real protein. A hydrophobicity
profile states approximately where in the protein sequence a protein
is hydrophobic and where in the sequence it is hydrophilic.

The hydrophobicity dispersion (psi) is dependent upon the size of the
block that is being used. This program can accomodate a dispersion block
range. It will average the psis of each block length within the range
in order to obtain the mean psi for a given polypeptide. 

Each real protein will be extracted from a MySQL database. The
scrambled polypeptides possessing similar hydrophobic dispersion
will be returned to a flat file; they will
not be inserted directly into a/the MySQL database. Extracting
them to a flat file will permit the user to do as she/he pleases
with the scrambled protein sequences. Additionally, it is generally
not wise to insert new ROWS into a MySQL database from a script due
to the risk of errors and database contamination.

=head2 REQUIREMENTS:

This program requires the DBD::mysql module (and the
library [DBD::mysql directory]). Although DBD and DBI are both
part of the core modules, the DBD:mysql driver is not and must be
downloaded from CPAN.

Finally, no protein unique ID may contain a space or tab.

#####################################################################
=head1 # ARGUMENTS ##################################################
#####################################################################

The hardcoded values below are the default values for each variable.
Additionally, below these are descriptions of each variable.

=head2 ARGUMENT VALUES

The following input variables are lexical and can be utilized
thoughout the entire program.

=cut

use warnings;
use strict;

print "Time of program execution:\n";
main::print_time();
print "\n";

my  ($min_s,
    $max_s,
    $max_trial_scrambled_polypeptides,
    $output_file_name,
    $cores_to_use,

    $mysql_extraction_query,

    $database_host_ip_address,
    $mysql_account_username,
    $mysql_account_password,
    $mysql_database_name,
    );

# Default Equation Values:
# ________________________

$min_s = 6;
$max_s = 6;
$max_trial_scrambled_polypeptides = 1000;
$output_file_name = 'OUTPUT_FILE_NAME.tab';
$cores_to_use = 32;

# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT UniqueID, ProteinSequence FROM TableName WHERE Designation = "CodingProtein"'; 

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$[min|max]_s = The minimum and maximum number of amino acids per block (see Irback [2000]
for more information on the dispersion algorithm). The program will calculate
multiple block lengths for each sequence. These block lengths will range from the
$min_s to the $max_s block length. 

$max_trial_scrambled_polypeptides = For each coding protein, we will generate
numerous scrambled polypeptides from which the program will choose the best
psi value. This scalar is a number that indicates how many scrambled polypeptides
the program will generate per coding protein.

$output_file_name = This is the name of the output file where the
program will place the scrambled polypeptide sequences.

$cores_to_use = This program can be run in parallel. This scalar is the number of
cores one would like this program to use.

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. If the user wishes to calculate the
dispersion of only certain protein sequences, the user may add a
"WHERE" statement to the MySQL command.  

$database_host_ip_address = This is the IP address of the computer hosting
the MySQL database. If the MySQL database is located on the local computer
that is executing this program, the IP address should be set as
"localhost" or "127.0.0.1". 

$mysql_account_username = This is the username required to access the MySQL
database. Importantly, this is NOT the username required to access the
remote computer. That is, this is the MySQL username, NOT the computer login
username.

$mysql_account_password = This is the password required to access the MySQL
database. Importantly, this is NOT the password required to access the
remote computer. That is, this is the MySQL password, NOT the computer login
password.

$mysql_database_name = This is the name of the MySQL database containing the
information to be extracted. 

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($database_connection_object,
    $database_extraction_object,
    $row_array,
    
    $protein_sequences_in_parallel,
    $protein_sequences_for_one_core,
    $core_number,
    $reader_file_handle,
    $writer_file_handle,
    $pid,
    %child_process_ids,
    %child_file_handles,
    
    $unique_id,
    $protein_sequence,
    $max_s_for_sequence,
    $min_s_for_sequence,
    
    $truncated_sequences,
    $block_length,
    $hydrophobicity_scores_array,
    $psi,
    @all_truncated_psi,
    %protein_psi_per_s,
    
    $scrambled_polypeptide,
    %scrambled_psi_per_s,
    
    $protein_mean_psi,
    $scrambled_mean_psi,
    $best_mean_psi,
    $best_scrambled_polypeptide,
    
    $output_line,
    @control_sequences,
    @con_seq_from_child,
    );

# We must first access the MySQL database.

use DBI;

print 'Connecting to host database = ' . "$database_host_ip_address\n";

$database_connection_object = DBI->connect("dbi:mysql:$mysql_database_name;$database_host_ip_address",
                                           $mysql_account_username, $mysql_account_password);

print "Database connected\n";

# We now pre-prepare our MySQL commands that will be used to both
# extract and insert data.

$database_extraction_object = $database_connection_object->prepare($mysql_extraction_query);

$database_extraction_object->execute();

print "Database extraction command executed\n";

$core_number = 1;

while ($row_array = $database_extraction_object->fetchrow_arrayref)
    {
    # The $row_array is a reference to an array, each element of which
    # is a value from a row of the database. Specifically, each element
    # contains the value returned by the "SELECT" statement in the
    # $mysql_extraction_query. In this case, the zeroth
    # element is the UniqueID and the first element is the
    # ProteinSequence.

    $unique_id = $row_array->[0];    
    $protein_sequence = $row_array->[1];

    undef $row_array;
    
    $protein_sequences_in_parallel->{$core_number}->{$unique_id} = $protein_sequence;
    
    $core_number++;
    
    if ($core_number > $cores_to_use) {$core_number = 1};
    };

$database_connection_object->disconnect;

no DBI;  

# Now that the data is extracted from the MySQL database, we
# can begin running the calculations in parallel.

use IO::Handle;
  
CORE_FOREACH_LOOP: foreach $core_number (1..$cores_to_use)    
    {
    pipe ($reader_file_handle, $writer_file_handle);        
    $writer_file_handle->IO::Handle::autoflush(1);

    # We will generate an independent data segment to save
    # RAM. Only retaining the particular data segment for each
    # child process is much better than having to retain the
    # entire $protein_sequences_in_parallel data structure in
    # all the child processes.
    
    $protein_sequences_for_one_core = $protein_sequences_in_parallel->{$core_number};

    $pid = fork(); # "pid" is the "process identification number"
    defined $pid or die "Fork process failed: $!"; 

    if ($pid)
        {
        # This is the parent process.
        
        close $writer_file_handle;
        
        undef $protein_sequences_for_one_core;

        print "PARENT PROCESS: CORE_NUM: $core_number\tPID: $pid\tREADER_FILE_HANDLE: $reader_file_handle\n";

        $child_process_ids{$core_number} = $pid;       
        $child_file_handles{$core_number} = $reader_file_handle;
        
        undef $reader_file_handle;
        }
    
    else
        {
        # This is the child process.
        
        close $reader_file_handle;

        undef $protein_sequences_in_parallel;  

        while (($unique_id, $protein_sequence) = each %$protein_sequences_for_one_core)
            {    
            print "Unique_ID: $unique_id\n";  
            
            $best_mean_psi = 99999999999999;
            
            # The maximum s that can be utilized for a given protein
            # ($max_s_for_sequence) must be shorter than half the length of the
            # protein; therefore, we will temporarily shorten $max_s if necessary
            # to accommodate for shorter proteins. We must have the $max_s variable
            # so that $max_s_for_sequence can be reset for each protein.
            
            $max_s_for_sequence = $max_s;
            $min_s_for_sequence = $min_s;
            
            if ( length($protein_sequence) / 2 < $max_s_for_sequence ) { $max_s_for_sequence = length($protein_sequence) / 2 };
            
            if ($max_s_for_sequence < $min_s_for_sequence) {$min_s_for_sequence = $max_s_for_sequence};
            
            ######################
            # Coding Protein Psi #
            ######################    
            
            # First we must calculate the psi for the coding protein. This
            # code is not located in the "Scrambled Polypeptide Psi" loop
            # because there is no need to keep recalculating it.
            
            # If we calculate a "block" length of only one amino acid, the psi
            # for the protein sequence and all the scrambled polypeptides will be
            # the same because this psi score is dependent solely on amino acid
            # content (with a block length of 1).
            
            foreach $block_length ($min_s_for_sequence..$max_s_for_sequence)
                {        
                $truncated_sequences = &truncate_sequence($protein_sequence, $block_length);
                
                foreach (@$truncated_sequences)
                    {            
                    $hydrophobicity_scores_array = &hydrophobicity_scores($_);
                    
                    $psi = &calculate_psi($hydrophobicity_scores_array, $block_length);
                    
                    push(@all_truncated_psi, $psi);
                    
                    undef $hydrophobicity_scores_array;
                    undef $psi;                
                    };
                
                $protein_psi_per_s{$block_length} = &sigma(\@all_truncated_psi) / scalar(@all_truncated_psi);
                
                undef @$truncated_sequences;
                undef @all_truncated_psi;
                };    
            
            # Now that we have the psi values for each s for the query protein
            # (%protein_psi_per_s), we can calculate the mean of the psi
            # values for each sequence and determine how close these values
            # are. Importantly, to compensate for the spreading regression
            # lines, we must transform the data on the y-axis from "psi" to
            # "psi/s". 
    
            while (($block_length, $psi) = each %protein_psi_per_s)
                {
                $protein_psi_per_s{$block_length} = $psi / $block_length;   
                };        
            
            $protein_mean_psi = &sigma( [values %protein_psi_per_s] ) / scalar(values %protein_psi_per_s);
        
            undef %protein_psi_per_s;
        
            #############################
            # Scrambled Polypeptide Psi #
            #############################
            
            # Foreach gene, we must generate X number of scrambled polypeptides;
            # where "X" is the $max_trial_scrambled_polypeptides. We must
            # calculate psi values for all of these (stored in the @psi_values)
            # array and the original gene.
            
            foreach ( 1..$max_trial_scrambled_polypeptides )
                {        
                $scrambled_polypeptide = &scrambler($protein_sequence);
                
                # The psi calculation does not depend upon the sequence length
                # (N) or the total hydrophobicity score (M). Therefore, we can
                # truncate different amino acids off the ends of the polypeptides
                # and still compare the psi values for each block length (s). 
                
                foreach $block_length ($min_s_for_sequence..$max_s_for_sequence)
                    {        
                    $truncated_sequences = &truncate_sequence($scrambled_polypeptide, $block_length);
                    
                    foreach (@$truncated_sequences)
                        {            
                        $hydrophobicity_scores_array = &hydrophobicity_scores($_);
                        
                        $psi = &calculate_psi($hydrophobicity_scores_array, $block_length);
                        
                        push(@all_truncated_psi, $psi);
                        
                        undef $hydrophobicity_scores_array;
                        undef $psi;                
                        };
                    
                    $scrambled_psi_per_s{$block_length} = &sigma(\@all_truncated_psi) / scalar(@all_truncated_psi);
                    
                    undef @$truncated_sequences;
                    undef @all_truncated_psi;
                    };
                
                while (($block_length, $psi) = each %scrambled_psi_per_s)
                    {
                    $scrambled_psi_per_s{$block_length} = $psi / $block_length;   
                    };        
                
                $scrambled_mean_psi = &sigma( [values %scrambled_psi_per_s] ) / scalar(values %scrambled_psi_per_s);        
                
                undef %scrambled_psi_per_s;
                
                # As this algorithm generates scrambled polypeptides and calculates the
                # psi value, the best psi value (and its associated scrambled polypeptide
                # sequence) will be stored.
                
                if ( abs($scrambled_mean_psi - $protein_mean_psi) < abs($best_mean_psi - $protein_mean_psi) )
                    {
                    # This scrambled polypeptide has a psi (close to a hydrophobicity
                    # profile) that is closer to the protein than our previous "best"
                    # scrambled polypeptide.
                    
                    $best_mean_psi = $scrambled_mean_psi;
                    $best_scrambled_polypeptide = $scrambled_polypeptide;
                    };
                };
    
            # Just as a check.
           
            if ($best_mean_psi == 99999999999999)
                {
                print "The best mean psi is still its default value!\n";
                print "BEST MEAN PSI: $best_mean_psi\n";
                print "UNIQUE ID: $unique_id\n";
                print "PROTEIN SEQUENCE: $protein_sequence\n";
                die();
                };
            
            $output_line = "$unique_id\t$protein_mean_psi\t$best_mean_psi\t$best_scrambled_polypeptide";

            push(@control_sequences, $output_line);
            
            undef $output_line;    
            }; # Close for the "while" loop.

        foreach (@control_sequences)
            {            
            # The child process will pause at the following print statement until the
            # parent is ready to accept the child's data (i.e., the $reader_file_handle
            # is open). Therefore, I do not have to worry about losing the @control_sequences
            # array just because the parent is not waiting for it.

            (print $writer_file_handle "$_\n") or (die "Failed to pass message to child: $!");

            undef $_;
            };

        undef @control_sequences;
        
        close $writer_file_handle;
       
        exit(0); # exit the child process    
        }; # Close for the "else" child process conditional      
    }; # Close for the Core Foreach Loop

# The data generated from each child process will now be transfered
# to the parent process. The parent process will wait for each
# child to finish; upon finishing the next child will transfer its
# data.

foreach $core_number (keys %child_process_ids) 
    {
    $pid = $child_process_ids{$core_number};    
    $reader_file_handle = $child_file_handles{$core_number};
    
    print "CHILD_PROCESS: CORE_NUM: $core_number\tPID: $pid\tREADER_FH: $reader_file_handle\n";
    
    push( @control_sequences, (@con_seq_from_child = <$reader_file_handle>) );
    
    waitpid($pid, 0);
    
    undef @con_seq_from_child;
    undef $reader_file_handle;
    undef $output_line;
    undef $pid;
    };
    
no IO::Handle;  

# Now that we have the $best_scrambled_polypeptide for the given gene
# (based upon the $best_mean_psi score), we can now print this information
# to a flat file.

open(SCRAMBLED_SEQ_OUTPUT, ">", $output_file_name) || die("Cannot open the $output_file_name file: $!\n");

foreach (@control_sequences) { print SCRAMBLED_SEQ_OUTPUT "$_\n" };

close SCRAMBLED_SEQ_OUTPUT; 

print "\nTime of program completion\n";
main::print_time();
print "\nProgram finished";

exit;

####################################################################################
# SUBROUTINES ######################################################################
####################################################################################

# HYDROPHOBICITY_SCORES

# The following subroutine will calculate the hydrophobicity scores
# for each amino acid. Input for the subroutine is a protein
# sequence. Output is an array of scores. Note that each index
# of the array will represent the position minus 1. For example,
# index 0 will actually contain the hydrophobicity score for the
# first amino acid.

sub hydrophobicity_scores
    {
    
    # Input Variables
    
    my 	($protein_sequence) = @_;
    
    # Subroutine Variables

    my  (%hydrophobicity_score_vector,
        @protein_sequence_array,
        @hydrophobicity_scores,
        );

    # The following hydrophobicity scores were obtained from
    # http://blanco.biomol.uci.edu/hydrophobicity_scales.html
    
    %hydrophobicity_score_vector =
        (
        'A' => -1, 
        'R' => -1, 
        'N' => -1,
        'D' => -1, 
        'C' => -1,
        'Q' => -1, 
        'E' => -1, 
        'G' => -1,
        'H' => -1, 
        'I' => 1,
        'L' => 1, 
        'K' => -1, 
        'M' => 1, 
        'F' => 1, 
        'P' => -1, 
        'S' => -1, 
        'T' => -1, 
        'W' => 1, 
        'Y' => -1, 
        'V' => 1,
        'X' => 0,
        'U' => -1
        );

    @protein_sequence_array = split('', $protein_sequence);
    
    foreach (@protein_sequence_array)
        {
        push( @hydrophobicity_scores, $hydrophobicity_score_vector{$_} );
        };
    
    return \@hydrophobicity_scores;

    };

########################################################################################

# PRINT_TIME

# The following subroutine will use the "localtime" builtin perl function
# to print the date and time. Note that this is the local time on the machine.
# If the machine is a remote login server in a different time zone, then the
# time will be incorrect according to the time of the user's local computer.

sub print_time
    {
    
    # There are no input variables	
    
    # Subroutine Variables
    
    my (@master_time_variables,
        @month_array,
        
        $seconds,
        $minutes,
        $hours,
        $month,
        $day_number,
        $year,
        
        $time,
        $date,
        ); 
    
    @master_time_variables = localtime(time);

    # The @master_time_variables array contains all the information one needs
    # to print the local time. The array indexes contain information as follows:
    
    # [0] => Seconds past the minute
    # [1] => Minutes past the hour
    # [2] => Hours past midnight
    # [3] => Day of the month
    # [4] => Months past the start of the year
    # [5] => Number of years since 1900
    # [6] => Number of days since the start of the week (Sunday)
    # [7] => Number of days since the start of the year
    # [8] => Whether or not daylight savings is active
    
    @month_array = qw/January February March April May June July August September October November December/;

    $seconds = $master_time_variables[0];
    $minutes = $master_time_variables[1];
    $hours = $master_time_variables[2];
    $month = $month_array[ $master_time_variables[4] ];
    $day_number = $master_time_variables[3];
    $year = 1900 + $master_time_variables[5];
    
    $time = "$hours". ':' . "$minutes" . ':' . "$seconds";
    $date = "$year $month $day_number";
    
    print "\nDATE: $date, TIME: $time\n\n";
    
    };

########################################################################################

# CALCULATE_PSI

# The following subroutine will calculate the psi value of a sequence. This
# psi value is similar to a hydrobobicity profile score. Specifically,
# it measures how bunched hydrophobic and hydrophilic amino acid residues
# are in a sequence. Psi is calculated using an equation from the following
# reference:
#
# Irback, A. & Sandelin, E. On hydrophobicity correlations in protein chains.
# BIOPHYSICAL JOURNAL. Nov 2000; 79: 2252-2258.
#
# All variables are the same as in the journal article.
#
# Input for this subroutine is an array containing all the hydrophobicity
# values for each amino acid (with index 0 being the first amino acid).
# Importantly, whatever hydrophobicity scale is used, hydrophobic residues
# must be positive and hydrophilic residues must be negative.
#
# Requirements: This subroutine operates under the assumption that the
# length of the sequence can be divided by the block size without a
# remainder.

sub calculate_psi
    {
    
    # Input Variables
    
    my 	($hydrophobicity_scores_array,
        $lowercase_s, # $block_size
        ) = @_;
    
    # Subroutine Variables

    my  ($capital_n,
        $capital_m,
        $capital_k,
        $psi,
        $psi_sub_k,
        $sigma_sub_k, 
        
        @hydrophobicity_scores_for_block, 
        @all_psi_sub_k, 
        );
    
    # We must first check to make sure our sequence can be divided evenly
    # by our block size ($lowercase_s).
    
    unless ( scalar(@hydrophobicity_scores_for_block) % $lowercase_s == 0)
        {
        print "The sequence is not divisible by the block size\n";
        print "BLOCK SIZE: $lowercase_s\n";
        print "SEQUENCE: @$hydrophobicity_scores_array\n";
        die "ENDING PROGRAM\n";     
        };
       
    # Now we can begin the calculation of $psi.    
    
    $capital_n = scalar(@$hydrophobicity_scores_array);
    
    $capital_m = &sigma($hydrophobicity_scores_array);
    
    $capital_k = ( $capital_n**2 - $capital_m**2 ) * ( 1 - ($lowercase_s / $capital_n) ) / ( $capital_n**2 - $capital_n );
    
    if ($capital_k == 0) {$capital_k = 0.0000000001};    

    # We will now split our sequence into blocks.
    
    foreach (@$hydrophobicity_scores_array)
        {
        push(@hydrophobicity_scores_for_block, $_);
        
        unless( scalar(@hydrophobicity_scores_for_block) == $lowercase_s ) {next};
        
        $sigma_sub_k = &sigma(\@hydrophobicity_scores_for_block);
        
        undef @hydrophobicity_scores_for_block;
        
        $psi_sub_k = ( 1 / $capital_k ) * (( $sigma_sub_k - ($lowercase_s * $capital_m / $capital_n) )**2);
        
        push(@all_psi_sub_k, $psi_sub_k);
        };        
        
    # We will now calculate the mean of all the $psi_sub_k scalars.
    
    $psi = ( $lowercase_s / $capital_n ) * &sigma(\@all_psi_sub_k);
    
    undef @all_psi_sub_k;
    
    return $psi;
    
    };
    
########################################################################################

# SIGMA

# The following subroutine function as a mathematical Greek sigma (capital
# sigma). That is, it will add all the values given in the input array. 

sub sigma
    {
    
    # Input Variables
    
    my 	($values_array) = @_;
    
    # Subroutine Variables

    my  ($sum,
        );
    
    $sum = 0;
    
    foreach (@$values_array)
        {
        $sum = $sum + $_;    
        };
    
    return $sum;
    
    };   
    
########################################################################################

# SCRAMBLER

# The following subroutine scrambles the input sequence. That is,
# it will randomly order the characters that are present in the
# input string.

sub scrambler
    {
    
    # Input Variables
    
    my 	($original_sequence) = @_;
    
    # Subroutine Variables
    
    my  (@original_seq_array,
        $remaining_aa_number,
        $index,
        @scrambled_seq_array, 
        $scrambled_sequence,
        );
    
    @original_seq_array = split('', $original_sequence);
    
    # We can now begin scrambling the sequence.
    
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
        # new scrambled array
        
        push( @scrambled_seq_array, $original_seq_array[$index] );
        
        splice(@original_seq_array, $index, 1);
        };
    
    undef @original_seq_array;
    
    $scrambled_sequence = join('', @scrambled_seq_array);
    
    undef @scrambled_seq_array;

    return $scrambled_sequence;
    
    };

########################################################################################  

# TRUNCATE_SEQUENCE

# The following subroutine will truncate a sequence so that it is
# of a length that is divisible by the given $block_length. The subroutine
# will then output several sequences of the appropriate length
# (located in the @truncated_sequences array). Each truncated sequence is
# derived by a sliding window based upon the length that the truncated
# sequence needs to be. Therefore, the number of entries in the output
# array will depend upon how many times the sliding window has to "slide".

# For example, if my sequence is of length 22 and my block length is 5.
# I can generate four blocks, each with a length of 5, out of the sequence.
# There will be two characters remaining at the end of the sequence that
# are not placed into a block. The truncated sequence will be of length 20
# (4 * 5 = 20). The "window" of length 20 will begin at position 1 (and
# extend to position 20); this will be the first truncated sequence. The
# "window" will then slide +1 so that it will now range from 2 to 21; this
# sequence will now be the second truncated sequence. The final truncated
# sequence will extend from position 3 to position 22. These three truncated
# sequences will be returned in the form of the @truncated_sequences array.

sub truncate_sequence
    {
    
    # Input Variables
    
    my 	($original_sequence, $block_length) = @_;
    
    # Subroutine Variables
    
    my  ($remainder,
        $truncated_sequence,
        @truncated_sequences,
        );
    
    $remainder = length($original_sequence) % $block_length;

    foreach (0..$remainder)
        {        
        $truncated_sequence = substr( $original_sequence, $_, length($original_sequence) - $remainder );
        
        push(@truncated_sequences, $truncated_sequence);
        
        undef $truncated_sequence;
        };

    return \@truncated_sequences;
    
    };
    
########################################################################################

