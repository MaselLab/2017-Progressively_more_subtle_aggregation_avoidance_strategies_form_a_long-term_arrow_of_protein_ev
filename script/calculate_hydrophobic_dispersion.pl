#!/usr/local/bin/perl

#################################################################################
=head1 # CALCULATE HYDROPHOBIC DISPERSION #######################################
#################################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will calculate the hydrophobicity dispersion (as calculated
by Irback [2000]) of all polypeptide sequences in a MySQL table. Essentially,
the hydrophobicity dispersion score is the variance divided by the mean:

    psi = VAR(hydrophobicity_per_block) / MEAN(hydrophobicity_per_block)

The hydrophobicity dispersion is dependent upon the size of the
block that is being used. This program can accomodate a dispersion block
range. It will average the psis of each block length within the range
in order to obtain the mean psi for a given polypeptide. 

The hydrophobic dispersion scores will be placed back into the database
after they are calculated. 

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
        
    $mysql_extraction_query,
    $mysql_insertion_statement,

    $database_host_ip_address,
    $mysql_account_username,
    $mysql_account_password,
    $mysql_database_name,
    );

# Default Equation Values:
# ________________________

$min_s = 6;
$max_s = 6;

# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT PrimaryKey, ProteinSequence FROM TableName'; 
$mysql_insertion_statement = 'UPDATE TableName SET HydrophobicDispersionScore = ? WHERE PrimaryKey = ?';

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

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. If the user wishes to calculate the
dispersion of only certain protein sequences, the user may add a
"WHERE" statement to the MySQL command.  

$mysql_insertion_statement = This is the MySQL command that inserts the
various hydrophobic dispersion scores back into the MySQL database. For the
most part, this command should be input exactly as it would be input into a
MySQL command line interface. However, the command requires 2 placeholders
(question marks). The first placeholder is for the hydrophobic score and the second
placeholder is the unique ID so that the program knows where to insert the
dispersion information.

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
information to be extracted. This is also the database into which the newly
generated values will be placed (i.e., the one containing the output information
from the bioinformatics program). 

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($database_connection_object,
    $database_extraction_object,
    $database_insertion_object,
    $row_array,
    
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
    
    $protein_mean_psi,
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
$database_insertion_object = $database_connection_object->prepare($mysql_insertion_statement);

$database_extraction_object->execute();

print "Database extraction command executed\n";

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

    # The maximum s that can be utilized for a given protein
    # ($max_s_for_sequence) must be shorter than half the length of the
    # protein; therefore, we will temporarily shorten $max_s if necessary
    # to accommodate for shorter proteins. We must have the $max_s variable
    # so that $max_s_for_sequence can be reset for each protein.
    
    $max_s_for_sequence = $max_s;
    $min_s_for_sequence = $min_s;
    
    if ( length($protein_sequence) / 2 < $max_s_for_sequence ) { $max_s_for_sequence = length($protein_sequence) / 2 };
    
    if ($max_s_for_sequence < $min_s_for_sequence) {$min_s_for_sequence = $max_s_for_sequence};
    
    # We must now calculate the psi for the input polypeptide. If we
    # calculate a "block" length of only one amino acid, the psi for
    # the polypeptide sequence will be the same because this psi score
    # is dependent solely on amino acid content (with a block length of 1).

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
    # values for the polypeptide sequence. Importantly, to compensate
    # for the spreading regression lines, we must transform the data
    # on the y-axis from "psi" to "psi/s". 
    
    while (($block_length, $psi) = each %protein_psi_per_s)
        {
        $protein_psi_per_s{$block_length} = $psi / $block_length;   
        };        
    
    $protein_mean_psi = &sigma( [values %protein_psi_per_s] ) / scalar(values %protein_psi_per_s);

    undef %protein_psi_per_s;

    # Now that we have the psi value for the query polypeptide sequence, we can
    # insert this information into the database. 
    
    $database_insertion_object->execute($protein_mean_psi, $unique_id);
    
    undef $max_s_for_sequence;
    undef $min_s_for_sequence;
    undef $protein_mean_psi;
    undef $unique_id;
    undef $protein_sequence;
    };

$database_connection_object->disconnect;

no DBI;  

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

