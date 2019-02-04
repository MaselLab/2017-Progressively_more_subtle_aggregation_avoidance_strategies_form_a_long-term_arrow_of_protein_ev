#!/usr/local/bin/perl

##################################################################
=head1 # CALCULATE AGGREGATION ###################################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program identifies oligopeptides with a high aggregation
propensity from polypeptide sequences located inside a MySQL
database. This program is similar to a shell script in that it
performs several small, independent tasks. The first task extracts
data from the MySQL database one row at a time. For each row, the
TANGO or the WALTZ program is executed on the row's polypeptide sequence
in order to generate the aggregation propensity for each amino acid
of that sequence. Finally, the generated aggregation-prone
oligopeptide sequences will be placed back into the database. Note
that this program extracts only individual rows from MySQL because
this makes it easier for the algorithm to return the detected
aggregation-prone oligopeptides to their respective entity in the
database.

Finally, this script should be executed on the computer running
the bioinformatics program; the MySQL database does not need to
be located on this same computer.

=head2 REQUIREMENTS:

The "WALTZ" and/or "TANGO" programs.

Additionally, this program requires the DBD::mysql module (and the
library [DBD::mysql directory]). Although DBD and DBI are both
part of the core modules, the DBD:mysql driver is not and must be
downloaded from CPAN.

Finally, no protein unique ID may contain a space or tab.

#####################################################################
=head1 # ARGUMENTS ##################################################
#####################################################################

The hardcoded values below are the default values for each variable.
Additionally, below there are descriptions of each variable.

=head2 ARGUMENT VALUES

The following input variables are package variables and can be
utilized thoughout the entire CALCULATE_AGGREGATION program.

=cut

use warnings;
use strict;

package Input_Variables;   

use vars qw
    ($aggregation_program
    $execute_tango_program_command
    $execute_waltz_program_command
    $tango_agg_threshold
    $waltz_agg_threshold

    $mysql_extraction_query
    $mysql_insertion_statement
    
    $database_host_ip_address
    $mysql_account_username
    $mysql_account_password
    $mysql_database_name
    );

# Default Aggregation Program:
# ____________________________

$aggregation_program = 'tango';
$execute_tango_program_command = '/home/$USER/bin/tango/tango_x86_64';
$execute_waltz_program_command = '/home/$USER/bin/waltz.pl';
$tango_agg_threshold = 5.0;
$waltz_agg_threshold = 92.0; # These must be 77 (high sensitivity), 92 (good overall performance; default), or 98 (high specificity)    

# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT PrimaryKey, ProteinSequence FROM TableName'; 
$mysql_insertion_statement = 'UPDATE TableName SET ProteinSequence = ?, OrigDataAminoAcidAggPropScores = ?, NumAggProneRegions = ?, AminoAcidNumInAPRs = ? WHERE PrimaryKey = ?';

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$aggregation_program = After the command, the name of the aggregation
program to be used must be stated. The current options for the
aggregation programs are:

    tango (default)
    waltz

$execute_[tango|waltz]_program_command = This is the Linux path and command for
TANGO or WALTZ. The entire command is not necessary because it is too complex and
is generated inside a subroutine. Additional variables can be added as
additional aggregation programs are added.

$[tango|waltz]_agg_threshold = This is the amino acid aggregation score
that must be surpassed in order to declare that amino acid as "aggregation-prone".

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. If the user wishes to calculate the
aggregation propensity of only certain protein sequences, the user may add a
"WHERE" statement to the MySQL command. 

$mysql_insertion_statement = This is the MySQL command that inserts the various
aggregation propensity scores back into the MySQL database. For the most part,
this command should be input exactly as it would be input into a MySQL command
line interface. However, the command requires 5 placeholders (question marks),
The first placeholder re-inserts the protein sequence back into the database.
This is necessary because this program serves a secondary function to "clean up"
the protein sequence by removing all non-letter characters from the protein
sequence and then returning the sequence with these removed. The second placeholder
is for the string of aggregation values per amino acid, the third is the number
of aggregation prone regions (APRs), and the fourth is for the number of amino
acids present in these APRs. The fifth placeholder is the unique ID so that
the program knows where to insert the aggregation information.

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
generated aggregation information will be placed (i.e., the one containing the
output information from the bioinformatics program). 

=cut

package main;

print "Time of program execution:\n";
main::print_time();
print "\n";

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($database_connection_object,
    $database_extraction_object,
    $database_insertion_object,
    $row_array,
    
    $unique_id,
    $protein_sequence_with_x,
    $protein_sequence,
    
    $aggregation_propensity_scores,
    $aggregation_threshold,
    $aggregation_propensity_string,
    
    $num_of_agg_prone_regions,
    $num_of_aa_in_agg_prone_regions,    
    );

# We must first access the MySQL database.

use DBI;

print 'Connecting to host database = ' . "$Input_Variables::database_host_ip_address\n";

$database_connection_object = DBI->connect("dbi:mysql:$Input_Variables::mysql_database_name;$Input_Variables::database_host_ip_address",
                                           $Input_Variables::mysql_account_username, $Input_Variables::mysql_account_password);

print "Database connected\n";

# We now pre-prepare our MySQL commands that will be used to both
# extract and insert data.

$database_extraction_object = $database_connection_object->prepare($Input_Variables::mysql_extraction_query);
$database_insertion_object = $database_connection_object->prepare($Input_Variables::mysql_insertion_statement);

$database_extraction_object->execute();

print "Database extraction command executed\n";

while ($row_array = $database_extraction_object->fetchrow_arrayref)
    {
    chomp @$row_array;
    
    # The $row_array is a reference to an array, each element of which
    # is a value from a row of the database. Specifically, each element
    # contains the value returned by the "SELECT" statement in the
    # $mysql_extraction_query.In this case, the zeroth
    # element is the PrimaryKey and the first element is the
    # ProteinSequence.
    
    $unique_id = $row_array->[0];    
    $protein_sequence_with_x = $row_array->[1];    

    undef $row_array;
  
    # We must first check the integrity of the protein sequence.
    # We will later return this sequence back to the database
    # in order to facilitate the use of the checked sequence.
    
    $protein_sequence_with_x =~ s/[^a-z]//ig; 
    
    # Unfortunately, the aggregation programs are unable to read an
    # unknown amino acid ("X") or selenocystein ("U"). Therefore, we must
    # remove all "X"s and "U"s from the protein sequence before executing
    # these programs. These unknown amino acids (or selenocysteins) will be
    # inserted into the program's output. Additionally, in order to conserve
    # the integrity of the sequence, we will check to see if any
    # sequence has more than three consecutive unknown amino acids (or
    # selenocysteins). 
   
    if ($protein_sequence_with_x =~ m/(X|U)(X|U)(X|U)/ig)
        {
        print "\nWARNING: The following protein sequence contains three or more\n";
        print "consecutive unknown amino acids or selenocysteins:\n";
        print "Unique ID: $unique_id\n";
        print "Sequence: $protein_sequence_with_x\n";
        
        next;
        };
  
    $protein_sequence = $protein_sequence_with_x;
    
    $protein_sequence =~ s/(X|U)//gi;

    # This subroutine will input the MySQL data into any of our
    # optional aggregation programs (i.e., TANGO, WALTZ, etc.), 
    # and then store the output aggregation scores for each
    # amino acid in the @$aggregation_propensity_scores array. 

    if ($Input_Variables::aggregation_program =~ m/^waltz$/i)
        {
        $aggregation_threshold = $Input_Variables::waltz_agg_threshold;
        
        $aggregation_propensity_scores = &execute_waltz($unique_id, $protein_sequence);     
        }
        
    else # TANGO will execute by default.
        {
        $aggregation_threshold = $Input_Variables::tango_agg_threshold;
        
        $aggregation_propensity_scores = &execute_tango($unique_id, $protein_sequence);    
        };
   
    # We will now insert the X's and U's into the
    # @$aggregation_propensity_scores array.
    
    $aggregation_propensity_scores = &insert_x($aggregation_propensity_scores, $protein_sequence_with_x);
    
    # Now that we know all the amino acid scores (or have
    # abbreviations in place for those that we do not know),
    # we can now calculate any output scores we wish.
    
    ($num_of_agg_prone_regions, $num_of_aa_in_agg_prone_regions) = &calc_num_of_aprs_and_aa_in_aprs($aggregation_propensity_scores, $aggregation_threshold);   
    
    # We will now conglomerate the output aggregation amino acid scores
    # into a string of comma-separated scores. We will then
    # output all our information to the database.
    
    $aggregation_propensity_string = join(',', @$aggregation_propensity_scores);
    
    undef @$aggregation_propensity_scores;
    
    # We must now insert the $agg_propensity_string, the $num_of_agg_prone_regions, and
    # the $num_of_aa_in_agg_prone_regions scalars into the database for this specific row.

    $database_insertion_object->execute
        (
        $protein_sequence_with_x,
        $aggregation_propensity_string,
        $num_of_agg_prone_regions,
        $num_of_aa_in_agg_prone_regions,
        $unique_id,
        );
    
    undef $aggregation_propensity_string;
    undef $num_of_agg_prone_regions;
    undef $num_of_aa_in_agg_prone_regions;
    undef $unique_id;
    undef $protein_sequence;
    undef $protein_sequence_with_x;
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

# EXECUTE_WALTZ

# The following subroutine will generate a FASTA file containing
# the protein sequence to be used as input for WALTZ. The output
# from WALTZ will be retrieved from STDOUT. The output will be
# saved to the @waltz_output_lines array. Finally, the WALTZ output
# will be curated into the @aggregation_propensity_scores array, which
# lists the aggregation scores for each individual amino acid.

sub execute_waltz
    {
    
    # Input Variables
    
    my 	($unique_id,
        $protein_sequence,
        ) = @_;
    
    # Subroutine Variables
    
    my  ($ph,
        $aggregation_threshold,

        $waltz_input_file,
        $execute_waltz_command,
        @waltz_output_lines,
        @aggregation_propensity_scores,
        );

    # Below are the WALTZ input parameters.

    $ph = 7.0;
    $aggregation_threshold = $Input_Variables::waltz_agg_threshold; 
    
    # WALTZ can only read a sequence from a FASTA format. Therefore,
    # we must generate a FASTA file for the current protein sequence.    
    
    $waltz_input_file = './temp_fasta_input_for_waltz.fasta';
  
    open(FASTA_INPUT_FOR_WALTZ, ">", $waltz_input_file) || die("Cannot open the $waltz_input_file file: $!\n");

    print FASTA_INPUT_FOR_WALTZ '>' . "$unique_id\n$protein_sequence\n\n";
    
    close FASTA_INPUT_FOR_WALTZ;
  
    # We must now generate the WALTZ command.
    
    $execute_waltz_command = "$Input_Variables::execute_waltz_program_command $waltz_input_file $aggregation_threshold tango $ph";

    # We will now execute the WALTZ command.
   
    @waltz_output_lines = `$execute_waltz_command` or die "Cannot input the temporary FASTA file into WALTZ";

    !system "rm $waltz_input_file" or die ("Cannot delete $waltz_input_file: $!\n");

    chomp @waltz_output_lines;
   
    # We will now remove the individual amino acid scores from the WALTZ output
    # and place them into an array.   
    
    # The @waltz_output_lines array should only contain one line.
    
    @aggregation_propensity_scores = split("\t", $waltz_output_lines[0]);
    
    undef @waltz_output_lines;
    
    # The first element of the array is the unique ID of the protein.
    
    shift(@aggregation_propensity_scores);
    
    return (\@aggregation_propensity_scores);
    
    };
    
########################################################################################  
    
# INSERT_X

# The following subroutine inserts the "X" and "U"
# amino acids into the aggregation scores calculated by WALTZ.

sub insert_x
    {
    
    # Input Variables
    
    my 	($aggregation_propensity_scores, $protein_sequence_with_x) = @_;
    
    # Subroutine Variables
    
    my  (@scores_with_x,
        );    

    foreach ( split('', $protein_sequence_with_x) )
        {
        if ( $_ =~ m/(X|U)/i )
            {
            push(@scores_with_x, $_);    
            }
        
        else
            {
            push( @scores_with_x, shift(@$aggregation_propensity_scores) );
            };
        };

    undef @$aggregation_propensity_scores;
    
    return \@scores_with_x;
    
    };
    
########################################################################################

# CALC_NUM_OF_APRS_AND_AA_IN_APRS

# The following subroutine will accept an array containing
# the amino acid aggregation scores and determine the
# number of aggregation prone regions (APRs) present as well as
# the number of amino acids within these APRs. An APR is
# defined as an oligopeptide that is at least five amino
# acids in length (pentapeptide) because this is the minimum length
# required for an aggregating beta-strand. This subroutine outputs
# the number of APRs and the total number of amino acids within
# the APRs for a single protein sequence.

sub calc_num_of_aprs_and_aa_in_aprs
    {
    
    # Input Variables
    
    my 	($aggregation_propensity_scores, $aggregation_threshold) = @_;
    
    # Subroutine Variables
    
    my  ($num_of_aa_in_agg_prone_regions,
        $num_of_agg_prone_regions,

        @single_apr,        
        );
   
    $num_of_aa_in_agg_prone_regions = 0;
    $num_of_agg_prone_regions = 0;
    undef @single_apr;
    
    # We must add a really small number to the end of the
    # @$aggregation_propensity_scores array so that the algorithm
    # will count the final APR.
    
    push(@$aggregation_propensity_scores, 0.0);
    
    foreach (@$aggregation_propensity_scores)
        {
        if ( ($_ =~ m/^\d/) and ($_ >= $aggregation_threshold) )
            {
            push(@single_apr, $_);    
            }
        
        elsif ($_ =~ m/^[a-z]$/i)
            {
            # A single "X" or "U" is considered to be a null value.
            # That is, they do not count for or against the number
            # of aggregation-prone regions unless it is located within
            # an APR. 
            
            push(@single_apr, $_);    
            }
        
        elsif ( ($_ =~ m/^\d/) and ($_ < $aggregation_threshold) )
            {
            # Regardless of the status of the @single_apr array, this marks
            # the end of it.

            unless ( @single_apr ) {next};
            
            # If the @single_apr array is defined, we must remove any unnumbered
            # scores from the front or back of the array. We do not want these
            # to contribute towards the length of the APR. At most we can have
            # two unknown amino acids on the ends (sequences containing three
            # or more consecutive amino acids were removed early in the proram).
            # The hidden "defined" commands are placed repeatedly because the status
            # (i.e., "defined" or "undefined") of the @single_apr array can chnage
            # after any of the commands.
            
            if ( @single_apr and ($single_apr[0] =~ m/[a-z]/i) ) { shift(@single_apr) };
            if ( @single_apr and ($single_apr[0] =~ m/[a-z]/i) ) { shift(@single_apr) };

            if ( @single_apr and ($single_apr[-1] =~ m/[a-z]/i) ) { pop(@single_apr) };
            if ( @single_apr and ($single_apr[-1] =~ m/[a-z]/i) ) { pop(@single_apr) };
            
            if ( @single_apr and (scalar(@single_apr) >= 5) )
                {
                $num_of_aa_in_agg_prone_regions = $num_of_aa_in_agg_prone_regions + scalar(@single_apr);    
                
                $num_of_agg_prone_regions++;
                };
                
            undef @single_apr;
            };
        };

    pop(@$aggregation_propensity_scores);
        
    return ($num_of_agg_prone_regions, $num_of_aa_in_agg_prone_regions);
    
    };    
    
########################################################################################

# EXECUTE_TANGO

# TANGO does not accept a standard FASTA file as input. Instead,
# the information is input into TANGO on the command line (as STDIN).
# Furthermore, TANGO does not generate the output to STDOUT; instead,
# TANGO produces a file that contains the output. Because the
# I/O for TANGO is a bit complex, we must utilize the following
# subroutine. 

# The following subroutine will input the protein sequence to
# TANGO on the command line (STDIN) and retrieve the output from
# the generated file. The output will be saved to the
# @tango_output_lines array. Finally, the TANGO output will be
# curated into the @aggregation_propensity_scores array, which
# lists the aggregation scores for each individual amino acid.

sub execute_tango
    {
    
    # Input Variables
    
    my 	($unique_id,
        $protein_sequence,
        ) = @_;
    
    # Subroutine Variables
    
    my  ($ph,
        $temp,
        $ionic_strength,
        $concentration,
        
        $input_command_line,
        $output_file_name,
        @tango_output_lines,
        
        @aggregation_propensity_scores,
        );
    
    # Below are the TANGO input parameters.

    $ph = 7.0;
    $temp = 298.15; # Kelvin
    $ionic_strength = 0.02; # M
    $concentration = 1.0; # M
    
    # We must now generate the TANGO command. Notice that the TANGO command
    # redirects the TANGO STDOUT to a file called "redirection.tmp". This is
    # because TANGO does not allow the user to ignore (not print) a message
    # to STDOUT. Over the course of the entire database, TANGO will print
    # to STDOUT several thousand times. This file is a way to prevent TANGO
    # from printing to STDOUT.
    
    $input_command_line = "$Input_Variables::execute_tango_program_command $unique_id ct=\"N\" nt=\"N\" ph=\"$ph\" te=\"$temp\" io=\"$ionic_strength\" seq=\"$protein_sequence\" \> next_redirection.tmp";
    
    $output_file_name = "$unique_id" . '.txt';
    
    # We will now execute Tango.
    
    !system "$input_command_line" or die ("Cannot execute Tango for $unique_id: $!\n");
    
    open(INPUT_FROM_TANGO, "<", $output_file_name) || die("Cannot open the $output_file_name file: $!\n");

    @tango_output_lines = <INPUT_FROM_TANGO>;
    
    chomp @tango_output_lines;
    
    close INPUT_FROM_TANGO;
    
    !system "rm $output_file_name" or die ("Cannot delete $output_file_name: $!\n");
    !system 'rm next_redirection.tmp' or die ("Cannot delete next_redirection.tmp: $!\n"); 
    
    undef $input_command_line;
    undef $output_file_name;
    
    # We will now remove the individual amino acid scores from the TANGO output
    # and place them into an array.  

    foreach (@tango_output_lines)
        {
        if ( $_ =~ m/(\d+\.\d+)\t+\d+\.\d+$/g )
            {
            # In the above regular expresion, "$1" is the aggregation propensity.

            push(@aggregation_propensity_scores, $1);           
            };
        };
    
    undef @tango_output_lines;
        
    return (\@aggregation_propensity_scores);
    
    };
    
########################################################################################  
