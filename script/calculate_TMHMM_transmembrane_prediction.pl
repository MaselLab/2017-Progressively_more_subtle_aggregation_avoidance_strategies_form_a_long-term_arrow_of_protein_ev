#!/usr/local/bin/perl

##################################################################
=head1 # CALCULATE TMHMM TRANSMEMBRANE PREDICTION ################
##################################################################

=head2 AUTHOR: Scott Foy, Ph.D.

This program will calculate the transmembrane information 
of polypeptide sequences located inside a MySQL database. This
program is similar to a shell script in that it performs several
small, independent tasks. The first task extracts data from the
MySQL database one row at a time. For each row, the TMHMM program
is executed on the row's polypeptide sequence in order to generate
the transmembrane information of that sequence. Note that, although
TMHMM will accept multiple input sequences, we extract only individual
rows from MySQL to keep to program simple. Finally, the generated
transmembrane information will be placed back into the database.

The transmembrane information that is currently returned to the
database is a binary value specifying whether the given polypeptide
sequence is a transmembrane protein (i.e., a value of "yes" or "no").
However, this program can be modified to return to the database
any information that can be extracted from the output of the
TMHMM program.

Finally, this script should be executed on the computer running
the TMHMM program; the MySQL database does not need to
be located on this same computer.

=head2 REQUIREMENTS:

The "TMHMM" program. 

Additionally, this program requires the DBD::mysql module (and the
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
thoughout the entire CALCULATE TMHMM TRANSMEMBRANE PREDICTION program.

=cut

use warnings;
use strict;

print "Time of program execution:\n";
main::print_time();
print "\n";

my  ($mysql_extraction_query,
    $mysql_insertion_statement,
    
    $database_host_ip_address,
    $mysql_account_username,
    $mysql_account_password,
    $mysql_database_name,
    );

# Default Input MySQL Commands:
# ____________________________

$mysql_extraction_query = 'SELECT PrimaryKey, ProteinSequence FROM TableName'; 
$mysql_insertion_statement = 'UPDATE TableName SET TMHMMTransmembraneBinary = ? WHERE PrimaryKey = ?';

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. If the user wishes to calculate the
transmembrane status of only certain protein sequences, the user may add a
"WHERE" statement to the MySQL command. 

$mysql_insertion_statement = This is the MySQL command that inserts the
various transmembrane scores back into the MySQL database. For the most part, this
command should be input exactly as it would be input into a MySQL command
line interface. However, the command requires 2 placeholders (question marks).
The first placeholder is a binary indicating whether the protein is a transmembrane
protein. The second placeholder is the unique ID so that the program knows where to
insert the transmembrane information.

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
generated transmembrane information will be placed (i.e., the one containing the
output information from TMHMM). 

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
    $transmembrane_binary,
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
    # $mysql_extraction_query.In this case, the zeroth
    # element is the PrimaryKey and the first element is the
    # ProteinSequence.
    
    $unique_id = $row_array->[0];    
    $protein_sequence = $row_array->[1];

    undef $row_array;   

    $transmembrane_binary = &transmembrane_score($unique_id, $protein_sequence);

    # We must now insert the $transmembrane_binary score
    # into the database for this specific row.

    $database_insertion_object->execute($transmembrane_binary, $unique_id);

    undef $transmembrane_binary;
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

# TRANSMEMBRANE_SCORE

# The following subroutine will execute and utilize the output from the TMHMM
# program to calculate whether a protein is a transmembrane protein
# (i.e., the $transmembrane_binary). Any protein with more than 18
# transmembrane helical amino acids is considered to be a transmembrane
# protein (this threshold comes from the TMHMM journal article). The output
# $transmembrane_binary can have two possible values: "yes" or "no", to
# indicate if the protein is a transmembrane protein.

sub transmembrane_score
    {
    
    # Input Variables
    
    my 	($unique_id, $protein_sequence) = @_;
    
    # Subroutine Variables
    
    my  (@tmhmm_output_lines,
        $tmhmm_execution_command,
        $number_of_tmhs,
        $number_of_aa_in_tmhs,
        $transmembrane_binary,
        );
    
    # TMHMM accepts a FASTA file via STDIN. Instead of generating a FASTA file,
    # we will just feed a FASTA format directly into the TMHMM executable.
    
    # In BASH, "printf" is a more consistent and relable version of "echo".

    $tmhmm_execution_command = 'printf "' . "\>$unique_id\n$protein_sequence" . '" | tmhmm';

    @tmhmm_output_lines = `$tmhmm_execution_command` or die "Cannot input the temporary FASTA format directly into TMHMM";

    chomp @tmhmm_output_lines;

    # The @tmhmm_output_lines contains the long output for the TMHMM program.
    # We will now extract numerous pieces of information, even if we do
    # not need all of them. Note that "tmhs" = transmembrane helices.
    
    $tmhmm_output_lines[1] =~ m/TMHs:\s+(\d*\.*\d*)$/;
    $number_of_tmhs = $1;

    $tmhmm_output_lines[2] =~ m/TMHs:\s+(\d*\.*\d*)$/;
    $number_of_aa_in_tmhs = $1;
    
    undef @tmhmm_output_lines;
    
    # We will now determine the $transmembrane_binary based upon the number of
    # amino acids in all the helices ($number_of_aa_in_tmhs).
    
    undef $transmembrane_binary;
    
    if ( $number_of_aa_in_tmhs > 18 ) { $transmembrane_binary = 'yes' }
    
    elsif ( $number_of_aa_in_tmhs <= 18 ) { $transmembrane_binary = 'no' }
    
    else { die( "The output for the number of amino acids in the helices is not a number. Number: $number_of_aa_in_tmhs\n" ) };
    
    print "UID: $unique_id,\tNumber_of_TM_AAs: $number_of_aa_in_tmhs,\tBinary: $transmembrane_binary\n";
    
    undef $number_of_tmhs;
    undef $number_of_aa_in_tmhs;
    
    return $transmembrane_binary;
    
    };

########################################################################################
