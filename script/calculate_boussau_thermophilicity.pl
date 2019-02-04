#!/home/sfoy/perl-5.14.2/perl -w

use warnings;
use strict;

#################################################################################
=head1 # CALCULATE BOUSSAU THERMOPHILICITY ######################################
#################################################################################

This program will calculate a thermophilicity score for a given
protein sequence. The algorithm to calculate the thermophilicity
score is in Boussau (2008). This algorithm will calculate the thermophilicity
score (M) and then divide it by length (N). 

Each polypeptide sequence will be extracted from a MySQL database. The
script will return the Boussau thermophilicity scores to the database. 

=head2 REQUIREMENTS:

This program requires the DBD::mysql module (and the
library [DBD::mysql directory]). Although DBD and DBI are both
part of the core modules, the DBD:mysql driver is not and must be
downloaded from CPAN.

#####################################################################
=head1 # ARGUMENTS ##################################################
#####################################################################

The hardcoded values below are the default values for each variable.
Additionally, below these are descriptions of each variable.

=head2 ARGUMENT VALUES

The following input variables are lexical and can be utilized
thoughout the entire CALCULATE BOUSSAU THERMOPHILICITY program.

=cut

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
$mysql_insertion_statement = 'UPDATE TableName SET BoussauThermophilicityPercent = ? WHERE PrimaryKey = ?';

# Default Database Information:
# ____________________________

$database_host_ip_address = '127.0.0.1'; # "localhost"
$mysql_account_username = ''; 
$mysql_account_password = ''; 
$mysql_database_name = ''; 

=head2 ARGUMENT DESCRIPTIONS

$mysql_extraction_query = This is the MySQL command that extracts data from
the MySQL database. This command should be input exactly as it would be
input into a MySQL command line interface. 

$mysql_insertion_statement = This is the MySQL command that inserts the ISD
score back into the MySQL database. For the most part, this command should
be input exactly as it would be input into a MySQL command line interface.
However, the command requires 2 placeholders (question marks). The first
placeholder is for the thermophilicity value, the second placeholder is the
unique ID so that the program knows where to insert the thermophilicity information.

$database_host_ip_address = This is the IP address of the computer hosting
the MySQL database. If the MySQL database is located on the local computer
that is executing this program, the IP address should be set as
"localhost". 

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
generated thermophilicity values will be placed (i.e., the one containing the
output information from the bioinformatics program). 

=cut

#####################################################################
# PRIMARY PROGRAM ###################################################
#####################################################################

my  ($database_connection_object,
    $database_extraction_object,
    $database_insertion_object,
    $row_array,
    
    $sequence_uid,
    $protein_sequence,
    $thermophilicity_score,
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
    
    $sequence_uid = $row_array->[0];
    $protein_sequence = $row_array->[1];
        
    undef $row_array;
    
    $thermophilicity_score = &calculate_thermophilicity($protein_sequence);
    
    # We must now insert the $protein_sequence into the database for this specific row.
    
    $database_insertion_object->execute($thermophilicity_score, $sequence_uid);
    
    undef $protein_sequence;
    undef $sequence_uid;
    };

$database_connection_object->disconnect;

no DBI;  

print "\nProgram finished";

exit;

####################################################################################
# SUBROUTINES ######################################################################
####################################################################################

# CALCULATE_THERMOPHILICITY

# The following subroutine will calculate the thermophilicity scores
# for each amino acid. Input for the subroutine is a protein
# sequence. Output is a thermophilicity score divided by the number
# of amino acids in the protein. 

sub calculate_thermophilicity
    {
    
    # Input Variables
    
    my 	
        ($protein_sequence) = @_;
    
    # Subroutine Variables

    my  (%thermophilicity_score_vector,
        @protein_sequence_array,
        $thermophilicity_score,
        );

    %thermophilicity_score_vector =
        (
        'A' =>  0, 
        'R' =>  1, 
        'N' =>  0,
        'D' =>  0, 
        'C' =>  0,
        'Q' =>  0, 
        'E' =>  1, 
        'G' =>  0,
        'H' =>  0, 
        'I' => 1,
        'L' => 1, 
        'K' =>  0, 
        'M' => 0, 
        'F' => 0, 
        'P' => 0, 
        'S' =>  0, 
        'T' =>  0, 
        'W' => 1, 
        'Y' =>  1, 
        'V' => 1,
        'X' => 0,
        'U' => 0
        );

    $thermophilicity_score = 0;
        
    @protein_sequence_array = split('', $protein_sequence);
    
    foreach (@protein_sequence_array)
        {
        $thermophilicity_score = $thermophilicity_score + $thermophilicity_score_vector{$_};
        };
    
    undef @protein_sequence_array;
    
    $thermophilicity_score = $thermophilicity_score / length($protein_sequence);
    
    return $thermophilicity_score;

    };

########################################################################################
