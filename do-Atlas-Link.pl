#!/usr/bin/perl
$|++;

use strict;
use FindBin;
use lib "$FindBin::Bin/lib";
use Atlas::vEXAMer;
use Getopt::Long;


=head1 NAME
            
do-A-Link.pl
            
=head1 OPTIONS
    
-f A configure file.
    
-l library file .
    
-t contig file

-a (option) A agp file, if you want to do de-novo scaffolding or gap-filling in the first step, you need define this file, otherwise, no

-o output folder

-h help
=head1 SYNOPSIS

perl do-A-Link.pl (-options)

=cut


my $configure_file; my $library_file; my $contig_file; my $agp_file; my $output;my $help;
GetOptions (
            "configure_file|f=s"=> \$configure_file,
                "library_file|l=s"=> \$library_file,
                "contig_file|t=s"=> \$contig_file,
            "output_folder|o=s"=> \$output,
            "agp_file|a=s"=> \$agp_file,
            "help|h"=>\$help
            );
if ($help  or !$configure_file or !$library_file or !$contig_file  ){
    exec ("perldoc $0");
    exit;
}

if ($output){
    unless (-d $output){
            mkdir $output;
    }
}
else{
    $output='.';
}

my $vEXAMer_obj;
if ($agp_file){
        print STDERR "a AGP file was given\n";
        $vEXAMer_obj=new Atlas::vEXAMer (-configure_file=>$configure_file,
                     -library_file=>$library_file,
                     -contig_file=>$contig_file,
                     -agp_file=>$agp_file
                    );
}
else{
        $vEXAMer_obj=new Atlas::vEXAMer (-configure_file=>$configure_file,
                     -library_file=>$library_file,
                     -contig_file=>$contig_file
                    );

}
$vEXAMer_obj->do_scaffolding_step_by_step();                                                
