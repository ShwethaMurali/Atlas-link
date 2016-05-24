package Atlas::vEXAMer; 

use strict;
use FindBin;
use lib "$FindBin::Bin/lib";
#use lib "$FindBin::Bin/Lib/Graph_94/lib";
#use lib "$FindBin::Bin/Lib/XML-DOM-1-44/lib";
#use lib "$FindBin::Bin/Lib/XML-RegExp-0-03/lib";
#use lib "$FindBin::Bin/Lib/XML-Parser-2-40/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi/";
#use lib "$FindBin::Bin/Lib/Statistics-Descriptive-3-0200/lib/";
use Atlas::EXAMer;
use XML::DOM;
use Graph::Directed;
use Atlas::AGP;
#use lib "$FindBin::Bin/Lib/Statistics-Descriptive-Weighted-0.4/lib/";
use Statistics::Descriptive::Weighted;
our @ISA = qw(Exporter);
use 5.008008;

require Exporter;


# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration       use vEXAMer ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(
        
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
        
);

our $VERSION = '0.01';


# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.




sub new{
    my ($caller,@args) = @_; # some parameter to initiate a genome
    my $class = ref($caller) || $caller;
    my $self = {};
    bless $self, $class;
    $self->_initialize(@args );
    return $self;
}

sub _initialize{
    my ($self,@args)=@_;
    @{$self->{possible_configure_keys}}=qw{remove_redundant_reads dominant_ratio min_links deviate_factor treat_as_no_link_if_less_than min_unit_to_form_cluster distance force_elongate_min_link overrule_strength_for_reinitiate  treat_minus_gap_size_as remove_node_if_shorter_than force_inclusion_if_unit_size_more_than remove_node_bridges_more_than excessive_bridges_dev_factor  excessive_mate_pairs_for_a_edge_dev_factor favor_density min_singleton_contig_length};
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase key
    
    if  ($param{'-agp_file'}){
        $self->{'agp_file'}=$param{'-agp_file'};
    } 
   
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    $year += 1900; $mon+=1; $hour++;
    $self->{'local_time'}=join ('-', ($year, $mon, $mday,$hour,$min,$sec)) ;
    my $log_dir;
    unless ($param{'-outfile_dir'}){
        $log_dir = '.';
    }
    else{
        $log_dir=$param{'-outfile_dir'};
    }
    $log_dir=~s /\/+$//;
    $self->{'outdir'}=$log_dir;
    my $log_file=$log_dir."/".$self->{'local_time'}.".vEXAMer.log";
    open (LOG, ">$log_file")||die "can not write to log file\n";
    #LOG->autoflush(1);
    die "no configure file defined\n" unless (-f $param{'-configure_file'});
    $self->_read_configure_file($param{'-configure_file'});
    print LOG "configure file: $param{'-configure_file'}\n";


    foreach my $step (keys %{$self->{'configure'}}){
        $self->{'configure'}->{$step}->{'log_file'}=*LOG;
    }

	open (LIB, "$param{'-library_file'}")||die "no library file are defined\n";
    print LOG "Library file: $param{'-library_file'}\n";
    $self->_read_library_file(*LIB);

    if ($param{'-evidence_file'}){
        print LOG "evidence file: $param{'-evidence_file'}\n";
	    $self->_read_evidence_file($param{'evidence_file'}) ;
    }

    if ($param{'-contig_file'}){
        open (CONTIG, "$param{'-contig_file'}") || die "can not open contig file\n" ;
        print STDERR "read contig file\n";
        print LOG "contig file: $param{'-contig_file'}\n";
        $self->_read_contig_file(*CONTIG);
        print STDERR "done with reading contig file\n";
    }

	$self->_calculate_stepwise_sd_of_library_size();

}

sub _initiate_by_AGP_file{
    my ($self, $file, $step, $type)=@_;
    my $AGP_obj=new AGP (-infile=>$file);
	my $contig_length=$AGP_obj->get_contig_length_hash();
	foreach my $i (keys %$contig_length){
		$self->{'contig_length'}->{$i}=$contig_length->{$i};
	}
    my $hash=$AGP_obj->get_scaffold_hash();
    my $graph=Graph::Directed->new ;
    my $super_scaffold_hash; my $scaffold_superscaffold_hash;
    my $count=0;
    foreach my $key (keys %$hash){
        next if split (/-/, $hash->{$key}->{'path'}) == 1 and !($step == 1 and $type == 2);  # I do not want singleton to be included in the scaffold graph and data structures. They should in link edge graph only.
        $count++;
        my $last_unit;
        foreach my $unit (split (/-/, $hash->{$key}->{'path'}) ){
            $scaffold_superscaffold_hash->{$unit}=$count;
            if ($last_unit){
                $graph->add_edge($last_unit, $unit) if $type==1;
            }
	        $last_unit=$unit;
        }
		die "path is not defined for $count\n" if length($hash->{$key}->{'path'}) <=0;
        if ($hash->{$key}->{'dir'} =~ /^-+$/){ # found in honeybee project, some scaffold contain contig with all - direction, it should be translate into + since my scaffolder default all + for this case
            #$hash->{$key}->{'dir'} =~ tr/-/+/;
			($hash->{$key}->{'path'}, $hash->{$key}->{'gap'}, $hash->{$key}->{'dir'} )=reverse_DNA_string($hash->{$key}->{'path'}, $hash->{$key}->{'gap'}, $hash->{$key}->{'dir'});
        }
		$super_scaffold_hash->{$count}->{'main_path'}=$hash->{$key}->{'path'};
        $super_scaffold_hash->{$count}->{'direct'}=$hash->{$key}->{'dir'};
        $super_scaffold_hash->{$count}->{'gap'}=$hash->{$key}->{'gap'}; 
    }
    $self->{'scaffold_hash'}->{0}=$super_scaffold_hash if $step==1 and $type ==2 ;
    $self->{'scaffold_superscaffold'}->{0}=$scaffold_superscaffold_hash if $step==1 and $type==2;    
    return ($graph, $super_scaffold_hash, $scaffold_superscaffold_hash);
}  



sub do_scaffolding_step_by_step{
	my $self=shift;
    my $EXAMer_obj;
	foreach my $step (sort {$a<=>$b} keys %{$self->{'configure'}}){
        print STDERR  "###Starting step $step, type is $self->{'configure'}->{$step}->{'type'} \n";
        my $logfh=$self->{'configure'}->{$step}->{'log_file'};
        if ($self->{'configure'}->{$step}->{'type'} == 1 ){
            print $logfh "##############\n#############\nStarting step $step: Superscaffolding + Gap filling\n########\n########\n";
        }
        else{
            print $logfh "##############\n#############\nStarting step $step: Superscaffolding only\n########\n########\n";
        }
        if ($step == 1){
            $self->{'configure'}->{$step}->{'contig_length'}=$self->{'contig_length'};
            if ($self->{'configure'}->{$step}->{'type'} == 1 ){
                if ($self->{'agp_file'}){
                    print STDERR "Initializing vEXAMer with a AGP file for type 1 scaffolding at step 1: $self->{'agp_file'}\n";
                    print $logfh "Initializing vEXAMer with a AGP file for type 1 scaffolding at step 1: $self->{'agp_file'}\n";
                    ($self->{'configure'}->{1}->{'graph_to_build'},  $self->{'configure'}->{1}->{'super_scaffold'}, $self->{'configure'}->{1}->{'scaffold_superscaffold'})=$self->_initiate_by_AGP_file($self->{'agp_file'}, $step, $self->{'configure'}->{$step}->{'type'});
                }
            }
            if ($self->{'configure'}->{$step}->{'type'} == 2 ){
                die "no given scaffold structure in AGP file, did you give AGP file using -a option?" if not exists $self->{'agp_file'};
                print STDERR "Initializing and reformat vEXAMer with a AGP file for type 2 scaffolding at step 1: $self->{'agp_file'}\n";
                print $logfh "Initializing and reformat vEXAMer with a AGP file for type 2 scaffolding at step 1: $self->{'agp_file'}\n";
                $self->_initiate_by_AGP_file($self->{'agp_file'}, $step, $self->{'configure'}->{$step}->{'type'});
                $self->reformat_according_to_formed_scaffold($step);
            }
        }
        else{
            if ($self->{'configure'}->{$step}->{'type'} == 1 ){
                ($self->{'configure'}->{$step}->{'graph_to_build'},  $self->{'configure'}->{$step}->{'super_scaffold'}, $self->{'configure'}->{$step}->{'scaffold_superscaffold'})=$self->_initiate_by_AGP_file(  $self->_get_step_AGP($step-1), $step, $self->{'configure'}->{$step}->{'type'} );
                $self->{'configure'}->{$step}->{'contig_length'}=$self->{'contig_length'};
            }
            else{
                $self->reformat_according_to_formed_scaffold($step);
            }
        }
        $EXAMer_obj=new Atlas::EXAMer (-configure_hash=>$self->{'configure'}->{$step},
                            -link_infor_hash=>$self->{'link_infor_hash'}->{$step},
                            -step_sd_size=>$self->{'step_sd_size'}->{$step},
                            -lib_hash=>$self->{'lib'},
                            -lib_prefix_len=>$self->{'lib_configure'}->{'prefix_length'},
							-reformat_map=>$self->{'decode_path_hash'}->{$step}
                      );
        foreach my $key (sort keys %{$self->{'configure'}->{$step}}){
            print $logfh "$key=$self->{'configure'}->{$step}->{$key}\n" if exists {map {$_=>1} @{$self->{possible_configure_keys}} }->{$key};
            
        }
        $EXAMer_obj->do_superscaffolding();
        $self->{'scaffold_hash'}->{$step}=$EXAMer_obj->get_super_scaffold_hash();
        $self->{'scaffold_superscaffold'}->{$step}=$EXAMer_obj->get_scaffold_superscaffold_hash();
	    my $out_file=$self->{'outdir'}."/".$self->{'local_time'}.".vEXAMer.step$step.AGP";
        $self->{step_AGP}->{$step}=$out_file;
	    open (OUT, ">$out_file")||die "can not open to write\n";
        my $superscaffold_hash=$self->{'scaffold_hash'}->{$step};
        if (exists $self->{'decode_path_hash'}->{$step}){
            my $scaffold_superscaffold_hash;
            ($superscaffold_hash, $scaffold_superscaffold_hash)=$self->decode_scaffold_hash($step); 
            $self->{'scaffold_hash'}->{$step}=$superscaffold_hash;
            $self->{'scaffold_superscaffold'}->{$step}=$scaffold_superscaffold_hash;
        }
        $self->print_superscaffold_AGP(*OUT, $superscaffold_hash, $self->{'contig_length'}, $self->{'configure'}->{1}->{'treat_minus_gap_size_as'});
		my $AGP_obj=new Atlas::AGP(-infile=>$out_file);
		my $AGP_stat=$AGP_obj->get_scaffold_stat_no_gap();
		close OUT ;
        $self->print_scaffold_stat($logfh, $AGP_stat);
	}
}

sub get_original_path_for_reformatted_contig_at_step{
	my $self=shift;
	my $contig=shift || die "no contig defined\n";
	my $step=shift || die "no step defined\n";
	return $self->{'decode_path_hash'}->{$step}->{$contig};
}

sub _get_step_AGP{
    my $self=shift;
    my $step=shift;
    return $self->{step_AGP}->{$step};
}

sub decode_scaffold_hash{
    my ($self, $step)= @_;
    my $superscaffold_hash=$self->{'scaffold_hash'}->{$step};
    my $decode_hash=$self->{'decode_path_hash'}->{$step};
    my $new_superscaffold_hash;
    my $new_scaffold_superscaffold_hash;
	my $logfh=$self->{'configure'}->{$step}->{'log_file'};
    foreach my $i (keys %$superscaffold_hash){
        my ($new_path, $new_dir, $new_gap);
        my @old_path=split /-/, $superscaffold_hash->{$i}->{'main_path'};
        my @old_gap=split /-/, $superscaffold_hash->{$i}->{'gap'};
        my @old_direct=split //, $superscaffold_hash->{$i}->{'direct'};
        foreach my $j (0..$#old_path){
            my $replace_path_string=$decode_hash->{$old_path[$j]} || die "can not find replace_path_string for $old_path[$j] \n";
            my @path_for_replace=split(/-/, $replace_path_string);
            my $replace_dir_string=$self->{'scaffold_hash'}->{$step-1}->{$self->{'scaffold_superscaffold'}->{$step-1}->{$path_for_replace[0]} }->{'direct'} ||  die "can not find replace_dir_string for $old_path[$j] whose path is $replace_path_string\n" ;
            my $replace_gap_string=$self->{'scaffold_hash'}->{$step-1}->{$self->{'scaffold_superscaffold'}->{$step-1}->{$path_for_replace[0]} }->{'gap'} ;
            die "can not find replace_gap_string for $old_path[$j] whose path is $replace_path_string\n" if !defined $replace_gap_string and @path_for_replace >1;
			die "direction $old_direct[$j] in decode for $old_path[$j] : $replace_path_string is not legal.\n" if $old_direct[$j] !~ /-|\+/;
            if ($old_direct[$j] eq '-'){
                ($replace_path_string, $replace_gap_string, $replace_dir_string)=reverse_DNA_string($replace_path_string, $replace_gap_string, $replace_dir_string);    
            }
            $new_path.='-'.$replace_path_string;
            $new_dir.=$replace_dir_string;
            $new_gap.='-'.$replace_gap_string if @path_for_replace >1; 
            #a bug caught on 12/17/09,should add the below line
            $new_gap.='-'.$old_gap[$j] if @old_gap>0 and exists $old_gap[$j];
            ########
        }
        $new_path=~ s/^-*//;
        $new_gap=~ s/^-*//;
        $new_superscaffold_hash->{$i}->{'main_path'}=$new_path;
        $new_superscaffold_hash->{$i}->{'gap'}=$new_gap;
        $new_superscaffold_hash->{$i}->{'direct'}=$new_dir;
        foreach my $unit (split (/-/, $new_path) ){
            $new_scaffold_superscaffold_hash->{$unit}=$i;
        }
    }
    return ($new_superscaffold_hash, $new_scaffold_superscaffold_hash) ;
}

sub reformat_according_to_formed_scaffold{
	my ($self, $step)=@_;
	my $scaffold_superscaffold_hash=$self->{'scaffold_superscaffold'}->{$step-1};
	my $super_scaffold_hash=$self->{'scaffold_hash'}->{$step-1};
    $self->{'configure'}->{$step}->{'all_edge_graph'}=Graph::Directed->new;
    my $old_contig_length_hash=$self->{'contig_length'};
    delete $self->{'configure'}->{$step-1}->{'contig_length'} if exists $self->{'configure'}->{$step-1}->{'contig_length'};
    my $new_contig_length_hash;
    my $scaffold_reformatcontig;
    foreach my $i (keys %{$super_scaffold_hash}){
        my $reformat_contig= (scalar keys %{$self->{'decode_path_hash'}->{$step}} ) + 1 ;
        $self->{'decode_path_hash'}->{$step}->{$reformat_contig}=$super_scaffold_hash->{$i}->{'main_path'};
		#print "contig is ", $super_scaffold_hash->{$i}->{'main_path'}, "\n";
        $scaffold_reformatcontig->{$i}->{reformatcontig}=$reformat_contig;
        #$scaffold_reformatcontig->{$i}->{ori_path}=$super_scaffold_hash->{$i}->{'main_path'};
        $self->{'configure'}->{$step}->{'all_edge_graph'}->add_vertex($reformat_contig);
        my @gap=split (/-/, $super_scaffold_hash->{$i}->{'gap'} );
        my @path=split (/-/,$super_scaffold_hash->{$i}->{'main_path'});
        my $length=$old_contig_length_hash->{$path[0]};
        my $treat_minus_gap_size_as=$self->{'configure'}->{$step}->{'treat_minus_gap_size_as'} || die "treat gap zero\n";
        foreach my $i ( 1..$#path ) {
            if ( $gap[$i-1] =~/^\{.*\}$/){
                $length+=$treat_minus_gap_size_as;
            }
            else{
                $length+=int($gap[$i-1]);
            }
            $length+=$old_contig_length_hash->{$path[$i]};
        }
        $new_contig_length_hash->{$reformat_contig}=$length;
		#print "length of $super_scaffold_hash->{$i}->{'main_path'} is $length\n";
    }
    $self->{'configure'}->{$step}->{'contig_length'}=$new_contig_length_hash;
    
    my $link_infor_hash=$self->{'link_infor_hash'}->{$step};
	foreach my $source (keys %{$link_infor_hash}){
		if ($source eq 'mate_pair'){
			foreach my $contig_pair (sort keys %{$link_infor_hash->{'mate_pair'}}){
				my ($contig1, $contig2)=split(/=/, $contig_pair);
                #die "no scaffold found for contig1:$contig1 or contig2:$contig2\n" if not exists $scaffold_superscaffold_hash->{$contig1} or not exists $scaffold_superscaffold_hash->{$contig2};               
				if (not exists $scaffold_superscaffold_hash->{$contig1} or not exists $scaffold_superscaffold_hash->{$contig2}){ #some contigs are not thrown away in last step. Should throw away in this step too.
					delete $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$contig_pair};
                    next;
                }
                if ($scaffold_superscaffold_hash->{$contig1} == $scaffold_superscaffold_hash->{$contig2}){
                    delete $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$contig_pair};
                    next;
                }
                my $reformat_contig1=$scaffold_reformatcontig->{$scaffold_superscaffold_hash->{$contig1}}->{reformatcontig};  
                my $reformat_contig2=$scaffold_reformatcontig->{$scaffold_superscaffold_hash->{$contig2}}->{reformatcontig};
                die "how can reformat contig 1 $reformat_contig1 same as reformat contig 2 $reformat_contig2\n" if $reformat_contig1==$reformat_contig2;
				my $reformat_contig_pair=join ("=",  sort ($reformat_contig1, $reformat_contig2)) ;
                foreach my $template (keys %{$link_infor_hash->{'mate_pair'}->{$contig_pair}}){
                    foreach my $chem ('F3', 'R3'){
                        my $old_contig=$link_infor_hash->{'mate_pair'}->{$contig_pair}->{$template}->{$chem}->{'contig'}||die "can not find contig for a internal hash\n";
                        my $flag=0;
                        $flag++ if (exists {map { $_ => 1 } split (/-/, $self->{'decode_path_hash'}->{$step}->{$reformat_contig1}) }->{$old_contig}  );
                        $flag=$flag+2 if (exists {map { $_ => 1 } split (/-/, $self->{'decode_path_hash'}->{$step}->{$reformat_contig2}) } ->{$old_contig}  );
                        die "flag is is $flag $old_contig is in both path1 $self->{'decode_path_hash'}->{$step}->{$reformat_contig1} and path2 $self->{'decode_path_hash'}->{$step}->{$reformat_contig2} or neither\n" if $flag!=1 and $flag!=2;
                        my $reformat_contig=$flag==1 ? $reformat_contig1 : $reformat_contig2;
                        my $old_contig_dir=&direct_of_a_unit_in_path($old_contig, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$old_contig}}->{'main_path'}, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$old_contig}}->{'direct'});
                        $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$reformat_contig_pair}->{$template}->{$chem}->{'contig'}=$reformat_contig;
                        foreach my $start_end ('start', 'end'){
							#print "for $template $chem $start_end before is $contig_pair $link_infor_hash->{'mate_pair'}->{$contig_pair}->{$template}->{$chem}->{$start_end} ";
                            $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$reformat_contig_pair}->{$template}->{$chem}->{$start_end}=$self->convert_old_coordinate_2_new($link_infor_hash->{'mate_pair'}->{$contig_pair}->{$template}->{$chem}->{$start_end}, $old_contig_dir, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$old_contig} }, $old_contig, $old_contig_length_hash, $step );
							#print "After is $reformat_contig_pair $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$reformat_contig_pair}->{$template}->{$chem}->{$start_end} \n";
                        }
                    }
                }  
                delete $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$contig_pair};
			}
		}
        else{
            print STDERR "Need to test the code for non mate pair part\n";
            foreach my $lib (keys %{$link_infor_hash->{'non_mate_pair'}}){
                foreach my $contig_pair (sort keys %{$link_infor_hash->{'non_mate_pair'}->{$lib}}){
                    my ($contig1, $contig2)=split(/=/, $contig_pair);
                    die "no scaffold found for contig1:$contig1 or contig2:$contig2\n" if (not exists $scaffold_superscaffold_hash->{$contig1}) or (not exists $scaffold_superscaffold_hash->{$contig2}) ;                
                    if ($scaffold_superscaffold_hash->{$contig1} == $scaffold_superscaffold_hash->{$contig2}){
                        delete $self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$contig_pair};
                        next;
                    }
                    my $reformat_contig1=$scaffold_reformatcontig->{$scaffold_superscaffold_hash->{$contig1}}->{reformatcontig};
                    my $reformat_contig2=$scaffold_reformatcontig->{$scaffold_superscaffold_hash->{$contig2}}->{reformatcontig};
                    my $reformat_contig_pair=join ("=",  sort ($reformat_contig1, $reformat_contig2)) ;
                    $self->{'link_infor_hash'}->{$step}->{'non_mate_pair'}->{$lib}->{$reformat_contig_pair}=$self->convert_old_to_new_for_nonMatePair($link_infor_hash->{'non_mate_pair'}->{$lib}->{$contig_pair}, $contig1, $contig2 , $step,$old_contig_length_hash);
                    delete  $self->{'link_infor_hash'}->{$step}->{'non_mate_pair'}->{$lib}->{$contig_pair};
                }
            }
        }
	}
}

sub convert_old_to_new_for_nonMatePair{
    my ($self, $hash, $contig1, $contig2, $step, $old_contig_length_hash)=@_;
    my $new_hash;
    my $super_scaffold_hash=$self->{'scaffold_hash'}->{$step-1};
    my $scaffold_superscaffold_hash=$self->{'scaffold_superscaffold'}->{$step-1};
    
    my @main_path1=split(/-/, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$contig1}}->{'main_path'} ) ;
    my @main_path2=split(/-/, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$contig2}}->{'main_path'} );
    my @dir1=split(//, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$contig1}}->{'dir'}) ;
    my @dir2=split(//, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$contig2}}->{'dir'}) ;
    my @gap1=split(/-/, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$contig1}}->{'gap'} ) ;
    my @gap2=split(/-/, $super_scaffold_hash->{$scaffold_superscaffold_hash->{$contig2}}->{'gap'} ) ;
    my $new_dir_contig1; my $new_dir_contig2;
    my $flag1=0; my $flag2=1;
    my $new_gap;
    my $treat_minus_gap_size_as=$self->{'configure'}->{$step}->{'treat_minus_gap_size_as'} || die "treat gap zero\n";
    foreach my $i ( 0..$#main_path1 ){
        if ($main_path1[$i] eq $contig1){
            $new_dir_contig1=$dir1[$i];
            $flag1=1;
        }
        if (defined $hash->{'gap'} and $flag1 ){
            if ( $gap1[$i] =~/^\{.*\}$/){
                $new_gap=-$treat_minus_gap_size_as;
            }
            else{
                $new_gap=-$gap1[$i];
            }
            $new_gap-=$old_contig_length_hash->{$main_path1[$i]};
        }
    }
    foreach my $i ( 0..$#main_path2 ){
        if ($main_path2[$i] eq $contig2){
            $new_dir_contig2=$dir2[$i];
            $flag2=0;
        }
        if (defined $hash->{'gap'} and $flag2 ){
            if ( $gap2[$i] =~/^\{.*\}$/){
                $new_gap=-$treat_minus_gap_size_as;
            }
            else{
                $new_gap=-$gap2[$i];
            }
            $new_gap-=$old_contig_length_hash->{$main_path2[$i]};
        }
    }
    $new_hash->{'gap'}=$new_gap;
    $new_hash->{'pair_dir'}=$new_dir_contig1.$new_dir_contig2;
    $new_hash->{'weight'}=$hash->{'weight'};
    return $new_hash;
}

sub convert_old_coordinate_2_new{
    my ($self, $coord, $dir, $hash, $old_contig, $old_contig_length_hash, $level)=@_;
    die "no input for function\n" if (!defined($coord) or !defined($dir) or !defined($hash));
    my $treat_minus_gap_size_as=$self->{'configure'}->{$level}->{'treat_minus_gap_size_as'} || die "treat gap zero\n"; 
    my @path=split(/-/, $hash->{'main_path'});
    my @gap=split(/-/, $hash->{'gap'});
    my $offset=0;
    foreach my $i ( 0..$#path ) {
        if ($path[$i] eq $old_contig ){
            if ($dir eq '+' ){
                if ($coord >= 0 ){
                    return $coord+$offset;
                }    
                else{
                    return $coord-$offset;
                }
            }
            else{
                my $rev_coor=$old_contig_length_hash->{$path[$i]} - $coord + 1 ;
                #die "for $path[$i] in superscaffold $hash->{'main_path'}, its orig coordinate is $coord, length of of contig is $old_contig_length_hash->{$path[$i]},  rev coor is less than zero: $rev_coor\n" if $rev_coor<=0;
                if ($coord >= 0 ){
                    return 0-( abs($coord) + $offset ) ;
                }    
                else{
                    return abs($coord) + $offset;                
                }
            }
        }
        else{
            $offset+=$old_contig_length_hash->{$path[$i]};
            if ( $gap[$i] =~/^\{.*\}$/){
                $offset+=$treat_minus_gap_size_as;
            }
            else{
                $offset+=int($gap[$i]);
            }
        }
    }
}
		
sub _read_configure_file{
    my ($self, $file)=@_;
    my $parser = new XML::DOM::Parser;
    my $doc = $parser->parsefile ($file);
    my $nodes = $doc->getElementsByTagName ("step");
    my $n = $nodes->getLength;
    for (my $i = 0; $i < $n; $i++){
        my $node = $nodes->item ($i);
        my $level = $node->getAttribute ("level") ||die "NO level defined for some step\n";
        die "in configure file level is not correctly defined for step tag, should be an integer from 1 to n step in order\n" if $level !~ /^\d+$/ or $level <1 or $level != 1+ $i;
        my $type = $node->getAttribute ("type") ||die "NO scaffolding type defined for some link in evidence file\n";
        die "type is not defined as 1 or 2\n" if $type !=1 and $type != 2;
        $self->{'configure'}->{$level}->{'type'}=$type;
        for my $kid ($node->getChildNodes) {
            next if $kid->getNodeName() eq '#text';
            if ($kid->getNodeName() ne "mate_pair" and $kid->getNodeName() ne "non_mate_pair"){ #configure options
                my $key=$kid->getNodeName(); my $value=&traverse($kid);
                die "tag $key not recognized by vEXAMer\n" if not exists {map {$_=>1} @{$self->{possible_configure_keys}} }->{$key};
				#if ($key  eq 'remove_redundant_reads' and $value ne '454' and $value !~ /Solid/i){
				#	die "$key can not defined as $value, should be anyone of 454 or Solid\n";
				#}
                $self->{'configure'}->{$level}->{$key}  = $value;
            }
            elsif ($kid->getNodeName() eq "mate_pair"){
                for my $kid_m ($kid->getChildNodes) {
                    if ($kid_m->getNodeName() eq "lib"){
                        push (@{$self->{'step_matelib'}->{$level}->{'lib'}}, &traverse($kid_m) );
                    }
                    if ($kid_m->getNodeName() eq "lib_size_limit"){
                        my $min= $kid_m->getAttribute ("min"); die "min size is not defined in lib_size_limit tag\n" if !defined $min;
                        my $max= $kid_m->getAttribute ("max"); die "min size is not defined in lib_size_limit tag\n" if !defined $max;
                        $self->{'step_matelib'}->{$level}->{'mate_pair_min_lib_size'} = $min;
                        $self->{'step_matelib'}->{$level}->{'mate_pair_max_lib_size'} = $max;             
                    }
                }
            }
            elsif ($kid->getNodeName() eq "non_mate_pair"){
                for my $kid_m ($kid->getChildNodes) {
                    if ($kid_m->getNodeName() eq "lib"){
                        my $lib=&traverse($kid_m) ;
                        push (@{$self->{'step_nonmatelib'}->{$level}->{'lib'} }, $lib );
                    }
                }
            }
        }
        die "when use type 2 scaffoding at step 1, you should give a initial AGP file use -a option\n" if ($level == 1 and $type ==2 and not exists $self->{'agp_file'});
        #give default value if undef
        $self->{'configure'}->{$level}->{'min_links'}=2 if not exists $self->{'configure'}->{$level}->{'min_links'};
        $self->{'configure'}->{$level}->{'treat_as_no_link_if_less_than'}=1 if not exists$self->{'configure'}->{$level}->{'treat_as_no_link_if_less_than'};
        $self->{'configure'}->{$level}->{'min_unit_to_form_cluster'}=1 if not exists $self->{'configure'}->{$level}->{'min_unit_to_form_cluster'};
        $self->{'configure'}->{$level}->{'distance'}=1000 if not exists $self->{'configure'}->{$level}->{'distance'};
        $self->{'configure'}->{$level}->{'force_elongate_min_link'}=0 if not exists $self->{'configure'}->{$level}->{'force_elongate_min_link'};
        $self->{'configure'}->{$level}->{'remove_redundant_reads'}='454' if not exists $self->{'configure'}->{$level}->{'remove_redundant_reads'};
        $self->{'configure'}->{$level}->{'deviate_factor'} = 15 if not exists $self->{'configure'}->{$level}->{'deviate_factor'};
        $self->{'configure'}->{$level}->{'dominant_ratio'} = 0.6 if not exists $self->{'configure'}->{$level}->{'dominant_ratio'};
        $self->{'configure'}->{$level}->{'treat_minus_gap_size_as'}=50 if not exists $self->{'configure'}->{$level}->{'treat_minus_gap_size_as'}; #all minus gap will be treated as 50 bp
        $self->{'configure'}->{$level}->{'overrule_strength_for_reinitiate'}=1 if not exists $self->{'configure'}->{$level}->{'overrule_strength_for_reinitiate'}; # if the current_edge's_weight is bigger than this_factorXcurrent_edge's_weight, it can re-initiate the possible extending or merging evaluation. 
		$self->{'configure'}->{$level}->{'favor_density'}=0 if not exists $self->{'configure'}->{$level}->{'favor_density'};
	}
}

sub _read_library_file{
    my ($self, $fh)=@_;
    my $read_lib;
    while (my $line =<$fh>){
        chomp $line;
        if ($line=~ /^prefix_length/){
            my @arr=split(/\s+/, $line);
            $self->{'lib_configure'}->{'prefix_length'}=$arr[1] ||die"no correct prefix length defined\n";
        }
        if ($line=~ /^F_suffix/){
            my @arr=split(/\s+/, $line);
            $self->{'lib_configure'}->{'fwd_suffix'}=$arr[1] || die "no fwd suffix define\n";
        }
        if ($line=~ /^R_suffix/){
            my @arr=split(/\s+/, $line);
            $self->{'lib_configure'}->{'rev_suffix'}=$arr[1] || die "no rev suffix define\n";
        }
        if ($line=~ /^library/){
            my @arr=split(/\s+/, $line);
            $self->{'lib'}->{$arr[1]}->{'ave_size'}=$arr[2];
            $self->{'lib'}->{$arr[1]}->{'sd_size'}=$arr[3];
        }
    }
	die "fwd suffix $self->{'lib_configure'}->{'fwd_suffix'} is same as rev_suffix $self->{'lib_configure'}->{'rev_suffix'} \n" if $self->{'lib_configure'}->{'fwd_suffix'} eq $self->{'lib_configure'}->{'rev_suffix'};
}

sub _read_evidence_file{
    my ($self, $file)=@_;
	my $non_mate_pairs_hash;
	my $parser = new XML::DOM::Parser;
	my $doc = $parser->parsefile ("$file");
	my $nodes = $doc->getElementsByTagName ("LINK");
	my $n = $nodes->getLength;
	for (my $i = 0; $i < $n; $i++){
		my ($contig_pair, $contig1, $contig2, $dir, $threePrime_contig, $gap, $weight);
		my $node = $nodes->item ($i);
		my $lib = $node->getAttribute ("name") ||die "NO ID defined for some link in evidence file\n";
		for my $kid ($node->getChildNodes) {
			if ($kid->getNodeName() eq "contig1"){
				$contig1 = &traverse($kid);
			}
			if ($kid->getNodeName() eq 'contig2'){
				$contig2 = &traverse($kid);
			}
			if ($kid->getNodeName() eq 'dir'){
				$dir= &traverse($kid);
			}
			if ($kid->getNodeName() eq 'gap'){
                $gap= &traverse($kid);
            }
			if ($kid->getNodeName() eq 'weight'){
                $weight= &traverse($kid);
            }
		}
		die "contig information missing in evidence file of LINK ID $lib\n" if !$contig1 or !$contig2;
		$contig_pair=join ('=', sort ($contig1, $contig2));
		die "die information missing for misdefine when dir is $dir\nCorrect dir should be -- ++ -+ or +-\n" if !$dir or ($dir=~ /[^+-]/) or (length($dir) !=2);
		$non_mate_pairs_hash->{$lib}->{$contig_pair}->{'gap'}=$gap if defined $gap;
		$non_mate_pairs_hash->{$lib}->{$contig_pair}->{'pair_dir'}=$dir;
		$non_mate_pairs_hash->{$lib}->{$contig_pair}->{'weight'}=$weight;
    }
	foreach my $lib (keys %$non_mate_pairs_hash){
		foreach my $step (keys %{$self->{'step_nonmatelib'}}){
			if (grep (/^$lib$/, @{$self->{'step_nonmatelib'}->{$step}->{'lib'}}) ){
				$self->{'link_infor_hash'}->{$step}->{'non_mate_pair'}->{$lib}= $non_mate_pairs_hash->{$lib} ;
			}
		}
	}
}


sub _read_contig_file{
	my ($self, $fh)=@_;
	my $contig;
	my $template_reads;
	my $template_lib;
    while (<$fh>){
		chomp;
		if ($_ =~ /^##(\S+)\s+\d+\s+(\d+)/){
			$contig=$1;
			$self->{'contig_length'}->{$contig}=$2 if not exists $self->{'contig_length'}->{$contig}; #contig length has been initiated when reading AGP. So, if length different here, use AGP's length
		}
		elsif ($_ =~ /^#(\S+)\(/){
			my $id=$1;
			my $len=$self->{'lib_configure'}->{'prefix_length'};
			$id =~ /^(\S{$len})/ || die "read id $id is not expected regexp\n";
			my $lib = $1;
			die "library $lib of the read $id is not found in library\n" if not exists $self->{'lib'}->{$lib}->{'ave_size'};
			my $suffix_f=$self->{'lib_configure'}->{'fwd_suffix'} || die "not define fwd suffix\n";
			my $suffix_r=$self->{'lib_configure'}->{'rev_suffix'} || die "not define rev suffix\n";
			my $f_r;
			my $template;
			if ($id=~ /^(\S+?)$suffix_f$/){
				$f_r= 'f';
				$template=$1;
			}
			elsif ($id=~ /^(\S+?)$suffix_r$/){
				$f_r= 'r' ;
				$template=$1;
			}
			else{
				die"read name format wrong, no fwd or rev suffix found when id is $id\n"
			}
			$template_lib->{$template}=$lib;
			$self->{'lib'}->{$lib}->{'count'}++;
			my $hash;
			$hash->{'f_r'}=$f_r;
			$hash->{'contig'}=$contig;
			$hash->{'dir'}=$_ =~ /\[RC\]/i ? "-" : "+";
			$_ =~ /<(\d+)\s+(\d+)>/ ||die "$_ is not of pattern expected\n";
			$hash->{'start'} =$1;
			$hash->{'end'} =$2;
			my $uniformized_hash=$self->_convert_input_mate_pair_format($hash);
			my @keys=keys %$uniformized_hash;
			$template_reads->{$template}->{$keys[0]}= $uniformized_hash->{$keys[0]};
		}
		else{
			die "contig file format not recogized for $_\n";
		}
    }
	foreach my $template (keys %$template_reads){
		if (keys %{$template_reads->{$template}} !=2){
			#warn "$template do not have either F or R reads or not just have F R reads, it only has ", join (' ', keys %{$template_reads->{$template}}) , "\n";
			next;
		}
		else{
			next if $template_reads->{$template}->{'F3'}->{'contig'} eq $template_reads->{$template}->{'R3'}->{'contig'}; # add to be able to deal with some reads in same contig 2/1/10
			#print "template is $template f3 contig is ", $template_reads->{$template}->{'F3'}->{'contig'}, " and r3 contig is ", $template_reads->{$template}->{'R3'}->{'contig'}, "\n" if $template_reads->{$template}->{'F3'}->{'contig'} =~/contig00021/ or $template_reads->{$template}->{'R3'}->{'contig'}=~ /contig00021/;
			my $contig_pair=join('=', sort( $template_reads->{$template}->{'F3'}->{'contig'}, $template_reads->{$template}->{'R3'}->{'contig'}) );
			my $steps=$self->_get_which_steps_of_the_mate_pair_lib($template_lib->{$template})  ;
			die "how can steps is empty \n" if @$steps==0;
			foreach my $step (@$steps){
				$self->{'link_infor_hash'}->{$step}->{'mate_pair'}->{$contig_pair}->{$template} = $template_reads->{$template};
            }
		}
	}
}

sub _get_which_steps_of_the_mate_pair_lib{
	my ($self, $lib)=@_;
	my $size=$self->{'lib'}->{$lib}->{'ave_size'};	
	my @steps;
	foreach my $step (keys %{$self->{'step_matelib'}}){
		if (defined ($self->{'step_matelib'}->{$step}->{'mate_pair_max_lib_size'} ) ){
			if ($size >=$self->{'step_matelib'}->{$step}->{'mate_pair_min_lib_size'} and $size <$self->{'step_matelib'}->{$step}->{'mate_pair_max_lib_size'}) {
				push (@steps, $step);
				push (@{$self->{'step_matelib'}->{$step}->{'lib'}}, $lib) if (not exists {map { $_ => 1 } @{$self->{'step_matelib'}->{$step}->{'lib'}} }->{$lib});
			}
		}
        elsif (exists {map { $_ => 1 } @{$self->{'step_matelib'}->{$step}->{'lib'}} }->{$lib}   ){
			push (@steps, $step) ;
		}
	}
	return &unique(\@steps);
}


sub _convert_input_mate_pair_format{ 
    # contig file is defined similarly as BAMBUS use. It is sanger format, so need to convert to solid format, which is my EXAMer want.
	my ($self, $hash)=@_;
	my @keys=keys %$hash;
	die "should only has five keys\n" if @keys !=5 ;
	my $dir =  $hash->{'dir'} ;
	my $convert;
	my $F3_or_R3 ;
	if ($hash->{'f_r'} eq 'r'){
		$dir =  $hash->{'dir'} eq '-' ? "+" : "-";
	}
	$F3_or_R3 = $hash->{'f_r'} eq 'r' ? 'F3' : 'R3';
	$convert->{$F3_or_R3 }->{'contig'}=$hash->{'contig'};
	if ($dir eq '+'){
		$convert->{$F3_or_R3}->{'start'}=$hash->{'start'};
		$convert->{$F3_or_R3}->{'end'}=$hash->{'end'};
	}
	else{
		$convert->{$F3_or_R3}->{'start'}=$dir.$hash->{'start'};
		$convert->{$F3_or_R3}->{'end'}=$dir.$hash->{'end'};
	}
	#print "hash start is $hash->{'start'} and hash end is $hash->{'end'} ; start is $convert->{$F3_or_R3}->{'start'} and end is $convert->{$F3_or_R3}->{'end'}\n";
	return $convert;
}


sub _calculate_stepwise_sd_of_library_size{
	my $self=shift;
	foreach my $step (sort {$a<=>$b} keys %{$self->{'configure'}}){
        my $logfh=$self->{'configure'}->{$step}->{'log_file'};
		my @sd_sizes;
		my @weights;
		my @libs_at_this_step=@{$self->{'step_matelib'}->{$step}->{'lib'}} if exists $self->{'step_matelib'}->{$step}->{'lib'};
		print $logfh  "Step $step contains libraries: ",join (' ', @libs_at_this_step), "\n" if @libs_at_this_step >0 ;
        print $logfh  "Step $step is non-mate-pair-only\n" unless  @libs_at_this_step;
		foreach my $lib (@libs_at_this_step){
			if (defined $self->{'lib'}->{$lib}->{'sd_size'} and $self->{'lib'}->{$lib}->{'count'} ){
				push (@sd_sizes, $self->{'lib'}->{$lib}->{'sd_size'}) ;
				push (@weights, $self->{'lib'}->{$lib}->{'count'}) ;
				print $logfh "lib $lib meansize is $self->{'lib'}->{$lib}->{'ave_size'}; sdsize is $self->{'lib'}->{$lib}->{'sd_size'} and count is $self->{'lib'}->{$lib}->{'count'}\n";
			}
		}
        if (@libs_at_this_step){
		    my $stat  = Statistics::Descriptive::Weighted::Full->new();
		    $stat->add_data(\@sd_sizes,\@weights);
		    $self->{'step_sd_size'}->{$step}=$stat->mean();
		    print  $logfh "Step $step weighted mean of standard deviation of all libraries is $self->{'step_sd_size'}->{$step}\n";
        }
	}
}


sub traverse{
  my($node)= @_;
  if ($node->getNodeType == ELEMENT_NODE)
  {
      foreach my $child ($node->getChildNodes())
     {
       return traverse($child);
     }
  }
  elsif ($node->getNodeType() == TEXT_NODE)
  {
    my $data=$node->getData;
    return $data;
  }
}

sub get_link_info_hash_for_step{
	my $self=shift;
	my $step=shift;
	my $type=shift;
	return $self->{'link_infor_hash'}->{$step}->{$type};
}

sub check_circle_scaffold_of_step_A_using_links_of_step_B{
	my ($self, $stepA, $stepB)=@_;
	my $EXAMer_obj=new EXAMer (-configure_hash=>$self->{'configure'}->{$stepB},
                            -link_infor_hash=>$self->{'link_infor_hash'}->{$stepB},
                            -step_sd_size=>$self->{'step_sd_size'}->{$stepB},
                            -lib_hash=>$self->{'lib'},
                            -lib_prefix_len=>$self->{'lib_configure'}->{'prefix_length'},
                      );
	$EXAMer_obj->Add_weighted_edges_to_Edge_set();
	my $edge_graph=$EXAMer_obj->get_edge_graph;
	foreach my $keys (keys %{$self->{'scaffold_hash'}->{$stepA}}){
		my @path=split(/-/, $self->{'scaffold_hash'}->{$stepA}->{'main_path'}) ;
		next if @path<3; #important!!!. If two contigs forms a circular genome, then I can not detect them because I assume two neighboring contig can only has one link orientation. I do not allow this case 1->2 and 2<-1. However, if 1->2->3->1, then, I report circular;
		if($edge_graph->has_edge($path[0], $path[-1])  or $edge_graph->has_edge($path[-1], $path[0]) ){
			print STDERR $self->{'scaffold_hash'}->{$stepA}->{'main_path'}, " could be a circlar genome\n";
		}
	}
}

sub print_superscaffold_AGP{
    my ($self, $fh, $hash, $contig_length, $treat_minus_gap_size_as)=@_;
    if (ref $hash ne 'HASH' or ref $contig_length ne 'HASH'){
        print STDERR  "scaffold or contig length hash not ready\n" ;
        return ;
    }
    $treat_minus_gap_size_as || die "treat gap zero\n";
    my $scaffold_idx=1;
    foreach my $i (sort {my $al=$hash->{$a}->{'main_path'};
                         my $bl=$hash->{$b}->{'main_path'};
                         split(/-/, $bl) <=>split(/-/, $al) } keys %$hash){
        my @dir=split(//, $hash->{$i}->{'direct'});
        my @gap=split(/-/, $hash->{$i}->{'gap'});
        my @main_path=split (/-/, $hash->{$i}->{'main_path'} );
        die "When try to print AGP, found the $scaffold_idx scaffold, main path is $hash->{$i}->{'main_path'}, dir is $hash->{$i}->{'direct'} and gap is $hash->{$i}->{'gap'}, number not match\n" if (@main_path != @dir or @main_path != @gap+1 ) ;
        my $offset=0;
        my $unit_idx=1;
        foreach my $unit (@main_path) {
            print $fh "Scaffold$scaffold_idx\t", $offset+1, "\t", $offset+$contig_length->{$unit},"\t",$unit_idx,"\t", "W", "\t",$unit,"\t","1\t", $contig_length->{$unit},"\t",shift @dir,"\n";
            $unit_idx++;
            $offset+=$contig_length->{$unit};
            if (my $gap_length=shift @gap){
                if ($gap_length =~ /\{\S+\}/){
                    $gap_length=$treat_minus_gap_size_as;
                }
                $gap_length=int $gap_length; die "gap length is print function found is $gap_length, less than 0\n" if $gap_length<0;
				$gap_length=$treat_minus_gap_size_as if $gap_length<=$treat_minus_gap_size_as;
                #print $fh "Scaffold$scaffold_idx\t",$offset+1, "\t", $gap_length+$offset, "\t$unit_idx\tN\t$gap_length\tFragment\tyes\n" if $gap_length>0;
				print $fh "Scaffold$scaffold_idx\t",$offset+1, "\t", $gap_length+$offset, "\t$unit_idx\tN\t$gap_length\tfragment\tyes\n";
                $unit_idx++;
                $offset+=$gap_length;
            }
        }
        $scaffold_idx++;
    }	
	close $fh;
}

sub print_scaffold_stat{
	my ($self, $fh, $hash)=@_;
    foreach my $key (sort keys %$hash){
        print $fh $key,"\t", $hash->{$key},"\n";
    }
}

1;

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 vEXAMer 

Atlas::vEXAMer - stands for versatile EXtending And Merging scaffolder
This is the upgraded version of my pervious Super_Scaffolder, compare to previous version, it can
1, hierarchical use mate pairs from user designiated lirabries in building scaffold(step by step, any libraries can be used at any step, however, it is intuitive to use small insert library before large ones).
2, use other supporting evidence in building scaffold
3, can upgrade existing scaffold structure  by giving linking information. And two options can be choosen 1.Aggressive upgrade 2. Conserve updgrade
        aggressive upgrade, type=1 Superscaffolding plus gap filling. Singleton contigs or contigs from one scaffold can inserted into the middle of another scaffold if enough information support this interculation. To do this, a AGP f
ile need to be given and all the mate pairs coordinate should be on contigs. Original order of contigs in a original scaffold will not change except new contigs from other scaffold can be inserted into them.
        conserved upgrade, type=2  Superscaffolding only. Mate pairs coordinate should given on scaffold. Scaffold are regarded as contig and a de novo scaffolding will be carried out.


=head1 SYNOPSIS

  use vEXAMer;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for EXAMer, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

Jixin Deng, E<lt>jdeng@bcm.edu<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Jixin Deng

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut

