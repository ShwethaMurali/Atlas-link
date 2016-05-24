package Atlas::EXAMer;

use strict;
use FindBin;
use lib "$FindBin::Bin/lib";
use lib "$FindBin::Bin/Lib/Graph_94/lib";
use Algorithm::ClusterPoints;
use Graph::Directed;
use List::Util qw[min max];
require Exporter;

our @ISA = qw(Exporter);
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration       use vEXAMer ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK

our @EXPORT = qw(direct_of_a_unit_in_path reverse_DNA_string completement_path unique print_superscaffold_AGP_export);

our %EXPORT_TAGS = ( 'all' => [ qw(
        
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our $VERSION = '0.01';


# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.



sub new{
    my ($caller,@args ) = @_; # some parameter to initiate a genome
    my $class = ref($caller) || $caller;
    my $self = {};
    bless $self, $class;
    $self->_initialize(@args );
    return $self;
}

sub _initialize{
	my ($self,@args)=@_;
	my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase key
	foreach my $key (sort keys %{$param{'-configure_hash'}}){
		$self->{$key}=$param{'-configure_hash'}->{$key};
	}
	unless (exists $self->{'graph_to_build'} ){
		$self->{'graph_to_build'} = Graph::Directed->new;
    }
    else{
        $self->{'indx'}=scalar keys %{$self->{'super_scaffold'}};
    }
	unless (exists $self->{'all_edge_graph'} ){
		$self->{'all_edge_graph'} = Graph::Directed->new;
	}
    foreach my $unit (keys %{$self->{'contig_length'}}){
		die "unit is not defined \n" if length($unit)==0;
        $self->{'all_edge_graph'}->add_vertex($unit);
    }
	$self->{'sd_lib_size'}=$param{'-step_sd_size'};
	$self->{'lib'}=$param{'-lib_hash'} ||die "no lib information for EXAMer\n";
	$self->{'lib_prefix_len'}=$param{'-lib_prefix_len'} ||die "no lib prefix len information for EXAMer\n";
    $self->{'link_infor_hash'}=$param{'-link_infor_hash'} ||die "no link infor hash for EXAMer\n";
	$self->{'reformat_map'}=$param{'-reformat_map'} if exists $param{'-reformat_map'};
    $self->{'updated_flag'}=0;
    my $logfh=$self->{'log_file'};
    if (exists $self->{'super_scaffold'} ){
        print $logfh "importing pre-existed ", scalar keys %{$self->{'super_scaffold'}}, " scaffolds\n";
    }
    foreach my $key (sort keys %$self){
        if (ref($self->{$key}) eq 'SCALAR'){
            print $logfh  "$key = $self->{$key} \n";
        }
    }
    print $logfh "step specific standard deviation of library size for this step is $self->{'sd_lib_size'}\n";
	$self->{'stat'}->{'Matepair1 N_of_input_matepairs'}=0;
	$self->{'stat'}->{'Matepair1 effective_clone_cover_span'}=0;
	$self->{'stat'}->{'Matepair2 N_of_input_matepairs_aft_remove_redd_reads'}=0;
	$self->{'stat'}->{'Matepair2 effective_clone_cover_span'}=0;
	$self->{'stat'}->{'Matepair3.1 N_of_failed_dominant_rule'}=0;
	$self->{'stat'}->{'Matepair3.2 N_of_removed_minorities'}=0;
	$self->{'stat'}->{'Matepair3.3 N_of_less_than_min_link'}=0;
	$self->{'stat'}->{'Matepair4 N_of_matepairs_pass_filter'}=0;
	$self->{'stat'}->{'Matepair4.1 N_assoicated_removed_edges'}=0;
	$self->{'stat'}->{'Matepair4.2 N_assoicated_excessive_insertsize_edges'}=0;
	$self->{'stat'}->{'Matepair5 N_of_matepairs_in_use'}=0;
	$self->{'stat'}->{'Matepair6 N_of_matepairs_Successed'}=0;
}

sub do_superscaffolding{
    my $self=shift;
    my $logfh=$self->{'log_file'};
    $self->Add_weighted_edges_to_Edge_set();
    print STDERR "finish adding edge\n";
	print STDERR "excessive_bridges_dev_factor is ", $self->{'excessive_bridges_dev_factor'}, "\n";
	$self->_remove_noisy_edges if $self->{'excessive_bridges_dev_factor'}>0 ;
	print "done with remove noisy edge\n";
    $self->_remove_multiple_briging_nodes() if exists $self->{'remove_node_bridges_more_than'} ;
	$self->_remove_short_length_nodes() if exists $self->{'remove_node_if_shorter_than'}; 
	###### Should first do _remove_multiple_briging_nodes then _remove_short_length_nodes
	#because otherwise, if you first remove short length nodes, some multple bridging nodes will longer multiple bridging because they are linking to short contigs;
	########################
    if ( defined $self->{'excessive_mate_pairs_for_a_edge_dev_factor'}){        
        print $logfh "will remove edges with ", int($self->{'clustered_mate_pair_count_hash'}->{'mean'}+$self->{'excessive_mate_pairs_for_a_edge_dev_factor'}*$self->{'clustered_mate_pair_count_hash'}->{'standard_deviation'}),"mate pairs or more in its largest cluster\n";
    }
    my $loop=0;
    while ($self->{'updated_flag'} or $loop==0){
        $loop++;
        print STDERR "recursive loop $loop\n";
		$self->{'stat'}->{'1 N_of_bridges_remov_due_excessive_matepair'}=0;
		$self->{'stat'}->{'2 N_of_edges_excessive_inst_size'}=0;
		$self->{'stat'}->{'3 N_of_bridges_for_EXAMing'}=0;
		$self->{'stat'}->{'4 N_of_bridge_for_denovo_join'}=0;
		$self->{'stat'}->{'5 N_of_bridges_for_extend'}=0;
		$self->{'stat'}->{'6 N_of_bridges_for_merge'}=0;
        $self->_Build_superscaffold_graph($loop);
        $self->print_graph_usage_stat($logfh);
    }
    $self->_pick_up_singletons_into_hash();
}

sub print_graph_usage_stat{
    my $self=shift;
    my $fh=shift;
    foreach my $key (sort keys %{$self->{'stat'}}){
        print $fh "$key $self->{'stat'}->{$key}\n" if $fh;
        print $fh "$key $self->{'stat'}->{$key}\n" unless $fh;
    }
}

sub _pick_up_singletons_into_hash{
    my $self=shift;
    foreach my $vertex ($self->{'all_edge_graph'}->vertices){
        unless ($self->{'scaffold_superscaffold'}->{$vertex}){
			die "vertex is $vertex, not defined, but picked up in _pick_up_singletons_into_hash\n" if !$vertex;
			if (exists $self->{'min_singleton_contig_length'}){
				next if $self->{'contig_length'}->{$vertex} < $self->{'min_singleton_contig_length'} ;
			}
            my $current_superscaffold=$self->{'indx'}+1;
            $self->{'indx'}++;
            $self->{'super_scaffold'}->{$current_superscaffold}->{'direct'}='+';
            $self->{'super_scaffold'}->{$current_superscaffold}->{'main_path'}=$vertex;
            $self->{'scaffold_superscaffold'}->{$vertex}=$current_superscaffold;
        }
    }
}


sub Add_weighted_edges_to_Edge_set{
	my $self=shift;
    my $logfh=$self->{'log_file'};
	if (exists $self->{'link_infor_hash'}->{'mate_pair'}){
		print STDERR "bundling edges for mate pairs\n";
		$self->_bundle_edges_from_different_sources('mate_pair', $self->{'link_infor_hash'}->{'mate_pair'});
		if ($self->{'edge_count'} != 0 ){
        	$self->{'clustered_mate_pair_count_hash'}->{'standard_deviation'}= sqrt( $self->{'clustered_mate_pair_count_hash'}->{'var'} / ( $self->{'edge_count'} -1 ) );
        	print $logfh "For mate pair, bundled ", scalar $self->{'all_edge_graph'}->vertices  , " total vetecies and ", scalar $self->{'all_edge_graph'}->edges , " edges\nMean mate pairs for a edge is $self->{'clustered_mate_pair_count_hash'}->{'mean'} sdev is $self->{'clustered_mate_pair_count_hash'}->{'standard_deviation'}\n";
		}	
		else{
			print $logfh "For mate pair, No valid edge are bundled\n";
		}
	}
	if (exists $self->{'link_infor_hash'}->{'non_mate_pair'}){
        print STDERR "bundling edges for non mate pairs\n";
        foreach my $library (keys %{ $self->{'link_infor_hash'}->{'non_mate_pair'} }){
            $self->_bundle_edges_from_different_sources('non_mate_pair',$self->{'link_infor_hash'}->{'non_mate_pair'}->{$library}->{'hash'});
        }
        print $logfh "After bundling non  mate pair, bundled ", scalar $self->{'all_edge_graph'}->vertices  , " total vetecies and ", scalar $self->{'all_edge_graph'}->edges , " edges\n";
	}
}

sub _remove_short_length_nodes{
    my $self=shift;
    my $logfh=$self->{'log_file'};
	foreach my $contig (keys %{$self->{'contig_length'}}){
		if ($self->{'contig_length'}->{$contig} < $self->{'min_length_to_be_node'}){
			print $logfh "delete contig $contig for short length\n";
			$self->{'all_edge_graph'}->delete_vertex($contig);
			push (@{$self->{'removed_nodes'}}, $contig);
		}
	}
}

sub _remove_multiple_briging_nodes{
	my $self=shift;
    my $logfh=$self->{'log_file'};
	foreach my $vertex ($self->{'all_edge_graph'}->vertices){
		if ($self->{'force_inclusion_if_unit_size_more_than'}){
			if ($self->{'force_inclusion_if_unit_size_more_than'} < $self->{'contig_length'}->{$vertex} ){
				print $logfh "always include $vertex for large length no matter how many bridges linking it\n";
				next;
			}
		}
		if (scalar $self->{'all_edge_graph'}->predecessors($vertex) + scalar $self->{'all_edge_graph'}->successors($vertex) > $self->{'remove_node_bridges_more_than'} ){
			print $logfh "$vertex has predessor links count ",  $self->{'all_edge_graph'}->predecessors($vertex), " and successor link count ", $self->{'all_edge_graph'}->successors($vertex), " vertexes, is regarded as repetitive nodes, will be removed from graph\n";
			$self->{'all_edge_graph'}->delete_vertex($vertex);
			push (@{$self->{'removed_nodes'}}, $vertex);
		}
	}
}

sub _remove_noisy_edges_alternative{
	my $self=shift;
	my $logfh=$self->{'log_file'};
	my $deviate_factor=$self->{'deviate_factor'};
	my $mean=0; my $var=0; my $total=0;
	my $vertex_effective_ps;
	foreach my $vertex ($self->{'all_edge_graph'}->vertices){
		foreach my $p ($self->{'all_edge_graph'}->predecessors($vertex) ){
			my @edge=($p, $vertex);
			my $wt=$self->{'all_edge_graph'}->get_edge_weight(@edge);
			my $gap_size=$self->get_gap_size_for_a_edge(\@edge);
			if ($self->get_weight_element_at($wt, 2) >= $self->{'min_links'} and $gap_size > (0 - $deviate_factor * $self->{'sd_lib_size'})   ) {
				$total++;
				push (@{$vertex_effective_ps->{$vertex}->{p}} , $p );
			}
		}
		($mean, $var)=$self->renew_stat( scalar @{$vertex_effective_ps->{$vertex}->{p}}, $mean, $var, $total) if exists $vertex_effective_ps->{$vertex}->{p} ;
		foreach my $s ($self->{'all_edge_graph'}->successors($vertex) ){
			my @edge=($vertex, $s);
			my $wt=$self->{'all_edge_graph'}->get_edge_weight(@edge);
			my $gap_size=$self->get_gap_size_for_a_edge(\@edge);
            if ($self->get_weight_element_at($wt, 2) >= $self->{'min_links'} and $gap_size > (0 - $deviate_factor * $self->{'sd_lib_size'})  ) {
                $total++;
                push (@{$vertex_effective_ps->{$vertex}->{s}} , $s );
            }
        }
		($mean, $var)=$self->renew_stat(scalar @{$vertex_effective_ps->{$vertex}->{s}}, $mean, $var, $total) if exists $vertex_effective_ps->{$vertex}->{s} ;
	}
	my $std_dev=sqrt ($var/($total-1) ) ;
	print $logfh "mean edges with min $self->{'min_links'}  links on each side of contig is $mean and standard dev is $std_dev\n";
	my $edges_to_delete;
	foreach my $vertex ($self->{'all_edge_graph'}->vertices){
		if (exists $vertex_effective_ps->{$vertex}->{p} and @{$vertex_effective_ps->{$vertex}->{p}} > $mean + $self->{'excessive_bridges_dev_factor'} * $std_dev){
			print $logfh "predecessors edges associated with $vertex need to be removed\n";
			foreach my $p (@{$vertex_effective_ps->{$vertex}->{p}}){
				print $logfh "remove $p->$vertex\n";
				my $edge=join ('-', ($p, $vertex) ) ;
                $edges_to_delete ->{$edge} =1 ;
			}
		}
		if (exists $vertex_effective_ps->{$vertex}->{s} and @{$vertex_effective_ps->{$vertex}->{s}} > $mean + $self->{'excessive_bridges_dev_factor'} * $std_dev){
            print $logfh "successor edges associated with $vertex need to be removed\n";
            foreach my $s (@{$vertex_effective_ps->{$vertex}->{s}}){
                print $logfh "remove $vertex->$s\n";
				my $edge=join ('-',( $vertex , $s) );
                $edges_to_delete ->{$edge} =1 ;
            }
        }
	}
	foreach my $edge(keys %$edges_to_delete){
		my @edge=split(/-/, $edge);
		my $weight=$self->{'all_edge_graph'}->get_edge_weight(@edge);
		$self->{'all_edge_graph'}->delete_edge(@edge) ;
		$self->{'stat'}->{'Matepair4.1 N_assoicated_removed_edges'} += $self->get_weight_element_at($weight, 2);
	}
}


sub _remove_noisy_edges{
    my $self=shift;
	print "go to remove noisy edge\n";
    my $logfh=$self->{'log_file'};
    my $mean=0; my $var=0; my $total=0;
    foreach my $vertex ($self->{'all_edge_graph'}->vertices){
        my $edge_p=$self->{'all_edge_graph'}->predecessors($vertex);
        $total+=$edge_p;
        ($mean, $var)=$self->renew_stat($edge_p, $mean, $var, $total) if $edge_p>0;
        my $edge_s= $self->{'all_edge_graph'}->successors($vertex);
        $total+= $edge_s;
        ($mean, $var)=$self->renew_stat($edge_s, $mean, $var, $total) if $edge_s>0;
		print "current mean $mean and var $var\n";
    }
    my $std_dev=sqrt ($var/($total-1) ) ;
    print $logfh "mean edges on each side of contig is $mean and standard dev is $std_dev\n";
	print "mean edges on each side of contig is $mean and standard dev is $std_dev\n";

    my $edges_to_delete;
    foreach my $vertex ($self->{'all_edge_graph'}->vertices){
        my $edge_p=$self->{'all_edge_graph'}->predecessors($vertex);
        if ($edge_p > $mean + $self->{'excessive_bridges_dev_factor'} * $std_dev){
            print $logfh "predecessors edges associated with $vertex need to be removed\n";
            foreach my $p ($self->{'all_edge_graph'}->predecessors($vertex)){
                print $logfh "remove $p->$vertex\n";
                my $edge=join ('-', ($p, $vertex) ) ;
                $edges_to_delete ->{$edge} =1 ;
            }
        }
        my $edge_s= $self->{'all_edge_graph'}->successors($vertex);
        if ($edge_s > $mean + $self->{'excessive_bridges_dev_factor'} * $std_dev){            
			print $logfh "successors edges associated with $vertex need to be removed\n";
            foreach my $s ($self->{'all_edge_graph'}->successors($vertex)){
                print $logfh "remove $vertex->$s\n";
                my $edge=join ('-',( $vertex , $s) );
                $edges_to_delete ->{$edge} =1 ;
            }
        }
    }
    foreach my $edge(keys %$edges_to_delete){
        my @edge=split(/-/, $edge);
        my $weight=$self->{'all_edge_graph'}->get_edge_weight(@edge);
        $self->{'all_edge_graph'}->delete_edge(@edge) ;
        $self->{'stat'}->{'Matepair4.1 N_assoicated_removed_edges'} += $self->get_weight_element_at($weight, 2);
    }
}



sub _Build_superscaffold_graph{
    my $self=shift;
    my $loop=shift;
    my $logfh=$self->{'log_file'};
    $self->{'updated_flag'}=0;
    foreach my $edge (sort {my @a=split(/_/, $self->{'all_edge_graph'}->get_edge_weight(@$a) ); my @b=split(/_/, $self->{'all_edge_graph'}->get_edge_weight(@$b) ) ;  $b[0]<=>$a[0] || $b[1]<=>$a[1] || $b[2]<=>$a[2] || $b[3]<=>$a[3] || $a[4]<=>$b[4] || $b[5]<=>$a[5]  }  $self->{'all_edge_graph'}->edges){
        my @wt=split (/_/, $self->{'all_edge_graph'}->get_edge_weight(@$edge) );
        last if $wt[0] =~ /^-1/; 
        $self->_update_superscaffold_graph($edge, $loop);
    }
}

sub get_weight_element_at{
    my ($self, $w, $at)=@_;
    my @wt=split(/_/, $w);
    die "weigth is $w, not correct format" if (@wt <3 );
    return $wt[$at];
}

sub reformat_contig_look_up{
	my $self=shift; 
	my $contig=shift;
	if (ref($self->{'reformat_map'}) eq 'HASH'){
		return $self->{'reformat_map'}->{$contig};
	}
	return;
}

sub get_gap_size_for_a_edge{
	my ($self, $edge)=@_;
	die "not valid edge of @$edge\n" if ref ($edge) ne 'ARRAY' or  @$edge !=2;
	my $key=join ("=", ($edge->[0], $edge->[1] ) ) ;
	my $gap_size;
	if (exists $self->{'linking_hash'}->{$key }->{'cluster_gap'}){
        $gap_size = $self->{'linking_hash'}->{$key }->{'cluster_gap'};
    }
    elsif (exists $self->{'linking_hash'}->{$key }->{'gap'}){
        my @gap_arr=@{$self->{'linking_hash'}->{$key }->{'gap'}};
        for my $i (0..$#gap_arr){
             $gap_size+=$gap_arr[$i];
        }
        $gap_size/=@gap_arr;
    }
	return $gap_size;
}

sub _update_superscaffold_graph{
	my ($self, $edge , $loop)=@_;
    my $logfh=$self->{'log_file'};
	my $graph_to_build=$self->{'graph_to_build'};
	my $deviate_factor=$self->{'deviate_factor'};
	my $key=join ("=", ($edge->[0], $edge->[1] ) ) ;
    my $case=$self->{'linking_hash'}->{$key}->{'case'} || die "internal error, no orientation case found for a edge ", $edge->[0], " and ", $edge->[1], "\n" ;
    my $f3_contig=$edge->[1];
    my $r3_contig= $edge->[0];
	my $gap_size=$self->get_gap_size_for_a_edge($edge);
    my $print_case;
    if ($case ==1){
        $print_case='++';
    }
    elsif ($case == 3){
        $print_case='-+';    
    }
    elsif ($case == 6){
        $print_case='+-';    
    }
    else{
        die "orientation case is $case, internal error\n";
    }
	my $wt=$self->{'all_edge_graph'}->get_edge_weight(@$edge) ;
	my $total_bundled_mate_pairs=$self->get_weight_element_at($wt, 2);
	my $total_clustered_bundled_mate_pairs=$self->get_weight_element_at($wt, 1);
	if ( defined $self->{'excessive_mate_pairs_for_a_edge_dev_factor'} and $total_clustered_bundled_mate_pairs > $self->{'clustered_mate_pair_count_hash'}->{'mean'}+$self->{'excessive_mate_pairs_for_a_edge_dev_factor'}*$self->{'clustered_mate_pair_count_hash'}->{'standard_deviation'} ){
		print $logfh "&&&&&&&&&&&&&&&&&&&&\n>", join ("->", @$edge), "  weight=", $wt, "  gap=$gap_size orientation=$print_case loop=$loop status=excessive links\n"  ;
		if (ref($self->{'reformat_map'}) eq 'HASH'){
			print $logfh $edge->[0],": ", $self->reformat_contig_look_up($edge->[0]),"\n", $edge->[1],": ", $self->reformat_contig_look_up($edge->[1]),"\n" ;
		}
		print $logfh "Total clustered mate pairs $total_clustered_bundled_mate_pairs is too excessive compared to mean $self->{'clustered_mate_pair_count_hash'}->{'mean'} and stdev $self->{'clustered_mate_pair_count_hash'}->{'standard_deviation'}\n";
		$self->{'all_edge_graph'}->delete_edge(@$edge);
		$self->{'stat'}->{'1 N_of_bridges_remov_due_excessive_matepair'}++;
		$self->{'stat'}->{'Matepair4.1 N_assoicated_removed_edges'} += $self->get_weight_element_at($wt, 2);
		return ;
	}
	my $non_matepair_only_flag=0;
	$non_matepair_only_flag=1 unless exists $self->{'linking_hash'}->{$key }->{'cluster_gap'};	
	unless ( $gap_size > (0 - $deviate_factor * $self->{'sd_lib_size'})  or $non_matepair_only_flag ){
		print $logfh "&&&&&&&&&&&&&&&&&&&&\n>", join ("->", @$edge), "  weight=", $wt, "  gap=$gap_size orientation=$print_case loop=$loop status=excessive insert_size\n"  ;
		if (ref($self->{'reformat_map'}) eq 'HASH'){ 
			print $logfh $edge->[0],": ", $self->reformat_contig_look_up($edge->[0]),"\n", $edge->[1],": ", $self->reformat_contig_look_up($edge->[1]),"\n" ;           
        }
        print $logfh  "gap size $gap_size suggest excessive insert size or redundance sequence in linking two contigs\n" ;
        $self->{'stat'}->{'2 N_of_edges_excessive_inst_size'} ++ ;
		$self->{'stat'}->{'Matepair4.2 N_assoicated_excessive_insertsize_edges'}+=$self->get_weight_element_at($wt, 2);
        $self->{'all_edge_graph'}->delete_edge(@$edge);
        return;
    }
	else{
		$self->{'stat'}->{'3 N_of_bridges_for_EXAMing'}++;
		$self->{'stat'}->{'Matepair5 N_of_matepairs_in_use'}+=$total_bundled_mate_pairs;
		print $logfh "&&&&&&&&&&&&&&&&&&&&\n>", join ("->", @$edge), "  weight=", $wt, "  gap=$gap_size orientation=$print_case loop=$loop\n"  ;
		if (ref($self->{'reformat_map'}) eq 'HASH'){        
			print $logfh $edge->[0],": ", $self->reformat_contig_look_up($edge->[0]),"\n", $edge->[1],": ", $self->reformat_contig_look_up($edge->[1]),"\n" ;                
        }
		unless ($graph_to_build->has_vertex($edge->[0]) or $graph_to_build->has_vertex($edge->[1]) ){
			print  $logfh "new edges, de novo join accepted \n", join ("->", @$edge), " status=de novo accept\n" ;
			$self->{'stat'}->{'Matepair6 N_of_matepairs_Successed'}+=$total_bundled_mate_pairs;
            $self->add_new_to_superscaffold_hash($gap_size, $case, $f3_contig, $r3_contig);
            $graph_to_build->add_edge($r3_contig, $f3_contig);
            $self->{'all_edge_graph'}->delete_edge(@$edge);
            $self->{'stat'}->{'4 N_of_bridge_for_denovo_join'} ++ ;
            $self->{'updated_flag'}=1;
		}
        elsif ($graph_to_build->has_vertex($edge->[0]) and $graph_to_build->has_vertex($edge->[1]) ){
            my $f3_superscaffold=$self->{'scaffold_superscaffold'}->{$f3_contig} || die "no superscaffold hash for f3 $f3_contig\n";
            my $r3_superscaffold=$self->{'scaffold_superscaffold'}->{$r3_contig} || die "no superscaffold hash for r3 $r3_contig\n";
            my $re_initiate=0;
            if ( exists $self->{'rejected_merge_path'}->{join ('+', sort ( $r3_superscaffold, $f3_superscaffold) ) }  and $self->{'used_edges_in_Hamiltonian'}->{$edge}  ){
                print $logfh "Edge used in Hamiltonian already! rejected for the same proposed merge\n", join ("->", @$edge), " status=merge reject\n";
                return 0;
            }
            elsif ( exists $self->{'rejected_merge_path'}->{join ('+', sort ( $r3_superscaffold, $f3_superscaffold) ) } ){
                my $wa=$self->{'rejected_merge_path'}->{join ('+', sort ( $r3_superscaffold, $f3_superscaffold) ) };
                my $wb=$self->{'all_edge_graph'}->get_edge_weight(@$edge);
                if ($self->get_weight_element_at($wa, 1)-$self->get_weight_element_at($wb, 1) >= $self->{'overrule_strength_for_reinitiate'}  * $self->get_weight_element_at($wb,1) ){
                    print $logfh "Edge not used in Hamiltonian! Rejected for the same proposed merge by weight $wa. Re-initiate failed, lack of overrule power\n" , join ("->", @$edge), " status=merge reject\n";
                    return 0;
                }
                else{
                    $re_initiate=1;
                    print $logfh "re-initiate EXAMer for merge!\n";
                }
            }
            my $path_f3 = $self->{'super_scaffold'}->{$f3_superscaffold}->{'main_path'} || die "no path for f3 $f3_contig\n"; 
            my $path_r3 = $self->{'super_scaffold'}->{$r3_superscaffold}->{'main_path'} || die "no path for r3 $r3_contig\n"; 
            my $dir_f3 = $self->{'super_scaffold'}->{$f3_superscaffold}->{'direct'} || die "no direct for f3 $f3_contig\n";
            my $dir_r3 = $self->{'super_scaffold'}->{$r3_superscaffold}->{'direct'} || die "no direct for r3 $r3_contig\n";
            my $gap_f3 = $self->{'super_scaffold'}->{$f3_superscaffold}->{'gap'} ;
            my $gap_r3 = $self->{'super_scaffold'}->{$r3_superscaffold}->{'gap'} ;
            defined $gap_f3 || die "gap is $gap_f3 no gap for f3 $f3_contig with path $path_f3\n";
            defined $gap_r3 || die "gap is $gap_r3 no gap for r3 $r3_contig with path $path_r3\n";
			if ($f3_superscaffold != $r3_superscaffold){
				if ($case == 1){
                    if (&direct_of_a_unit_in_path($f3_contig, $path_f3, $dir_f3) eq '-'){                                                
						($path_f3, $gap_f3, $dir_f3)=&reverse_DNA_string($path_f3, $gap_f3, $dir_f3);
                    }
                    if (&direct_of_a_unit_in_path($r3_contig, $path_r3, $dir_r3) eq '-'){
                        ($path_r3, $gap_r3, $dir_r3)=&reverse_DNA_string($path_r3, $gap_r3, $dir_r3);
                    }
                }
                elsif ($case == 3){
                    if (&direct_of_a_unit_in_path($f3_contig, $path_f3, $dir_f3) eq '-'){
                        ($path_f3, $gap_f3, $dir_f3)=&reverse_DNA_string($path_f3, $gap_f3, $dir_f3);
                    }
                    if (&direct_of_a_unit_in_path($r3_contig, $path_r3, $dir_r3) eq '+'){
                        ($path_r3, $gap_r3, $dir_r3)=&reverse_DNA_string($path_r3, $gap_r3, $dir_r3);
                    }
                }
                elsif ($case == 6){
                    if (&direct_of_a_unit_in_path($f3_contig, $path_f3, $dir_f3) eq '+'){
                        ($path_f3, $gap_f3, $dir_f3)=&reverse_DNA_string($path_f3, $gap_f3, $dir_f3);
                    }
                    if (&direct_of_a_unit_in_path($r3_contig, $path_r3, $dir_r3) eq '-'){
                        ($path_r3, $gap_r3, $dir_r3)=&reverse_DNA_string($path_r3, $gap_r3, $dir_r3);
                    }
                }
				print  $logfh "go to scaffold vs scaffold branch\nscaffold: $path_r3\nscaffold: $path_f3\n" ;

				my ($suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging) = $self -> most_possible_path_if_merging_two_supperscaffold($path_f3, $path_r3, $gap_f3, $gap_r3, $dir_f3, $dir_r3, $f3_contig, $r3_contig, $gap_size );
				if ($suggested_path_if_merging){
					$self->update_superscaffold_hash('merge', $f3_contig, $r3_contig, $suggested_path_if_merging, $suggested_gap_if_merging, $suggested_dir_if_merging );
					$self->update_build_graph_edge($path_f3, $path_r3, $suggested_path_if_merging);
					$self->{'all_edge_graph'}->delete_edge(@$edge);
                    $self->{'updated_flag'}=1;
					print $logfh "New Path:$suggested_path_if_merging\n", join ("->", @$edge), " status=merge accept\n" ;
                    $self->{'stat'}->{'6 N_of_bridges_for_merge'} ++;
					$self->{'stat'}->{'Matepair6 N_of_matepairs_Successed'}+=$total_bundled_mate_pairs;
				}
				else{
					$self->{'rejected_merge_path'}->{join ('+', sort ($r3_superscaffold, $f3_superscaffold) )} = $self->{'all_edge_graph'}->get_edge_weight(@$edge) if $re_initiate == 0;
                    print $logfh  join ("->", @$edge), " status=merge reject\n" ;
                    $self->{'all_edge_graph'}->delete_edge(@$edge) if $re_initiate == 1;
				}
			}
			else{
				print $logfh "two edges are in same scaffold\n" ;
			}
		}
		else { # should exist one and only one vertex now
			my $new_node =  $graph_to_build->has_vertex($f3_contig) ? $r3_contig : $f3_contig;
			my $existed_node = $graph_to_build->has_vertex($f3_contig) ? $f3_contig : $r3_contig;  
            my $existed_node_superscaffold=$self->{'scaffold_superscaffold'}->{$existed_node};
            my $re_initiate=0;
            if ( exists $self->{'rejected_merge_path'}->{join ('+', ( $existed_node_superscaffold, $new_node) ) }  and $self->{'used_edges_in_Hamiltonian'}->{$edge}  ){
                print $logfh "Edge used in Hamiltonian already! rejected for the same proposed extend\n", join ("->", @$edge), " status=extend reject\n";
                return 0;
            }
            if ( exists $self->{'rejected_merge_path'}->{join ('+',  ( $existed_node_superscaffold, $new_node) ) }){
                my $wa=$self->{'rejected_merge_path'}->{join ('+', ($existed_node_superscaffold, $new_node) ) };
                my $wb=$self->{'all_edge_graph'}->get_edge_weight(@$edge);
                if ($self->get_weight_element_at($wa, 1)-$self->get_weight_element_at($wb, 1) > $self->{'overrule_strength_for_reinitiate'}  * $self->get_weight_element_at($wb,1) ){
                    print $logfh "Edge not used in Hamiltonian! Rejected for the same proposed extend by weight $wa. Re-initiate failed, lack of overrule power\n",  join ("->", @$edge), " status=extend reject\n";
                    return 0;
                }
                else{
                    $re_initiate=1;
                    print $logfh "re-initiate EXAMer for extend!\n";
                }
            }
			die "internal error\n" if (!$new_node or !$existed_node or ($new_node eq $existed_node));
			my ($path_f3, $path_r3, $gap_f3, $gap_r3, $dir_f3, $dir_r3);
            my $existed_path;
			if ($existed_node eq $f3_contig){
			 	$path_f3=$self->{'super_scaffold'}->{$self->{'scaffold_superscaffold'}->{$existed_node}  }->{'main_path'};
			 	$gap_f3=$self->{'super_scaffold'}->{$self->{'scaffold_superscaffold'}->{$existed_node}  }->{'gap'};			 
			 	$dir_f3=$self->{'super_scaffold'}->{$self->{'scaffold_superscaffold'}->{$existed_node}  }->{'direct'};
				$path_r3=$new_node;
                $existed_path=$path_f3;
			}
			else{
				$path_r3=$self->{'super_scaffold'}->{$self->{'scaffold_superscaffold'}->{$existed_node}  }->{'main_path'};
                $gap_r3=$self->{'super_scaffold'}->{$self->{'scaffold_superscaffold'}->{$existed_node}  }->{'gap'};                               
                $dir_r3=$self->{'super_scaffold'}->{$self->{'scaffold_superscaffold'}->{$existed_node}  }->{'direct'};
				$path_f3=$new_node;
                $existed_path=$path_r3;
			}
            print  $logfh "go to contig vs scaffold branch\ncontig:$new_node\nscaffold:$existed_path\n" ;

			#my $suggested_path_if_merging ;
	        #my $suggested_dir_if_merging;
	        #my $suggested_gap_if_merging;

	        if ($case == 1){
				if (length($dir_f3) ){
					if (&direct_of_a_unit_in_path($f3_contig, $path_f3, $dir_f3) eq '-'){
						($path_f3, $gap_f3, $dir_f3)=&reverse_DNA_string($path_f3, $gap_f3, $dir_f3);
                    }
				}
				if (length($dir_r3)) {
                    if (&direct_of_a_unit_in_path($r3_contig, $path_r3, $dir_r3) eq '-'){
						($path_r3, $gap_r3, $dir_r3)=&reverse_DNA_string($path_r3, $gap_r3, $dir_r3);
                    }
				}
				$dir_f3='+' unless $dir_f3;
				$dir_r3='+' unless $dir_r3;
	        }
	        elsif ($case == 3){
				if (length($dir_f3)){
                    if (&direct_of_a_unit_in_path($f3_contig, $path_f3, $dir_f3) eq '-'){
						($path_f3, $gap_f3, $dir_f3)=&reverse_DNA_string($path_f3, $gap_f3, $dir_f3);
                    }
                }
                if (length($dir_r3)){
                    if (&direct_of_a_unit_in_path($r3_contig, $path_r3, $dir_r3) eq '+'){
						($path_r3, $gap_r3, $dir_r3)=&reverse_DNA_string($path_r3, $gap_r3, $dir_r3);
                    }
                }
				$dir_r3='-' unless $dir_r3;
				$dir_f3='+' unless $dir_f3;
	        }
	        elsif ($case == 6){
				if (length($dir_f3)){
                    if (&direct_of_a_unit_in_path($f3_contig, $path_f3, $dir_f3) eq '+'){
						($path_f3, $gap_f3, $dir_f3)=&reverse_DNA_string($path_f3, $gap_f3, $dir_f3);
                    }
                }
                if (length($dir_r3)){
                    if (&direct_of_a_unit_in_path($r3_contig, $path_r3, $dir_r3) eq '-'){
						($path_r3, $gap_r3, $dir_r3)=&reverse_DNA_string($path_r3, $gap_r3, $dir_r3);
                    }
                }
				$dir_f3='-' unless $dir_f3;
				$dir_r3='+' unless $dir_r3;
	        }
			my ($suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging) = $self -> most_possible_path_if_merging_two_supperscaffold($path_f3, $path_r3, $gap_f3, $gap_r3, $dir_f3, $dir_r3, $f3_contig, $r3_contig, $gap_size );
            if ($suggested_path_if_merging){
				$self->{'all_edge_graph'}->delete_edge(@$edge);
				print $logfh "New Path:$suggested_path_if_merging\n", join ("->", @$edge), " status=extend accept\n";
				$self->{'stat'}->{'Matepair6 N_of_matepairs_Successed'}+=$total_bundled_mate_pairs;
                $self->update_superscaffold_hash('extend', $f3_contig, $r3_contig, $suggested_path_if_merging, $suggested_gap_if_merging, $suggested_dir_if_merging );
                $self->update_build_graph_edge($path_f3, $path_r3, $suggested_path_if_merging);
                $self->{'stat'}->{'5 N_of_bridges_for_extend'}++;
                $self->{'updated_flag'}=1;
			}
			else{
                $self->{'rejected_merge_path'}->{join ('+',  ( $existed_node_superscaffold, $new_node) ) }=$self->{'all_edge_graph'}->get_edge_weight(@$edge) if $re_initiate==0;
                $self->{'all_edge_graph'}->delete_edge(@$edge) if $re_initiate==1;
				print $logfh  "Extending rejected\n",  join ("->", @$edge), " status=extend reject\n";
			}
		}
	}
}

sub update_build_graph_edge{
        my ($self, $path_f3, $path_r3, $suggested_path_if_merging )=@_;
        my $graph_to_build=$self->{'graph_to_build'};
        my @arr_f3=split(/-/, $path_f3); my @arr_r3=split(/-/, $path_r3);
        foreach my $i (0..$#arr_f3){
			if ($i >=1){
					$graph_to_build->delete_edge($arr_f3[$i-1], $arr_f3[$i]);
			}
        }
        foreach my $i (0..$#arr_r3){
			if ($i >=1){
					$graph_to_build->delete_edge($arr_r3[$i-1], $arr_r3[$i]);
			}
        }
        my @suggested_path_if_merging=split(/-/,$suggested_path_if_merging );
        foreach my $i (0..$#suggested_path_if_merging){
			if ($i >=1){
					$graph_to_build->add_edge($suggested_path_if_merging[$i-1], $suggested_path_if_merging[$i]);
			}
        }
}

sub add_new_to_superscaffold_hash{
        my ($self, $gap_size, $case, $f3_contig, $r3_contig )=@_;
        my $current_superscaffold=$self->{'indx'}+1;
        $self->{'indx'}++;
        my $gap=special_format($gap_size);
        $self->{'super_scaffold'}->{$current_superscaffold}->{'direct'}='++';
        $self->{'super_scaffold'}->{$current_superscaffold}->{'gap'}=$gap;
        $self->{'super_scaffold'}->{$current_superscaffold}->{'main_path'}=$r3_contig.'-'.$f3_contig;
        $self->{'scaffold_superscaffold'}->{$f3_contig}=$current_superscaffold;
        $self->{'scaffold_superscaffold'}->{$r3_contig}=$current_superscaffold;
        if ($case == 1){
        }
        elsif ($case == 3){
                $self->{'super_scaffold'}->{$current_superscaffold}->{'direct'}='-+';
        }
        elsif ($case == 6){
                $self->{'super_scaffold'}->{$current_superscaffold}->{'direct'}='+-';
        }
        else{
            die "no else\n"
        }
}

sub update_superscaffold_hash{
        my ($self, $what, $f3_contig, $r3_contig , $suggested_path_if_merging, $suggested_gap_if_merging, $suggested_dir_if_merging)=@_;
        my $current_superscaffold=$self->{'indx'}+1;
        $self->{'indx'}++;
        if ($what eq 'extend' ){
            my $graph_to_build=$self->{'graph_to_build'};
            my $existed_node = $graph_to_build->has_vertex($f3_contig) ? $f3_contig : $r3_contig;
            my $new_node= $graph_to_build->has_vertex($f3_contig) ? $r3_contig : $f3_contig;
		    my $old_superscaffold_idx=$self->{'scaffold_superscaffold'}->{$existed_node};
            $self->{'scaffold_superscaffold'}->{$new_node}=$current_superscaffold;
		    foreach my $contig (split(/-/, $self->{'super_scaffold'}->{$old_superscaffold_idx}->{'main_path'}) ){
                $self->{'scaffold_superscaffold'}->{$contig}=$current_superscaffold;
            }
		    delete $self->{'super_scaffold'}->{$old_superscaffold_idx};
        }
        elsif ($what eq 'merge'){
            my $f3_superscaffold=$self->{'scaffold_superscaffold'}->{$f3_contig};
            my $r3_superscaffold=$self->{'scaffold_superscaffold'}->{$r3_contig};
            foreach my $contig (split(/-/, $self->{'super_scaffold'}->{$f3_superscaffold}->{'main_path'}) ){
                    $self->{'scaffold_superscaffold'}->{$contig}=$current_superscaffold;
            }
            foreach my $contig (split(/-/, $self->{'super_scaffold'}->{$r3_superscaffold}->{'main_path'}) ){
                    $self->{'scaffold_superscaffold'}->{$contig}=$current_superscaffold;
            }
		    delete $self->{'super_scaffold'}->{$f3_superscaffold};
            delete $self->{'super_scaffold'}->{$r3_superscaffold};
        }
        if ( ($what eq 'extend')  or ($what eq 'merge') ){
            $self->{'super_scaffold'}->{$current_superscaffold}->{'main_path'}=$suggested_path_if_merging;
            $self->{'super_scaffold'}->{$current_superscaffold}->{'gap'}=$suggested_gap_if_merging ;
            $self->{'super_scaffold'}->{$current_superscaffold}->{'direct'}=$suggested_dir_if_merging;
        }
}



sub most_possible_path_if_merging_two_supperscaffold{
	my ($self, $path_f3, $path_r3, $gap_f3, $gap_r3, $dir_f3, $dir_r3, $f3_contig, $r3_contig, $gap_size )=@_;
    my $logfh=$self->{'log_file'};
	my ($suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging);
	my @arr_f3=split(/-/, $path_f3); my @arr_r3=split(/-/, $path_r3);
	my $end_joining=0;
	if ($arr_f3[0] eq $f3_contig and $arr_r3[-1] eq $r3_contig){
		print $logfh  "Easy end extending: $path_r3 to $path_f3\n";
		$end_joining=1; 
	}
	else{
        print $logfh "Go to Hamiltonian algorithm for judge\n";
		my ($what, $path,$dir,$gap)=$self->is_possible_to_merge_two_paths($path_f3, $path_r3, $gap_f3, $gap_r3, $dir_f3, $dir_r3, $f3_contig, $r3_contig, $gap_size);
		if ($what eq 'END_JOINING' ){  
			$end_joining=1;
            $gap_size=$path; #the $path actually stored the updated new gap size now
            print $logfh "Forced end extending\n";
		}
		elsif ($what eq 'INTERVENE'){
			$suggested_path_if_merging=$path;
			$suggested_dir_if_merging=$dir;
			$suggested_gap_if_merging=$gap;		
			return ($suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging);
		}	
	}
	if ($end_joining){
                $suggested_path_if_merging=$path_r3."-".$path_f3;
                $suggested_dir_if_merging=$dir_r3.$dir_f3;
                $suggested_gap_if_merging=$gap_r3."-".special_format($gap_size)."-".$gap_f3;
                $suggested_gap_if_merging=~ s/^-+//; $suggested_gap_if_merging =~ s/-+$//;
                return ($suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging);
        }
        else{
                return 0;
        }
}


sub is_possible_to_merge_two_paths{
	my ($self, $path_f3, $path_r3, $gap_f3, $gap_r3, $dir_f3, $dir_r3, $f3_contig, $r3_contig, $gap_size )=@_;
    my $logfh=$self->{'log_file'};
	my ($suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging);
	my @arr_f3=split(/-/, $path_f3); my @arr_r3=split(/-/, $path_r3);
	my @dir_f3=split(//, $dir_f3); my @dir_r3=split(//, $dir_r3);
    my @gap_f3=split(/-/, $gap_f3); my @gap_r3=split(/-/, $gap_r3);
	my $Merge_graph=Graph::Directed->new;
	$Merge_graph->add_edge($r3_contig, $f3_contig);
    my @in_house_edge=($r3_contig, $f3_contig);
	my $hash_f3; my $hash_r3;
	my @valid_edges;
	my $dir_hash;
	my $edge_f3r3;
	my $gap_hash;
	my $flag=1;
	
	##############################################################################
	#To improve performance, need to have better algorithm to detect cycle in DAG
	# Name a scalar to remember the restrict point. Any link if pass that point will
	# cause cycle. E.g
	#                 A1 -> A2 -> A3-> A4-> A5-> A6-> A7-> A8-> A9
	#                             V
	#	     B1-> B2-> B3-> B4-> B5-> B6-> B7-> B8-> B9-> B10     
	# A3 -> B5 has a edge, suggesting and link from B5 B6 B7 B8 B9 B10 to A1 A2 A3 are not allowed. 
	# Then let B_to_A_restrict_point=(4, 2). index of B5 and A3.
	#############################################################################
    
	my @r3_to_f3_restrict_point ;
	my @f3_to_r3_restrict_point ;
	foreach my $i (0..$#arr_f3){
		$hash_f3->{$arr_f3[$i]}=$i;
		$Merge_graph->add_edge($arr_f3[$i-1], $arr_f3[$i]) if $i>0;
		$dir_hash->{$arr_f3[$i]}=$dir_f3[$i] ;
		$gap_hash->{$arr_f3[$i-1].'-'.$arr_f3[$i]}=$gap_f3[$i-1] if $i>0;
		foreach my $j (0..$#arr_r3){
			if ($flag){
				$hash_r3->{$arr_r3[$j]}=$j ;
				$Merge_graph->add_edge($arr_r3[$j-1], $arr_r3[$j]) if $j>0 ;
				$dir_hash->{$arr_r3[$j]}=$dir_r3[$j] ;
				$gap_hash->{$arr_r3[$j-1].'-'.$arr_r3[$j]}=$gap_r3[$j-1] if $j>0;
			}
            next if ($arr_r3[$j] eq $r3_contig and $arr_f3[$i] eq $f3_contig);
            my @edge;
            if ($self->{'all_edge_graph'}->has_edge($arr_r3[$j],$arr_f3[$i] ) and $self->{'all_edge_graph'}->has_edge($arr_f3[$i],$arr_r3[$j] ) ){
                die "how can you have link edge $arr_r3[$j] -> $arr_f3[$i] as well as $arr_f3[$i] -> $arr_r3[$j] \n";
            }
		    elsif ($self->{'all_edge_graph'}->has_edge($arr_r3[$j],$arr_f3[$i] )){
			    @edge=($arr_r3[$j],$arr_f3[$i] );
		    }
		    elsif ($self->{'all_edge_graph'}->has_edge($arr_f3[$i],$arr_r3[$j])){
			    @edge=($arr_f3[$i],$arr_r3[$j]);
		    }
            
            if (@edge){
                my $edge=join ("=", (@edge) ) ;
				my $case=$self->{'linking_hash'}->{$edge}->{'case'} ||die "can not get case for $edge\n" ;
				if ( $dir_f3[$hash_f3->{$arr_f3[$i]}] ne $dir_r3[$hash_r3->{$arr_r3[$j]}] and $case==1  ){
                    next;
                }
                if ( $dir_f3[$hash_f3->{$arr_f3[$i]}] eq $dir_r3[$hash_r3->{$arr_r3[$j]}] and ($case==3 or $case == 6) ){
                    next;
                }

                my $gap_size=0;                    
                if (exists $self->{'linking_hash'}->{$edge }->{'cluster_gap'}){
                    $gap_size = $self->{'linking_hash'}->{$edge }->{'cluster_gap'};
                }
                elsif (exists $self->{'linking_hash'}->{$edge }->{'gap'}){
                    $gap_size = $self->{'linking_hash'}->{$edge }->{'gap'};
                }
                
                my $non_matepair_only_flag=0;
                $non_matepair_only_flag=1 unless exists $self->{'linking_hash'}->{$edge }->{'cluster_gap'};
                my $deviate_factor=$self->{'deviate_factor'};
                if ($gap_size <= (0 - $deviate_factor * $self->{'sd_lib_size'})  or $non_matepair_only_flag ){
                        next;
                }
                push (@valid_edges, \@edge);
                $edge_f3r3->{\@edge}->{'f3'}=$arr_f3[$i];
                $edge_f3r3->{\@edge}->{'r3'}=$arr_r3[$j];
			}
		}
        $flag=0;
    }
	@f3_to_r3_restrict_point=($hash_f3->{$f3_contig}, $hash_r3->{$r3_contig} );
	
	my $new_edge_f3_is_on_old_f3_path;
	my $find_end_join_link_flag;
	foreach my $edge(sort {my @a=split(/_/, $self->{'all_edge_graph'}->get_edge_weight(@$a) ); my @b=split(/_/, $self->{'all_edge_graph'}->get_edge_weight(@$b) ) ;  $b[0]<=>$a[0] || $b[1]<=>$a[1] || $b[2]<=>$a[2] || $b[3]<=>$a[3] || $a[4]<=>$b[4] || $b[5]<=>$a[5]  } @valid_edges){
		my $new_link=join ("=", @$edge )  ;
		my $case=$self->{'linking_hash'}->{ $new_link }->{'case'} || die "can not find case for @$edge\n";
        my $newlink_f3_contig=$edge->[1];        
		my @valid_new_link_edge;
        my $subgraph_flag;
		if ($case == 1  or $case==3){
            if (exists $hash_f3->{$newlink_f3_contig}){
                if ($dir_f3[$hash_f3->{$newlink_f3_contig}] eq '+'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'r3'}, $edge_f3r3->{$edge}->{'f3'});
                    $subgraph_flag=1;
                }
                elsif ($dir_f3[$hash_f3->{$newlink_f3_contig}] eq '-'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'f3'}, $edge_f3r3->{$edge}->{'r3'});
                    $subgraph_flag=2;
                }
                else {
                    die "direction not expected mark-1\n" ;
                }
				$new_edge_f3_is_on_old_f3_path=1;
            }
            elsif (exists $hash_r3->{$newlink_f3_contig}){
                if ($dir_r3[$hash_r3->{$newlink_f3_contig}] eq '+'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'f3'}, $edge_f3r3->{$edge}->{'r3'});
					$subgraph_flag=2;
                }
                elsif ($dir_r3[$hash_r3->{$newlink_f3_contig}] eq '-'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'r3'}, $edge_f3r3->{$edge}->{'f3'});
                    $subgraph_flag=1;
                }
                else {
                    die "direction not expected mark-2 \n" ;
                }
				$new_edge_f3_is_on_old_f3_path=0
            }
            else{
                die "$newlink_f3_contig is not found in dir hash 1\n";
            }
		}
		elsif ($case == 6){
            if (exists $hash_f3->{$newlink_f3_contig}){
                if ($dir_f3[$hash_f3->{$newlink_f3_contig}] eq '-'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'r3'}, $edge_f3r3->{$edge}->{'f3'});
                    $subgraph_flag=1;
                }
                elsif ($dir_f3[$hash_f3->{$newlink_f3_contig}] eq '+'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'f3'}, $edge_f3r3->{$edge}->{'r3'});
                    $subgraph_flag=2;
                }
                else {
                    die "direction not expected mark-3 \n" ;
                }
				$new_edge_f3_is_on_old_f3_path=1;
            }
            elsif (exists $hash_r3->{$newlink_f3_contig}){
                if ($dir_r3[$hash_r3->{$newlink_f3_contig}] eq '-'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'f3'}, $edge_f3r3->{$edge}->{'r3'});
                    $subgraph_flag=2;
                }
                elsif ($dir_r3[$hash_r3->{$newlink_f3_contig}] eq '+'  ){
                    @valid_new_link_edge=($edge_f3r3->{$edge}->{'r3'}, $edge_f3r3->{$edge}->{'f3'});
                    $subgraph_flag=1;
                }
                else {
                    die "direction not expected mark-4 \n" ;
                }
				$new_edge_f3_is_on_old_f3_path=0;
            }
            else{
                die "$newlink_f3_contig is not found in dir hash 2\n";
            }
        }

        $self->{'used_edges_in_Hamiltonian'}->{$edge}=1;
	
		my $re_f3=$edge->[1]; my $re_r3=$edge->[0];
		if ($new_edge_f3_is_on_old_f3_path == 0){
			$re_f3=$edge->[0]; $re_r3=$edge->[1];
		}
		if ($subgraph_flag==1 ){
			if (@r3_to_f3_restrict_point and $hash_r3->{$re_r3} >= $r3_to_f3_restrict_point[0] and $hash_f3->{$re_f3} <= $r3_to_f3_restrict_point[1]){
				#print $logfh "@valid_new_link_edge forms cycle\n";
				next;
			}
			else{
				if ($hash_r3->{$re_r3} == @arr_r3-1 and $hash_f3->{$re_f3} == 0){ 
					$find_end_join_link_flag=1;
					$suggested_path_if_merging=$path_r3."-".$path_f3;
					print $logfh "find end extending path during all inter edge checking\n";
					last; 
				}
				else{
					$Merge_graph->add_edge(@valid_new_link_edge);
				}
			}
			if ( $hash_f3->{$re_f3} <= $f3_to_r3_restrict_point[0] and  $hash_r3->{$re_r3} >= $f3_to_r3_restrict_point[1]   ){
				@f3_to_r3_restrict_point=($hash_f3->{$re_f3}, $hash_r3->{$re_r3});
			}
		}	
		else{
			if ($hash_r3->{$re_r3} <= $f3_to_r3_restrict_point[1] and $hash_f3->{$re_f3} >= $f3_to_r3_restrict_point[0] ){
				#print $logfh "@valid_new_link_edge forms cycle\n";
				next;
			}
			else{
				$Merge_graph->add_edge(@valid_new_link_edge);
			}
            if (@r3_to_f3_restrict_point==0 or  $hash_r3->{$re_r3} <= $r3_to_f3_restrict_point[0] and $hash_f3->{$re_f3} >= $r3_to_f3_restrict_point[1]   ){
                @r3_to_f3_restrict_point=($hash_r3->{$re_r3}, $hash_f3->{$re_f3});
            }
		}
    }

	unless ($find_end_join_link_flag){
		my @source=$Merge_graph->source_vertices;
		my @sink=$Merge_graph->sink_vertices;
		if (@source !=1 or @sink !=1 ){
			print $logfh "Hamiltonian graph has not exactly one end node, source node contain", join (' ', @source), " and sink node contain ", join(' ', @sink), "\n";
			my ($endjoin_yes, $newgap)=$self->LAST_CHANCE_END_JOINGING($path_r3, $path_f3, $gap_size, $r3_contig, $f3_contig) ;
			if ($endjoin_yes){
				return ('END_JOINING', $newgap);
			}
			return (0);
		}
		else{
			my @path=$Merge_graph->topological_sort(empty_if_cyclic => 1); 

			if (@path==0  ){
				print $logfh "No topological sort path find\n"  if (@path==0);
				my ($endjoin_yes, $newgap)=$self->LAST_CHANCE_END_JOINGING($path_r3, $path_f3, $gap_size, $r3_contig, $f3_contig) ;
				if ($endjoin_yes){
					return ('END_JOINING', $newgap);
				}
				return (0);
			}
			elsif ($Merge_graph->has_path(@path)){
				$suggested_path_if_merging=join ('-', @path);
				print  $logfh "Found a good Hamiltonian path immediately after topological sort\n";
			}
			else{
				my $predessor_with_longest_path_to;
				my $longest_length_to;
				foreach my $v (@path){
					foreach my $e_from_v ($Merge_graph-> edges_from ($v) ){
						my $sink=$e_from_v->[-1];
						my $temp_l=$longest_length_to->{$v}+1; 
						if ($longest_length_to->{$sink}<=$temp_l){
							$longest_length_to->{$sink}=$temp_l;
							$predessor_with_longest_path_to->{$sink}=$v;
						}
					}
				}
				if ($longest_length_to->{$sink[0]} < @path-1){
					print $logfh "longest path is $longest_length_to->{$path[-1]}, less than ", @path-1," while last node in topological sort is $sink[0]\n";
					my ($endjoin_yes, $newgap)=$self->LAST_CHANCE_END_JOINGING($path_r3, $path_f3, $gap_size, $r3_contig, $f3_contig) ;
					if ($endjoin_yes){
						return ('END_JOINING', $newgap);
					}
					return (0);
				}
				elsif ($longest_length_to->{$sink[0]} > @path-1){
					die "BUG!! longest path is $longest_length_to->{$path[-1]}, and \@path-1 is ", @path-1,"\n";
				}
				my @path_retrieve = $sink[0];
				while ($predessor_with_longest_path_to->{$path_retrieve[0]} ) {
					unshift (@path_retrieve, $predessor_with_longest_path_to->{$path[0]});
				}
				die "path retrieve @path_retrieve is not equal to @path\n" if @path_retrieve != @path;
				$suggested_path_if_merging=join ('-', @path_retrieve);
				print  $logfh "Found a good Hamiltonian path by dg\n";
			}
		}
	}
	my @suggested_path_if_merging=split(/-/, $suggested_path_if_merging);
	die "redundant elements in path of $suggested_path_if_merging\n" if @suggested_path_if_merging != @{unique(\@suggested_path_if_merging)};
	for my $i (0..$#suggested_path_if_merging){
		$suggested_dir_if_merging.=$dir_hash->{$suggested_path_if_merging[$i]};
		if ($i>=1){
			if ( defined ($gap_hash->{$suggested_path_if_merging[$i-1].'-'.$suggested_path_if_merging[$i]}) ){
				$suggested_gap_if_merging.='-'.$gap_hash->{$suggested_path_if_merging[$i-1].'-'.$suggested_path_if_merging[$i]};
			}
			else{
				my @edge1 =($suggested_path_if_merging[$i-1], $suggested_path_if_merging[$i]);
				my @edge2 =($suggested_path_if_merging[$i], $suggested_path_if_merging[$i-1]);
				if ($self->{'all_edge_graph'} -> has_edge(@edge1)){
					my $gap_size=0;
					if (exists $self->{'linking_hash'}->{join ("=", @edge1) }->{'cluster_gap'}){
						$gap_size = $self->{'linking_hash'}->{join ("=", @edge1) }->{'cluster_gap'};
					}
					elsif (exists $self->{'linking_hash'}->{join ("=", @edge1) }->{'gap'}){
						$gap_size = $self->{'linking_hash'}->{join ("=", @edge1) }->{'gap'};
					}
					$suggested_gap_if_merging.='-'.special_format($gap_size) ;
				}
				elsif ($self->{'all_edge_graph'} -> has_edge(@edge2)){
					my $gap_size=0;
					if (exists $self->{'linking_hash'}->{join ("=", @edge2) }->{'cluster_gap'}){
						$gap_size = $self->{'linking_hash'}->{ join ("=", @edge2) }->{'cluster_gap'};
					}
					elsif (exists $self->{'linking_hash'}->{join ("=", @edge2) }->{'gap'}){
						$gap_size = $self->{'linking_hash'}->{join ("=", @edge2) }->{'gap'};
					}
					$suggested_gap_if_merging.='-'.special_format($gap_size) ;
				}
				else{
					die "no edge @edge1 or @edge2 exist, internal error\n";
				}
			}
		}		
	}
	$suggested_gap_if_merging=~ s/^-//;
	die "gap is $suggested_gap_if_merging, Pattern not right\n" if $suggested_gap_if_merging=~ /--/ or $suggested_gap_if_merging=~ /^-/ or $suggested_gap_if_merging=~ /-$/;
	die "something wrong where intervened path is $suggested_path_if_merging, dir is $suggested_dir_if_merging and gap is $suggested_gap_if_merging \n", scalar split(/-/, $suggested_path_if_merging), scalar split(/-/, $suggested_gap_if_merging),scalar split(//,$suggested_dir_if_merging) ,"\n" 
	if ( (scalar split(/-/, $suggested_path_if_merging) != (scalar split(/-/, $suggested_gap_if_merging) +1) ) or 
	( scalar split(/-/, $suggested_path_if_merging) != scalar split(//,$suggested_dir_if_merging) )  );
	return ('INTERVENE',$suggested_path_if_merging, $suggested_dir_if_merging, $suggested_gap_if_merging) ;
}


sub LAST_CHANCE_END_JOINGING{           
        my ($self, $path_r3, $path_f3, $gap_size, $r3_contig, $f3_contig)=@_;
        my $logfh=$self->{'log_file'};
        my @edge=sort ($r3_contig, $f3_contig); 
        my $suggested_path_if_merging;          
        if ( $self->{'force_elongate_min_link'} > 0  and  $self->_total_bridges_for_this_edge(\@edge) >= $self->{'force_elongate_min_link'}){
                print $logfh  "Going to force end join if gap size is satisfied for matepair edge. Definitely join for non mate pair only edge\n";
                my $deviate_factor=$self->{'deviate_factor'};
                $suggested_path_if_merging=$path_r3."-".$path_f3;
                my $gap_size_between_superscaffolds=$gap_size;
                my @arr=split (/-/, $suggested_path_if_merging ) ;
                my ($idx_r3_contig, $idx_f3_contig);
                for my $i (0..$#arr){
                        $idx_r3_contig=$i if $arr[$i] eq $r3_contig;
                        $idx_f3_contig=$i if $arr[$i] eq $f3_contig;
                }
                die "in LAST_CHANCE_END_JOINGING function, idx r3 is bigger than or equal idx f3, not possible! the path is $suggested_path_if_merging and f3 contig is $f3_contig and r3 contig is $r3_contig\n" if $idx_r3_contig >= $idx_f3_contig;
                if ($idx_f3_contig-$idx_r3_contig != 1 ){
                        for my $i ($idx_r3_contig+1..$idx_f3_contig-1){
                                die "no scaffold length for $arr[$i]\n" unless $self->{'contig_length'}->{$arr[$i]};
                                $gap_size_between_superscaffolds-=$self->{'contig_length'}->{$arr[$i]};
                        }
                }
                my $non_matepair_only_flag=0;    
                $non_matepair_only_flag=1 unless exists $self->{'linking_hash'}->{join('=', ($r3_contig, $f3_contig) ) }->{'cluster_gap'};
                print $logfh "Change gap size from $gap_size to $gap_size_between_superscaffolds\n" unless $non_matepair_only_flag;
                if ($gap_size_between_superscaffolds > (0 - $deviate_factor * $self->{'sd_lib_size'})  or $non_matepair_only_flag ){
                        return (1, $gap_size_between_superscaffolds) unless $non_matepair_only_flag;
                        return (1, 0) if $non_matepair_only_flag;
                }
                else{
                        return (0);
                }
        }
        else{
                return (0);
        }
}



sub _bigger_weigth_of_two{
    my ($self, $weight1, $weight2)=@_;
    my @w1=split(/_/, $weight1);
    my @w2=split(/_/, $weight2);
	die "w1 is $weight1 and w2 is $weight2\n" if (@w1 <3 or @w2 <3);
	return -1 if $weight1 eq $weight2;
    return 1 if $w1[0]<$w2[0];
    return 1 if $w1[1]<$w2[1] and $w1[0]==$w2[0];
    return 1 if $w1[2]<$w2[2] and $w1[0]==$w2[0] and $w1[1]==$w2[1];
    return 1 if $w1[3]<$w2[3] and $w1[0]==$w2[0] and $w1[1]==$w2[1] and $w1[2]==$w2[2];
    return 1 if $w1[4]>$w2[4] and $w1[0]==$w2[0] and $w1[1]==$w2[1] and $w1[2]==$w2[2] and $w1[3]==$w2[3];
    return 1 if $w1[5]<$w2[5] and $w1[0]==$w2[0] and $w1[1]==$w2[1] and $w1[2]==$w2[2] and $w1[3]==$w2[3] and $w1[4]==$w2[4];
    return 2;
}


sub _bundle_edges_from_different_sources{
	my ($self, $link_type, $hash)=@_;
    my $logfh=$self->{'log_file'};
	die "no hash or role defined in function\n" if !$hash or !$link_type ;
	foreach my $contig_pair (keys %$hash){
		my ($contig1, $contig2)=split(/=/, $contig_pair);
        #if (grep (/^$contig1$/, @{$self->{'removed_nodes'}}) and grep (/^$contig2$/, @{$self->{'removed_nodes'}}) ){
		if ( exists {map {$_=> 1} @{$self->{'removed_nodes'} } }->{$contig1} or exists {map {$_=>1} @{$self->{'removed_nodes'} } }->{$contig2} ){
            next;
        }
		my $current_contig_pair_hash=$hash->{$contig_pair};
		my ($weight, $linking_hash, $f3_contig)=$self->_get_weight_and_linking_stat_infor($contig_pair, $link_type, $current_contig_pair_hash);
		#print "$contig_pair weight=$weight\n" ;
        my $r3_contig= $contig1 eq $f3_contig ? $contig2:$contig1;
		if ($weight!=0){
			if ($self->{'all_edge_graph'}->has_edge($contig1, $contig2) or $self->{'all_edge_graph'}->has_edge($contig2, $contig1) ){
                die "$contig1, $contig2 has two different orientation edges, not bundles correctly before\n" if $self->{'all_edge_graph'}->has_edge($contig1, $contig2) and $self->{'all_edge_graph'}->has_edge($contig2, $contig1);
                my $exist_edge_f3_contig ;
                my $exist_edge_weight;
                if ($self->{'all_edge_graph'}->has_edge($contig1, $contig2) ){
                    $exist_edge_f3_contig=$contig2;
                    $exist_edge_weight=$self->{'all_edge_graph'}->get_edge_weight($contig1, $contig2);
                }
                else{
                    $exist_edge_f3_contig=$contig1; 
                    $exist_edge_weight=$self->{'all_edge_graph'}->get_edge_weight($contig2, $contig1);
                }
                $self->{'all_edge_graph'}->delete_edge($contig1, $contig2) if $self->{'all_edge_graph'}->has_edge($contig1, $contig2);
                $self->{'all_edge_graph'}->delete_edge($contig2, $contig1) if $self->{'all_edge_graph'}->has_edge($contig2, $contig1);
				if ($link_type eq 'mate_pair'){
					$self->{'all_edge_graph'}->add_weighted_edge($r3_contig, $f3_contig, $weight) ;
                    $self->_update_linking_hash($r3_contig, $f3_contig, $linking_hash);
				}
				else{
                    my $case=$linking_hash->{'case'};
					if ( ($self->{'linking_hash'}->{$contig_pair}->{'case'} != 1 and  $case !=1) or ($self->{'linking_hash'}->{$contig_pair}->{'case'} == 1 and  $case ==1 and $f3_contig eq $exist_edge_f3_contig ) ){
						my $new_weight=$self->_bundle_weight($exist_edge_weight, $weight); #bundle  $weight of existing edge.
						$self->{'all_edge_graph'}->add_weighted_edge($r3_contig, $f3_contig, $new_weight); #update weight	
						push (@{$self->{'linking_hash'}->{$contig_pair}->{'gap'}}, $linking_hash->{'gap'}) if exists $linking_hash->{'gap'};
					}
					else{
                        my $which=$self->_bigger_weigth_of_two($exist_edge_weight, $weight);
                        print $logfh "New non mate pair edge conflict with exist edge orientation and has same weigth. Existing edge $r3_contig->$f3_contig\n" if $which == -1;
						if ($which == 2){
                            print $logfh "Choose current non mate pair weight and orientation case\n";
							$self->{'all_edge_graph'}->add_weighted_edge($r3_contig, $f3_contig, $weight);
							$self->_update_linking_hash($r3_contig, $f3_contig, $linking_hash);
						}
                        else{
                            print $logfh "Keep exist weight and orientation case\n";
                            $self->{'all_edge_graph'}->add_weighted_edge($r3_contig, $f3_contig, $weight);
                        }
					}
				}
			}
			else{
				$self->{'all_edge_graph'}->add_weighted_edge($r3_contig, $f3_contig, $weight) ;
				$self->_update_linking_hash($r3_contig, $f3_contig, $linking_hash);
			}
		}
	}
}

sub _update_linking_hash{
	my ($self, $r3_contig, $f3_contig, $linking_hash)=@_;
    my $contig_pair=join ("=", ($r3_contig, $f3_contig) );
	$self->{'linking_hash'}->{$contig_pair}->{'case'}=$linking_hash->{'case'};
	$self->{'linking_hash'}->{$contig_pair}->{'cluster_gap'}=$linking_hash->{'cluster_gap'} if exists $linking_hash->{'cluster_gap'};
    $self->{'linking_hash'}->{$contig_pair}->{'gap'}=$linking_hash->{'gap'} if exists $linking_hash->{'gap'};
}


sub _bundle_weight{
	my ($self, $w1, $w2)=@_;
	my @arr1=split(/_/, $w1);
	die "weight1 $w1 is not well defined\n" if @arr1<2 or ($arr1[0] != 1 and $arr1[0] != -1) or $arr1[1] !~ /^[\d0\.]+$/;
	my @arr2=split(/_/, $w2);	
	die "weight2 $w2 is not well defined\n" if @arr2<2 or ($arr2[0] != 1 and $arr2[0] != -1) or $arr2[1] !~ /^[\d0\.]+$/;
	die "The weight suggest to bundle two sets of mate pair library, this is not supposed to be processed by this function\n" if @arr1>2 and @arr2>2;
	my @arr;
	$arr[1]=$arr1[1]+$arr2[1];
    if ($arr[1] >= $self->{'min_links'}){
        $arr[0]=1;
    }
    else{
        $arr[0]=-1;
    }
	die "this function is for bundling matepair link to non-matepair link, buth the two weight received is from two mate pair links\n" if @arr1 > 2 and @arr2 > 2;
	if ( @arr1 > 2 ){
		foreach my $i (2..$#arr1){
			$arr[$i]=$arr1[$i];
		}
		$arr[2]=$arr[2]+$arr2[1];
	}
	elsif (@arr2 > 2 ){
        foreach my $i (2..$#arr2){
                $arr[$i]=$arr2[$i];
        }
		$arr[2]=$arr[2]+$arr1[1];
    }
	else{
		$arr[2]=$arr[1];
		$arr[3]=1;
		$arr[4]=1;
		$arr[5]=-2; #edges with this weight are bundling bridges only from non mate library.
	}
	return (join ('_', @arr));
}

sub _get_weight_and_linking_stat_infor{
	my ($self, $contig_pair, $link_type, $hash)  = @_;
	if ($link_type eq 'mate_pair'){
		return $self->_get_weight_and_linking_stat_infor_for_mate_pair($hash);
	}
	else{
		return $self->_get_weight_and_linking_stat_infor_for_non_mate_pair($contig_pair,$hash);
	}
}

sub _get_weight_and_linking_stat_infor_for_non_mate_pair{
	my ($self, $contig_pair, $non_mate_pairs_hash)=@_;
	my $weight=$non_mate_pairs_hash->{'weight'} || die "no weight are defined for non mate library\n"; #The weight for nonmate pairs means this link is comparable to how many mate pairs. 
	$weight="1_".$weight;
	my ($gap,$cases, $f3_on);
	$gap= $non_mate_pairs_hash->{'gap'} if exists $non_mate_pairs_hash->{'gap'};
    my ($contig1, $contig2)=sort split(/=/, $contig_pair);
    $f3_on=$contig2;
    $cases=1;
	if ($non_mate_pairs_hash->{'pair_dir'} eq '++' ){
	}
	elsif($non_mate_pairs_hash->{'pair_dir'} eq '--'){
        $cases=1;
		$f3_on=$contig1;
    }
	elsif ($non_mate_pairs_hash->{'pair_dir'} eq '-+' ){
		$cases = 3;
	}
	elsif ($non_mate_pairs_hash->{'pair_dir'} eq '+-' ){
        $cases = 6;
    }
	else{
		die "not correct pair contig direction defined\n";
	}

	my $linking_hash;
    $linking_hash->{'case'}=$cases;
	$linking_hash->{'gap'}=$gap if $gap;
	return ($weight, $linking_hash, $f3_on);
}

sub _get_weight_and_linking_stat_infor_for_mate_pair{
	my ($self, $mate_pairs_hash)=@_;
	my ($same_orientation_pair,$different_orientation_pair_case1,$different_orientation_pair_case2) = $self->remove_redundant_reads_and_classify_different_orientation_between_two_contigs($mate_pairs_hash);
	my $cases=0;
	$cases=scalar keys %$same_orientation_pair if $same_orientation_pair;
	$cases+=3 if $different_orientation_pair_case1;
	$cases+=6 if $different_orientation_pair_case2;

=pod
	my $tt;
	foreach my $key (keys %$same_orientation_pair){
		$tt+=keys %{$same_orientation_pair->{$key}};
		print "***************\nsame keys is $key\n";
	}
	foreach my $key (keys %$different_orientation_pair_case1){
        $tt+=keys %{$different_orientation_pair_case1->{$key}};
		print "diff1 keys is $key\n";
    }
	foreach my $key (keys %$different_orientation_pair_case2){
        $tt+=keys %{$different_orientation_pair_case2->{$key}};
		print "diff2 keys is $key\n";
    }
	print "Now tt is $tt case is $cases\n*******************\n";

=cut

	my $weight=0; 
	my ($l1, $l2, $r, $c, $d);
	my $f3_contig;
	my $hash;
	$r=1;
	
	my $index;
	if ($cases == 2){ #CASE1 and CASE2 conflict
		my @keys= keys %$same_orientation_pair; die "die case2 \n" if @keys !=2;
		($index,  $r)=$self->_get_dominant_case($same_orientation_pair->{$keys[0]}, $same_orientation_pair->{$keys[1]});	
        return $weight if $index ==0;
		$cases=1;
		if ($index == 1){
			delete $same_orientation_pair->{$keys[1]};
		}
		elsif ($index == 2){
			delete $same_orientation_pair->{$keys[0]};
                }
		else{
			die "no else possible in case2\n";
		}
		die "how possible\n" if keys %$same_orientation_pair !=1;
		
    }
    elsif ($cases == 4){ #CASE[12] CASE3 conflict
        my @key1=keys %$same_orientation_pair; die "die case4 same oritentaion problem\n" if @key1 !=1;
        my @key2=keys %$different_orientation_pair_case1; die "die case4 different orientation problem\n" if @key2 !=1;
        
        ($index,  $r)=$self->_get_dominant_case($same_orientation_pair->{$key1[0]}, $different_orientation_pair_case1->{$key2[0]});
        return $weight if $index ==0;
        if ($index == 1){
            $cases=1;
        }
        elsif ($index == 2){
            $cases=3;
        }
        else{
            die "no else possible in case4\n";
        }
    }
    elsif ($cases == 5){ #CASE1 and CASE2 and CASE3 conflict
        my @keys= keys %$same_orientation_pair; die "die case5 \n" if @keys !=2;
        my @key2=keys %$different_orientation_pair_case1; die "die case5 different orientation problem\n" if @key2 !=1;
        ($index,  $r)=$self->_get_dominant_case($same_orientation_pair->{$keys[0]}, $same_orientation_pair->{$keys[1]}, $different_orientation_pair_case1->{$key2[0]} );
        return $weight if $index ==0;
        if ($index == 1){
            $cases=1;
            delete $same_orientation_pair->{$keys[1]};
        }
        elsif ($index == 2){
            $cases=1;
            delete $same_orientation_pair->{$keys[0]};
        }
        elsif ($index == 3){
            $cases=3;
        }
        else{
            die "no else possible in case5\n";
        }
    }
    elsif ($cases == 7){ #CASE[12] and CASE4 conflict
        my @key1=keys %$same_orientation_pair; die "die case7 same oritentaion problem\n" if @key1 !=1;
        my @key2=keys %$different_orientation_pair_case2; die "die case7 different orientation problem\n" if @key2 !=1;
        ($index, $r)=$self->_get_dominant_case($same_orientation_pair->{$key1[0]}, $different_orientation_pair_case2->{$key2[0]});
        return $weight if $index ==0;
        if ($index == 1){
            $cases=1;
        }
        elsif ($index == 2){
            $cases=6;
        }
        else{
            die "no else possible in case7\n";
        }
    }
    elsif ($cases == 8){ #CASE1 and CASE2 and CASE4 conflict
        my @keys= keys %$same_orientation_pair; die "die case8 \n" if @keys !=2;
        my @key2=keys %$different_orientation_pair_case2; die "die case8 different orientation problem\n" if @key2 !=1;
        ($index,  $r)=$self->_get_dominant_case($same_orientation_pair->{$keys[0]}, $same_orientation_pair->{$keys[1]}, $different_orientation_pair_case2->{$key2[0]} );
        return $weight if $index ==0;
        if ($index == 1){
            $cases=1;
            delete $same_orientation_pair->{$keys[1]};
        }
        elsif ($index == 2){
            $cases=1;
            delete $same_orientation_pair->{$keys[0]};
        }
        elsif ($index == 3){
            $cases=6;
        }
        else{
            die "no else possible in case8\n";
        }

    }
    elsif ($cases == 9){ #CASE3 and CASE4 conflict
        my @key1=keys %$different_orientation_pair_case1; die "die case9 same oritentaion problem\n" if @key1 !=1;
        my @key2=keys %$different_orientation_pair_case2; die "die case9 different orientation problem\n" if @key2 !=1;
        ($index, $r)=$self->_get_dominant_case($different_orientation_pair_case1->{$key1[0]}, $different_orientation_pair_case2->{$key2[0]});
        return $weight if $index ==0;
        if ($index == 1){
            $cases=3;
        }
        elsif ($index == 2){
            $cases=6;
        }
        else{
            die "no else possible in case 9\n";
        }
    }
    elsif ($cases == 10){ #CASE[12] and CASE3 and CASE4 conflict
        my @keys= keys %$same_orientation_pair; die "die case8 \n" if @keys !=1;
        my @key1=keys %$different_orientation_pair_case1; die "die case10 different orientation problem\n" if @key1 !=1;
        my @key2=keys %$different_orientation_pair_case2; die "die case10 different orientation problem\n" if @key2 !=1;
        ($index, $r)=$self->_get_dominant_case($same_orientation_pair->{$keys[0]}, $different_orientation_pair_case1->{$key1[0]}, $different_orientation_pair_case2->{$key2[0]} );
        return $weight if $index ==0;
        if ($index == 1){
            $cases=1;
        }
        elsif ($index == 2){
            $cases=3;
        }
        elsif ($index == 3){
            $cases=6;
        }
        else{
            die "no else possible in case10\n";
        }
    }
    elsif ($cases == 11){ #CASE1 and CASE2 and CASE3 and CASE4 conflict
        my @keys= keys %$same_orientation_pair; die "die case11 \n" if @keys !=2;
        my @key1=keys %$different_orientation_pair_case1; die "die case11 different orientation problem\n" if @key1 !=1;
        my @key2=keys %$different_orientation_pair_case2; die "die case11 different orientation problem\n" if @key2 !=1;
        ($index,  $r)=$self->_get_dominant_case($same_orientation_pair->{$keys[0]}, $same_orientation_pair->{$keys[1]}, $different_orientation_pair_case1->{$key1[0]}, $different_orientation_pair_case2->{$key2[0]} );
        return $weight if $index ==0;
        if ($index == 1){
            $cases=1;
            delete $same_orientation_pair->{$keys[1]};
        }
        elsif ($index == 2){
            $cases=1;
            delete $same_orientation_pair->{$keys[0]};
        }
        elsif ($index == 3){
                $cases=3;
        }
        elsif ($index == 4){
            $cases=6;
        }
        else{
            die "no else possible in case11\n";
        }
    }
    else{
        die "cases is $cases something wrong\n" if $cases != 1  and $cases != 3 and $cases != 6;
    }

	if ($cases == 1 ){ #CASE[12]
		$hash=$same_orientation_pair;
	}
	elsif ($cases == 3 ){#CASE3
		$hash=$different_orientation_pair_case1;
	}
	elsif ($cases == 6 ){#CASE4
		$hash=$different_orientation_pair_case2;
	}
	die "still ambiguous\n" if keys %$hash>1;
	($f3_contig)=keys %$hash;
	############################################
	if (keys %{$hash->{$f3_contig}} < $self->{'treat_as_no_link_if_less_than'}) {
		$self->{'stat'}->{'Matepair3.3 N_of_less_than_min_link'}+=keys %{$hash->{$f3_contig}} ;
		return $weight;
	}
	##########################################
	$l2=scalar keys %{$hash->{$f3_contig}};
	my $gap_of_largest_cluster;
    ($l1, $c, $d, $gap_of_largest_cluster)=$self->new_do_clustering($hash->{$f3_contig}, $cases, $f3_contig ) ;
	#print "max cluster size is $l1, no of links is $l2, ratio is $r, cluster is $c, densentiy is $d, cutoff is $cutoff\n";
	if ($self->{'favor_density'}){
		 $weight=$l1."_".$l2."_".$r."_".$c."_".eval {1-$d/$self->{'distance'}};
	}
	else{
    	$weight=$l1."_".$l2."_".$r."_".$c."_".$d/$self->{'distance'};
	}
    if ($l2  <  $self->{'min_links'} ){
        $weight="-1_".$weight;  #attach the weight a negative value for edge which do not have enough linking
		$self->{'stat'}->{'Matepair3.3 N_of_less_than_min_link'}+=$l2;
    }
    else{
        $weight="1_".$weight;  
		$self->{'stat'}->{'Matepair4 N_of_matepairs_pass_filter'}+=$l2;
    }
    $self->{'edge_count'}++;
	($self->{'clustered_mate_pair_count_hash'}->{'mean'}, $self->{'clustered_mate_pair_count_hash'}->{'var'})=$self->renew_stat($l1, $self->{'clustered_mate_pair_count_hash'}->{'mean'}, $self->{'clustered_mate_pair_count_hash'}->{'var'}, $self->{'edge_count'});
    #my $old_mean=$self->{'clustered_mate_pair_count_hash'}->{'mean'};
    #$self->{'clustered_mate_pair_count_hash'}->{'mean'}+= ($l1 - $old_mean)/$self->{'edge_count'};
    #$self->{'clustered_mate_pair_count_hash'}->{'var'}+= ($l1 - $old_mean) * ($l1- $self->{'clustered_mate_pair_count_hash'}->{'mean'});
	my $linking_hash;
    $linking_hash->{'cluster_gap'}=$gap_of_largest_cluster;
    $linking_hash->{'case'}=$cases;
	$self->check_matepair_stat();
    return ($weight, $linking_hash, $f3_contig);
}

sub renew_stat{
	my ($self, $new_ele, $mean, $var, $total_ele) =@_;
	my $old_mean=$mean;
	$mean+=($new_ele - $old_mean)/$total_ele;
	$var+=($new_ele - $old_mean)*($new_ele-$mean);
	return ($mean, $var);
}

sub check_matepair_stat{
	my $self=shift;
	if ($self->{'stat'}->{'Matepair2 N_of_input_matepairs_aft_remove_redd_reads'} != $self->{'stat'}->{'Matepair3.1 N_of_failed_dominant_rule'} + $self->{'stat'}->{'Matepair3.2 N_of_removed_minorities'} + $self->{'stat'}->{'Matepair3.3 N_of_less_than_min_link'} + $self->{'stat'}->{'Matepair4 N_of_matepairs_pass_filter'} ){
		print STDERR "Mate Pair Accounting failed\n";
	
		foreach my $i (sort keys %{$self->{'stat'}} ){
            print "$i ", $self->{'stat'}->{$i}, "\n" if $i =~ /Matepair/;
       }
	   die "Sorry, code may contain bug\n";
	}
}


sub remove_redundant_reads_and_classify_different_orientation_between_two_contigs{
	my $self=shift;
	my $mate_pairs_hash=shift;
	my $same_orientation_pair;
	my $different_orientation_pair_case1;
	my $different_orientation_pair_case2;
	########################################################################################################################
	#Must pay enough attention to reads mapping orientation, To connect two contigs:       A----------------------B (contig1)      to C---------------------------------D (contig2), there are four possibilities, ABCD CDAB CDBA DCAB
	#CASE1 (ABCD or DCBA): a mate pair of F3 + and  R3 + and a mate pair of  F3 - and R3 - are both suggesting a same linking direction if F3+ at contig2 and F3- at contig1; will record in $same_orientation_pair->{$contig2} 
	#Note, for CASE1, F3+ and R3+ mate pair and F3- and R3- mate pair are not suggesting same linking direction if F3+ and F3- are in same contig.
	#CASE2 (CDAB or BADC): therefore, is similar as CASE1, the only different is F3+ is at contig1
	#CASE3 (DCAB or BACD): F3 + and R3 -  are suggesting same linking direction no matter which contig F3 and R3 end located. The linking direction is different from above case if the same two scaffolds are concerned; I will record in $different_orientation_pair_case1 
	#CASE4 (ABDC or CDBA) F3 - and R3 + are suggesting same linking direction but different from the above cases if for same two  scaffold; I will record in $different_orientation_pair_case2
	#########################################################################################################################
	my $unique_reads_infor_hash;
	$self->{'stat'}->{'Matepair1 N_of_input_matepairs'}+=scalar keys %$mate_pairs_hash;
	foreach my $id (keys %$mate_pairs_hash){
		 $self->{'stat'}->{'Matepair1 effective_clone_cover_span'}+=$self->span_of_a_mate_pair($id);
		my $hash;
		$hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'}=$mate_pairs_hash->{$id}->{'F3'}->{'start'};
		$hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'}=$mate_pairs_hash->{$id}->{'F3'}->{'end'};
		$hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'}=$mate_pairs_hash->{$id}->{'R3'}->{'start'};
		$hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'}=$mate_pairs_hash->{$id}->{'R3'}->{'end'};
		#if ($mate_pairs_hash->{$id}->{'F3'}->{'contig'} =~ /contig00021/ or $mate_pairs_hash->{$id}->{'R3'}->{'contig'} =~ /contig00021/){ 
		#	print "id is $id;", $mate_pairs_hash->{$id}->{'F3'}->{'contig'}, " ", $mate_pairs_hash->{$id}->{'R3'}->{'contig'},"\n";
		#}
		if ($self->{'remove_redundant_reads'} =~ /Solid/i ){
			my $flag=0;
			foreach my $i (keys %$unique_reads_infor_hash){
				if ($hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'}
				and $hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'} 
				and $hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'} 
				and $hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'}
				){
					$flag=1;
					last;
				}
			}
			next if $flag;
			my $key=1;
			$key+=scalar keys %$unique_reads_infor_hash;
			$unique_reads_infor_hash->{$key}=$hash;
		}
		elsif ($self->{'remove_redundant_reads'} eq '454' ){
			my $flag=0;
            foreach my $i (keys %$unique_reads_infor_hash){
                if ( ($hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'}
                and $hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'} )
                and ($hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'}
                or $hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'} == $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'})
                ){
                    $flag=1;
					my $key=1;
            		$key+=scalar keys %$unique_reads_infor_hash;
            		$unique_reads_infor_hash->{$key}=$hash;
					#print "$hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'}  $hash->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'} $hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'}  $hash->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'}  redd to $i:  $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'start'} $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{'end'} $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'start'} $unique_reads_infor_hash->{$i}->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{'end'}\n";
                    last;
                }
            }
            next if $flag;
			my $key=1;
            $key+=scalar keys %$unique_reads_infor_hash;
            $unique_reads_infor_hash->{$key}=$hash;
		}
		elsif (exists $self->{'remove_redundant_reads'}  ){
			#die "remove_redundant_reads option not correctly defined\n";
			#do nothing
		}
		$self->{'stat'}->{'Matepair2 N_of_input_matepairs_aft_remove_redd_reads'}++;
		$self->{'stat'}->{'Matepair2 effective_clone_cover_span'}+=$self->span_of_a_mate_pair($id);
		if ($mate_pairs_hash->{$id}->{'F3'}->{'start'} > 0 and $mate_pairs_hash->{$id}->{'R3'}->{'start'} > 0 ){
			$same_orientation_pair->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{$id}=$hash;
		}
		elsif ($mate_pairs_hash->{$id}->{'F3'}->{'start'} < 0 and $mate_pairs_hash->{$id}->{'R3'}->{'start'} < 0 )  {
            $same_orientation_pair->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{$id}=$self->switch_the_direction_sign_of_hash($hash);
		}
		elsif ($mate_pairs_hash->{$id}->{'F3'}->{'start'} > 0 and $mate_pairs_hash->{$id}->{'R3'}->{'start'} < 0  ) {
			if (exists $different_orientation_pair_case1->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}} ){
                $different_orientation_pair_case1->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{$id}=$self->switch_the_direction_sign_of_hash($hash);
			}
			else{
				$different_orientation_pair_case1->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{$id}=$hash;
			}
		}
		elsif ($mate_pairs_hash->{$id}->{'F3'}->{'start'} < 0 and $mate_pairs_hash->{$id}->{'R3'}->{'start'} > 0  ) {
			if (exists $different_orientation_pair_case2->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}} ){
                $different_orientation_pair_case2->{$mate_pairs_hash->{$id}->{'R3'}->{'contig'}}->{$id}=$self->switch_the_direction_sign_of_hash($hash);
			}
			else{
				$different_orientation_pair_case2->{$mate_pairs_hash->{$id}->{'F3'}->{'contig'}}->{$id}=$hash;
			}
		}
		else{
			die "reads coordinate should not be zero, check whether you did not meet this requirement\n";
		}
	}
	return ($same_orientation_pair,$different_orientation_pair_case1,$different_orientation_pair_case2);
}

sub switch_the_direction_sign_of_hash{
    my ($self, $hash)=@_;
    foreach my $contig (keys %$hash){
        $hash->{$contig}->{'start'}=0-$hash->{$contig}->{'start'};
        $hash->{$contig}->{'end'}=0-$hash->{$contig}->{'end'};
    }
    return $hash;
}




sub new_do_clustering{
	my ($self, $hash, $case, $f3_contig)=@_;
	my ($l1, $c, $d)=(0,0,0);
	my $cluster_obj=Algorithm::ClusterPoints->new (dimension => 2,
                                           radius => $self->{'distance'},
                                           minimum_size => $self->{'min_unit_to_form_cluster'},
                                            );
	my $ix_gap;
	my $ct=0;
	foreach my $read_id (sort {$a cmp $b} keys %$hash){
		my @contigs=sort keys %{$hash->{$read_id}};
		die "internal hash is wrong, mate $read_id has ", scalar @contigs, " contig positions, they are @contigs\n" if @contigs !=2;
        $cluster_obj->add_points(abs($hash->{$read_id}->{$contigs[0]}->{'start'})/2+abs($hash->{$read_id}->{$contigs[0]}->{'end'})/2, abs($hash->{$read_id}->{$contigs[1]}->{'start'})/2+abs($hash->{$read_id}->{$contigs[1]}->{'end'})/2 );
		my $lib=$self->_find_lib_of_a_template($read_id);
		my $gap_size=$self->calculate_gap_size_for_a_mate_pair($read_id, $hash->{$read_id}, $self->{'lib'}->{$lib}->{'ave_size'} , $f3_contig, $case);
		$ix_gap->{$ct}=$gap_size;
		$ct++;
    }
	my @clusters_ix = sort {scalar @$b<=>scalar @$a} $cluster_obj->clusters_ix;
	my @clusters;
	my $cluster_gap;
	$ct=0;
	for my $c (@clusters_ix) {
		my @cluster_coords;
		my $sum_gap_size;
     	for my $index (@$c) {
       		push (@cluster_coords, $cluster_obj->point_coords($index) );
			$sum_gap_size+=$ix_gap->{$index};
		}
		$cluster_gap->{$ct}=$sum_gap_size/@$c;
		push (@clusters, \@cluster_coords);
		$ct++;
	}
	
	$l1=scalar @{$clusters[0]}/2 if @clusters;
	$c=scalar @clusters;
	my $largest_cluster_ix=$self->_break_tie_of_largest_cluster_by_gap_size (\@clusters, $cluster_gap) ;
	$d=$self->calculate_density_of_a_cluster(@{$clusters[$largest_cluster_ix]});
	my $gap=$cluster_gap->{$largest_cluster_ix};		
	return ($l1, $c, $d, $gap);	
}

sub _find_lib_of_a_template{
	my ($self, $id)=@_;
	my $len=$self->{'lib_prefix_len'};
	#$id =~ /^([\w-]{$len})/ || die "template id $id is not expected regexp\n";
	$id =~ /^(\S{$len})/ || die "template id $id is not expected regexp\n";
    return $1;
}

sub span_of_a_mate_pair{
	my ($self, $id)=@_;
	my $lib=$self->_find_lib_of_a_template($id) || die "no lib for template $id\n" ;
	return $self->{'lib'}->{$lib}->{'ave_size'}
}


sub calculate_gap_size_for_a_mate_pair{
	my ($self, $read_id, $hash, $lib_ave_size, $f3_contig , $case)=@_;
	my @edge=keys %$hash;
    $f3_contig  || die "internal error, no f3 contig found for a edge" ;
    my $r3_contig= $edge[0] eq $f3_contig ? $edge[1] : $edge[0];
    
    my $contig_length=$self->{'contig_length'};
    $contig_length->{$r3_contig} || die "can not get scaffold length for $r3_contig,  the f3_contig is $f3_contig\n";
    $contig_length->{$f3_contig} || die "can not get scaffold length for $f3_contig, the r3_contig is $r3_contig\n";

    my $f3_tail=abs($hash->{$f3_contig}->{'start'}) < abs($hash->{$f3_contig}->{'end'}) ? abs($hash->{$f3_contig}->{'start'}) : abs($hash->{$f3_contig}->{'end'});
    my $r3_tail= abs($hash->{$r3_contig}->{'end'}) > abs($hash->{$r3_contig}->{'start'}) ? abs($hash->{$r3_contig}->{'end'}) : abs($hash->{$r3_contig}->{'start'}) ;
    $r3_tail= $contig_length->{$r3_contig} - $r3_tail ;
    
    if ($case == 1){
    }
    elsif ($case == 3){
        $r3_tail= abs( $hash->{$r3_contig}->{'end'} ) < abs( $hash->{$r3_contig}->{'start'} ) ? abs( $hash->{$r3_contig}->{'end'} ) : abs( $hash->{$r3_contig}->{'start'} ) ; 
    }
    elsif ($case == 6){
        $f3_tail= abs($hash->{$f3_contig}->{'end'}) > abs($hash->{$f3_contig}->{'start'}) ? abs($hash->{$f3_contig}->{'end'}) : abs($hash->{$f3_contig}->{'start'});
        $f3_tail=$contig_length->{$f3_contig}- $f3_tail ;
    }
    else{
        die "not possible\n";
    }
    warn "$read_id f3 tail $f3_tail is not positive for edge on $f3_contig case $case, contig length is $contig_length->{$f3_contig}\n" if  $f3_tail < 0;
    warn "$read_id r3 tail $r3_tail is not positive for edge on $r3_contig case $case, contig length is $contig_length->{$r3_contig}\n" if  $r3_tail < 0;
 	my $gap_size=$lib_ave_size-($f3_tail + $r3_tail)- (abs($hash->{$f3_contig}->{'start'} - $hash->{$f3_contig}->{'end'} ) ) - (abs($hash->{$r3_contig}->{'start'} - $hash->{$r3_contig}->{'end'} ) );
	return ($gap_size);
}

sub calculate_density_of_a_cluster{
	my ($self, @coords_arr)=@_;
	my @rearrange_arr;
	die "wrong coords arr given\n" if @coords_arr%2 !=0 or @coords_arr==0;
	while(@coords_arr != 0){
		my @ta=splice(@coords_arr, 0, 2);
		push (@rearrange_arr, \@ta);
	}
	my $dist_sum;
	my $c=0;
	for my $i (0..$#rearrange_arr){
		for my $j (1..$#rearrange_arr){
			next if $j <=$i;
			$c++;
			my $dist=sqrt(  ($rearrange_arr[$i]->[0] - $rearrange_arr[$j]->[0] )**2 + ($rearrange_arr[$i]->[1] - $rearrange_arr[$j]->[1])**2 ) ;
			$dist_sum+=$dist;
		}
	}
	if ($c != 0){
		return ($dist_sum/$c) ;
	}
	else{
		return ($self->{'distance'}) ;
	}
}

sub _break_tie_of_largest_cluster_by_gap_size{
	my ($self, $cluster, $cluster_gap)=@_;
	my @cluster=@$cluster;
	my $cluster_size=scalar @{$cluster[0]};
	my $largest_cluster_ix=0;
	for my $i (1..$#cluster){
		last if scalar @{$cluster[$i]} < scalar @{$cluster[$largest_cluster_ix]};
		$largest_cluster_ix=$i if $cluster_gap->{$i} > $cluster_gap->{$largest_cluster_ix};
	}
	return ($largest_cluster_ix);
}


sub _get_dominant_case{
	my ($self, @hashes)=@_;
	die "did not get enought information for judging dominant\n" if @hashes<2;
	my $dominant_ratio=$self->{'dominant_ratio'};
	my $c=0;
	my $case_cnt;
	my $total_mate_pairs=0;
	foreach my $hash (@hashes){
		$c++;
		$case_cnt->{$c}= scalar keys %$hash;
		$total_mate_pairs+=scalar keys %$hash;
	}
	my $max=max (values %$case_cnt);

	#print "total mate pair is $total_mate_pairs ;  and max is $max\n" ;
	#foreach my $case (keys %$case_cnt){
	#	print "case $case: $case_cnt->{$case}\n";
	#}

	#print "max is $max and values are ", join (" " , values %$case_cnt) ,"\n";
	if (grep (/^$max$/, values %$case_cnt) >1 ){
		#print "max is $max, values of casecnt is ", join (' ',  values %$case_cnt), "\n";
		$self->{'stat'}->{'Matepair3.1 N_of_failed_dominant_rule'}+=$total_mate_pairs;
		return 0;
	}
	my $case;
	my $sum;
	foreach my $i (keys %$case_cnt){
		$sum+=$case_cnt->{$i};
		if ($case_cnt->{$i} == $max ){
			$case=$i;
		}
	}
	if ($max/$sum > $dominant_ratio){
		#print "total mate pair is $total_mate_pairs ; sum is $sum and max is $max\n";
		$self->{'stat'}->{'Matepair3.2 N_of_removed_minorities'}+=($sum-$max) ;
		return ($case,  $max/$sum );
	}
	else{
		$self->{'stat'}->{'Matepair3.1 N_of_failed_dominant_rule'}+=$total_mate_pairs;
		return 0;
	}
}


sub _total_bridges_for_this_edge{        
        my ($self, $edge)=@_;
        my $weight=$self->{'all_edge_graph'}->get_edge_weight(@$edge);
        my @arr=split(/_/, $weight);
        return ($arr[2]);
}

sub _total_bridges_in_biggest_cluster_for_this_edge{
        my ($self, $edge)=@_;
        my $weight=$self->{'all_edge_graph'}->get_edge_weight(@$edge);
        my @arr=split(/_/, $weight);
        return ($arr[1]);
}

sub special_format{
        my $number=shift;
        if ($number < 0){
                my $v=abs ($number); $v='{'.$v.'}';
                return $v;
        }
        else{
                return $number;
        }
} 


sub print_superscaffold{
    my $self=shift;
    my $total_superscaffold=0;
    my $total_scaffold_involved=0;
	 my $logfh=$self->{'log_file'};
    foreach my $i (sort {my $al=$self->{'super_scaffold'}->{$a}->{'main_path'};
                         my $bl=$self->{'super_scaffold'}->{$b}->{'main_path'};
                         split(/-/, $bl) <=>split(/-/, $al) } keys %{$self->{'super_scaffold'}}  ){
                    my $ct=scalar split(//, $self->{'super_scaffold'}->{$i}->{'direct'});
                    print $logfh "superscaffold$i\t$ct\t$self->{'super_scaffold'}->{$i}->{'main_path'}\t$self->{'super_scaffold'}->{$i}->{'direct'}\t$self->{'super_scaffold'}->{$i}->{'gap'} \n";
                    $total_superscaffold++;
                    $total_scaffold_involved+=$ct;
            #}
    }
    print "total superscaffld has $total_superscaffold\n";
    print "total scaffold involved in is $total_scaffold_involved\n";
}

sub get_super_scaffold_hash{
	my $self=shift;
	return $self->{'super_scaffold'};
}

sub get_scaffold_superscaffold_hash{
	my $self=shift;
	return $self->{'scaffold_superscaffold'};
}

sub get_edge_graph{
	my $self=shift;
	return $self->{'all_edge_graph'} ;
}

sub get_built_graph{
	my $self=shift;
	return $self->{'graph_to_build'}
}

sub get_linking_hash{
	my $self=shift;
	return $self->{'linking_hash'};
}

sub get_removed_vertex{
    my $self=shift;
    return $self->{'removed_nodes'};
}

sub direct_of_a_unit_in_path {
        my ($unit, $path, $dir)=@_;
        my @path=split(/-/, $path);
        my @dir=split(//, $dir);
        my ($ind ) =grep $path[$_] eq $unit, 0..$#path;
        die "direction of $unit in $path found is $dir[$ind], not + or -\n" if $dir[$ind] ne '+' and $dir[$ind] ne '-';
        return $dir[$ind];
}
sub reverse_DNA_string{
        my ($path, $gap, $dir) = @_;
        $path=&completement_path($path);
        $gap=&completement_path($gap);
        $dir=join('', reverse(split(//,$dir) ) );
        $dir=~ tr/+-/-+/;
        return ($path, $gap, $dir);
}

sub completement_path{
        my ($path)=@_;
        return join ('-', reverse(split(/-/, $path) ) )   ;
}
sub unique {
        my $array=$_[0];
        my %hsh;
        undef @hsh{@$array};
        my @unique_array = keys %hsh;
        return (\@unique_array);
}


1;


__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Atlas::EXAMer - Perl extension for blah blah blah

=head1 SYNOPSIS

  use EXAMer;
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

