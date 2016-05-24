package Atlas::AGP;
use List::Util qw[sum];
use 5.008008;
use strict;

require Exporter;

our @ISA = qw(Exporter);

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
    $self->{'infile'}=$param{'-infile'}  ;
    print STDERR "infile in AGP  is $param{'-infile'} \n";
    open (IN, $self->{'infile'}) || die "can not open In file of AGP file\n";
    my $last_scaffold;
    my $last_is;
	while (<IN>){
		chomp;
		my @arr=split (/\t/, $_);
		#$arr[0]=~ /(Scaffold\d+)/i || die "pattern not expected in Scaffold agp file\n";
        $arr[0]=~ /(\S+)/;
		my $scaffold_id=$1;
		if ($arr[4] ne 'N' ){
            if (exists $self->{'contigs'}->{$arr[5]}){
                $self->{'redd_contig'}->{$arr[5]}++;
                $arr[5].=".redd".$self->{'redd_contig'}->{$arr[5]}; #some AGP file has one contig present twice or more in AGP, I do not like the feature, but I have to deal with.
            }
			$self->{'contigs'}->{$arr[5]}->{'scaffold'}=$scaffold_id;
			$self->{'contigs'}->{$arr[5]}->{'scaffold_start'}=$arr[1];
            $self->{'contigs'}->{$arr[5]}->{'scaffold_end'}=$arr[2];
			$self->{'contigs'}->{$arr[5]}->{'contig_start'}=$arr[6];
			$self->{'contigs'}->{$arr[5]}->{'contig_end'}=$arr[7];
			$self->{'contigs'}->{$arr[5]}->{'contig_direct'}=$arr[8];
			$self->{'scaffold_length'}->{$1}=$arr[2];
            $self->{'scaffold'}->{$scaffold_id}->{'path'}.="-".$arr[5];
            if ($arr[8] ==0 and $arr[8]  !~ /\+|-/ ){ # for honeybee AGP, the contig are of dir 0, it actually should be +
                $arr[8] = '+'; 
            }
            $self->{'scaffold'}->{$scaffold_id}->{'dir'}.=$arr[8];
            if ($last_is ne 'N' and $last_scaffold eq $scaffold_id){
                $self->{'scaffold'}->{$scaffold_id}->{'gap'}.="-".0;
                $self->{'scaffold'}->{$scaffold_id}->{'gap'}=~ s/^-//;
            }
            $self->{'scaffold'}->{$scaffold_id}->{'path'}=~ s/^-//;
            die "end $arr[2] less than  start $arr[1] for $_\n" if $arr[1]>$arr[2];
            $self->{'span_id'}->{$scaffold_id}->{$arr[1].'-'.$arr[2]}->{'contig'}=$arr[5];
            $self->{'span_id'}->{$scaffold_id}->{$arr[1].'-'.$arr[2]}->{'direct'}=$arr[8];
		}
		elsif ($arr[4] eq 'N' ){
			my $key = scalar keys %{$self->{'scaffold_gap'}->{$scaffold_id}};
			$self->{'scaffold_gap'}->{$scaffold_id}->{$key}->{'start'}=$arr[1];
			$self->{'scaffold_gap'}->{$scaffold_id}->{$key}->{'size'}=$arr[5];
            $self->{'scaffold'}->{$scaffold_id}->{'gap'}.="-".$arr[5];
            $self->{'scaffold'}->{$scaffold_id}->{'gap'}=~ s/^-//;
		}
        $last_scaffold=$scaffold_id;
        $last_is=$arr[4];
	}	
}                       

sub trim_AGP_based_NCBI_list{
	my ($self, $retain_span)=@_;
	my $new_retain_span;
	foreach my $contig (keys %$retain_span){
		#print "dealing with contig $contig\n";
		my $scaffold=$self->get_contig_scaffold($contig);
		my $ct=0;
		my $flag=0;
		my $flag2=0;
		my $last_end;
		my $last_index;
		my @path=split(/-/, $self->{'scaffold'}->{$scaffold}->{'path'});
		my $origal_units=@path;
		my @gap=split(/-/, $self->{'scaffold'}->{$scaffold}->{'gap'});
		my @dir=split(//, $self->{'scaffold'}->{$scaffold}->{'dir'});
		my ($index)=grep $path[$_] eq "$contig", 0..$#path;
		#print "contig $contig in path @path index is $index gap is @gap\n" if $scaffold eq 'Scaffold2150';
		die "can not find $contig in @path\n" if !defined $index;
		my $splicecontig=splice (@path, $index, 1);
		#print "path is @path\n;spliced $splicecontig for scaffold $scaffold when contig is $contig\n";
        my $contig_dir=splice (@dir, $index, 1);

		foreach my $span ( sort {my ($s1, $e1)=split(/-/, $a); my ($s2, $e2)=split(/-/, $b); $s1<=>$s2}   @{$retain_span->{$contig} }){
			$ct++;
			if ($span =~ /0-0/){
				#print "contig is $contig HH gap is @gap\n" if $scaffold eq 'Scaffold2150';
				#print "path is @path\n" if $scaffold eq 'Scaffold2150';
				push (@{$new_retain_span->{$contig}}, $span);
				if ( $self->is_singleton($scaffold) ){
					delete $self->{'scaffold'}->{$scaffold};
					delete $self->{'scaffold_length'}->{$scaffold};
					delete $self->{'span_id'}->{$scaffold};
					delete $self->{'scaffold_gap'}->{$scaffold};
					last;
				}
				else{
					print "scaffold $scaffold is not singleton\n";
					if ($index==0){
						shift @gap;
					}
					elsif ($index == $origal_units-1){
						pop @gap;
					}
					else{
						my ($gap_before_contig, $gap_aft_contig);
						$gap_before_contig=$gap[$index-1]  ;
						$gap_aft_contig=$gap[$index] ;
						my $new_gap=$gap_before_contig+($self->{'contigs'}->{$contig}->{'contig_end'}-$self->{'contigs'}->{$contig}->{'contig_start'}+1) + $gap_aft_contig ;
						splice (@gap, $index-1, 1, $new_gap);
						splice (@gap, $index, 1  );
					}
					print "path is @path and gap is @gap and dir is @dir\n" if $scaffold eq 'Scaffold2150';	
					$self->re_define_scaffold($scaffold, \@path, \@gap, \@dir);
				}
				delete $self->{'contigs'}->{$contig};
				last;
			}
			else{
				my ($ss, $se)=split(/-/, $span);
				if ($ss<=$self->{'contigs'}->{$contig}->{'contig_start'}  and $se >= $self->{'contigs'}->{$contig}->{'contig_end'} ){
					push (@{$new_retain_span->{$contig}}, $self->{'contigs'}->{$contig}->{'contig_start'}.'-'.$self->{'contigs'}->{$contig}->{'contig_end'});
					last;
				}
				elsif ($se<=$self->{'contigs'}->{$contig}->{'contig_start'}){
					#do nothing
				}
				elsif ($ss>=$self->{'contigs'}->{$contig}->{'contig_end'} ){
					#do nothing;
				}
				else{
					$flag=1;
					if ($ss<=$self->{'contigs'}->{$contig}->{'contig_start'} ){
						$ss=$self->{'contigs'}->{$contig}->{'contig_start'};
					}
					if ( $se >= $self->{'contigs'}->{$contig}->{'contig_end'} ){
						$se=$self->{'contigs'}->{$contig}->{'contig_end'};
					}
					push (@{$new_retain_span->{$contig}}, $ss.'-'.$se);
					my $new_id=$contig.'.'.$ct;
					if ($ct == 1){
						splice (@path, $index, 0, $new_id);
						splice (@dir, $index, 0, $contig_dir);
						#print "hehe index is $index and gap is @gap\n";
						if ($index==0){
							$flag2=1;
							#print "begin of scaffold! path is @path gap now is @gap\n";
							#my $gap=shift @gap;
							#$gap+=$ss-$self->{'contigs'}->{$contig}->{'contig_start'};
                            #unshift (@gap, $gap);
						}
						elsif ($index == $origal_units-1){
							my $gap=pop @gap;
							die "no gap here can be found, bug\n" if !defined $gap;
							#print "end of scaffold\n";
							$gap+=$ss-$self->{'contigs'}->{$contig}->{'contig_start'};
							push (@gap, $gap);
						}
						else{
							my $gap=splice (@gap, $index-1,1);
							die "Bug, no gap found here\n" if !defined $gap;
							$gap=$gap+ ($ss-$self->{'contigs'}->{$contig}->{'contig_start'}) ;
							splice (@gap, $index-1, 0 , $gap );
							$flag2=1;
						}
						$last_index=$index;
					}
					else{
						splice (@path, $last_index+1, 0, $new_id);
						#my $gap=splice (@gap, $last_index,1);
						my $gap+=$ss-$last_end-1;
						splice (@gap, $last_index, 0 , $gap);
						splice (@dir, $last_index+1, 0 , $contig_dir);
						$last_index++;
					}
					$self->{'contigs'}->{$new_id}->{'contig_start'}=1;
					$self->{'contigs'}->{$new_id}->{'contig_end'}=$se-$ss+1;
					$last_end=$se;
				}
			}
		}
		if ($flag){
			if (@gap and $flag2){
				my $gap=splice (@gap, $last_index,1);
				if (defined $gap){
            		$gap+=$self->{'contigs'}->{$contig}->{'contig_end'}-$last_end;
            		splice (@gap, $last_index, 0 , $gap);
				}
			}
			$self->re_define_scaffold($scaffold, \@path, \@gap, \@dir);
			delete $self->{'contigs'}->{$contig};
		}
	}
	return $new_retain_span;
} 

sub is_singleton{
	my ($self, $scaffold)=@_;
	return 0  if length($self->{'scaffold'}->{$scaffold}->{'gap'})>0 ;
	return 1;
}


sub re_define_scaffold{
	my ($self, $scaffold, $path, $gap, $dir)=@_;
	my @path=@$path; my @gap=@$gap; my @dir=@$dir;
	$self->{'scaffold'}->{$scaffold}->{'path'}=join('-', @path);
	$self->{'scaffold'}->{$scaffold}->{'gap'}=join('-', @gap) if @gap>0;
	undef $self->{'scaffold'}->{$scaffold}->{'gap'} if @gap==0;
	$self->{'scaffold'}->{$scaffold}->{'dir'}=join('', @dir);
	my $off_set=0;
	my $c=0;
	if ($scaffold eq 'Scaffold16652'){
		#print "path is @path gap is @gap dir is @dir\n";
		#print "gap larthn 0 lala" if @gap > 0;
		#print "gap now is $self->{'scaffold'}->{$scaffold}->{'gap'}\n";
	}

	for my $c (@path){
		my ($old_start, $old_end)=($self->{'contigs'}->{$c}->{'contig_start'}, $self->{'contigs'}->{$c}->{'contig_end'});
		die "not defined old_start or/and old_end for $c\n" if !defined $old_start or !defined $old_end;
		$self->{'contigs'}->{$c}->{'scaffold_start'}=$off_set+1 ;
		$self->{'contigs'}->{$c}->{'scaffold_end'}=$self->{'contigs'}->{$c}->{'scaffold_start'}+($old_end-$old_start);
		$off_set+=($self->{'contigs'}->{$c}->{'scaffold_end'}-$self->{'contigs'}->{$c}->{'scaffold_start'})  ;
		$off_set+=$gap[$c];
		$c++;
	}
	$self->{'scaffold_length'}->{$scaffold} = $off_set;
}

sub get_contig_and_pos_and_direct_given_scaffold_and_pos{
    my ($self, $scaf, $pos, $extension_allowed )=@_;
    my @spans=keys %{$self->{'span_id'}->{$scaf}} ;
    my $split_arr=split_arr(\@spans);
    my ($index, $arrp)=find_index($split_arr, abs($pos));
    print "index is $index and ", join (' ', @$arrp),"\n" if $scaf eq 'Scaffold29250';
    my $span;
    if ($index%2 !=1) {
        if ($arrp->[$index+1]-abs($pos) >=0 and $arrp->[$index+1]-abs($pos) <=$extension_allowed){
            $span=$arrp->[$index+1].'-'.$arrp->[$index+2] ;
        }
        elsif (abs($pos)- $arrp->[$index-1] >=0 and abs($pos)- $arrp->[$index-1] <=$extension_allowed){
            $span=$arrp->[$index-2].'-'.$arrp->[$index-1] ;
        }
        else{
            die "given pos $pos in the scaffold $scaf is out of boundary\nindex is $index; arrp from index -2 to index +1 are -2:",$arrp->[$index-2], "-1:", $arrp->[$index-1],'0:', $arrp->[$index], '+1:', $arrp->[$index+1],' ', " extension allowed is $extension_allowed\n";
        }
    }
    else{
        $span=$arrp->[$index-1].'-'.$arrp->[$index+1] ;
    }
    die "no span $span of $scaf are found for pos $pos where index is $index; arrp from index -2 to index +2 are -2:",$arrp->[$index-2], ", -1:", $arrp->[$index-1],', 0:', $arrp->[$index],', +1:', $arrp->[$index+1],' ', ', +2:', $arrp->[$index+2], "\n" unless $self->{'span_id'}->{$scaf}->{$span};
    my ($span_start, $span_end)=split (/-/, $span);
    #print "span is $span \n" if $scaf eq 'Scaffold29250';
    my $newpos=$pos-$span_start+1;
    if ($newpos<=0){
	    $newpos=1;
    }
    if ($newpos>$span_end-$span_start+1){
	    $newpos=$span_end-$span_start+1;
    }
    return ($self->{'span_id'}->{$scaf}->{$span}->{'contig'}, $newpos, $self->{'span_id'}->{$scaf}->{$span}->{'direct'} );
}
           
sub get_scaffold_and_pos_and_direct_given_contig_and_pos_and_direct{
	#This function is written to deal with wallaby old solid data when the reference is contig fasta where we want to link scaffolds which means should use scaffold fasta reference.
	my ($self, $contig, $pos, $direct)=@_;
	my $scaffold=$self->{'contigs'}->{$contig}->{'scaffold'};
	return (-1,0) if length($scaffold == 0);
	if ($pos >= $self->{'contigs'}->{$contig}->{'contig_start'} and $pos <= $self->{'contigs'}->{$contig}->{'contig_end'} ) {
		my $scaffold_pos= $self->{'contigs'}->{$contig}->{'scaffold_start'}+$pos - $self->{'contigs'}->{$contig}->{'contig_start'}  ;
		my $return_direct;
		if ($self->{'contigs'}->{$contig}->{'contig_direct'} eq '-'){
			$return_direct=$direct eq '+' ? '-': '+';
		}
		else{
			$return_direct=$direct;
		}
		$return_direct='' if $return_direct eq '+';
		return ($scaffold, $scaffold_pos, $return_direct);
	}
	else{
		return (-1, 0)
	}
}

sub _caculate_scaffold_stat{
	my $self=shift;
    my $total_nonsinglton_scaf=0;
    my $total_nonsinglton_cont=0;
	my $total_cont=0;
    if (keys %{$self->{'scaffold'}} == 0){
        print STDERR "no stat available\n";
        return ;
    }
	my $lenstat; my $totallen; 
	my $contig_length=$self->get_contig_length_hash();  
    foreach my $i (keys %{$self->{'scaffold'}} ){
        my $len_plus_gap;
		my $len_only;
        my @paths=split(/-/, $self->{'scaffold'}->{$i}->{'path'});
		$total_cont+=@paths;
        for my $j (@paths){
            $len_plus_gap+=$contig_length->{$j} ;
			$len_only+=$contig_length->{$j} ;
        	my @gaps=split(/-/, $self->{'scaffold'}->{$i}->{'gap'} );
			for my $g (@gaps){
				$len_plus_gap+=int($g);
			}
		}
		if ($len_only == 0){
			print STDERR "scaffold $i is of length 0\n";
			die;
		}
        $lenstat->{withgap}->{int $len_plus_gap}++;
		$lenstat->{nogap}->{int $len_only}++;
        $totallen->{withgap}+=int $len_plus_gap;
		$totallen->{nogap}+=int $len_only;
        if (@paths>1){
            $total_nonsinglton_scaf++;
            $total_nonsinglton_cont+=@paths;        
		}
    }
	foreach my $flag (keys %{$lenstat} ){
		my $switch=1; my $sum;
		foreach my $len (sort {$b <=> $a} keys %{$lenstat->{$flag}}) {
			$self->{'stat'}->{$flag}->{'maxlen'}=$len if (not exists $self->{'stat'}->{$flag}->{'maxlen'} ) or (exists $self->{'stat'}->{$flag}->{'maxlen'} and $len > $self->{'stat'}->{$flag}->{'maxlen'});
			$self->{'stat'}->{$flag}->{'minlen'}=$len if (not exists $self->{'stat'}->{$flag}->{'minlen'} ) or (exists $self->{'stat'}->{$flag}->{'minlen'} and $len < $self->{'stat'}->{$flag}->{'minlen'}); 
			$sum->{$flag} += $len*$lenstat->{$flag}->{$len};
			$self->{'stat'}->{$flag}->{'N50cnt'} += $lenstat->{$flag}->{$len} if ($switch);
			if ($sum->{$flag} >=  $totallen->{$flag}/2 && $switch) {
				$self->{'stat'}->{$flag}->{'N50'} = $len; $switch=0;
			}
		}
		$self->{'stat'}->{$flag}->{'scaffold_count'}=scalar keys %{$self->{'scaffold'}};
        $self->{'stat'}->{$flag}->{'contig_count'}=$total_cont;
		$self->{'stat'}->{$flag}->{'avelen'}=int ($totallen->{$flag}/$self->{'stat'}->{$flag}->{'scaffold_count'} ) ;
		$self->{'stat'}->{$flag}->{'total_len'}=$totallen->{$flag};
    	$self->{'stat'}->{$flag}->{'total_nonsinglton_scaffold'}=$total_nonsinglton_scaf;
    	$self->{'stat'}->{$flag}->{'total_nonsinglton_contig'}=$total_nonsinglton_cont;
	}
}

sub get_scaffold_stat_no_gap{
	my $self=shift;
	my $count_gap=shift;
	$self->_caculate_scaffold_stat() unless exists $self->{'stat'}->{nogap}; #If you call get_scaffold_stat multiple times in a code, I do not want to calculate stat as many times. I just do once.
	return $self->{'stat'}->{nogap} ;
}

sub get_scaffold_stat_with_gap{
	my $self=shift;
	$self->_caculate_scaffold_stat() unless exists $self->{'stat'}->{withgap}; #If you call get_scaffold_stat multiple times in a code, I do not want to calculate stat as many times. I just do once.
	return $self->{'stat'}->{withgap} 
}
	

sub get_scaffold_length_hash{ 
	my $self=shift;
	return $self->{'scaffold_length'};
}

sub get_contig_length_hash{
    my $self=shift;
    my $hash;
    foreach my $key (keys %{$self->{'contigs'}} ){
        $hash->{$key}=$self->{'contigs'}->{$key}->{'contig_end'}-$self->{'contigs'}->{$key}->{'contig_start'}+1;
    }
    return $hash;
}

sub get_scaffold_gap_info{
	my $self=shift;
	return $self->{'scaffold_gap'}
}

sub get_scaffold_hash{
    my $self=shift;
    return $self->{'scaffold'};
}

sub get_contig_scaffold{
    my $self=shift;
    my $contig=shift;
    return $self->{'contigs'}->{$contig}->{'scaffold'};
}

sub get_contig_scaffold_start{
    my $self=shift;
    my $contig=shift;
    return $self->{'contigs'}->{$contig}->{'scaffold_start'};
}

sub get_contig_scaffold_end{
    my $self=shift;
    my $contig=shift;
    return $self->{'contigs'}->{$contig}->{'scaffold_end'};
}

sub get_contig_scaffold_direct{
	 my $self=shift;
     my $contig=shift;
	 return $self->{'contigs'}->{$contig}->{'contig_direct'};
}

sub get_contig_contig_start{
    my $self=shift;
    my $contig=shift;
    return $self->{'contigs'}->{$contig}->{'contig_start'};
}

sub get_contig_contig_end{
    my $self=shift;
    my $contig=shift;
    return $self->{'contigs'}->{$contig}->{'contig_end'};
}

sub get_redd_contig{
    my $self=shift;
    my $hash;
    foreach my $i (keys %{$self->{'redd_contig'}} ){
	$hash->{$i}=1;
    }
    return $hash;
}

sub get_redd_contig_redd_times{
    my $self=shift; my $contig=shift;
    return $self->{'redd_contig'}->{$contig};
}

sub get_singleton_contigs{
	my $self=shift;
	my @contigs;
	foreach my $scaffold_id (keys %{$self->{'scaffold'}}){
		my @path=split (/-/, $self->{'scaffold'}->{$scaffold_id}->{'path'} ) ;
		push (@contigs, @path) if @path==1;
	}
	return \@contigs;
}

sub split_arr{
	my $arr=shift;
	my @new_arr;
	foreach my $i (@$arr){
			push (@new_arr, split('-', $i));
	}
	return \@new_arr;
}


sub find_index{
        my ($newarr, $pos)=@_;
        my @new_arr=sort {$a<=>$b} (@$newarr,$pos) ;
        my ($index) = grep $new_arr[$_] eq $pos, 0..$#new_arr;
        return ($index, \@new_arr);
}       
        
sub get_an_array_of_leaked_contig_length_and_Num_of_total_errors{
	#if scaffold is like contig1-contig2-contig5-contig6-contig7, then the leaked contig is contig3 and contig4, and the leaked length array is one element with value length(contig3)+length(contig4)
	#So, I need a reference genome to get the correct order map of all contigs in AGP. 
	#If the multiple chromosome organism, should use this function in context of one chromosome per time. If circle chromsome (or genome like bacterial, should give the option);
	my ($self, @args)=@_;
	my %param = @args;
	@param{ map { lc $_ } keys %param } = values %param; # lowercase key
	my $contig_order_map_file=$param{contig_order_map_file} || die "no contig_order_map_file given\n";
	my $circular_ref=0;
	$circular_ref=$param{circular} if defined $param{circular};
	#print "contig_order_map_file is $contig_order_map_file and cir is $circular_ref\n";
	#my $total_ref_len=$param{total_ref_len} if defined $param{total_ref_len}; #for cirular genome, it is required;
	#die "cicular reference genome has not total length information, I can not caculate correctly the leaked length for this case\n" if ($circular_ref) and !$total_ref_len;
	my $verbose=1;
	$verbose=$param{verbose} if defined $param{verbose};
	my $max_beared_leaked_contigs_length=0;
	$max_beared_leaked_contigs_length=$param{max_beared_leaked_contigs_length} if defined $param{max_beared_leaked_contigs_length};
	my @leaked_contigs_length;
	open (Map, $contig_order_map_file)||die "can not open contig order map file\n";
	my $c=0;
	my $contig_order;
	my $order_contig;
	my $contig_dir;
	my @mapped_contig;
	my $contigs=$self->get_contig_length_hash;
	while (<Map>){
		chomp;
		die "map file not correct format: $_\n" if !/\S+\s+\d+\s+\d+\s+\d+\s+[+-]$/;
		my @arr=split(/\s+/, $_);
		push (@mapped_contig, $arr[0]);
		if (exists $contigs->{$arr[0]}){ #only check contigs in the AGP obj. If there are contigs mapped but not in this AGP obj, I do not care about it.
			$c++;
			push (@{$contig_order->{$arr[0]}}, $c ) ;
			$order_contig->{$c}=$arr[0];
			push (@{$contig_dir->{$arr[0] }}, $arr[4]);
		}
	}
	@mapped_contig=@{unique (\@mapped_contig) };
	print "if between two contigs in a scaffold, there is a $max_beared_leaked_contigs_length base or less missing, I regarded as no error\n" if $verbose;
	print "total mapped contigs are ", scalar @mapped_contig, "; The largest mapped index is $c; total contigs are ", scalar keys %$contig_order, "\n" if $verbose;
	my $scaffold_feature = $self->get_scaffold_hash();
	my $order_errors =0;
	my $dir_errors =0;
	foreach my $c (sort {my $al=$scaffold_feature->{$a}->{'path'};
						 my $bl=$scaffold_feature->{$b}->{'path'};
						 split(/-/, $bl) <=>split(/-/, $al) }
		keys %$scaffold_feature){
		my @previous_contig_orders;
		print "**************scaffold $c************************\n" if $verbose;
		my @AGP_dir=split(//,$scaffold_feature->{$c}->{'dir'} ) ;
		my $valid_contig=0;
		my $temp_conflict_dir=0;
		foreach my $contig (split (/-/, $scaffold_feature->{$c}->{'path'})){
			my $dir=shift @AGP_dir;
			unless (exists $contig_order->{$contig}){ # If there are contigs in AGP but not mapped, I do not care about it either.
				next;
			}
			print "$contig ", join (' ',  @{$contig_order->{$contig}} ) , "\n" if $verbose;
			if (@previous_contig_orders){
				unless ($self->is_continuous_of_two_contigs_order(\@previous_contig_orders, $contig_order->{$contig}, $circular_ref, $order_contig) ){
					if ($temp_conflict_dir ==0){
						print "Scaffold no dir error\n" if $verbose;
					}
					elsif ($temp_conflict_dir == $valid_contig){
						print "Scaffold map to reverse completement of ref, nothing wrong\n" if $verbose;
					}
					else{
						my $min=$temp_conflict_dir < $valid_contig-$temp_conflict_dir ? $temp_conflict_dir : $valid_contig-$temp_conflict_dir;
						$dir_errors+=$min;
						print "Direction Error: $min\n";
					}
					$valid_contig=0;
					$temp_conflict_dir = 0;

					my $missed_contig_len=$self->calculate_the_shortest_contig_len_in_between(\@previous_contig_orders, $contig_order->{$contig} , $circular_ref, $order_contig ) ;
					if ($missed_contig_len<$max_beared_leaked_contigs_length){
						print "Warn: From @previous_contig_orders to ", join (' ', @{$contig_order->{$contig}} ), " Miss_Length $missed_contig_len Max_allowed $max_beared_leaked_contigs_length\n" if $verbose;
					}
					else{
						print "Order Error: From @previous_contig_orders to ", join (' ', @{$contig_order->{$contig}} ), " Miss_Length $missed_contig_len\n" if $verbose;
						$order_errors++;
					}
					push (@leaked_contigs_length,  $missed_contig_len);
				}
				#else{
				#	
				#}
			}
			$valid_contig++;
            if ( not exists {map{$_=>1} @{$contig_dir->{$contig}} }->{$dir} ){
                $temp_conflict_dir++   ;
                #print "AGP dir is $dir but contig dir is ",@{$contig_dir->{$contig}},"\n";
            }
			@previous_contig_orders=@{$contig_order->{$contig}}
		}
	}
	print "Total order errors $order_errors; Total dir errors $dir_errors\n" if $verbose;
	return (\@leaked_contigs_length,  $dir_errors) ;
}

sub get_errors_by_compare_to_reference{
    #If the multiple chromosome organism, should use this function in context of one chromosome per time. If circle chromsome (or genome like bacterial, should give the option);
    my ($self, @args)=@_;
    my %param = @args;
    @param{ map { lc $_ } keys %param } = values %param; # lowercase key
    my $contig_order_map_file=$param{contig_order_map_file} || die "no contig_order_map_file given\n";
    my $circular_ref=0;
    $circular_ref=$param{circular} if defined $param{circular};
    #print "contig_order_map_file is $contig_order_map_file and cir is $circular_ref\n";
    #my $total_ref_len=$param{total_ref_len} if defined $param{total_ref_len}; #for cirular genome, it is required;
    #die "cicular reference genome has not total length information, I can not caculate correctly the leaked length for this case\n" if ($circular_ref) and !$total_ref_len;
    my $verbose=1;
    $verbose=$param{verbose} if defined $param{verbose};
    my $max_beared_leaked_contigs_length=0;
    $max_beared_leaked_contigs_length=$param{max_beared_leaked_contigs_length} if defined $param{max_beared_leaked_contigs_length};
    my @leaked_contigs_length;
	####################################
	my $num_of_synteny_wrong_contigs = 0; # This is similar to Bambus paper # of errors in Table 1. Bambus paper also count the mapping dir errors by this var. I will put dir errors in another var.
	my $disrupt_by_error_cont=0;
	my @scaffold_who_contain_synteny_wrong_contig;
	my $no_synteny_error_total_len=0;
	my @all_continous_contigs_length;
	my @current__continous_contigs_length;
	my @current__continous_contigs;
	my @overall_length;
	#######################################
    open (Map, $contig_order_map_file)||die "can not open contig order map file\n";
    my $c=0;
    my $contig_order;
    my $order_contig;
    my $contig_dir;
    my @mapped_contig;	

    my $contigs=$self->get_contig_length_hash;
    while (<Map>){
        chomp;
        die "map file not correct format: $_\n" if !/\S+\s+\d+\s+\d+\s+\d+\s+[+-]$/;
        my @arr=split(/\s+/, $_);
        push (@mapped_contig, $arr[0]);
        if (exists $contigs->{$arr[0]}){ #only check contigs in the AGP obj. If there are contigs mapped but not in this AGP obj, I do not care about it.
            $c++;
            push (@{$contig_order->{$arr[0]}}, $c ) ;
            $order_contig->{$c}=$arr[0];
            push (@{$contig_dir->{$arr[0] }}, $arr[4]);
        }
    }
    @mapped_contig=@{unique (\@mapped_contig) };
    print "if between two contigs in a scaffold, there is a $max_beared_leaked_contigs_length base or less missing, I regarded as no error\n" if $verbose;
    print "total mapped contigs are ", scalar @mapped_contig, "; The largest mapped index is $c; total contigs are ", scalar keys %$contig_order, "\n" if $verbose;
    my $scaffold_feature = $self->get_scaffold_hash();

	foreach my $c (sort {my $al=$scaffold_feature->{$a}->{'path'};
                         my $bl=$scaffold_feature->{$b}->{'path'};
                         split(/-/, $bl) <=>split(/-/, $al) }
        keys %$scaffold_feature){
		print "**************scaffold $c************************\n" if $verbose;
		my $synteny_dir; #can be +, which suggest order is such as 1->3 (incremental) or -, which suggest order is such as 3->1 (decremental)		
		my $flag=0;
		my $scaffold_span_circle_genome_start_point=0;
		my $last_order;
		my $valid_contig=0;
		my $temp_synteny_wrong_contigs=0;
		my $valid_len=0;
		foreach my $contig (split (/-/, $scaffold_feature->{$c}->{'path'})){
            unless (exists $contig_order->{$contig}){ # If there are contigs in AGP but not mapped, I do not care about it either.
                next;
            }
			$valid_contig++;
			$valid_len+=$contigs->{$contig};
			print "$contig ", join (' ',  @{$contig_order->{$contig}} ) , "\n" if $verbose;
			if ($circular_ref and defined $last_order){
				$scaffold_span_circle_genome_start_point = $self->is_scaffold_span_circle_genome_start_point($last_order, $contig_order->{$contig}, $order_contig);
			}
			if (defined $synteny_dir){
				unless ($self->is_in_synteny($contig_order->{$contig}, $last_order, $synteny_dir, $scaffold_span_circle_genome_start_point )){
					$temp_synteny_wrong_contigs++;
					$flag=1;
					print "synteny wrong now\n" if $verbose;
				}
			}
			if ($last_order and not defined $synteny_dir){                
				$synteny_dir=$self->get_synteny_dir($last_order, $contig_order->{$contig} , $scaffold_span_circle_genome_start_point, $order_contig) ;
            }
			if ($last_order){
				if ($self->is_continuous_of_two_contigs_order($last_order, $contig_order->{$contig}, $circular_ref, $order_contig) ){
					push (@current__continous_contigs_length, $contigs->{$contig});
					push (@current__continous_contigs, $contig);
				}
				else{
					if (@current__continous_contigs_length == 1){
					}
					elsif (@current__continous_contigs_length > 1 ){
						my $sum=sum (@current__continous_contigs_length);
						push (@all_continous_contigs_length, $sum);
						print "continuous length now add $sum \n" if $verbose ;
						print "CCC @current__continous_contigs\n";
						$disrupt_by_error_cont++;
					}
					else{
						die "not possible to come here\n";
					}
					undef @current__continous_contigs_length;
					undef @current__continous_contigs;
					push (@current__continous_contigs_length, $contigs->{$contig});
					push (@current__continous_contigs, $contig);
				}
			}
			else{
				push (@current__continous_contigs_length, $contigs->{$contig});
				push (@current__continous_contigs, $contig);
			}
			$last_order=$contig_order->{$contig};
		}
		push (@overall_length , $valid_len) if $valid_len;
		if ($flag){
			push (@scaffold_who_contain_synteny_wrong_contig, $c);
			my $alternative_count_errors=$valid_contig - $temp_synteny_wrong_contigs -1; # the first contig need not to consider, so -1;
			$alternative_count_errors=1 if $alternative_count_errors == 0; # you should have at least one error
			my $current_synteny_errors=$temp_synteny_wrong_contigs < $alternative_count_errors ? $temp_synteny_wrong_contigs : $alternative_count_errors;
			$num_of_synteny_wrong_contigs +=$current_synteny_errors;
			print "synteny wrong now add $current_synteny_errors\n" if $verbose;
		}
		if (@current__continous_contigs_length == 1){
		}
		elsif (@current__continous_contigs_length > 1 ){
			my $sum=sum (@current__continous_contigs_length);
			push (@all_continous_contigs_length, $sum);
			print "continuous length now add $sum at end of scaffold\n" if $verbose ;
			print "CCC @current__continous_contigs at end of scaffold\n";
		}
		else{
			#die "not possible to come here at the end of scaffold\n";
		}
		undef @current__continous_contigs_length;
		undef @current__continous_contigs;
	}
	print "synteny error total count is $num_of_synteny_wrong_contigs\n";
    my $discontinue_errors =0;
    my $dir_errors =0;
    foreach my $c (sort {my $al=$scaffold_feature->{$a}->{'path'};
                         my $bl=$scaffold_feature->{$b}->{'path'};
                         split(/-/, $bl) <=>split(/-/, $al) }
        keys %$scaffold_feature){
		next if (exists {map{$_=>1} @scaffold_who_contain_synteny_wrong_contig }->{$c} );
        my @previous_contig_orders;
        print "^^^^^^scaffold $c^^^^^^^^\n" if $verbose;
        my @AGP_dir=split(//,$scaffold_feature->{$c}->{'dir'} ) ;
        my $valid_contig=0;
        my $temp_conflict_dir=0;
        foreach my $contig (split (/-/, $scaffold_feature->{$c}->{'path'})){
			unless (exists $contig_order->{$contig}){ # If there are contigs in AGP but not mapped, I do not care about it either.
				shift @AGP_dir;
                next;
            }
			$no_synteny_error_total_len+=$contigs->{$contig};
            my $dir=shift @AGP_dir;
            unless (exists $contig_order->{$contig}){ # If there are contigs in AGP but not mapped, I do not care about it either.
                next;
            }
            print "$contig ", join (' ',  @{$contig_order->{$contig}} ) , "\n" if $verbose;
            if (@previous_contig_orders){
                unless ($self->is_continuous_of_two_contigs_order(\@previous_contig_orders, $contig_order->{$contig}, $circular_ref, $order_contig) ){
                    if ($temp_conflict_dir ==0){
                        print "Scaffold no dir error\n" if $verbose;
                    }
                    elsif ($temp_conflict_dir == $valid_contig){
                        print "Scaffold map to reverse completement of ref, nothing wrong\n" if $verbose;
                    }
                    else{
                        my $min=$temp_conflict_dir < $valid_contig-$temp_conflict_dir ? $temp_conflict_dir : $valid_contig-$temp_conflict_dir;
                        $dir_errors+=$min;
                        print "Direction Error: $min\n";
                    }
                    $valid_contig=0;
                    $temp_conflict_dir = 0;

                    my $missed_contig_len=$self->calculate_the_shortest_contig_len_in_between(\@previous_contig_orders, $contig_order->{$contig} , $circular_ref, $order_contig ) ;
                    if ($missed_contig_len<$max_beared_leaked_contigs_length){
                        print "Warn: From @previous_contig_orders to ", join (' ', @{$contig_order->{$contig}} ), " Miss_Length $missed_contig_len Max_allowed $max_beared_leaked_contigs_length\n" if $verbose;
                    }
                    else{
                        print "discontinue Error: From @previous_contig_orders to ", join (' ', @{$contig_order->{$contig}} ), " Miss_Length $missed_contig_len\n" if $verbose;
                        $discontinue_errors++;
                    }
                    push (@leaked_contigs_length,  $missed_contig_len);
                }
                #else{
                #   
                #}
            }
            $valid_contig++;
            if ( not exists {map{$_=>1} @{$contig_dir->{$contig}} }->{$dir} ){
                $temp_conflict_dir++   ;
                #print "AGP dir is $dir but contig $contig  dir is ",@{$contig_dir->{$contig}},"\n";
            }
            @previous_contig_orders=@{$contig_order->{$contig}}
        }
    }
    print "Total discontinue errors $discontinue_errors; Total dir errors $dir_errors\n " if $verbose;
    return (\@leaked_contigs_length, \@all_continous_contigs_length, \@overall_length, $num_of_synteny_wrong_contigs, $dir_errors,  $no_synteny_error_total_len, $disrupt_by_error_cont ) ;
}

sub is_in_synteny{
    my ($self, $order, $last_order, $synteny_dir, $span_start_point) = @_;
    unless ($span_start_point){
        foreach my $i (@$last_order){
            foreach my $j (@$order){
                return 1 if $j > $i and $synteny_dir eq '+' ;
                return 1 if $j < $i and $synteny_dir eq '-' ;
            }
        }
        return 0;
    }
    else{
        print STDERR "Your tell me is it in synteny or not? input 'y' or 'n' and enter\n";
        my $input=<STDIN>;
        if ($input =~ /y/i){
         return 1
        }
        elsif ($input =~ /n/i){
            return 0
        }
        else{
            die "you did not give y or n, program exit, rerun the program\n";
        }
    }       
    return 0    
}       

sub is_scaffold_span_circle_genome_start_point{
	my ($self, $last_order, $order , $order_contig)=@_;
	foreach my $i (@$last_order){
    	foreach my $j (@$order){
			return 1 if abs($j - $i ) >= int (keys (%$order_contig) /2 );
		}
	}
	return 0;
}

sub get_synteny_dir{
	my ($self, $last_order, $order , $scaffold_span_circle_genome_start_point, $order_contig)=@_;
	my $flag;
	foreach my $i (@$last_order){
        foreach my $j (@$order){
			if ($flag eq '+' and $j < $i){
				if ($scaffold_span_circle_genome_start_point){
					print STDERR "you input the synteny dir by enter '+' or '-' \n";
					my $input = <STDIN>; chomp $input;
					die "not correct input\n" if $input ne '+' and $input ne '-';
					return $input;
				}
				else{
					die "reference mapping results suggest conflict synteny information between two contigs when at least on of them is mapped in two positions when last $order_contig->{$last_order->[0]} is ", join (' ', @$last_order), " and current $order_contig->{$order->[0]}  is " , join (' ', @$order), "\n";
				}
			}
			elsif ($flag eq '-' and $j > $i){
				if ($scaffold_span_circle_genome_start_point){                    
					print "you input the synteny dir by enter '+' or '-' \n";
                    my $input = <STDIN>;  chomp $input;                  
					die "not correct input\n" if $input ne '+' and $input ne '-';
					return $input;
                }                
				else{
                	die "reference mapping results suggest conflict synteny information between two contigs when at least on of them is mapped in two positions when last $order_contig->{$last_order->[0]} is ", join (' ', @$last_order), " and current $order_contig->{$order->[0]}  is " , join (' ', @$order), "\n";
				}
            }
            $flag ='+' if $j > $i;
			$flag ='-' if $j < $i;
        }
    }
    return $flag;
}


sub get_len_stat{
    my ($self, $arr)=@_;
    my $lenstat;
    my $totallen=0;
    for my $i (@$arr){
        $lenstat->{$i}++;
        $totallen+=$i;
    }
    my $switch=1; my $sum;
    my $hash;
    foreach my $len (sort {$b <=> $a} keys %$lenstat) {
        $hash->{'maxlen'}=$len if (not exists $hash->{'maxlen'} ) or (exists $hash->{'maxlen'} and $len > $hash->{'maxlen'});
        $hash->{'minlen'}=$len if (not exists $hash->{'minlen'} ) or (exists $hash->{'minlen'} and $len < $hash->{'minlen'});
        $sum += $len*$lenstat->{$len};
        $hash->{'N50cnt'} += $lenstat->{$len} if ($switch);
        if ($sum >=  $totallen/2 && $switch) {
            $hash->{'N50'} = $len; $switch=0;
        }
    }
    $hash->{'count'}=scalar @$arr;
	$hash->{'avelen'}=0;
    $hash->{'avelen'}=int ($totallen/@$arr ) if @$arr;
    $hash->{'totallen'}=$totallen;
    return $hash;
}


sub is_continuous_of_two_contigs_order{
    my ($self,$pre_orders, $now_orders, $circular_ref , $order_contig)=@_;
    foreach my $i (@$pre_orders){
        foreach my $j (@$now_orders){
            return 1 if  abs ($j - $i)  ==1 ;
			if ($circular_ref and (($i == keys %$order_contig and $j ==1) or ($j == keys %$order_contig and $i ==1) )   ){
				return 1;
			}
        }
    }
    return 0
}

sub calculate_the_shortest_contig_len_in_between{
    my ($self, $pre_orders, $now_orders, $circular_ref , $order_contig)=@_;
    my $min_length;
	my $contigs=$self->get_contig_length_hash;
    foreach my $i (@$pre_orders){
        foreach my $j (@$now_orders){
            my ($start_order, $end_order)=sort {$a<=>$b} ($i, $j);
            my $length;
            my @all_between_contigs;
            for my $k ($start_order+1..$end_order-1){
                push (@all_between_contigs, $order_contig->{$k});
            }
            foreach my $c (@{unique(\@all_between_contigs)} ){
                $length+=$contigs->{$c};
            }
            $min_length=$length if $min_length> $length or !defined $min_length;
			if ($circular_ref){
				my $length_alt;
            	my @all_between_contigs_alt;
				for my $k ($end_order+1..scalar keys %$order_contig ){
					push (@all_between_contigs_alt, $order_contig->{$k});
				}
				foreach my $k (1..$start_order-1){
					push (@all_between_contigs_alt, $order_contig->{$k});
                }
				foreach my $c (@{unique(\@all_between_contigs_alt)} ){
    	            $length_alt+=$contigs->{$c};
	            }
				$min_length=$length_alt if $min_length> $length_alt or !defined $min_length;
			}
        }
    }
    return $min_length;
}

sub unique {
    my $array=$_[0];
    my %hsh;
    undef @hsh{@$array};
    my @unique_array = keys %hsh;
    return (\@unique_array);
}

sub print_AGP{
    my ($self, $fh)=@_;
    my $scaffold_idx=1;
	my $hash=$self->{'scaffold'};
    foreach my $i (sort  keys %$hash){
        my @dir=split(//, $hash->{$i}->{'dir'});
        my @gap=split(/-/, $hash->{$i}->{'gap'});
        my @main_path=split (/-/, $hash->{$i}->{'path'} );
        die "When try to print AGP, found the  scaffold $i, main path is $hash->{$i}->{'path'}, dir is $hash->{$i}->{'dir'} and gap is $hash->{$i}->{'gap'}, number not match\n" if (@main_path != @dir or @main_path != @gap+1 ) ;
        my $offset=0;
        my $unit_idx=1;
        foreach my $unit (@main_path) {
			my $length=$self->{'contigs'}->{$unit}->{'contig_end'}-$self->{'contigs'}->{$unit}->{'contig_start'}+1;
			die "length is $length for $unit, how possible\n" if $length<=0;
            print $fh "$i\t", $offset+1, "\t", $offset+$length,"\t",$unit_idx,"\t", "W", "\t",$unit,"\t","1\t", $length,"\t",shift @dir,"\n";
            $unit_idx++;
            $offset+=$length;
            if (my $gap_length=shift @gap){
                $gap_length=int $gap_length; 
				if ($gap_length>0){
                	print $fh "$i\t",$offset+1, "\t", $gap_length+$offset, "\t$unit_idx\tN\t$gap_length\tfragment\tyes\n" ;
                	$unit_idx++;
				}
                $offset+=$gap_length;
            }
        }
        $scaffold_idx++;
    }
    close $fh;
}


1

__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Atlas::AGP - Perl extension for blah blah blah

=head1 SYNOPSIS

  use vEXAMer;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for vEXAMer, created by h2xs. It looks like the
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

Jixin Deng, E<lt>jdeng@localdomainE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2011 by Jixin Deng

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut

