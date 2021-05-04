#! /usr/bin/perl -w


package Read_PDB;
use Carp;
use File::Copy;
use strict 'vars';

my $includeUNK = 0;
	# if set, any unknown residues will be converted to UNK and still printed by the various printXXX fns.
	# if not set, these residues will simply be dropped.  Either way, they will be noted in the .trans file
my $convertedHETATM = 0;
    # how many HETATM records are added to the ATOM array due to recognized amino-acid name in MODRES?

# A_ for atom   section
# S_ for seqres section
# M_ for modres section
# C_ for conflict section 
# 
my ($H_caOnly, $H_engineered,
	
    @A_atom_index,    @A_atom_name,    @A_residue_name,  @A_chain,
    @A_residue_index, @A_X,   @A_Y,     @A_Z,
    @A_occupancy,     @A_tempFactor,    @A_theRest,       @A_index,     
    @A_bond_len,      @A_ins,           @A_alternate,   
    @A_was_HETATM,
    
    @S_residue_name,
    @S_chain,
    
    @M_residue_name, @M_chain, @M_chainBegin, @M_residue_index, @M_std_res_name,
	%M_chain_resindex, %M_translations,
    
    @C_pdb_res_name, @C_chain,   @C_pdb_index,     @C_seq_index, 
    @C_database,     @C_accession_number,          @C_conflict,
    @C_seq_res_name,
    
  );

#########################################################################################

sub new {
    my( $that, $pdbfile, $code, $chain )  = @_;
    
    # if $chain is only whitespace, all chains will be loaded.  Otherwise only the specified chain will be loaded.
    
    if( ! -e $pdbfile ) {
      return undef;  
    }
    
    # Initializes lists for constructing data members.
    undef $H_caOnly;
    undef $H_engineered;
    
    undef @A_atom_index;
    undef @A_atom_name;
    undef @A_residue_name;
    undef @A_chain;
    undef @A_residue_index;
    undef @A_X;
    undef @A_Y;
    undef @A_Z;
    undef @A_occupancy;
    undef @A_tempFactor;
    undef @A_theRest;
    undef @A_ins;
    undef @A_alternate;
    undef @A_index;
    undef @A_bond_len;
    undef @A_was_HETATM;
    
    undef @S_residue_name ;  
    undef @S_chain      ;  
    
    undef @M_residue_name;
    undef @M_chain;
    undef @M_chainBegin;
    undef @M_residue_index;
    undef @M_std_res_name;
    undef %M_chain_resindex;
    undef %M_translations;
    
    undef @C_pdb_res_name;
    undef @C_seq_res_name;
    undef @C_chain;
    undef @C_pdb_index;
    undef @C_seq_index;
    undef @C_database;
    undef @C_accession_number;
    undef @C_conflict;
    
    # Parse the PDB file.
    Parse_PDB($pdbfile, $chain);
    
    # Save information.
    my $this;
  
    $this = {
        "ATOM" => {
            "atom_index"    => [@A_atom_index],
            "atom_name"     => [@A_atom_name],
            "residue_name"  => [@A_residue_name],
            "chain"         => [@A_chain],
            "residue_index" => [@A_residue_index],
            "modified_index"=> [@A_index],
            "ins"           => [@A_ins],
            "alternate"     => [@A_alternate],
            "X_coordinate"  => [@A_X],
            "Y_coordinate"  => [@A_Y],
            "Z_coordinate"  => [@A_Z],
            "occupancy"     => [@A_occupancy],
            "tempFactor"    => [@A_tempFactor],
            "theRest"       => [@A_theRest],
            "bond_len"      => [@A_bond_len],
            "was_HETATM"    => [@A_was_HETATM],
            "sizeOf"        => $#A_atom_index,
        },
        
        "SEQRES" => {
            "residue_name"   =>[@S_residue_name],
            "chain"          =>[@S_chain],
            "sizeOf"         =>$#S_residue_name,
        },
        "MODRES" => {
            "residue_name"   =>[@M_residue_name],
            "chain"          =>[@M_chain],
            "chainBegin"     =>[@M_chainBegin],
            "residue_index"  =>[@M_residue_index],
            "std_res_name"   =>[@M_std_res_name],
            "chain_resindex" =>%M_chain_resindex,
            "translations"   =>%M_translations,
            "sizeOf"         =>$#M_residue_index,
        },
        "SEQADV" =>{
            "pdb_res_name"     =>[@C_pdb_res_name],
            "seq_res_name"     =>[@C_seq_res_name],
            "chain"            =>[@C_chain],
            "pdb_index"        =>[@C_pdb_index],
            "data_base"        =>[@C_database],
            "accession_number" =>[@C_accession_number],
            "seq_index"        =>[@C_seq_index],
            "conflict"         =>[@C_conflict],
            "sizeOf"           =>$#C_conflict,
        },
        "NAME"  => $pdbfile,
        "CODE"  => $code,
        "CHAIN" => $chain,
    };

	  bless $this, $that;
}

###################################################################################################################################
#Parsing of PDB file
#
sub Parse_PDB(){
    my( $pdbfile, $chain) = @_;
    my( $END_MODEL, $ModResFlag, $i, $LABEL, $SET_LABEL, $INCLUDE );
    
    my $modifiedIndex = 0;
        # new 9 feb 2010: used to uniquely number residue indices in the face of insertion codes
        # this is an attempt to solve the problem in which insertions appear in the middle of
        # chains with residue indices that are non-sequential and that create duplicates if
        # used unaltered.  The previous logic for dealing with this just added a value to the
        # residue indices seen, but this is not robust against insertions using a residue index
        # that is not adjacent to the previously seen index.
    my $lastRes = "";
    my $lastResNumeric = 0;
    
    $ModResFlag=0; $SET_LABEL=0; $LABEL="@";$INCLUDE=1;
    open(PDB, "<$pdbfile") || carp("Can't open $pdbfile. $!\n");
    my @PDB = <PDB>;
    close(PDB) || print "Can't close $pdbfile. $!<br>\n";
    
    $convertedHETATM = 0;
    # reset the HETATM count
    
    my %firstChainResidue;
    # track the first residue number encountered for each chain; this is useful
    # when applying the MODRES later
    
    # Gets only the specified chain.
    if ($chain =~ /\S/) { #chain==chain
        my @PDB_tmp;
        foreach (@PDB) {
           #if (/^ATOM.{17}$chain.*/ || /^HETATM.{15}$chain.*/ || /^SEQRES..(..).$chain.*/ || /^MODRES.(....).(...).$chain.*/ || /^MODEL.{5}$chain.*/ || /^ENDMDL{5}$chain.*/ || /^SEQADV.(....).(...).$chain.*/ ) {
            if (/^ATOM.{17}$chain.*/ || /^HETATM.{15}$chain.*/ || /^SEQRES..(..).$chain.*/ || /^MODRES.(....).(...).$chain.*/ || /^MODEL.{5}.*/ || /^ENDMDL.{5}.*/ || /^SEQADV.(....).(...).$chain.*/ ) {
				push(@PDB_tmp, $_) 
			}
        }
        @PDB = @PDB_tmp;
    }
    $END_MODEL=0;
    $i=0;
  
  foreach (@PDB) {
	if (/^MDLTYP.+CA ATOMS ONLY/ ) {
		# Note that this model only contains the backbone CA atoms
		# and is probably not appropriate for homology modelling via LOOPP
		$H_caOnly = 1;
	}
	if( /ENGINEERED: YES/ ) {
		$H_engineered = 1;
			# a not completely robust attempt at recognizing structures
			# that are engineered
	}
    if (/^MODEL/){
	  # we don't really do anything here...
    }
    if (/^ENDMDL/){
	  # but here we note the model has ended.  Set a flag such that we read no more ATOM or HETATM records.
	  # We could exit the loopp at this point, but we choose to remain in the loopp in case there are other
	  # records we may be interested in that are out of order.
      $END_MODEL=1;
    }
    #parse modified residue name
    if (/^MODRES.(....).(...).(.).(....)(.).(...)(.*)/){
      print "PARSING: MODRES resName=$2 chain=$3 resIdx=$4 stdName=$6\n";
      push(@M_residue_name,  $2);    # residue name
      push(@M_chain ,        $3);    # chain
      push(@M_residue_index, $4);    # residue index
      push(@M_std_res_name,  $6);    # standard residue name
	  $M_chain_resindex{ "$3_$4_$5" } = $6;
		# fast lookup when modifiy atom records later
        # $5 is the insertion code - you need it to uniquely identify, since residue indices
        # can be duplicated 
	  $M_translations{ $2 } = $6;
		# this lookup is used to supplement the "checkRes" functionality.  In particular we need this
		# to translate non-standard names found in the SEQRES section
      $ModResFlag=1;
		# IMPORTANT: note that we always are referring to the residue index in $4.  Later when we
		# use the residue index to look this up, we must use the resiude index as found in the
		# original PDB file, and not the "modified_index" which is the result of our renumbering scheme.
    }
    # Parse ATOM section.
    if (/^ATOM..(.....).(....)(.)(...).(.)(....)(.)...(.{8})(.{8})(.{8})(.{6})(.{6})(.*)/ && !$END_MODEL){      
      
        if       ($3 eq " ")                         {$SET_LABEL=0; $INCLUDE=1;}
        elsif    ($3 eq $LABEL && $SET_LABEL )       {$INCLUDE=1              ;} 
        elsif    ($3 ne $LABEL && $SET_LABEL )       {$INCLUDE=0              ;}
        elsif    (! $SET_LABEL && $3 ne " ")         {$SET_LABEL=1; $LABEL=$3 ; $INCLUDE=1}
          # the above appears to be logic to include only a given "alternate atom" entry in the
          # case that some atoms have multiple entries;  we include the first one seen, and ensure
          # that subsequent atoms included that have an alternate marker also choose the same one.
          # tfb notes nov 2009
	  
	my $occupancy = "  1.00";
		# since we only include 1 atom per set of alternates, we make the occupancy 1.0    

        if ( $INCLUDE ) {
       
            # deal with insertion codes.  insertion codes seem like a pdb-format hack, and are
            # often problematic.  e.g., see 2zcy, in which we find residues 102, 10A, and 103 sequentially.
            # An examination of SEQRES indicates this residue does in fact belong between 102 and 103, but
            # we can't just call it res 10, because a res 10 already exists, and this will cause problems
            # downstream (notably modeller).  Another example: 1eu3 - note the 2A residue index that starts
            # chain A, followed by residues 1,2,3... So we deal with this by always renumbering residue indices
            # starting at 1.  We also watch the sequence of residue indices to preserve gaps in the numbering
            # scheme which can be used in an attempt to identify atom "segments" for some alignment schemes.
            #
            if( $lastRes ne "$5_$6_$7" ) {
                $lastRes = "$5_$6_$7";
                my $diff = $6 - $lastResNumeric;
                $lastResNumeric = $6;
                $modifiedIndex++;
                if( $diff != 1 && $diff != 0 ) {
                    if( $diff < 1 ) {
                        $modifiedIndex++;
                    }
                    elsif( $diff > 1 ) {
                        $modifiedIndex += $diff - 1;
                    }
                        # preserve non-sequential residue numbering, required by getATOMSegments()
                        # note that in the case of $diff < 1, we just indicate a gap by a single inc;
                        # otherwise we indicate the size of the gap with multiple inc
                }
            }
            
            if( !$firstChainResidue{$5} ) {
                $modifiedIndex = 1;
                $firstChainResidue{$5} = $6;
                    # track the first residue number for each chain; used to help apply MODRES to SEQRES in applyModres
            }

            
            push(@A_alternate,     $3);     
            push(@A_atom_index,    $1);     
            push(@A_atom_name,     $2);     
            push(@A_residue_name,  $4);    
            push(@A_chain,         $5);    
            push(@A_residue_index, $6);     
            push(@A_ins,           $7);     
            push(@A_X,             $8);           
            push(@A_Y,             $9);           
            push(@A_Z,             $10);               
            push(@A_occupancy,     $occupancy);    
            push(@A_tempFactor,    $12);   
            push(@A_theRest,       $13);     
            push(@A_index,         $modifiedIndex);

            push( @A_was_HETATM, 0 );
        }
        else {
            #print " $1 3=$3 I=$INCLUDE S=$SET_LABEL L=$LABEL <br>";
        }
    }
    elsif (/^HETATM(.....).(....)(.)(...).(.)(....)(.)...(.{8})(.{8})(.{8})(.{6})(.{6})(.*)/ && !$END_MODEL){
		
		my $res = $M_chain_resindex{ "$5_$6_$7" };
		if( $res ne "" && checkRes( $res ) ne "UNK" ) {
			# here we accept HETATM records as ATOM only if they appear in the MODRES section AND
			# correspond to a recognized amino-acid name.
            $convertedHETATM++;

            if       ($3 eq " ")                         {$SET_LABEL=0; $INCLUDE=1;}
            elsif    ($3 eq $LABEL && $SET_LABEL )       {$INCLUDE=1              ;} 
            elsif    ($3 ne $LABEL && $SET_LABEL )       {$INCLUDE=0              ;}
            elsif    (! $SET_LABEL && $3 ne " ")         {$SET_LABEL=1; $LABEL=$3 ; $INCLUDE=1}
              # the above appears to be logic to include only a given "alternate atom" entry in the
              # case that some atoms have multiple entries;  we include the first one seen, and ensure
              # that subsequent atoms included that have an alternate marker also choose the same one.
              # tfb notes nov 2009
	      
	      
	    my $occupancy = "  1.00";
		# since we only include 1 atom per set of alternates, we make the occupancy 1.0    

            if ($INCLUDE ){
                # see long comment in ATOM section parse above regarding residue renumbering
                if( $lastRes ne "$5_$6_$7" ) {
                    # if this is not the same chain-residue-insertioncode...(it needs a unique res index)
                    $lastRes = "$5_$6_$7";
                    my $diff = $6 - $lastResNumeric;
                    $lastResNumeric = $6;
                    $modifiedIndex++;
                    if( $diff != 1 && $diff != 0 ) {
                        if( $diff < 1 ) {
                            $modifiedIndex++;
                        }
                        elsif( $diff > 1 ) {
                            $modifiedIndex += $diff - 1;
                        }
                            # preserve non-sequential residue numbering, required by getATOMSegments()
                            # note that in the case of $diff < 1, we just indicate a gap by a single inc;
                            # otherwise we indicate the size of the gap with multiple inc
                    }
                }
                
                if( !$firstChainResidue{$5} ) {
                    $modifiedIndex = 1;
	                $firstChainResidue{$5} = $6;
	                    # track the first residue number for each chain; used to help apply MODRES to SEQRES in applyModres
                }



                push(@A_alternate,     $3);     
                push(@A_atom_index,    $1);     
                push(@A_atom_name,     $2);     
                push(@A_residue_name,  $4);    
                push(@A_chain,         $5);    
                push(@A_residue_index, $6);     
                push(@A_ins,           $7);     
                push(@A_X,             $8);           
                push(@A_Y,             $9);           
                push(@A_Z,             $10);               
                push(@A_occupancy,     $occupancy);    
                push(@A_tempFactor,    $12);   
                push(@A_theRest,       $13);     
                push(@A_index,         $modifiedIndex);

                push( @A_was_HETATM, 1 );
                    # I am marking that this record was a HETATM record, but we are placing it into the ATOM array
                    # and it will get printed as an ATOM record in our cleaned pdb file.
                
              
            }
            else {
              #print " $1 3=$3 I=$INCLUDE S=$SET_LABEL L=$LABEL <br>";
            }
        }
    }
		
    # parse reseq section
    if (/^SEQRES..(..).(.).(....)..(...).(...).(...).(...).(...).(...).(...).(...).(...).(...).(...).(...).(...)(.*)/){ 
       push(@S_chain,        $2);
       push(@S_residue_name, $4);
       push(@S_chain,        $2);
       push(@S_residue_name, $5);
       push(@S_chain,        $2);
       push(@S_residue_name, $6);
       push(@S_chain,        $2);
       push(@S_residue_name, $7);
       push(@S_chain,        $2);
       push(@S_residue_name, $8);
       push(@S_chain,        $2);
       push(@S_residue_name, $9);
       push(@S_chain,        $2);
       push(@S_residue_name, $10);
       push(@S_chain,        $2);
       push(@S_residue_name, $11);
       push(@S_chain,        $2);
       push(@S_residue_name, $12);
       push(@S_chain,        $2);
       push(@S_residue_name, $13);
       push(@S_chain,        $2);
       push(@S_residue_name, $14);
       push(@S_chain,        $2);
       push(@S_residue_name, $15);
       push(@S_chain,        $2);
       push(@S_residue_name, $16);
	   
	   # Now get rid of any of those that were just blank
	   my $lastSeqRes;
	   do {
			$lastSeqRes = $S_residue_name[$#S_residue_name];
			if( $lastSeqRes eq "   " ) {
				pop @S_residue_name;
				pop @S_chain;
			}
	   } while ( $lastSeqRes eq "   " );
    }
	# Defined conflict between RESSEQ and ATOM section
	if(/^SEQADV.(....).(...).(.).(....)..(....).(.{8}).(...).(....).(.{19})/){
		 push(@C_pdb_res_name,    $2);
		 push(@C_chain,           $3);
		 push(@C_pdb_index,       $4);
		 push(@C_database,        $5);
		 push(@C_accession_number,$6);
		 push(@C_seq_res_name,    $7);
		 push(@C_seq_index,       $8);
		 push(@C_conflict,        $9);
	}	  
  }
  
  foreach( @M_chain ) {
    push( @M_chainBegin, $firstChainResidue{$_} );
  }
}

###################################################################################################################################
sub getATOMSegments( $$$$ ) {
    my( $this, $chain, $atomStart, $atomEnd ) = @_;
    
    # return reference to array which holds the indices of $this->{ATOM} entries
    # that define beginning of segments.  A segment is a series of ATOM records
    # for which the residue number is sequential.
    my @segments;
    my @currentSegment;
    my $lastIndex  = -99999;
    foreach( 0..$this->{ATOM}{sizeOf} ) {
		my $res = $this->{ATOM}{residue_name}[$_];
        $res = checkRes($res); 
        
		if( $this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chain &&
            $this->{ATOM}{modified_index}[$_] >= $atomStart && $this->{ATOM}{modified_index}[$_] <= $atomEnd ) {
			my $curIndex = $this->{ATOM}{modified_index}[$_];
			if( $curIndex != $lastIndex + 1 ) {
				if( scalar @currentSegment > 0 ) {
					my @copy = @currentSegment;
					push @segments, \@copy;
				}
				@currentSegment = ();
			}
			push @currentSegment, { modified_index => $curIndex,
									residue_name   => $this->{ATOM}{residue_name}[$_] };
			$lastIndex = $curIndex;
        }
    }
    if( scalar @currentSegment > 0 ) {
		push @segments, \@currentSegment;
    }
    return \@segments;
}


###################################################################################################################################

sub alignATOMSegsToSeqres( $$$$ ) {
	my ( $this, $chain, $atomStart, $atomEnd ) = @_;
	
	my @alignment;
	
	my $segsRef = $this->getATOMSegments( $chain, $atomStart, $atomEnd );
	my @segs    = @$segsRef;
	my $seqResRef = $this->getSeqresSeqForChain( $chain, 1 );
	my $seqres;
	map{ $seqres .= amino3to1( $_ ); } @$seqResRef;

	my $seqresOffset = 0;
	my $alignIndex = 0;
	foreach( @segs ) {
		my @seg = @$_;
		my $segment;
		map { $segment .= amino3to1( checkRes( $_->{residue_name} ) ); } @seg;
		my $offset = index( substr( $seqres, $seqresOffset), $segment );
		if( $offset < 0 ) {
			print "ERROR alignATOMSegsToSeqres: \'$segment\' not found in SEQRES.\n";
			print "  (The PDB file ATOM records are inconsistent with the SEQRES.)\n";
			return undef;
		}
		
		# insert the gaps since the previous aligned segment
		if( $alignIndex > 0 ) {
			# if alignIndex is 0, we don't add gaps since these are gaps on the front
			# of the atom sequence that we don't care about.
			for( my $i=0; $i<$offset; $i++ ) {
				my $lastResIndex = $alignment[ $#alignment ]{ atom_resindex };
				push @alignment, {	"atom_residue"      => 'GAP',
									"atom_resindex"	    => $lastResIndex + 1,
                                    "seqres_residue"    => $seqResRef->[ $seqresOffset + $i ],
                                    "seqres_resindex"   => $seqresOffset + $i };
				$alignIndex++;
			}
		}
        $seqresOffset += $offset;
            # now pointing at portion of seqres that aligned to present atom segment
		
		# insert the current aligned segment
		for( my $i=0; $i < length( $segment ); $i++ ) {
			push @alignment, {	"atom_residue"      => $seg[ $i ]->{ residue_name },
								"atom_resindex"	    => $seg[ $i ]->{ modified_index },
								"seqres_residue"    => $seqResRef->[ $seqresOffset + $i ],
								"seqres_resindex"   => $seqresOffset + $i + 1 };
			$alignIndex++;
		}
		$seqresOffset += length( $segment );
	}
	
	return \@alignment;
}

sub convertSegsAlignToBL2SeqFormat( $$ ) {
    # this is an adaptor function: much of this module was written with an alignment format I wrote
    # when useing the bl2seq, but I've created a new format for the segment-based alignment,
    # and rather than editing every routien in the module, we'll just convert to the used format
    # here, and save the tedious work for later.  ;)
    
    my ( $this, $alignRef ) = @_;
    my @align = @$alignRef;
    
    my $gaps;
    my @query;
    map { push @query, amino3to1( $_->{ atom_residue } ); if( $query[$#query] eq '-' ) { $gaps++; } } @align;
    my $queryStart = $align[0]->{ atom_resindex };
    my $queryEnd   = $align[ $#align ]->{ atom_resindex };
    
    my @subject;
    map { push @subject, amino3to1( $_->{ seqres_residue } ); } @align;
    my $subjectStart = $align[0]->{ seqres_resindex };
    my $subjectEnd   = $align[ $#align ]->{ seqres_resindex };

    my %alignment;
   
	$alignment{query} = \@query;
		# holds array of residues in the query, including '-' chars for gaps
	$alignment{queryStart} = $queryStart;
	$alignment{queryEnd} = $queryEnd;
		# the starting and ending residue of alignment, from original query sequence	
	$alignment{subject} = \@subject;
		# holds array of residues in the subject
	$alignment{subjectStart} = $subjectStart;
	$alignment{subjectEnd} = $subjectEnd;
		# the starting and ending residue of alignment, from original subject sequence
	$alignment{gaps} = $gaps;
    $alignment{identities} = (scalar @query) - $gaps;


	return ( \%alignment, $gaps );
}

sub printAtomToSeqresAlignment( $$$ ) {
    my( $this, $alignRef, $fhandle ) = @_;
    my @alignment = @$alignRef;
    
    
    if( ! $fhandle ) {
        $fhandle = STDOUT;
    }

    printf $fhandle "First,Last ATOM residue index   : %d, %d\n", $alignRef->[0]->{ atom_resindex }, $alignRef->[$#alignment]->{ atom_resindex };
    printf $fhandle "First,Last SEQRES residue index : %d, %d\n", $alignRef->[0]->{ seqres_resindex }, $alignRef->[$#alignment]->{ seqres_resindex };
    
    my $totalLen = 0;
    while( $totalLen <= $#alignment ) {
        for( my $i=0; $i < 60; $i++ ) {
            if( $i + $totalLen > $#alignment ) { last; }
            print $fhandle amino3to1( $alignment[ $i + $totalLen ]->{ atom_residue } );
        }
        print $fhandle "\n";
        for( my $i=0; $i < 60; $i++ ) {
            if( $i + $totalLen > $#alignment ) { last; }
            print $fhandle amino3to1( $alignment[ $i + $totalLen ]->{ seqres_residue } );
        }
        print $fhandle "\n\n";
        $totalLen += 60;
    }
}
###################################################################################################################################


###################################################################################################################################
sub logConvertedHETATM($$$) {
    # print the entries that were converted from HETATM to ATOM due to a MODRES record
    # that indicates a nonstandard residue modified from a recognized amino-acid name.
    # (These entries were converted during parsePDB -- see HETATM processing there)
    # The common example is MSE -> MET, in which the MSE atoms were given as HETATM in
    # the original PDB, but we added them to the ATOM array as they were seen because
    # we saw them in the MODRES section (see parsePDB for MODRES/HETATM records)
    my( $this, $code, $name ) = @_;
    my $count = 0;
    if( $convertedHETATM > 0 ) {
      	my %uniqueResChain = {};
    	open HET, ">$name";
    	for( my $i=0; $i<=$this->{ATOM}{sizeOf} && $count < $convertedHETATM; $i++ ) {
            if( $this->{ATOM}{was_HETATM}[$i] ) {
                print HET "$code $this->{ATOM}{atom_index}[$i] $this->{ATOM}{atom_name}[$i] $this->{ATOM}{residue_name}[$i] $this->{ATOM}{residue_index}[$i] $this->{ATOM}{chain}[$i]\n";
                $count++;        
            }
        }
        close HET;
    }
}

###################################################################################################################################
sub logResidueNameTranslations( $$$$ ) {
	# check the residue names of the atoms and print a list of residues using "nonstandard" amino acid names
	# so we'll know what we will have translated when we clean these up.
    # This should be called *before* applyModRes since the latter will change the residue names to standard
    # for any that appeared in the MODRES section.
   	
    # Note: this only prints unique translations per chain, so it may be that multiple MSE->MET translations
	# occur in a given chain, but only one entry will be printed for brevity.

	my ($this, $code, $chain, $name) = @_;
	my %uniqueResChain = {};
        # hash to track which residues we've already noted in this chain: see "Note" above
	open TR, ">$name";
	my $count = 0;
	my $dropped;
	for( my $i=0; $i<=$this->{ATOM}{sizeOf}; $i++ ) {
		my $res = $this->{ATOM}{residue_name}[ $i ];
		
		#my $modres = $this->{ATOM}{chain_resindex}{ "$this->{ATOM}{chain}[$i]" . "_$this->{ATOM}{residue_index}[$i]" . "_$this->{ATOM}{ins}[$i]" };
		my $modres = $M_chain_resindex{ "$this->{ATOM}{chain}[$i]" . "_$this->{ATOM}{residue_index}[$i]" . "_$this->{ATOM}{ins}[$i]" };
			# note: we really should use the $this reference to this hash to support a true object-use but somehow I'm not
			# getting the reference to the hash setup correctly.
		
		# and then check any translation due to 'checkRes'
		my $finalRes;
		if( $modres ne "" ) {
		  $finalRes = checkRes( $modres );
		}
		else {
		  $finalRes = checkRes( $res );
		}
		
		$dropped = "";
		if( $finalRes eq "UNK" && ! $includeUNK ) {
			$dropped = "DROPPED!";
		}
		
		# print the translation, if any
		if( $res ne $finalRes ) {
			my $resChain = $res . $finalRes . $this->{ATOM}{chain}[$i] . $this->{ATOM}{residue_index}[$i];
			if( ! $uniqueResChain{ $resChain } ) {
                print TR "$code [$chain] : $this->{ATOM}{chain}[$i] $this->{ATOM}{residue_index}[$i] $res ";
                if( $modres ne "" ) {
                    print TR "MODRES--> $modres ";
                    if( $modres ne $finalRes ) {
                        print TR "STD--> $finalRes ";
                    }
                }
                elsif( $res ne $finalRes ) {
                    print TR "STD--> $finalRes ";
                }
                if( $this->{ATOM}{was_HETATM}[$i] ) {
                    print TR "(was HETATM)"
                }
                print TR " $dropped\n";
                  # when the various files that are generated from the ATOM records are printed, the UNK residues
                  # will be dropped if includeUNK is not set.
                $uniqueResChain{ $resChain } = 1;
                $count++;
			}
		}
	}

	foreach (0 .. $this->{SEQRES}{sizeOf}) {
		my $res = $this->{SEQRES}{residue_name}[$_];
		my $finalRes = checkRes($res);
		if( $finalRes ne $res ) {
			$dropped = ( $finalRes eq "UNK" && ! $includeUNK ) ? "DROPPED!" : "";
			print TR "$code [$chain] : $this->{SEQRES}{chain}[$_] SEQRES residue $_ $res --> $finalRes  $dropped\n";
			$count++;
		}
	}
	
	close TR;
	if( $count == 0 ) {
		unlink( $name );
	}
}

###################################################################################################################################
# the atom names as they should appear in the standard residue.  The main example we
# have is the MSE residue, in which a selenium atom has replaced the standard sulfur
# of MET.  In applyModRes below, the MSE will get converted to MET, and at the same
# time we will translate the atom names as appropriatiate.  Note that you can add
# any number of atom translations per residue name.
# NOTE: we try to make the resulting atom align with the other atoms by knowing that
# the pdb format allows 4 spaces for an atom name.
# NOTE: I have added a 2-element array for this translation: the first element is the
# name to use in the atom_name field; the second is used in the atom element field -
# these may often be the same, but need not be (as in the case of the many unique
# names given to the carbon atoms in a residue, but which all have element C)
my %modAtomNames = ( 'MSE' => { 'SE'=>[' S  ', 'S' ]  },
						# this says that in any MSE residues, for any atom named SE, we should
						# changed the name to S (first S), and give it element name S (second S)
						# Will we need to extend this to explicity name what the element previously was?
                   );

sub applyModRes($) {
  my($this) = @_;
  # this used to be done in the ParsePDB function, but I have factored it out
  # because I want to be able to call logResidueNameTranslations before the MODRES
  # are applied (and in that function note which ones will be changed by MODRES)

  # Here we alter residue names from the modified name to the standard name
  # as specified in the MODRES records of the pdb, parsed earlier.
  for( my $i=0; $i<=$this->{ATOM}{sizeOf}; $i++ ) {
    #my $newRes = $this->{MODRES}{chain_resindex}{ "$this->{ATOM}{chain}[$i]" . "_$this->{ATOM}{residue_index}[$i]"  . "_$this->{ATOM}{ins}[$i]"};
	my $newRes = $M_chain_resindex{ "$this->{ATOM}{chain}[$i]" . "_$this->{ATOM}{residue_index}[$i]" . "_$this->{ATOM}{ins}[$i]" };
			# note: we really should use the $this reference to this hash to support a true object-use but somehow I'm not
			# getting the reference to the hash setup correctly.
		
    if( $newRes ne "" ) {
	    print "ATOM: Mapping chain $this->{ATOM}{chain}[$i] residue $this->{ATOM}{residue_index}[$i] from $this->{ATOM}{residue_name}[$i] to $newRes\n";
		if( $modAtomNames{ $this->{ATOM}{residue_name}[$i] } ) {
			my $oldName = trim( $this->{ATOM}{atom_name}[$i] );
			my $newName = $modAtomNames{ $this->{ATOM}{residue_name}[$i] }->{ $oldName }[0];
			if( $newName ne "" ) {
				print "   atom $this->{ATOM}{atom_name}[$i] changed to $newName\n";
				$this->{ATOM}{atom_name}[$i] = $newName;
				# also try to replace the atom element name in the end of the record:
				my $newElement = $modAtomNames{ $this->{ATOM}{residue_name}[$i] }->{ $oldName }[1];
				$this->{ATOM}{theRest}[$i] =~ s/$oldName/$newElement/;
			}
		}
	    $this->{ATOM}{residue_name}[ $i ] = $newRes;
    } 
  }
  
    # Also map the residue names for the SEQRES records;
    for( my $i=0; $i<=$this->{MODRES}{sizeOf}; $i++) {
        my $chain  = $this->{MODRES}{chain}[$i];            
        my $chainResIndex  = $this->{MODRES}{residue_index}[$i];            
        my $resMod = $this->{MODRES}{residue_name}[$i];            
        my $resStd = $this->{MODRES}{std_res_name}[$i];
        my $chainBegin = $this->{MODRES}{chainBegin}[$i];
            # chainBegin is the first residue index seen for ATOM records, and as tells us something about
            # how to apply MODRES to the seqres sequence: normally we think of SEQRES as starting at 1,
            # even if the first ATOM residue seen is some higher number ( because some ATOM residues are
            # missing ) but it can also happen that we see ATOM residue < 1, in which case the SEQRES
            # sequence also needs to be considered as starting below 1: so if this chainBegin is < 1,
            # we need to adjust our addressing of the seqres below.
			# Note: this still leaves the possibility that the first index seen was > 1, in which
			# case we cannot really hope to apply the change here without first doing an alignment
			# of ATOM->SEQRES to discover the mapping.  This could be done, but will not be at this
			# late stage of CASP preparation.  However, it is not a critical error, since the name
			# translation will occur in any case due to the MODRES translation information we make
			# use of in checkRes.
        my $offset = 1;
            # normally the seqres is assumed to start at index 1
        if( $chainBegin < 1 ) {
            $offset = $chainBegin;
        }
        
        my $chainIdx = 0;
        my $lastChain = "";
        my $found = 0;
        my $skipCount = 0;
        for( my $j=0; (! $found) && $j<=$this->{SEQRES}{sizeOf}; $j++ ) {
            my $curChain = $this->{SEQRES}{chain}[$j];
            $chainIdx++;
            if( $curChain ne $lastChain ) {
                $chainIdx  = $offset;
                $lastChain = $curChain;
                $skipCount = 0;
            }
            my $curRes = $this->{SEQRES}{residue_name}[$j];
            if( checkRes( $curRes ) eq 'UNK' ) {
                $skipCount++;
                    # skipcount is kind of a hack to deal with HETATM residues that are listed in the SEQRES,
                    # but did not appear in a MODRES record (if they had been in a MODRES, we would have added
                    # the HETATM as an ATOM record, and that residue index would have affected the $offset above).
                    # An example of this can be seen in 1CCR - the 0 residue is ACE - acetylation - it is not a MODRES
                    # and so not included in the ATOM residue sequence, but it does appear in the SEQRES record
                    # and therefore affects the index addressing of the residues.
                next;
            }

            if( $curChain eq $chain && $chainResIndex == ($chainIdx - $skipCount) ) {
                if( $curRes ne $resMod ) {
                    print "WARNING: (applyModRes - SEQRES) expected to find $resMod from MODRES chain $chain index $chainResIndex but found $curRes.\n";
                    print "          The residue numbering for this PDB file is probably not consistent between ATOM records and SEQRES records.\n";

                }
                else {
                    print "SEQRES: Mapping chain $curChain residue $chainIdx from $curRes to $resStd\n";
                    $this->{SEQRES}{residue_name}[$j] = $resStd;
                        # change the residue name per the MODRES record
                }
                $found = 1;
            }
        }
    }
}

###################################################################################################################################
sub PrintPdbLength($$){ 
    my($this, $chainL) = @_;
    
    my $n=0;
    my $res;
    foreach (0 .. $this->{ATOM}{sizeOf}) {
        $res = $this->{ATOM}{residue_name}[$_];
        $res = checkRes($res);
        if($this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chainL){
            $n++;
        }
    }
    return $n;
}

###################################################################################################################################
sub PrintPdbSeqresSeq($$){ 
    my($this, $fname, $chainL, $pdbc) = @_;
    
    my $FHANDLE; # File handle changer.
    if (!$fname =~ /\S/) {
        $FHANDLE = 'STDOUT';
    }
    else {
        open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
        $FHANDLE = 'OUT';
    }
    
    my $n=0;
    my $res;
    my $res1;
    print $FHANDLE ">$pdbc.seqres\n";
    foreach (0 .. $this->{SEQRES}{sizeOf}) {
        $res = $this->{SEQRES}{residue_name}[$_];
        $res = checkRes($res);
        if( ($includeUNK || $res ne "UNK") && $this->{SEQRES}{chain}[$_] eq $chainL){
            $n++;
            $res1 = amino3to1($res);
            print $FHANDLE "$res1";
            if($n%60==0){print $FHANDLE "\n";}
        }
    }
    close($FHANDLE);
}

###################################################################################################################################
sub PrintPdbAtomAll($$) { 
    my($this, $fname) = @_;
    my $FHANDLE; # File handle changer.
    if (!$fname =~ /\S/) {
        $FHANDLE = 'STDOUT';
    }
    else {
        open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
        $FHANDLE = 'OUT';
    }
    
    my $res;
    my $lastResidue=999999;
    
    # setup which residues will get printed: we will only print residues that
    # include a CA atom.
    #
    my @resHasCA;
    # array with entry per residue: set to 0 for res with no CA atoms, 1 for valid res w/CA atoms
    # note this is 1-based 
    my $currResIdx=0;
    # the current residue index (as found in 'modified_index')
    
    foreach (0 .. $this->{ATOM}{sizeOf} ) {
        $res = $this->{ATOM}{residue_name}[$_];
        $res = checkRes($res);
        
        # each time we find a new residue number, mark that residue 0 in resHasCA array,
        # it will get marked to 1 only if a CA atom is found.  This 
        if($this->{ATOM}{modified_index}[$_] != $lastResidue  && ($includeUNK || $res ne "UNK")){
            $currResIdx++;
            $lastResidue=$this->{ATOM}{modified_index}[$_];
            $resHasCA[$currResIdx]=0;
        }
        if($this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK")){
            $resHasCA[$currResIdx]=1;
        }
    }
    
    printf $FHANDLE "HEADER\n";
    
    $currResIdx=0;
    $lastResidue=999999;
    my $start=1;
    my $last_chain=" ";
    my $chain;
    my $flag=0;
    foreach (0 .. $this->{ATOM}{sizeOf} ) {
        $res = $this->{ATOM}{residue_name}[$_];
        $chain =$this->{ATOM}{chain}[$_]; 
        if ($start) { $start=0; $last_chain=$chain;}
        
        # print CTER if we found a new chain and our flag was set indicating we
        # have seen residues belonging to last_chain
        if ($last_chain ne $chain) {
            if ($flag>0) { printf $FHANDLE "CTER\n"; }
            $last_chain=$chain;
            $flag=0;
        }
        $res = checkRes($res);
        if($this->{ATOM}{modified_index}[$_] != $lastResidue  && ($includeUNK || $res ne "UNK")){
            $currResIdx++;
            $lastResidue=$this->{ATOM}{modified_index}[$_];
        }    
        if (($includeUNK || $res ne "UNK") && $resHasCA[$currResIdx]==1 ) {
            printf $FHANDLE "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
            $this->{ATOM}{atom_index}[$_],
            $this->{ATOM}{atom_name}[$_],
            $res,
            $this->{ATOM}{chain}[$_],
            $this->{ATOM}{modified_index}[$_],
            $this->{ATOM}{X_coordinate}[$_],
            $this->{ATOM}{Y_coordinate}[$_],
            $this->{ATOM}{Z_coordinate}[$_],
            $this->{ATOM}{occupancy}[$_],
            $this->{ATOM}{tempFactor}[$_],
            $this->{ATOM}{theRest}[$_];
            $flag++;
        }
    }
    if($flag>0){printf $FHANDLE "%-80s\n", 'CTER';}
    printf $FHANDLE "%-80s\n", 'END';
    
    if ($FHANDLE eq 'OUT') {
    close(OUT) || carp("Can't close $fname. $!\n");
    }
}

###################################################################################################################################
sub PrintPdbAtom($$){ 
	
    my($this, $fname, $chainL, $subjectStartRes, $subjectEndRes, $queryLength) = @_;
    
    my $FHANDLE; # File handle changer.
    if (!$fname =~ /\S/) {
        $FHANDLE = 'STDOUT';
    }
    else {
        open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
        $FHANDLE = 'OUT';
    }
    printf $FHANDLE "HEADER\n";
    
    my $res;
    my $lastResidue=999999;
    
    # setup which residues will get printed: we will only print residues that
    # include a CA atom.
    #
    my @resHasCA;
        # array with entry per residue: set to 0 for res with no CA atoms, 1 for valid res w/CA atoms
        # note this is 1-based 
    my $currResIdx=0;
        # counts unique residue indices we've run into
    
    foreach (0 .. $this->{ATOM}{sizeOf} ) {
		if( $this->{ATOM}{chain}[$_] eq $chainL ) {
			$res = $this->{ATOM}{residue_name}[$_];
			$res = checkRes($res);
			
			if( $includeUNK || $res ne "UNK" ) {
			
				# each time we find a new residue number, mark that residue 0 in resHasCA array,
				# it will get marked to 1 only if a CA atom is found.
				if($this->{ATOM}{modified_index}[$_] != $lastResidue ) {
					$currResIdx++;
					$lastResidue=$this->{ATOM}{modified_index}[$_];
					$resHasCA[$currResIdx]=0;
				}
				if($this->{ATOM}{atom_name}[$_] eq " CA " ){
					$resHasCA[$currResIdx]=1;
				}
			}
		}
    }
    
    
    #
    # get the starting and ending residue number we should include based on the alignment information
    # we were passed
    #
    my( $totalResCount, $resStartActual, $resEndActual, $resCountActual ) =
        getResidueRangeFromAlignmentAndQueryLength( $this, $chainL, $subjectStartRes, $subjectEndRes, $queryLength );
        
    $currResIdx=0;
    $lastResidue=999999;
    foreach (0 .. $this->{ATOM}{sizeOf}) {
		if( $this->{ATOM}{chain}[$_] eq $chainL ) {
			$res = $this->{ATOM}{residue_name}[$_];
			$res = checkRes($res);
			
			if( $includeUNK || $res ne "UNK" ) {
				if( $this->{ATOM}{modified_index}[$_] != $lastResidue ) {
					$currResIdx++;
					$lastResidue=$this->{ATOM}{modified_index}[$_];
				}    
				if ( $this->{ATOM}{modified_index}[$_]>=$resStartActual && $this->{ATOM}{modified_index}[$_]<=$resEndActual && $resHasCA[$currResIdx]==1 ) {
					printf $FHANDLE "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
					$this->{ATOM}{atom_index}[$_],
					$this->{ATOM}{atom_name}[$_],
					$res,
					$this->{ATOM}{chain}[$_],
					$this->{ATOM}{modified_index}[$_],
					$this->{ATOM}{X_coordinate}[$_],
					$this->{ATOM}{Y_coordinate}[$_],
					$this->{ATOM}{Z_coordinate}[$_],
					$this->{ATOM}{occupancy}[$_],
					$this->{ATOM}{tempFactor}[$_],
					$this->{ATOM}{theRest}[$_];
				}
			}
		}
    }
    
    printf $FHANDLE "%-80s\n", 'CTER';
    printf $FHANDLE "%-80s\n", 'END';
    
    if ($FHANDLE eq 'OUT') {
        close(OUT) || carp("Can't close $fname. $!\n");
    }
}

###################################################################################################################################
sub getResidueRangeFromAlignmentAndQueryLength( $$$$$ ) {
    # written for use by the routines which print the .seq, .log, .xyz for a residue range
    # given the subject start/end residues and queryLength of an alignment.  The goal of
    # the code this replaces, which appeared in each fn, seems to be to calculate a start and
    # end residue that will produce a segment which is as close to the query length
    # as possible - in contrast to just always printing the info for all of the residues.
    #
    # This is presumably to save computation on longer templates.
    
    # WARNING: this is not robust for dynamic database creation when we're trying
    # to clip a long structure to the query length becauseof (a) gaps in the residue numbering
    # (b) missing residues in the ATOM records and (c) it does not compensate for gaps in the
    # alignment that produced the $aliStart and $aliEnd arguments.
    #
    # (a) and (b) can be fixed by doing the ATOM->SEQRES alignment first, and using $aliStart and
    # $aliEnd to index into this alignment and retrieve the corresponding modified index to return.
    #
    # For now, I'm going instead to simply have dynamic entries always use the entire sequence as is
    # done when a database entry is created, since I know this works.  I'm not sure why this was
    # not always done, except that this can potentially take much longer for long sequences.
    
    my( $this, $chainL, $aliStart, $aliEnd, $queryLength ) = @_;
    
    my $totalResCount = 0;
    my $resStart = 0;
    my $resEnd = 0;
    my $resCount = 0;
        # the above values are all computed and returned from this subroutine
        
    my @CA_index;
        # index of the CA atoms for each residue for the passed chain
    foreach( 0 .. $this->{ATOM}{sizeOf} ) {
        my $res = $this->{ATOM}{residue_name}[$_];
        $res = checkRes($res);
        if( $this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chainL ) {
            $CA_index[$totalResCount]=$this->{ATOM}{modified_index}[$_];
            $totalResCount++;
        }
    }
    my @CA_index_sorted = sort {$a <=> $b} @CA_index;
        # this takes care of non-ascending residue numbers seen in some pdb files -- see pdb2f07.ent
    
    my ($indexA, $indexB);
    if( $totalResCount <= $queryLength ) {
        $indexA = 0;
        $indexB = $totalResCount - 1;
    }
    else {
        # in the case that our (subject) chain is longer than the query, some pains
        # are taken presumably to save computation downstream by making the subject
        # database entries (e.g. .log, .xyz, .2nd etc) only as long as the query,
        # but not longer.  The offsets have been computed which are equal length (within 1)
        # stretches applied to either end of the aligned section to lengthen the subject
        # in both directions where possible.
        my $extraLen = $queryLength - ( $aliEnd - $aliStart + 1 );
        my $offset1  = int($extraLen / 2);
        my $offset2  = $extraLen - $offset1;
            # offset1 and 2 are amounts to extend the start and end by, if possible.
			
		# note that the numeric params passed to us (aliStart, aliEnd) are 1-based indices
		# read from a blast alignment.  We want indexA and indexB to be 0-based since they
		# will index into the CA_index_sorted array.  
            
        if( $aliStart - $offset1 > 0 && $aliEnd + $offset2 <= $totalResCount ) {
            # lengthen subject in at each end
            $indexA = $aliStart - $offset1 - 1;
            $indexB = $aliEnd + $offset2 - 1;
        }
        elsif( $aliStart - $offset1 <= 0 ) {
            # start at the first residue, and apply all extra length to the end
            $indexA = 0;
            $indexB = $queryLength - 1;
        }
        else {
            # go to last residue, and apply all extra length to beginning
            $indexA = $totalResCount - $queryLength;
            $indexB = $totalResCount - 1;
        }
    }

    $resStart = $CA_index_sorted[ $indexA ];
    $resEnd   = $CA_index_sorted[ $indexB ];
    $resCount = $indexB - $indexA + 1;
    
    return ( $totalResCount, $resStart, $resEnd, $resCount );
}

###################################################################################################################################
sub PrintLOOPPseq(){ 
    my($this, $pdbcodeL, $chainL, $subjectStartRes, $subjectEndRes, $queryLength) = @_;
    
    my $fname;
    my $pdbCodeProcessed;
    
    if($chainL eq " "){
        $fname = "$pdbcodeL.seq";
        $pdbCodeProcessed = "$pdbcodeL";
    }
    else {
        $fname = "$pdbcodeL" . "_$chainL.seq";
        $pdbCodeProcessed = "$pdbcodeL" . "_$chainL";
    }
    
    my $FHANDLE; # File handle changer.
    if (!$fname =~ /\S/) {
        $FHANDLE = 'STDOUT';
    }
    else {
        open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
        $FHANDLE = 'OUT';
    }

    my( $totalResCount, $resStartActual, $resEndActual, $resCountActual ) =
        getResidueRangeFromAlignmentAndQueryLength( $this, $chainL, $subjectStartRes, $subjectEndRes, $queryLength );
        
    print $FHANDLE "$pdbCodeProcessed $resCountActual\n";

    foreach (0 .. $this->{ATOM}{sizeOf}) {
        my $res = $this->{ATOM}{residue_name}[$_];
        $res = checkRes($res);
        if($this->{ATOM}{atom_name}[$_] eq " CA " && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chainL && ($this->{ATOM}{modified_index}[$_]>=$resStartActual && $this->{ATOM}{modified_index}[$_]<=$resEndActual)){
            printf $FHANDLE "%s\n",$res;
        }
    }
    
    if ($FHANDLE eq 'OUT') {
        close(OUT) || carp("Can't close $fname. $!\n");
    }
}

###################################################################################################################################

sub PrintLOOPPlog() { 
    my( $this, $pdbcodeL, $chainL, $subjectStartRes, $subjectEndRes, $queryLength ) = @_;
    
    my $fname;
    my $pdbCodeProcessed;

    if($chainL eq " ") {
        $fname = "$pdbcodeL.log";
        $pdbCodeProcessed = "$pdbcodeL";
    }
    else {
        $fname = "$pdbcodeL" . "_$chainL.log";
        $pdbCodeProcessed = "$pdbcodeL" . "_$chainL";
    }

    my $FHANDLE; # File handle changer.
    if (!$fname =~ /\S/) {
        $FHANDLE = 'STDOUT';
    }
    else {
        open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
        $FHANDLE = 'OUT';
    }

    # old way: print log from directly from the ATOM records for the specified range
    
    #my( $totalResCount, $resStartActual, $resEndActual, $resCountActual ) =
    #    getResidueRangeFromAlignmentAndQueryLength( $this, $chainL, $subjectStartRes, $subjectEndRes, $queryLength );
    #    
    #print $FHANDLE "$pdbCodeProcessed $resCountActual\n";
    #
    #my $i=0;
    #foreach (0 .. $this->{ATOM}{sizeOf}) {
    #    my $res = $this->{ATOM}{residue_name}[$_];
    #    $res = checkRes($res);
    #    if($this->{ATOM}{atom_name}[$_] eq " CA " && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chainL && ($this->{ATOM}{modified_index}[$_]>=$resStartActual && $this->{ATOM}{modified_index}[$_]<=$resEndActual)){
    #        printf $FHANDLE "%d %d\n",$i,$this->{ATOM}{modified_index}[$_];
    #        $i++;
    #    }
    #}
    ## $i should equal $resCountActual
	
	# new way: print log from the alignment of ATOM to SEQRES, using the modified-index values stored there
    # for ATOM records that exists, or assigned there in the alignment routine for any residues that were
    # inserted from the SEQRES section.  Note that the alignment already  takes into account the
    # range of residues we are working with.
	
	my $queryRef     = $this->{ALIGNATOMTOSEQRES}{query};
	my @query 		 = @$queryRef;
    
	my $resCount = scalar @query;
	print $FHANDLE "$pdbCodeProcessed $resCount\n";

    my $i = 0;
    foreach( @query ) {
		printf $FHANDLE "%d %d\n", $i, $this->{ALIGNATOMTOSEQRES}{modifiedAtomIndices}->[ $i ];
        $i++;
    }

    if ($FHANDLE eq 'OUT') {
        close OUT;
    }
}

###################################################################################################################################
sub strucToLOOPP($) {
	# converts the output of dssp/sable to codes recognized by LOOPP
	my $struc = shift;
	my $code = 'X';
	if($struc eq "H")    { $code = "A"; }
	elsif($struc eq "E") { $code = "B";	}
	return $code;
}

sub PrintLOOPP2nd( $$$$ ) {
	# prints secondary structure information for each residue in the SEQ.  Note that this is now supplemented with
	# predicted structure for any residues that were missing from the ATOM records -- we have run SABLE on the SEQRES
	# sequence to get the predicted information.
	#
	# alignATOMToSEQRES and getSable... should have already been called.  The information printed is based specifically
	# on this alignment.
	#
	# It used to be that a separate script "makeSecond.pl" was called to do this work, but that seems unnecessary since
	# with the alignment we have all the information we need here.
	
	my( $this, $code, $chain, $fname ) = @_;
	
	open OUT, ">$fname" or die "Can't open $fname in printLOOPP2nd()\n";
	
	#my $queryRef     = $this->{ALIGNATOMTOSEQRES}{query};
	#my @query 		 = @$queryRef;
	my $subjectRef   = $this->{ALIGNATOMTOSEQRES}{subject};
	my @subject		 = @$subjectRef;
	my $subjectStart = $this->{ALIGNATOMTOSEQRES}{subjectStart};
		# query is the ATOM residue sequence, subject is the SEQRES residue sequence; we'll just read
		# the subject, which is identical to the query, except there are no gaps.
	my $resCount	 = scalar @subject;
	my $resNumBegin  = $this->{ALIGNATOMTOSEQRES}{modifiedAtomIndices}->[0];
		# this is the starting residue number original found in the atom records.  We're just going to
		# increment it since now the interior gaps in the ATOM records have been filled from the SEQRES
	my $dsspStrucRef = $this->{DSSP}{structure};
	my $sableStrucRef= $this->{SABLE}{structure};

	print OUT "$code" . "_$chain $resCount\n";

    my $i = 0;
    foreach( @subject ) {
		my $res      = amino1to3( $_ );
		my $resIndex = $this->{ALIGNATOMTOSEQRES}{modifiedAtomIndices}->[ $i ];
		my $struct   = $dsspStrucRef->[ $resIndex ];
			# dssp only has struct for existing ATOM residues, and not even all of those, so we index based
			# on the actual ATOM residue number as exists in our cleaned pdb (single chain) file which the
			# dssp output was created from.
		chomp( $struct );
		if( !$struct || $struct eq "" ) {
			$struct = $sableStrucRef->[ $subjectStart + $i - 1 ];
				# here we index into the sable results by residue number from SEQRES, -1 to convert
				# the $subjectStart into an array index (was a res index output from blast)
		}
		$struct = strucToLOOPP( $struct );
	    print OUT "$res $struct\n";
        $i++;
    }

	close OUT;
}

###################################################################################################################################

sub PrintLOOPPxyz(){ 
    my($this, $pdbcodeL, $chainL, $subjectStartRes, $subjectEndRes, $queryLength) = @_;

    my $fname;
    my $pdbCodeProcessed;
    
    if($chainL eq " "){
        $fname = "$pdbcodeL.xyz";
        $pdbCodeProcessed = "$pdbcodeL";
    }
    else{
        $fname = "$pdbcodeL" . "_$chainL.xyz";
        $pdbCodeProcessed = "$pdbcodeL" . "_$chainL";
    }
    
    my $FHANDLE; # File handle changer.
    if (!$fname =~ /\S/) {
        $FHANDLE = 'STDOUT';
    }
    else {
        open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
        $FHANDLE = 'OUT';
    }

    my( $totalResCount, $resStartActual, $resEndActual, $resCountActual ) =
        getResidueRangeFromAlignmentAndQueryLength( $this, $chainL, $subjectStartRes, $subjectEndRes, $queryLength );

    my $alignmentRef = $this->{ALIGNATOMTOSEQRES};
    my %alignment = %$alignmentRef;
    my $subjRef   = $alignment{subject};
    my @subject   = @$subjRef;
    my $length    = scalar @subject;
        # new: we will use information from the ATOM -> SEQRES alignment to insert dummy coordinates where we
        # know we will be inserting residues from SEQRES into the SEQ

    print $FHANDLE "$pdbCodeProcessed $length\n";
        
    my $ind=-99999;
    my $indf=0;
    my $nme="";
    my $chn="";
    my ( $xa, $ya, $za, $xb, $yb, $zb, $xg, $yg, $zg );
    my ( $aa, $bb, $gg );
    my $count = 0;
    foreach (0 .. $this->{ATOM}{sizeOf}) {
        if($this->{ATOM}{modified_index}[$_] != $ind && $indf>0 && ($includeUNK || $nme ne "UNK") && $chn eq $chainL && ($ind>=$resStartActual && $ind<=$resEndActual)) {
            # here we recognize that the current index is not the last found index, that we have indeed previously found an index (indf>0), we have the right chain
            # and the index (ind) is in range.  So print the coordinate information that was collected while scanning the atom entries for the previous residue index.
            if($aa==0){
                $xa=$ya=$za=999.0;
            }
            if($bb==0){
                $xb=$yb=$zb=999.0;
            }
            if($gg==0){
                $xg=$xa;
                $yg=$ya;
                $zg=$za;
            }
            else{
                $xg=$xg/$gg;
                $yg=$yg/$gg;
                $zg=$zg/$gg;
            }
            if($xa!=999.0){
                printf $FHANDLE "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$xg,$yg,$zg,$xa,$ya,$za,$xb,$yb,$zb;
                $count++;
                    # count the residues we are printing coordinate information for.  Check to see if we have come to a break in the ATOM
                    # structure according to our alignment with the SEQRES, and if so, insert dummy coordinates for the inserted residues.
                while( $alignment{query}[ $count ] eq '-' ) {
                        # [ 0 ] is the entry for the ATOM residue - it will be 'GAP' in cases where the residue was missing in the ATOM records.
                        # This is the case in which we want to insert dummy coords.
                        printf $FHANDLE "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",999.0,999.0,999.0,999.0,999.0,999.0,999.0,999.0,999.0;
                        $count++;
                }
            }
        }
        if($this->{ATOM}{modified_index}[$_] != $ind){
            # the "new residue found" case
            $ind=$this->{ATOM}{modified_index}[$_];
            $nme=$this->{ATOM}{residue_name}[$_];
            $nme=checkRes($nme);
            $chn=$this->{ATOM}{chain}[$_];
            $xa=$ya=$za=$xb=$yb=$zb=$xg=$yg=$zg=0.0;
            $aa=0;
            $bb=0;
            $gg=0;
        }
        if($this->{ATOM}{atom_name}[$_] eq " CA " && $this->{ATOM}{chain}[$_] eq $chainL){
            # found the CA atom for the current residue - set xa, ya, za, aa++
            # indf tracks number of unique residues found (with CA atom)
            $xa=$this->{ATOM}{X_coordinate}[$_];
            $ya=$this->{ATOM}{Y_coordinate}[$_];
            $za=$this->{ATOM}{Z_coordinate}[$_];
            $aa++;
            $indf++;
        }
        if($this->{ATOM}{atom_name}[$_] eq " CB " && $this->{ATOM}{chain}[$_] eq $chainL){
            # found the CB atom for the current residue - set xb, yb, zb, bb++
            $xb=$this->{ATOM}{X_coordinate}[$_];
            $yb=$this->{ATOM}{Y_coordinate}[$_];
            $zb=$this->{ATOM}{Z_coordinate}[$_];
            $bb++;
        }
        if(substr($this->{ATOM}{atom_name}[$_],1,1) ne "H" && $this->{ATOM}{atom_name}[$_] ne " OXT" && $this->{ATOM}{atom_name}[$_] ne " OX " && $this->{ATOM}{atom_name}[$_] ne " CA "
        && $this->{ATOM}{atom_name}[$_] ne " C  " && $this->{ATOM}{atom_name}[$_] ne " O  " && $this->{ATOM}{atom_name}[$_] ne " N  " && $this->{ATOM}{chain}[$_] eq $chainL)  {
            # summing coords other than those excluded for some kind of center of mass computation, gg++
            $xg=$xg+$this->{ATOM}{X_coordinate}[$_];
            $yg=$yg+$this->{ATOM}{Y_coordinate}[$_];
            $zg=$zg+$this->{ATOM}{Z_coordinate}[$_];
            $gg++;
        }    
    }
    
    if($indf>0 && ($includeUNK || $nme ne "UNK") && $chn eq $chainL && ($ind>=$resStartActual && $ind<=$resEndActual)){
        # a repeat of the print logic at top of foreach loopp: this one for the last residue we collected
        # coordinates for...
        if($aa==0){
            $xa=$ya=$za=999.0;
        }
        if($bb==0){
            $xb=$yb=$zb=999.0;
        }
        if($gg==0){
            $xg=$yg=$zg=999.0;
        }
        else{
            $xg=$xg/$gg;
            $yg=$yg/$gg;
            $zg=$zg/$gg;
        }
        printf $FHANDLE "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",$xg,$yg,$zg,$xa,$ya,$za,$xb,$yb,$zb;
        $count++;
            # count the residues we are printing coordinate information for.  Check to see if we have come to a break in the ATOM
            # structure according to our alignment with the SEQRES, and if so, insert dummy coordinates for the inserted residues.
            # Note that since this is the last ATOM record we were going to print, we should not find any more gaps in the case
            # that we filled breaks for "interior residues" only, but we place the while() below to handle the case of filling
            # a gap at the end of the ATOM residues if desired.
        while( $alignment{query}[ $count ] eq '-' ) {
                # [ 0 ] is the entry for the ATOM residue - it will be 'GAP' in cases where the residue was missing in the ATOM records.
                # This is the case in which we want to insert dummy coords.
                printf $FHANDLE "%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",999.0,999.0,999.0,999.0,999.0,999.0,999.0,999.0,999.0;
                $count++;
        }
    }
    
    if ($FHANDLE eq 'OUT') {
        close(OUT) || carp("Can't close $fname. $!\n");
    }
}

###################################################################################################################################

sub trim($) {
    # remove leading and trailing white-space characters
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub isCaOnly() {
  return $H_caOnly;
}

sub isEngineered() {
  return $H_engineered;
}

sub amino3to1($) {
  my($i) = @_;
  my $j="?";
  
  if($i eq "ALA"){ $j="A";}

  elsif($i eq "CYS"){ $j="C";}
  elsif($i eq "ABA"){ $j="C";}
  elsif($i eq "CSO"){ $j="C";}
  elsif($i eq "CSD"){ $j="C";}
  elsif($i eq "CYM"){ $j="C";}
  elsif($i eq "CME"){ $j="C";}
  elsif($i eq "CSX"){ $j="C";}
  elsif($i eq "CSE"){ $j="C";}
  
  elsif($i eq "ASP"){ $j="D";}
  elsif($i eq "IAS"){ $j="D";}
  elsif($i eq "ASX"){ $j="D";}
  
  elsif($i eq "GLU"){ $j="E";}
  elsif($i eq "GLX"){ $j="E";}
  
  elsif($i eq "PHE"){ $j="F";}
  
  elsif($i eq "GLY"){ $j="G";}
  elsif($i eq "GLZ"){ $j="G";}
  
  elsif($i eq "HIS"){ $j="H";}
  elsif($i eq "HIP"){ $j="H";}
  elsif($i eq "HID"){ $j="H";}
  elsif($i eq "HIE"){ $j="H";}
  elsif($i eq "HSD"){ $j="H";}
  elsif($i eq "HSE"){ $j="H";}
  elsif($i eq "HSP"){ $j="H";}
  elsif($i eq "DDE"){ $j="H";}
  
  elsif($i eq "ILE"){ $j="I";}
 
  elsif($i eq "LEU"){ $j="L";}

  elsif($i eq "LYS"){ $j="K";}
  elsif($i eq "MLY"){ $j="K";}

  elsif($i eq "MET"){ $j="M";}
  elsif($i eq "MSE"){ $j="M";}
  elsif($i eq "MHO"){ $j="M";}
  
  elsif($i eq "PRO"){ $j="P";}
  
  elsif($i eq "ARG"){ $j="R";}
  elsif($i eq "ASN"){ $j="N";}
  elsif($i eq "GLN"){ $j="Q";}

  elsif($i eq "SER"){ $j="S";}
  elsif($i eq "SEP"){ $j="S";}

  elsif($i eq "THR"){ $j="T";}
  elsif($i eq "TPO"){ $j="T";}

  elsif($i eq "TRP"){ $j="W";}
  elsif($i eq "TRN"){ $j="W";}
  
  elsif($i eq "TYR"){ $j="Y";}
  
  elsif($i eq "VAL"){ $j="V";}
  elsif($i eq "MVA"){ $j="V";}
  
  elsif($i eq "GAP"){$j="-";}
    
  return $j;
}

###################################################################################################################################
sub amino1to3($) {
  my($j)= @_;
  my $i="UNK";
  if( $j eq "A" ){ $i="ALA";}
  elsif( $j eq "C" ){ $i="CYS";}
  elsif( $j eq "D" ){ $i="ASP";}
  elsif( $j eq "E" ){ $i="GLU";}
  elsif( $j eq "F" ){ $i="PHE";}
  elsif( $j eq "G" ){ $i="GLY";}
  elsif( $j eq "H" ){ $i="HIS";}
  elsif( $j eq "I" ){ $i="ILE";}
  elsif( $j eq "K" ){ $i="LYS";}
  elsif( $j eq "L" ){ $i="LEU";}
  elsif( $j eq "M" ){ $i="MET";}
  elsif( $j eq "N" ){ $i="ASN";}
  elsif( $j eq "P" ){ $i="PRO";}
  elsif( $j eq "Q" ){ $i="GLN";}
  elsif( $j eq "R" ){ $i="ARG";}
  elsif( $j eq "S" ){ $i="SER";}
  elsif( $j eq "T" ){ $i="THR";}
  elsif( $j eq "V" ){ $i="VAL";}
  elsif( $j eq "W" ){ $i="TRP";}
  elsif( $j eq "Y" ){ $i="TYR";}
  elsif( $j eq "-" ){ $i="GAP";}
  return $i;
}

###################################################################################################################################
sub checkRes($){
  my ( $res)=@_;
  #my (%rule,@II,$r);
  #%rule= ( 
  #"ALA"=>1,"ARG"=>2,"ASN"=>3,"ASP"=>4,"CYS"=>5,
  #"GLN"=>6,"GLU"=>7,"GLY"=>8,"HIS"=>9,"ILE"=>10,
  #"LEU"=>11,"LYS"=>12,"MET"=>13,"PHE"=>14,"PRO"=>15,
  #"SER"=>16,"THR"=>17,"TRP"=>18,"TYR"=>19,"VAL"=>20,
  #"GLX"=>7,"ASX"=>3,"UNK"=>0,"PCA"=>15);
  # 
  # @II=("UNK","ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL");
  # $r= $II[$rule{$res}];
  # return $r;
  
  my $result = amino1to3( amino3to1( $res ) );
    # this has the effect of standardizing the residue name, or resulting in "UNK" for any unkown residue

  if( $result eq "UNK" ) {
    # try mapping via translations as found in the MODRES records;
    my $modres = $M_translations{ $res };
      # note: really we should be using $this{MODRES}{translations} passed to this sub if we want true obj orient
    if( $modres ) {
      $result = $modres;
      print "Translating via mapping found in MODRES record: $res -> $result\n";
    }
  }
  
  return $result;
}

###################################################################################################################################
sub getAtomSeqForChain($$$$) {
    my ( $this, $chain, $startRes, $endRes ) = @_;
    my @seq;
    my @modified_index;
    
    foreach (0 .. $this->{ATOM}{sizeOf} ) {
        my $res = $this->{ATOM}{residue_name}[$_];
        $res = checkRes($res);
        if( $this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chain &&
            $this->{ATOM}{modified_index}[$_] >= $startRes && $this->{ATOM}{modified_index}[$_] <= $endRes ) {
            push @seq, $res;
            push @modified_index, $this->{ATOM}{modified_index}[$_];
        }
    }

    return ( \@seq, \@modified_index );
}

###################################################################################################################################
sub getSeqresSeqForChain($$$) {
    my ( $this, $chain, $standardize ) = @_;
    my @seq;
    
    foreach ( 0 .. $this->{SEQRES}{sizeOf} ) {
        if( $this->{SEQRES}{chain}[$_] eq $chain ) {
            my $res =  $this->{SEQRES}{residue_name}[$_];
            if( $standardize ) {
                $res = checkRes( $res );
            }
            if( $res ne 'UNK' || $includeUNK ) {
                push @seq, $res;
            }
        }
    }
    
    return \@seq;
}

###################################################################################################################################
sub printFasta( $$$$ ) {
	my ($this, $name, $info, $seqRef ) = @_;
	my @seq = @$seqRef;
	
	open OUT, ">$name" or die "can't open $name";
	print OUT "$info\n";
	
	map { print OUT amino3to1( $_ ); } @seq;
	print OUT "\n";
}

###################################################################################################################################
#sub getAtomToSeqresAlginmentForChain($$$$$$) {
#    my ( $this, $chain, $atomSeqRef, $atomModifiedIndexRef, $seqresSeqRef, $interiorOnly ) = @_;
#    my @atomSeq = @$atomSeqRef;
#    my @atomModifiedIndex = @$atomModifiedIndexRef;
#    my @seqresSeq = @$seqresSeqRef;
#    my @alignment;
#    my $gapsExist = 0;
#    
#    # @atomSeq - residue sequence based on ATOM section, obtained via getAtomSeqForChain()
#    # @seqresSeq - residue sequence based on SEQRES section, obtained via getSeqresSeqForChain()
#    # $interiorOnly - if set, gaps at the ends of the atom sequence are ignored
#    #
#    # It is assumed that atomSeq is a subset of seqresSeq -- that is, atomSeq may contain gaps
#    # at the beginning, in the interior, and at the end, making the alignment process easier.
#    # The routine fails if this is not the case, and alignment will be returned empty.
#    
#    
#    my $atomSeqIndex = 0;
#    my $alignmentIndex = 0;
#    my $inInterior = 0;
#    my $atomModifiedResidueIndex = 0;
#        # we have the modified_index for any ATOM entries that exist.  We'll try assign indicies that
#        # fit with this for the inserted SEQRES residues by incrementing from the last known ATOM index
#        # that was found.
#    for( my $seqIdx=0; $seqIdx <= $#seqresSeq; $seqIdx++ ) {
#        my $seqResResidue = checkRes( $seqresSeq[ $seqIdx ] );
#        if( $atomSeqIndex <= $#atomSeq && $atomSeq[ $atomSeqIndex ] ne $seqResResidue ) {
#            # insert a gap for the alignment: atomSeq doesn't match here.  Whether or not
#            #  we actually insert a gap depends on if we are in the interior or want exterior gaps.
#            if( (! $interiorOnly) || $inInterior ) {
#                # insert a gap in the alignment
#                $atomModifiedResidueIndex++;
#                push @alignment, [ 'GAP', $seqResResidue, $atomModifiedResidueIndex, $seqIdx ];
#                    # within alignment, we can find the atom residue, and the sequence residue
#                $gapsExist++;
#            }
#            next;
#        }
#        if( $atomSeqIndex <= $#atomSeq ) {
#            # we are at a matched residue.  This also means from now on, we are in the interior
#            $atomModifiedResidueIndex = $atomModifiedIndex[ $atomSeqIndex ];
#            push @alignment, [ $atomSeq[ $atomSeqIndex ], $seqResResidue, $atomModifiedResidueIndex, $seqIdx ];
#                # these two are the same
#            $atomSeqIndex++;
#            $inInterior = 1;
#        }
#        else {
#            # this final case occurs when we have exhausted the $atomSeq, but there are still
#            # remaining $seqresSeq entries.  In this case, if we want exterior gaps, we
#            # insert a gap.
#            if( ! $interiorOnly ) {
#                $atomModifiedResidueIndex++;
#                push @alignment, [ 'GAP', $seqResResidue, $atomModifiedResidueIndex, $seqIdx ];
#                    # within alignment, we can find the atom residue, and the sequence residue
#                $gapsExist++;
#            }
#        }
#        
#    }
#    if( $atomSeqIndex != $#atomSeq + 1 ) {
#        print "ERROR: The alignment bewteen ATOM and SEQRES failed!  Please check your PDB file.\n";
#        print STDERR "ERROR: The alignment bewteen ATOM and SEQRES failed!  Please check your PDB file.\n";
#        $gapsExist = -$gapsExist;
#    }
#    return ( \@alignment, $gapsExist );
#}
#
####################################################################################################################################
#sub getAtomToSeqresAlginmentForChainUsingSegments($$$$$$) {
#    my ( $this, $chain, $atomSeqRef, $atomModifiedIndexRef, $seqresSeqRef, $interiorOnly ) = @_;
#    my @atomSeq = @$atomSeqRef;
#    my @atomModifiedIndex = @$atomModifiedIndexRef;
#    my @seqresSeq = @$seqresSeqRef;
#    my @alignment;
#    my $gapsExist = 0;
#    
#    # @atomSeq - residue sequence based on ATOM section, obtained via getAtomSeqForChain()
#    # @seqresSeq - residue sequence based on SEQRES section, obtained via getSeqresSeqForChain()
#    # $interiorOnly - if set, gaps at the ends of the atom sequence are ignored
#    #
#    # It is assumed that atomSeq is a subset of seqresSeq -- that is, atomSeq may contain gaps
#    # at the beginning, in the interior, and at the end, making the alignment process easier.
#    # The routine fails if this is not the case, and alignment will be returned empty.
#    
#    # for each segment -- one or more residues -- found in the ATOM records, locate that segment
#    # in the SEQRES sequence to produce an alignment.
#    
#    
#    
#    my $atomSeqIndex = 0;
#    my $alignmentIndex = 0;
#    my $inInterior = 0;
#    my $atomModifiedResidueIndex = 0;
#        # we have the modified_index for any ATOM entries that exist.  We'll try assign indicies that
#        # fit with this for the inserted SEQRES residues by incrementing from the last known ATOM index
#        # that was found.
#    for( my $seqIdx=0; $seqIdx <= $#seqresSeq; $seqIdx++ ) {
#        my $seqResResidue = checkRes( $seqresSeq[ $seqIdx ] );
#        if( $atomSeqIndex <= $#atomSeq && $atomSeq[ $atomSeqIndex ] ne $seqResResidue ) {
#            # insert a gap for the alignment: atomSeq doesn't match here.  Whether or not
#            #  we actually insert a gap depends on if we are in the interior or want exterior gaps.
#            if( (! $interiorOnly) || $inInterior ) {
#                # insert a gap in the alignment
#                $atomModifiedResidueIndex++;
#                push @alignment, [ 'GAP', $seqResResidue, $atomModifiedResidueIndex, $seqIdx ];
#                    # within alignment, we can find the atom residue, and the sequence residue
#                $gapsExist++;
#            }
#            next;
#        }
#        if( $atomSeqIndex <= $#atomSeq ) {
#            # we are at a matched residue.  This also means from now on, we are in the interior
#            $atomModifiedResidueIndex = $atomModifiedIndex[ $atomSeqIndex ];
#            push @alignment, [ $atomSeq[ $atomSeqIndex ], $seqResResidue, $atomModifiedResidueIndex, $seqIdx ];
#                # these two are the same
#            $atomSeqIndex++;
#            $inInterior = 1;
#        }
#        else {
#            # this final case occurs when we have exhausted the $atomSeq, but there are still
#            # remaining $seqresSeq entries.  In this case, if we want exterior gaps, we
#            # insert a gap.
#            if( ! $interiorOnly ) {
#                $atomModifiedResidueIndex++;
#                push @alignment, [ 'GAP', $seqResResidue, $atomModifiedResidueIndex, $seqIdx ];
#                    # within alignment, we can find the atom residue, and the sequence residue
#                $gapsExist++;
#            }
#        }
#        
#    }
#    if( $atomSeqIndex != $#atomSeq + 1 ) {
#        print "ERROR: The alignment bewteen ATOM and SEQRES failed!  Please check your PDB file.\n";
#        print STDERR "ERROR: The alignment bewteen ATOM and SEQRES failed!  Please check your PDB file.\n";
#        $gapsExist = -$gapsExist;
#    }
#    return ( \@alignment, $gapsExist );
#}

###################################################################################################################################
sub bl2seqAlign( $$$$ ) {
	my ( $this, $seqRef1, $seqRef2, $blastEXE ) = @_;
	
	my @query;
	my $queryStart = -1;
	my $queryEnd = -1;
	my @subject;
	my $subjectStart = -1;
	my $subjectEnd = -1;
	my $gaps = 0;
    my $identities = 0;
	
	$this->printFasta( "seq1", ">seq1_ATOM", $seqRef1 );
	$this->printFasta( "seq2", ">seq2_SEQRES", $seqRef2 );
	
	system( "$blastEXE/bl2seq -F F -p blastp -i seq1 -j seq2 -o seq1seq2" );
	
	#my $seq1Beg = -1;
	#my $seq2Beg = -1;
	open ALI, "<seq1seq2" or die "Can't open alignment seq1seq2 in bl2seqAlign\n";
	while( <ALI> ) {
   		if( /.+Identities =\s*(\d+).+/ ) {
			$identities = $1;
		}
		if( /.+Gaps =\s*(\d+).+/ ) {
			$gaps = $1;
		}
		elsif( /Query: (\d+)\s+(\S+)\s+(\d+)/ ) {
			if( $queryStart == -1 ) {
				$queryStart = $1;
			}
			$queryEnd = $3;
			push @query, split( //, $2 );
		}
		elsif( /Sbjct: (\d+)\s+(\S+)\s+(\d+)/ ) {
			if( $subjectStart == -1 ) {
				$subjectStart = $1;
			}
			$subjectEnd = $3;
			push @subject, split( //, $2 );
		}
        elsif( /Score =.+/ ) {
            if( $subjectEnd != -1 ) {
                last;
                # we only want to look at the first alignment given, not subsequent chunks
            }
        }
	}
	close ALI;
    
    #system( "rm -f seq1 seq2 seq1seq2" );
    # I'll append these to the SEQALI file later for diagnostics
		
	my %alignment;
	$alignment{query} = \@query;
		# holds array of residues in the query, including '-' chars for gaps
	$alignment{queryStart} = $queryStart;
	$alignment{queryEnd} = $queryEnd;
		# the starting and ending residue of alignment, from original query sequence	
	$alignment{subject} = \@subject;
		# holds array of residues in the subject
	$alignment{subjectStart} = $subjectStart;
	$alignment{subjectEnd} = $subjectEnd;
		# the starting and ending residue of alignment, from original subject sequence
	$alignment{gaps} = $gaps;
    $alignment{identities} = $identities;

	return ( \%alignment, $gaps );
}

###################################################################################################################################

sub renumberResiduesUsingGappedAlignment( $ ) {
    # make the residue indices for the ATOM records match the SEQRES, assuming the SEQRES
    # starts at 1, and leaving numeric gaps where ATOM records were missing.
    
    # this is not used presently; in it's current form, it's only really useful to
    # renumber the residues of the first chain.  The alignment we are using is only
    # a single chain...
    
    my( $this ) = @_;
    my $alignmentRef = $this->{ALIGNATOMTOSEQRES};
    my %alignment    = %$alignmentRef;
	my $qRef         = $this->{ALIGNATOMTOSEQRES}{query};
	my @query		 = @$qRef;
    
    my $resIndex = 1;
    my $atomRecordIndex = 0;
    foreach ( @query ) {
        if( $_ eq "-" ) {
            $resIndex++;
            next;
        }
        my $resI = $this->{ATOM}{modified_index}[ $atomRecordIndex ];
        while( $this->{ATOM}{modified_index}[ $atomRecordIndex ] == $resI ) {
            $this->{ATOM}{modified_index}[ $atomRecordIndex ] = $resIndex;
            $atomRecordIndex++;
        }
        $resIndex++;
    }
}

sub printLOOPPSeqGapModifications( $$$ ) {
    my( $this, $name, $chain ) = @_;
    my $alignmentRef = $this->{ALIGNATOMTOSEQRES};
    my %alignment    = %$alignmentRef;
    my $gaps         = $this->{ALIGNATOMTOSEQRES}{gaps};
	my $qRef         = $this->{ALIGNATOMTOSEQRES}{query};
	my @query		 = @$qRef;
	my $sRef         = $this->{ALIGNATOMTOSEQRES}{subject};
	my @subject		 = @$sRef;
	
    
    File::Copy::move( "$name.seq", "$name.seqold" );
    
    open SEQ, ">$name.seq";
    open SEQGAP, ">$name.seqgap";
    
    my $count = $#query + 1;
    print SEQ "$name $count\n";
    print SEQGAP "$name $count\n";
    
	my $res;
    my $i = 0;
	map { $res = amino1to3( $_ ); print SEQ "$res\n"; $i++; } @subject;
        # note that 
	map { $res = amino1to3( $_ ); print SEQGAP "$res\n"; } @query;
		# this version actually shows the word "GAP" where there were gaps filled
    
    close SEQ;
    close SEQGAP;
    
    # diagnostic: print an alignment for easy visual comparison
    open SEQALI, ">$name.seqali";
	my $subjectText = join "", @subject;
	my $queryText = join "", @query;
    if( index( $subjectText, "-" ) >= 0 ) {
		# a  gap should not be required in the subject, since the subject is SEQRES, and query was SEQ
        print SEQALI "ERROR in alignment.  Probably due to a MODRES record, the ATOM record contains a residue the SEQRES does not.\n";
        print SEQALI "The GAP in the SEQRES resulting from the alignment was filled with the ATOM records for that resiude.\n";
        print "GAP found in SEQRES after alignment to SEQ: \n  Probably due to a MODRES record, the ATOM record contains a residue the SEQRES does not.\n";
    }
	if( $alignment{errorMessage} ) {
		# some alignment error occured
		print SEQALI $alignment{errorMessage};
	}
    print SEQALI "$name gaps: $gaps\n";
    
    print SEQALI "ATOM   : ";
    my ($seqRef, $idxRef) = $this->getAtomSeqForChain( $chain, -9999, 990000 );
        # the large numbers just ensure we catch all residues from the chain
    map { print SEQALI amino3to1($_); } @$seqRef;
    print SEQALI "\n\n";
    
    print SEQALI "SEQRES : ";
    $seqRef = $this->getSeqresSeqForChain( $chain );
    map { print SEQALI amino3to1($_); } @$seqRef;
    print SEQALI "\n\n";
    
    print SEQALI "ALIGNED:\n         ";
	print SEQALI $queryText;
    print SEQALI "\n         ";
	print SEQALI $subjectText;
    print SEQALI "\n\n============================================================================\n\n";
    
    if( $this->{alignSegsToSeqres} ) {
        $this->printAtomToSeqresAlignment( $this->{alignSegsToSeqres}, SEQALI );
    }
    elsif( -e "seq1seq2" ) {
        open S, "<seq1";
        my @s = <S>;
        print SEQALI @s; print SEQALI "\n\n";
        close S;

        open S, "<seq2";
        @s = <S>;
        print SEQALI @s; print SEQALI "\n\n";
        close S;

        open S, "<seq1seq2";
        @s = <S>;
        print SEQALI @s; print SEQALI "\n\n";
        close S;
    }
    else {
        print SEQALI "\n\nERROR: no alignment found!\n\n";
    }
    
    close SEQALI;
}

###################################################################################################################################
sub alignATOMToSEQRES( $$$$$$ ) {
    my ( $this, $chain, $aliStart, $aliEnd, $queryLen, $blastExe ) = @_;

    my( $totalResCount, $resStartActual, $resEndActual, $resCountActual ) =
        getResidueRangeFromAlignmentAndQueryLength( $this, $chain, $aliStart, $aliEnd, $queryLen );

    my ( $atomSeqRef, $atomModifiedIndexRef ) = $this->getAtomSeqForChain( $chain, $resStartActual, $resEndActual );
   
    my $seqresSeqRef = $this->getSeqresSeqForChain( $chain, 1 );
	my ( $alignRef, $gaps ) = $this->bl2seqAlign( $atomSeqRef, $seqresSeqRef, $blastExe );
    
    if( $alignRef->{queryStart} != 1 || $alignRef->{queryEnd} != scalar @$atomSeqRef  ||
        $alignRef->{identities} != scalar @$atomSeqRef ) {
        # the bl2seq can fail for our purposes for a couple of reasons:
        # (a) the alignment did not cover the entire range, usually because there are
        #     really large gaps or gaps near the end which would lower the alignment score.
        # (b) it aligned some non-identical residues
        #
        # so try the segment-based alignment, and if it succeeds, use it instead.
        my $segalignRef = $this->alignATOMSegsToSeqres( $chain, $resStartActual, $resEndActual );
        if( $segalignRef ) {
            #$this->printAtomToSeqresAlignment( $segalignRef );
                # debug
            unlink( seq1 );
            unlink( seq2 );
            unlink( seq1seq2 );
            ($alignRef, $gaps ) = $this->convertSegsAlignToBL2SeqFormat( $segalignRef );
            $this->{alignSegsToSeqres} = $segalignRef;
        }
		else {
			# we're in trouble.  The segment-based alignment failed, usually because there is an
			# inconsistency between the SEQRES and ATOM records.  Make a note of the
			# problem to be reported in the seqres file.
			my $alignError = "ERROR in alignment:\n";
			if( $alignRef->{queryStart} != 1 || $alignRef->{queryEnd} != scalar @$atomSeqRef ) {
				$alignError .= " (1) bl2seq failed because the alignment did not cover the entire query.\n";
			}
			if( $alignRef->{identities} != scalar @$atomSeqRef ) {
				$alignError .= " (2) bl2seq failed because some non-identical residues were aligned.\n";
			}
			$alignError .= " (3) alignATOMSegsToSeqres() failed because some ATOM residues are different from SEQRES.\n";
			$alignError .= "\nConsider not using this chain for the LOOPP database, or fixing the original pdb file.\n\n";
			$alignRef->{errorMessage} = $alignError;
		}
    }
    
    $alignRef->{modifiedAtomIndices} = $atomModifiedIndexRef;
    $this->{ALIGNATOMTOSEQRES} = $alignRef;
   
    return $gaps;
}

###################################################################################################################################
sub getSableResultsForChain( $$$$$$ ) {
    # load the predicted structure and solvent-exposed area per residue from sable results
    my ( $this, $pdbcode, $chain, $sableCache, $sableExe, $convertLoopp ) = @_;
	my $name = $pdbcode . "_$chain";
    
    # for now we will count on the fact that these cached results for all of our database already exist.
    my $sableResults = "$sableCache/$name.2nd.srf";
	
	print "\n------------ Getting SABLE results from $sableResults\n";

    if( ! -e $sableResults ) {
		my $input = "$name.seqres";
		system( "$sableExe $input");
		system( "$convertLoopp $name" );
		if( -e "$name.2nd.srf" ) {
			File::Copy::move( "OUT_SABLE_RES", "$sableCache/$name.2nd" );
			File::Copy::move( "OUT_SABLE_graph", "$sableCache/$name.2nd.graph" );
			File::Copy::move( "$name.2nd.srf", "$sableCache/$name.2nd.srf" );
			File::Copy::move( "$name.2nd.trans", "$sableCache/$name.2nd.trans" );
			chmod 0660, "$sableCache/$name.2nd";
			chmod 0660, "$sableCache/$name.2nd.graph";
			chmod 0660, "$sableCache/$name.2nd.srf";
			chmod 0660, "$sableCache/$name.2nd.trans";
			system( "rm -rf OUT_SABLE");
		}
    }
    
    my ( @predictedStructure, @predictedRSA );
	my $seqChainStart = 0;
	for( 0..$this->{SEQRES}{sizeOf} ) {
		if( $this->{SEQRES}{chain}[$_] eq $chain ) {
			$seqChainStart = $_;
			last;
		}
	}

    if( -e  $sableResults ) {
        open SABLE, "<$sableResults";
        my @sableLines = <SABLE>;
		print $sableLines[1];
		print $sableLines[2];
		print $sableLines[4];
        close SABLE;
		
		my $length = length( $sableLines[1] ) - 1;
			# -1 is to ignore the newline at end
        
        # read the sableText, matching with the seqres sequence we already have.  If it does not
        # match, this is an error.
            # it is *possible* that the seqres we print via this information is shorter if there were unknown residues
            # encountered, but this should not happen, and does not in the 2010 database.
        my ( $res, $struc, $rsa );
        for( my $i=0; $i<$length; $i++ ) {
            $res    = substr( $sableLines[1], $i, 1 );
            $struc  = substr( $sableLines[2], $i, 1 );
            $rsa    = substr( $sableLines[4], $i, 1 );
            my $seqRes = amino3to1( $this->{SEQRES}{residue_name}[ $seqChainStart + $i ] );
            if( $res ne $seqRes ) {
                print "ERROR: (getSableResultsForChain): the seqres sequences differs at index $i.  Expected $seqRes, found $res ($struc, $rsa)\n";
                print STDERR "ERROR: (getSableResultsForChain): the seqres sequences differs at index $i.  Expected $seqRes, found $res ($struc, $rsa)\n";
                return 0;
            }
            push @predictedStructure, $struc;
            push @predictedRSA, $rsa;
				# so these arrays ultimately contain one structure prediction per residue found in the seqres
        }
    }
    else {
        print "ERROR: could not find or generated sable results $sableResults\n";
        return 0;
    }

    $this->{SABLE}{structure} = \@predictedStructure;
    $this->{SABLE}{rsa} = \@predictedRSA;
    return 1;
}

###################################################################################################################################
# following values from http://www.imb-jena.de/IMAGE_AA.html
my %aminoAcidSurfaceArea = (
  A => 115.0,
  R => 225.0,
  D => 150.0,
  N => 160.0,
  C => 135.0,
  E => 190.0,
  Q => 180.0,
  G => 75.0,
  H => 195.0,
  I => 175.0,
  L => 170.0,
  K => 200.0,
  M => 185.0,
  F => 210.0,
  P => 145.0,
  S => 115.0,
  T => 140.0,
  W => 255.0,
  Y => 230.0,
  V => 155.0
);
sub getDsspResultsForChain( $$$$$ ) {
    my ( $this, $pdbcode, $chain, $dsspCache, $dsspExe ) = @_;
    my $dsspFile = "$pdbcode" . "_$chain.dssp";
	
	print "\n------------ Getting dssp results from $dsspFile.\n";
    
    if( ! -e "$dsspCache/$dsspFile" ) {
        my $entFile = "pdb" . "$pdbcode" . "_$chain.ent.new";
        if( -e $entFile ) {
            system( "$dsspExe $entFile $dsspFile 2> $dsspFile.err" );
			chmod 0660, $dsspFile; 
			File::Copy::move( "$dsspFile", $dsspCache );
			if( -e "$dsspFile.err" && ( -s "$dsspFile.err" < 90 ) ) {
				unlink( "$dsspFile.err" )
			}
			else {
				chmod 0660, "$dsspFile.err";
				File::Copy::move( "$dsspFile.err", $dsspCache );
			}
        }
        else {
            print "Can't generate a .dssp file without a .ent.new file!\n";
            return 0;
        }
    }
    if( -e "$dsspCache/$dsspFile" ) {
        
        open IN, "<$dsspCache/$dsspFile";
        my @dssplines = <IN>;
        close IN;
        
        open OUT, ">$dsspFile";
            # the file with dssp results supplemented with sable predictions is the new .dssp file
            # the original dssp file will be now called .dssp.old if there is a difference
		
		my @dsspResidueNames = ();
		my @dsspStructTypes  = ();
		my @dsspRSAVals      = ();
			# these will be filled with info from the dssp file and placed into this->{DSSP}{<name>}
		
        my $lastResNum;
        my $foundBreaks=0;
        my $alignmentRef = $this->{ALIGNATOMTOSEQRES};
        my %alignment = %$alignmentRef;
        my $alignIndex = 0;
        my $dsspIndex = -9999;
		my $seqresAlignStart = $alignment{subjectStart};
		my $seqresEnd = $alignment{subjectEnd};    
			# where the alignment begins in the SEQRES
        for( my $i=0; $i <= $#dssplines; $i++  ) {
            my $dsspline = $dssplines[ $i ];
		#	print " $dsspline";
            if( substr( $dsspline, 13, 1 ) eq '!' ) {
				# This may be a break that corresponds to a missing residue(s) in the ATOM records of
				# the PDB file, in which case we have predicted secondary-structure information from
				# SABLE that we can use to fill in this break.
				#
				# HOWEVER, it also happens ( ~1% of chains ) that dssp will insert break residues (!)
				# when it finds the C-N distance is greater than 2.5A even when according to the PDB
				# there is no residue(s) missing at this location.  In this case we cannot fill the break.

				$foundBreaks = 1;
				
                my $breakAtResidue = $lastResNum + 1;
					# breakAtResidue is the residue from the cleaned pdb file we output, and thus this
					# residue number matches the residue numbers we have in the alignment.  If there is a
					# GAP indicator at this position in the alignment, we'll fill the break in the dssp file,
					# otherwise we won't even write the break to the new dssp file, keeping it in sync with
					# our other files like .log and .2nd
					
                while( $alignIndex < $seqresEnd && $breakAtResidue >= $alignment{modifiedAtomIndices}->[ $alignIndex ] ) {
                        # find the location where the gap should be for this break, if it exists
                    $alignIndex++;
                }
				
				# note the condition - if we're in a gap of the alignment, we'll fill with info from sable,
				# otherwise we'll not print anything for this break.
                while( $alignIndex < $seqresEnd && $alignment{query}->[ $alignIndex ] eq '-') {
					my $seqresIndex     = $seqresAlignStart + $alignIndex - 1;
						# -1 converts from a 1-based residue number to a 0-based array index
					my $seqresResidue   = $alignment{subject}->[ $alignIndex ];
                    my $sablePredStruct = $this->{SABLE}{structure}->[ $seqresIndex ];
                    my $degreeOfBurial  = $this->{SABLE}{rsa}->[ $seqresIndex ];
                    my $aminoArea       = $aminoAcidSurfaceArea{ $seqresResidue };
                    my $sablePredRSA    = $degreeOfBurial / 9.0 * $aminoArea;
					my $resIndex        = $lastResNum + 1;
                    printf OUT "%5d %4d %s %s  %s  0 00     0   0  %3d     00,-0.0    -0,-0.0    -0,-0.0    -0,-0.0   0.000 000.0  00.0 -00.0 -00.0   00.0   00.0   00.0\n",
                        $dsspIndex, $resIndex, $chain, $seqresResidue, $sablePredStruct, int($sablePredRSA);
         #           printf "*%5d %4d %s %s  %s  0 00     0   0  %3d     00,-0.0    -0,-0.0    -0,-0.0    -0,-0.0   0.000 000.0  00.0 -00.0 -00.0   00.0   00.0   00.0\n",
         #              $dsspIndex, $resIndex, $chain, $seqresResidue, $sablePredStruct, int($sablePredRSA);

                    $dsspIndex++;
                    $alignIndex++;
					$lastResNum++;
                }
            }
            else {
                $lastResNum = substr( $dsspline, 6, 4 );
				my $output = $dsspline;
				if( $dsspIndex > 0 ) {
				   #   22   22 B V  H ><5S+     0   0   30     -4,-3.0     3,-0.9    -5,-0.2    -1,-0.2   0.945 114.6  46.2 -63.5 -46.6    9.4  -19.9   52.7
					$dsspline =~ /^\s*\d+\s+(\d+)\s+\S+\s+(\S+)..(.)..................\s*(\d+).+/;
						# we can always count on the first few fields, but after that we need to count columns to get at the struct and rsa
					my $dsspResIndex = $1;
						# dsspResIndex is the residue index from our cleaned PDB file (the one with a single chain)
					my $dsspResName = $2;
					my $dsspStruct  = $3;
					my $dsspRSA = $4;
					$dsspResidueNames[ $dsspResIndex ] = $dsspResName;
					$dsspStructTypes[ $dsspResIndex ]  = $dsspStruct;
					$dsspRSAVals[ $dsspResIndex ]      = $dsspRSA;
					
					$output = sprintf "%5d ", $dsspIndex;
					$output .= substr( $dsspline, 6 );
				}
				print OUT $output;
		#		print "*$output";
				if( substr( $dsspline, 2, 10 ) eq '#  RESIDUE' ) {
					$dsspIndex = 0;
						# the start of actual residue entries
				}
                $dsspIndex++;
            }
            
            
        }
        close OUT;
		
		chmod 0660, $dsspFile;
		$this->{DSSP}{resName} = \@dsspResidueNames;
		$this->{DSSP}{structure} = \@dsspStructTypes;
		$this->{DSSP}{rsa} = \@dsspRSAVals;
		
		if( $foundBreaks == 1 ) {
			File::Copy::copy( "$dsspCache/$dsspFile", "$dsspFile.old" );
				# if there were breaks we'll have the .dssp.old file parallel to make this obvious
		}
        return 1;
    }
    return 0;
}

###################################################################################################################################

############################################################################################
############################################################################################
# The following routines may not be used in the current loopp incarnation and are subject to
# removal.
############################################################################################
############################################################################################

sub PrintPdbAtomSeq($$){ 
  my($this, $fname, $chainL, $pdbc) = @_;

  my $FHANDLE; # File handle changer.
  if (!$fname =~ /\S/) {
    
    $FHANDLE = 'STDOUT';
  } else {
    open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
    $FHANDLE = 'OUT';
  }

  my $n=0;
  my $res;
  my $res1;
  print $FHANDLE ">$pdbc.atom\n";
  foreach (0 .. $this->{ATOM}{sizeOf}) {
   $res = $this->{ATOM}{residue_name}[$_];
   $res = checkRes($res);
   if($this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chainL){
    $n++;
    $res1 = amino3to1($res);
    print $FHANDLE "$res1";
    if($n%60==0){print $FHANDLE "\n";}
    }
   }
  close($FHANDLE);
 } 

#######################################################################################

sub PrintPdbAtomSeqPart($$){ 
  my($this, $fname, $chainL, $pdbc, $n1, $n2) = @_;

  my $FHANDLE; # File handle changer.
  if (!$fname =~ /\S/) {
    
    $FHANDLE = 'STDOUT';
  } else {
    open(OUT, ">$fname") || carp("Can't open $fname. $!\n");
    $FHANDLE = 'OUT';
  }

  my $nn=0;
  my $n=0;
  my $res;
  my $res1;
  print $FHANDLE ">$pdbc.atom\n";
  foreach (0 .. $this->{ATOM}{sizeOf}) {
   $res = $this->{ATOM}{residue_name}[$_];
   $res = checkRes($res);
   if($this->{ATOM}{atom_name}[$_] eq " CA "  && ($includeUNK || $res ne "UNK") && $this->{ATOM}{chain}[$_] eq $chainL){
    $nn++;
    if($nn>=$n1 && $nn<=$n2){
     $n++;
     $res1 = amino3to1($res);
     print $FHANDLE "$res1";
     if($n%60==0){print $FHANDLE "\n";}
     }
    }
   }
  close($FHANDLE);
 }

#######################################################################################

sub get_backbone($$$$){
   # map index residue of atombackbone to atom index
   my ($this,$atom,$backBone)=@_;
   my (%Seen,$i,$residue_index,$get_atom);

   
   for ($i=0; $i<=$this->{ATOM}{sizeOf}; $i++) {
	   $get_atom=$this->{ATOM}{atom_name}[$i];
	   $get_atom=~s/\s//g;
           if ($get_atom eq $atom ) {
                if (!$Seen{$this->{ATOM}{residue_index}[$i]}++) {
                        $residue_index=$this->{ATOM}{residue_index}[$i];
			%{$backBone}->{$residue_index}=$i;
                  }
           }
   }
}

#######################################################################################

sub define_Aseq_segments($$$$$$$){  
  my ($this, $A_C, $A_N, $A_CA, $segments,  $A_seq,$A_seq_index,$A_index,$R_seq)=@_;
  my (%Seen, $atom_index1, $atom_index2, $bond_len, $accepted, $index, $residue_index,$i,$n_A,$residue_index1);
  $index=0;
  
  
  atomSeq($this,$A_seq,$A_seq_index,$A_index); #get Atom section residues
  resSeq($this,$R_seq);                        #get SEQRES residues
  $n_A=$#$A_seq_index;
 
  
  # calculate segments in atom section
  for ($i=0; $i<$n_A-1; $i++) {
	  
    $accepted=0;
    $residue_index =@{$A_seq_index}->[$i];
    $residue_index1=@{$A_seq_index}->[$i+1];
    
    #print "$i =>$residue_index" . " AC=" . %{$A_C}->{$residue_index} . "  AN=" .  %{$A_N}->{$residue_index1} ."<br>";
    
    if (exists %{$A_C}->{$residue_index} && exists %{$A_N}->{$residue_index1} ) {
	    $atom_index1 = %{$A_C}->{$residue_index};
	    $atom_index2 = %{$A_N}->{$residue_index1};
	    $bond_len=dist($atom_index1,$atom_index2);
	    $accepted=$bond_len<2.25;
    }
    elsif (exists %{$A_CA}->{$residue_index} && exists %{$A_N}->{$residue_index1} ) {
	    $atom_index1 = %{$A_CA}->{$residue_index};
	    $atom_index2 = %{$A_N}->{$residue_index1}; 
	    $bond_len=dist($atom_index1,$atom_index2);
	    $accepted=$bond_len<7.29;
            print "CA-N $bond_len <br>";
    }
    elsif (exists %{$A_C}->{$residue_index} && exists %{$A_CA}->{$residue_index1} ) {
	    $atom_index1 = %{$A_C}->{$residue_index};
	    $atom_index2 = %{$A_CA}->{$residue_index1};
	    $bond_len=dist($atom_index1,$atom_index2);
	    $accepted=$bond_len<7.29;
            print "C-CA $bond_len <br>";
    }
    elsif (@{$A_seq_index}->[$i+1] - @{$A_seq_index}->[$i]==1){
	    $accepted=1;
            print "pdb  $i <br>";
    }
    
    if ($accepted) {push(@{$segments}, $index++);}
    else           {push(@{$segments}, 0);} 
    
   }
   
}

#######################################################################################

sub align_Aseq2Rseq($$$$$){
   my ($this, $R_seq, $A_seq, $segments,$A_index, $foundCoord, $n_coordDomain, $fname,)=@_;
   my ($r,$a,$res,$a0,$r0,$FHANDLE,$n_A,$n_R,$n_seg,$segs,$i_seg,$atom,$i,$j1,$j2);
   my $debug=0;
   my $count_seg=0;
   
   if (!$fname =~ /\S/) {
      # If $fname argument is omitted, print out to STDOUT.
      $FHANDLE = 'STDOUT';
   } 
   else {
      open(OUT, ">>$fname") || print ("Can't open $fname. $!<br>");
      $FHANDLE = 'OUT';
   }
  
  #print join("<br>",@{$A_seq});
  
  $n_A=$#$A_seq+1;
  $n_R=$#$R_seq+1;
  $n_seg=$#$segments+1;

  $segs=0;
  for($i=0;$i<$n_seg;$i++) {
	  if (${segments}->[$i]==0) {
		  $segs++;
	  }
  }
  $n_coordDomain=$segs;
  
  if ($debug){
     print "DEBUG IS ON <br> <br>";
     print "n_Aseq=$n_A n_Rseq=$n_R n_seg=$segs <br>";
  }

 

  $r=$r0=0; $a=$a0=0; $i_seg=0;
  while ($a0 < $n_A && $r0<$n_R){

     if ($debug) {
	     print "Start iteration with r=$r r0=$r0 a=$a a0=$a0 <br>";
     }

     
     while ((@{$A_seq}->[$a0] ne @{$R_seq}->[$r0]) && $r0 <$n_R) {
       if (0) {
	print "$r0 Atom<->SeqRes:......<font color='red'> @{$A_seq}->[$a0]</font> <==> <font color='blue'> @{$R_seq}->[$r0] </font><br>"; 
       }
       @{$foundCoord}->[$r0]=0; $r0++;  
     }
     $r=$r0;
     if ($debug) {
	     print "Recognize beging of segment in SEQRES  r=$r<br>";
     }


     
     while ((@{$A_seq}->[$a0] eq @{$R_seq}->[$r0]) && $r0 <$n_R) {
       if (0) {
	print "$r0 Atom<->SeqRes:......<font color='red'> @{$A_seq}->[$a0]</font> <==> <font color='blue'> @{$R_seq}->[$r0] </font><br>"; 
       }
       @{$foundCoord}->[$r0]=1; $a0++; $r0++;
     }
 	
     if ($debug) {
	     print "Recognize segment in SEQRES and ATOM  r0=$r0 a0=$a0 <br>";
     }
     
     if ($debug) {
      my $cond = @{$segments}->[$a0] == 0 || $a0 == $n_A || $r0 == $n_R;
      if ($cond){
	print "cond = $cond <br>";
        print "segment at $a0 is  @{$segments}->[$a0] <br>"; 
        print "SEQRES Starts at r=$r and ends at r0=$r0 <br>";
        print "ATOM   startsa at a=$a and ends at a0=$a0 <br>"; 
        for (my$Ares=$a, my$Sres=$r; $Ares<$a0;   $Ares++, $Sres++) {
	     print "$Ares:......<font color='red'> @{$A_seq}->[$Ares]</font> <==> <font color='blue'> @{$R_seq}->[$Sres] </font><br>";
       }
      }
     }
     

     if (@{$segments}->[$a0]==0 || $a0 == $n_A || $r0 == $n_R) {
	    $count_seg++;
	    for ($res=$a; $res<$a0; $res++) {
	      $j1 = @{$A_index}->[$res];
	      $j2 = @{$A_index}->[$res+1];
	      for($atom =$j1; $atom<$j2; $atom++){
		  
                  printf $FHANDLE "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
                  $this->{ATOM}{atom_index}[$atom],
	          $this->{ATOM}{atom_name}[$atom],
	          $this->{ATOM}{residue_name}[$atom],
       	          $this->{ATOM}{chain}[$atom],
                  $this->{ATOM}{modified_index}[$atom],
                  $this->{ATOM}{X_coordinate}[$atom],
                  $this->{ATOM}{Y_coordinate}[$atom],
	          $this->{ATOM}{Z_coordinate}[$atom],
	          $this->{ATOM}{occupancy}[$atom],
	          $this->{ATOM}{tempFactor}[$atom],
	          $this->{ATOM}{theRest}[$atom];
	      }
           }
	   $a=$a0; $r=$r0; $i_seg++;
     }
     else {  $r++; $r0=$r; $a0=$a }
  }

  if ($a0 < $n_A ) {
	  print "<font color='red'><H1> ERROR </H1></FONT>";
	  print "file=$this->{NAME}   chain=$this->{CHAIN} a0=$a0 n_A=$n_A<br>";
	  print "Failed to match seqres with atom section. Please check your PDB file<br>";
	  exit;
  }

  printf $FHANDLE "%-80s\n", 'END';

  if ($FHANDLE eq 'OUT') {
     close(OUT) || carp("Can't close $fname. $!\n");
  }

  if ($i_seg < $segs) {return 0;}
  return 1;
}

####################################################################################
sub atomSeq($$$$$){
  my($this, $seq, $res_index,$atom_index) = @_;
  my( %seen,$i,$j);

  
 
  $i=0;$j=0; %seen=();
  foreach (0 .. $this->{ATOM}{sizeOf}) {
	#print "$i $this->{ATOM}{residue_name}[$_] <br>";
      if (!$seen{$this->{ATOM}{modified_index}[$_]}++) {
	#print "$i $this->{ATOM}{residue_name}[$_] <br>";
	push(@{$seq},   $this->{ATOM}{residue_name}[$_]);
	push(@{$res_index}, $this->{ATOM}{residue_index}[$_]);
	push(@{$atom_index},$i);
	$j++;
      }
      $i++;
  }
  push(@{$atom_index},$this->{ATOM}{sizeOf}+1);
  #print "j=$j  $this->{ATOM}{sizeOf} <br>";
  
}

######################################################################################

sub resSeq($$$){
  my($this, $seq) = @_;
  my(@seq, %seen);

  foreach (0 .. $this->{SEQRES}{sizeOf}) {
	if ($this->{SEQRES}{residue_name}[$_] !~ /\s/ ){
	   push(@{$seq},$this->{SEQRES}{residue_name}[$_]);
        }
  }
}

######################################################################################

sub dist(){
	my ($p1,$p2)=@_;
    my ($dx,$dy,$dz,$r2);
	
	$dx = $A_X[$p1] - $A_X[$p2];
	$dy = $A_Y[$p1] - $A_Y[$p2];
    $dz = $A_Z[$p1] - $A_Z[$p2];

	$r2 = $dx*$dx + $dy*$dy + $dz*$dz;
	
	return ($r2);
}

#######################################################################################

sub printBackBonePDB(){
 my($this, $fname) = @_;

  my $FHANDLE; # File handle changer.
  if (!$fname =~ /\S/) {
    # If $fname argument is omitted, print out to STDOUT.
    $FHANDLE = 'STDOUT';
  } else {
    open(OUT, ">>$fname") || carp("Can't open $fname. $!\n");
    $FHANDLE = 'OUT';
  }

  foreach (0 .. $this->{ATOM}{sizeOf}) {
    if ($this->{ATOM}{atom_name}[$_] eq "C"  ||
	$this->{ATOM}{atom_name}[$_] eq "CA" ||
        $this->{ATOM}{atom_name}[$_] eq "N"  ||
        $this->{ATOM}{atom_name}[$_] eq "O" 
       ){
           printf $FHANDLE "ATOM  %5d %-4s %-3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f%s\n",
           $this->{ATOM}{atom_index}[$_],
	   $this->{ATOM}{atom_name}[$_],
	   $this->{ATOM}{residue_name}[$_],
	   $this->{ATOM}{chain}[$_],
           $this->{ATOM}{modified_index}[$_],
           $this->{ATOM}{X_coordinate}[$_],
           $this->{ATOM}{Y_coordinate}[$_],
	   $this->{ATOM}{Z_coordinate}[$_],
	   $this->{ATOM}{occupancy}[$_],
	   $this->{ATOM}{tempFactor}[$_],
	   $this->{ATOM}{theRest}[$_];
      }
  }
  printf $FHANDLE "%-80s\n", 'END';

  if ($FHANDLE eq 'OUT') {
    close(OUT) || carp("Can't close $fname. $!\n");
  }
}

#######################################################################################
1;


