#!/usr/bin/perl 
use lib '/org/centers/clsb/shruthi/loopp_scripts/';

# This script is currently not used directly by LOOPP.  I created it for use
# by persons wanting to "clean" pdb files in the way that they are processed
# by loopp when creating our databsae.  The output of this script is the "cleaned"
# version of the pdb for the given chain (or for all chains if no chain given)
# and corresponds to the format as found in the pdb_clean folder of loopp.
# (these files typically called <pdbXXXX.ent.new>

if( $#ARGV != 1 ) {
	print "\n  Usage: perl cleanpdb.pl <pdbcode> <folder with full pdb files>";
	print "\n     Ex: perl cleanpdb.pl 2PLS   /data/pdbs/     -- process all chains for 2PLS";
	print "\n     Ex: perl cleanpdb.pl 2PLS_A /data/pdbs/     -- process only chain A for 2PLS\n\n";
	exit;
}
	

$name=$ARGV[0];
$PDBdir = $ARGV[1]; 

use Read_PDB1_updateDB;
use File::Spec::Functions;

$pdbc  = substr ($name,0,4); 
$chain = " ";
my $fnameBase = "pdb$pdbc.ent";

if(substr ($name,4,1) eq "_"){
	$chain = substr( $name, 5, 1 );
	$fnameBase = "pdb$pdbc" . "_$chain.ent";
}

my $pdbFile= catfile( $PDBdir, "pdb$pdbc.ent" );
if( ! -e $pdbFile ) {
	print "$pdbFile does not exist, trying lowercase...\n";
	my $lc_pdbc = lc($pdbc);
	$pdbFile = catfile( $PDBdir, "pdb$lc_pdbc.ent" );
	if( ! -e $pdbFile ) {
		die "$pdbFile does not exist either!  Exiting.\n\n";
	}
}


$pdb = new Read_PDB( $pdbFile, $pdbc, $chain );
	# alloc and parse the pdb file
	
$pdb->logResidueNameTranslations( $pdbc, $chain, "$fnameBase.trans" );
	# dianostic: a file showing any translations to residue name from initial PDB entries to ours,
	# including those performed due to MODRES records, as well as standardization of amino-acid names.
	# This should be called *before* applyModRes, since once applied, the original residue name is lost.
	# Note that this only prints unique translations per chain, so it may be that multiple MSE->MET translations
	# occur in a given chain, but only one entry will be printed for brevity.

$pdb->applyModRes();
	# apply the MODRES records, which cause certain specific residues to be renamed (e.g. MSE->MET)
	
$pdb->logConvertedHETATM( $pdbc, "$fnameBase.hetatm" );
	# diagnostic: print any entries that *were* HETATM records, but that we converted to ATOM records
	# because we recognized the name of an amino-acid as the HETATM (e.g. MSE->MET)
	# This can be called before or after appplyModRes, with slightly different results.  If called
	# before, the records will indicate the name of the residue before any translation (e.g. MSE)
	# If called after, the records will indicate final residue name, after translation (e.g. MET)

if( $chain eq " " ) {
	$pdb->PrintPdbAtomAll( "$fnameBase.new" );
	$chain = "(all chains)";
}
else {
	$pdb->PrintPdbAtom( "$fnameBase.new", $chain, -9999, 990000, 999999 );
}

print "\n * $fnameBase.new contains the processed ATOMs for the specified chain(s).\n";
if( -e "$fnameBase.trans" ) {
	print " * $fnameBase.trans logged the name translations applied to some residues in this PDB.\n";
}
if( -e "$fnameBase.hetatm" ) {
	print " * $fnameBase.hetatm logged the residues that were HETATM in the PDB, but have been included in the .ent.new output.\n";
}
