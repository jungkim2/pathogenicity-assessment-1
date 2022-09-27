## process ClinVar file (eg. clinvar_20220507.vcf)

$clinvar_file = "input/clinvar_20220507.vcf"; ## clinvar location

## process and store clinVar information with star anootations
print "location\tvariant\tid\tCLNREVSTAT\tCLNSIG\n";

open(FIL,$clinvar_file) || die ("Can't read file clinvar file '$clinvar_file'\n");
while (<FIL>) {

	chomp;
	next if ($_=~/^\#/); ## skip comment lines
	$_ =~ s/^chr//;
  my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info) = split "\t"; ## retrive column values

	## annotate one star
	# search CLNREVSTAT field
	if ($info =~ /(criteria_provided,_single_submitter)/)
  {
		my $CLNREVSTAT = $1;

		# retrieve CLNSIG value
		if($info =~/CLNSIG\=(\w+)\;/)
		{
			print "chr".$chr.":".$pos."\t".$ref.":".$alt,"\t";
			print $id,"\t".$CLNREVSTAT."\tOne Star\t",$1,"\n";
		}
	}

	## annotate two stars
	elsif ($info =~ /(criteria_provided,_multiple_submitters)/ )
	{
		my $CLNREVSTAT = $1;
		if($info =~/CLNSIG\=(\w+)\;/)
		{
			print "chr".$chr.":".$pos."\t".$ref.":".$alt,"\t";
			print $id,"\t".$CLNREVSTAT."\tTwo Stars\t",$1,"\n";
		}
	}

	## annotate three stars
	elsif ($info =~ /(reviewed_by_expert_panel)/ )
  {
		my $CLNREVSTAT = $1;

		if($info =~/CLNSIG\=(\w+)\;/)
		{
			print "chr".$chr.":".$pos."\t".$ref.":".$alt,"\t";
			print $id,"\t".$CLNREVSTAT."\tThree Stars\t",$1,"\n";
		}
	}

	## annotate four stars
	elsif ($info =~ /(practice_guideline)/ )
  {
		my $CLNREVSTAT = $1;

		if($info =~/CLNSIG\=(\w+)\;/)
		{
			print "chr".$chr.":".$pos."\t".$ref.":".$alt,"\t";
			print $id,"\t".$CLNREVSTAT."\tFour Stars\t",$1,"\n";

		}
	}

	## annotate calls that must be resolved
	elsif ($info =~ /(criteria_provided,_conflicting_interpretations)/ )
	{
		my $CLNREVSTAT = $1;
		if($info =~/CLNSIG\=(\w+)\;/)
		{
			print "chr".$chr.":".$pos."\t".$ref.":".$alt,"\t";
			print $id,"\t".$CLNREVSTAT."\tNeeds Resolving\t",$1,"\n";
		}
	}

	## annotate calls that have no stars
  else {
		my $CLNREVSTAT = "";
		if($info=~/CLNREVSTAT\=(\w+)\;/)
		{
			$CLNREVSTAT = $1;
		}
			if($info =~/CLNSIG\=(\w+)\;/)
			{
				print "chr".$chr.":".$pos."\t".$ref.":".$alt,"\t";
				print $id,"\t".$CLNREVSTAT."\tNo Stars\t",$1,"\n";
			}
  	}
}
