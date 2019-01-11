use strict;
use Data::Dumper;
use lib '/home/dh2880/.bin';
use ARGH;
use POSIX qw(strftime);

###### NEED TO OPTIMISE AND THOROUGH TESTSING!?!? BUT RIGHT NOW IMMEDIATE POST-SEQ QC INCL CONTAMINATION AND EXPERIMENTID ARE THE PRIORISTIES

###### do we use ATAV for forced calls/ped generation?!?
###### should we use a higher threshold for cryptic_relatedness?!?

###### currently ped file creation checks ped table - put ped file location into qc table and use that instead,
###### ethrnicity also uses ped table and should instead use qc genotyping_rate > 0 and relatedness should just check combined list itself 
###### once back log of ped is done update - HOWEVER, may be better to go with ATAV for ped file?!?
###### should speed this up!?!?
###### should speed this up!?!?
###### should speed this up!?!?
###### ONLY CHECK NEW SAMPLES EXCEPT FOR ONCE A WEEK FULL UPDATE!?! - i.e. AVOID ALL BUT 'Not Checked'
###### DON'T STORE THE STRINGS?!?

##### http://people.virginia.edu/~wc9c/KING/manual.html
#####   >0.354              duplicate/MZ twin, 
#####   [0.177, 0.354],     1st-degree, 
#####   [0.0884, 0.177]     2nd-degree, and 
#####   [0.0442, 0.0884]    3rd-degree relationships

###### wtf are dups doing in there?!?
# grep 'n/a ' /nfs/fastq_temp/dsth/relatedness/combined.ped -v > x && mv x /nfs/fastq_temp/dsth/relatedness/combined.ped
# sort -u /nfs/fastq_temp/dsth/relatedness/combined.ped > x && mv x /nfs/fastq_temp/dsth/relatedness/combined.ped

################## EITHER USE 0.125 - i.e. 3rd-gen?!?, 1.5 is almost utterly arbitrary - BUT - REALLY DOES CLEAN IT UP?!?
################## EITHER USE 0.125 - i.e. 3rd-gen?!?, 1.5 is almost utterly arbitrary - BUT - REALLY DOES CLEAN IT UP?!?
################## EITHER USE 0.125 - i.e. 3rd-gen?!?, 1.5 is almost utterly arbitrary - BUT - REALLY DOES CLEAN IT UP?!?
################## try using one of the other ped files - or both?!?
################## try using one of the other ped files - or both?!?
################## try using one of the other ped files - or both?!?
################## NEED GENUINE 'CRYPTIC' RELATIONSHIPS - i.e. LOOK AT VALUES FOR GREAT-GRANDPARENTS AND GO LOWER FOR THRESHOLD?!?
#
#
################## do we try running 10K with each of the other ped files?!?!?
################## do we try running 10K with each of the other ped files?!?!?
################## do we try running 10K with each of the other ped files?!?!?
#
#
my$cryptic_low=0.125;
my$cryptic_med=0.15;
my$dups_low=0.3; # bit low?!?

# sort -u /nfs/fastq_temp/dsth/relatedness/combined.ped | grep 'n/a ' > x && mv x
# /nfs/fastq_temp/dsth/relatedness/combined.ped
my$run_full=0;
my$pk=0; # my$pk=1;

sub parse_w {
    ###### we don't actually bother putting in phi (known kinship)
    my$l=shift;
    ### snp number?!?
    my@p=split("\t",$l);
    return (split(q{_},$p[1]),split(q{_},$p[2]),@p[0,8]);
}

sub parse_b {
    ###### clealry, between family so aren't supposed ot be related...?!?
    my$l=shift;
    ### snp number?!?
    my@p=split("\t",$l);
    # print Dumper \@p;
    return (split(q{_},$p[1]),$p[0],split(q{_},$p[3]),$p[2],$p[7]);
}

sub check_for_either {
    my($p,$pp1,$pp2,$l)=@_;
    for my$r (@{$l}) {
      return 1 if( $p->{$pp1}->{FamilyRelationProband} eq $r || $p->{$pp2}->{FamilyRelationProband} eq $r );
    }
    return 0;
}

sub check_pair {
    my($p,$pp1,$pp2,$r1,$r2)=@_;
    if(!defined($r2) && ( $p->{$pp1}->{FamilyRelationProband} eq $r1 && $p->{$pp2}->{FamilyRelationProband} eq $r1) ) {
        return 1;
    }elsif(
      ( $p->{$pp1}->{FamilyRelationProband} eq $r1 && $p->{$pp2}->{FamilyRelationProband} eq $r2 )
      || ( $p->{$pp2}->{FamilyRelationProband} eq $r1 && $p->{$pp1}->{FamilyRelationProband} eq $r2 ) ) {
        return 1;
    }
    return 0;
}

sub rel_prob {
    my($PROBS,$pp,$pp1,$pp2,$s1,$s2,$kin,$l)=@_;
    push(@{$PROBS->{$pp1}{$s1}},[ $pp->{$pp1}{FamilyRelationProband}.'_'.$pp->{$pp2}{FamilyRelationProband}, $s2,$kin,$pp2,$l]);
    push(@{$PROBS->{$pp2}{$s2}},[ $pp->{$pp2}{FamilyRelationProband}.'_'.$pp->{$pp1}{FamilyRelationProband}, $s1,$kin,$pp1,$l]);
}

######## for now we use create_ped=1 but once backlog is dealt with merge ped&predict and key off of genotyping_rate=null
######## then just pull all with genotyping_rate>0.0 except customcapture and run here
# my@list=&ARGH::mq("select * from dragen_sample_metadata d "
my@list=&ARGH::mq("select PercentContamination,relatedness_summary,relatedness_report,FamilyRelationProband,RepConFamMem,SelfDeclGender,FamilyID,d.sample_type,upper(d.sample_type) lazy,d.pseudo_prepid,experiment_id,st.sample_id,BroadPhenotype,is_external,pt.chgvid,d.sample_name,is_merged,pt.exomekit from dragen_sample_metadata d "
  ." join dragen_ped_status p on d.pseudo_prepid=p.pseudo_prepid "
  ." join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid "
  ." join prepT pt on d.pseudo_prepid=pt.p_prepid "
  ." join Experiment ex on pt.experiment_id=ex.id "
  ." join SampleT st on ex.sample_id=st.sample_id "
  ." where ( d.sample_type = 'Exome' or d.sample_type = 'Genome_As_Fake_Exome' ) and create_ped = 1 ");
  # ." where is_merged in (40,100) and ( pt.sample_type = 'Exome' or pt.sample_type = 'Genome_As_Fake_Exome' ) and create_ped = 1 and Species = 'Human'");
 
my%pp;
for my$p (@list) { ##### wasteful...
    $pp{$p->{pseudo_prepid}}=$p;
}

if($run_full){

# print Dumper \%pp;
# print Dumper \@list;
my$datestring = strftime("%Y_%b_%e_%H-%M", localtime());
# print $datestring;
my$cp='/nfs/fastq_temp/dsth/relatedness/combined.ped';
my$old_cp=substr($cp,0,length($cp)-4)."_${datestring}.ped";
if(!-e $old_cp) {
    system(qq{cp $cp $old_cp}) 
}else{
    die qq{there's already $old_cp\n};
}

open(my$f,'<',$cp)||die;
my%appended;
# sort -u /nfs/fastq_temp/dsth/relatedness/combined.ped!?!
while(my$l=<$f>){
    chomp($l);
    # should we check length?!?
    # Family  Subject  Father  Mother  Sex Phenotype
    my($fam,$expt,$pai,$mae,$genero,$phenotype)=split(q{ },$l);
    # print qq{family=$fam\nexpt=$expt\n}; # $pai\n$mae\n$genero\n$phenotype\n};
    my($chgvid,$pp)=split('_',$expt);
    # print qq{chgvid=$chgvid\npseudo_prepid=$pp\n}; # $pai\n$mae\n$genero\n$phenotype\n};
    # sort -u /nfs/fastq_temp/dsth/relatedness/combined.ped > x && mv x /nfs/fastq_temp/dsth/relatedness/combined.ped
    die qq{attempt sort -u else there's something nasty going on?!? : $fam, $expt, $chgvid, $pp} if(exists$appended{$pp});
    $appended{$pp}=[$chgvid,$fam,($chgvid eq $fam?1:0)];
}
# print Dumper \%appended;

for my$eligible (@list) {
    if(exists$appended{$eligible->{pseudo_prepid}}){
        # print qq{already added\n};
        #### do some checking
        # print Dumper $eligible; 
    }else{
        my$ped=q{/nfs/fastq_temp/dsth/genotyping/}.$eligible->{sample_name}.q{.}.$eligible->{pseudo_prepid}.q{.}.$eligible->{sample_type}.q{/}
          .$eligible->{sample_name}.q{.}.$eligible->{pseudo_prepid}.q{.ped};
        if(!-e $ped) {
            print qq{ped file $ped is missing\n};
            next;
        }
        chomp(my$wcw=`wc -w $ped | awk '{print \$1}'`);
        die qq{there's a problem with the ped file $ped\n} if($wcw!=25686);
        print qq{need to add $ped ($wcw)\n};
        print Dumper $eligible; 
        #### generate a new version
        system(qq{cat $ped >> $cp}) && die;
    }
}

# print Dumper \@list;
} #### run_full
if($pk){
my$cmd=q{/nfs/goldstein/software/PLINK/PLINK_1.90_3.38/plink --make-bed --file  /nfs/fastq_temp/dsth/relatedness/combined --out tmp};
system($cmd) && die;
$cmd=q{/nfs/goldstein/software/king_relatedness/king -b tmp.bed --kinship --related --degree 3 --prefix output};
system($cmd) && die;
}
# clearly remove any previous files if byte size is same as combined and exit completely...?!?
# clearly remove any previous files if byte size is same as combined and exit completely...?!?
# clearly remove any previous files if byte size is same as combined and exit completely...?!?

my%PROBS;
my%TOTAL;

######## for anything with ANY issue do way more checking that usual!?!
######## for anything with ANY issue do way more checking that usual!?!
######## for anything with ANY issue do way more checking that usual!?!

{ # really is virtually nothing to do here - just flag related pairs that haven't be reported as related via 'familyid' 
my$w_fm_pw='output.kin0';
open(my$wfh,'<',$w_fm_pw)||die;
while(my$l=<$wfh>){
    next if($.==1);
    chomp($l);
    my($chgvid1,$pp1,$fam1,$chgvid2,$pp2,$fam2,$kin)=&parse_b($l);

    die if($pp{$pp1}{chgvid} eq $pp{$pp2}{chgvid});

    # if($pp{$pp1}{FamilyID} ne 'N/A' && $pp{$pp2}{FamilyID} ne 'N/A') {
        # print Dumper $pp{$pp1};
        # print Dumper $pp{$pp2};
        # die;
    # }

    die $l if(!exists$pp{$pp1}||!exists$pp{$pp2});
    die $l if($pp{$pp1}{chgvid} ne $chgvid1);
    die $l if($pp{$pp2}{chgvid} ne $chgvid2);

    my$s1=$chgvid1.':'.$pp{$pp1}{FamilyID}.':'.$pp{$pp1}{sample_type}.':'.$pp{$pp1}{exomekit};
    my$s2=$chgvid2.':'.$pp{$pp2}{FamilyID}.':'.$pp{$pp2}{sample_type}.':'.$pp{$pp2}{exomekit};

    ##### should check not seen before?!?
    $TOTAL{$pp1}=$s1;
    $TOTAL{$pp2}=$s2;

    if($kin > $dups_low ) {

        #### v. silly - should always be 1<->1 at this stage?!?
        push(@{$PROBS{$pp1}{$s1}},['Unregistered_Duplicate',$s2,$kin,$pp2]);
        push(@{$PROBS{$pp2}{$s2}},['Unregistered_Duplicate',$s1,$kin,$pp1]);

        # print qq{appears to be a duplicate - despite family difference\n};
        # print qq{$chgvid1, $pp1, $fam1\n};
        # print qq{$chgvid2, $pp2, $fam2\n};
        # print qq{$kin\n};
        # print Dumper $pp{$pp1};
        # print Dumper $pp{$pp2};
       
    }elsif($kin > $cryptic_med ) {

        push(@{$PROBS{$pp1}{$s1}},['Cryptic_Relatedness',$s2,$kin,$pp2]);
        push(@{$PROBS{$pp2}{$s2}},['Cryptic_Relatedness',$s1,$kin,$pp1]);

    }elsif($kin > $cryptic_low ) {

        push(@{$PROBS{$pp1}{$s1}},['Cryptic_Relatedness_Low',$s2,$kin,$pp2]);
        push(@{$PROBS{$pp2}{$s2}},['Cryptic_Relatedness_Low',$s1,$kin,$pp1]);

        # print qq{appears to be related - despite family difference\n};
        # print qq{$chgvid1, $pp1, $fam1\n};
        # print qq{$chgvid2, $pp2, $fam2\n};
        # print qq{$kin\n};
        # print Dumper $pp{$pp1};
        # print Dumper $pp{$pp2};
    }
}
}

{

####### ignoring RepConFamMem : i.e. would otherwise check if both are 'repconfam' then issue 'repconfam' error if kinship is > 0.1?!? - i.e. they are NEVER supposed to be related!?!vb
####### we aren't allowing unknown relationships...?!?

###### pre-load all relevant info as eligible step?!?

###### within SAME family pairwise relations
my$w_fm_pw='output.kin';
open(my$wfh,'<',$w_fm_pw)||die;
while(my$l=<$wfh>){
    next if($.==1);
    chomp($l);
    my($chgvid1,$pp1,$chgvid2,$pp2,$fam,$kin)=&parse_w($l);
    # print qq{\nhave $chgvid1, $pp1, $chgvid2, $pp2, $fam, $kin\n};
    die $l if(!exists$pp{$pp1}||!exists$pp{$pp2});
    die $l if($pp{$pp1}{chgvid} ne $chgvid1);
    die $l if($pp{$pp2}{chgvid} ne $chgvid2);
    die if($pp{$pp1}{FamilyID} ne $pp{$pp2}{FamilyID});

    # next if(&check_pair(\%pp,$pp1,$pp2,'N/A'));

    my$s1=$chgvid1.':'.$pp{$pp1}{FamilyID}.':'.$pp{$pp1}{sample_type}.':'.$pp{$pp1}{exomekit};
    my$s2=$chgvid2.':'.$pp{$pp2}{FamilyID}.':'.$pp{$pp2}{sample_type}.':'.$pp{$pp2}{exomekit};

    $TOTAL{$pp1}=$s1;
    $TOTAL{$pp2}=$s2;

    if ($pp{$pp1}{sample_id}==$pp{$pp2}{sample_id}) {

        die if ($pp{$pp1}{FamilyRelationProband} ne $pp{$pp2}{FamilyRelationProband});

        # why the f' is there an upper-limit
        if( $kin < 0.45 || $kin > 0.5) {                
            # print qq{Withing Family Relationship error : doesn't appear to be the same sample\n};       
            push(@{$PROBS{$pp1}{$s1}},['Sample_Mismatch',$s2,$kin,$pp2]);
            push(@{$PROBS{$pp2}{$s2}},['Sample_Mismatch',$s1,$kin,$pp1]);
        }
        next;

    ######## after same sample...
    }elsif(&check_pair(\%pp,$pp1,$pp2,'Parent')) {

        # print Dumper $pp{$pp1};
        # kprint Dumper $pp{$pp2};
        ##### sanity check?!?
        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin>$cryptic_low);           next;

    }elsif ( ### allow for proband-proband here as it's clearly been entered incorrectly with different sample_id!?!?!
        ($pp{$pp1}{FamilyRelationProband} eq $pp{$pp2}{FamilyRelationProband}) || &check_pair(\%pp,$pp1,$pp2,'Proband','N/A')
    ) {

        if($kin>$dups_low) { ####### bothered giving same family but changed it's chgvid?!?
            # print qq{Withing Family Relationship error : this may be a sample duplicate\n};
            push(@{$PROBS{$pp1}{$s1}},['Unregistered_Duplicate',$s2,$kin,$pp2]);
            push(@{$PROBS{$pp2}{$s2}},['Unregistered_Duplicate',$s1,$kin,$pp1]);
        }elsif($kin<0.1) {
            # print qq{the strange thing below\n};
            &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__);
        }
        next;

    # }elsif(&check_pair(\%pp,$pp1,$pp2,'Proband')) {

    ######## this 'might' be excessive but it doesn't make sense not to apply half-sibling/sibling range?!?
    }elsif(
        &check_pair(\%pp,$pp1,$pp2,'Parent','Sibling') ||
        #### this should 'really' be Aunt-Uncle 0.088388348 0.176776695
        &check_pair(\%pp,$pp1,$pp2,'Sibling','Child') ||
        #### this should Grandparent
        &check_pair(\%pp,$pp1,$pp2,'Parent','Child')
    ) {

        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin < 0.088388348 || $kin > 0.353553391 );       next;
        # &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin < 0.088388348 || $kin > 0.176776695 );       next;
    }elsif(&check_pair(\%pp,$pp1,$pp2,'Parent','Sibling')) {

        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin < 0.088388348 || $kin > 0.353553391 );       next;

    ######## again 'might' be excessive but it doesn't make sense not to apply sibling
    }elsif(&check_pair(\%pp,$pp1,$pp2,'Sibling')) {

        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin < 0.176776695 || $kin > 0.353553391 );       next;
    
    ##### just doing this for sanity purposes atm?!? don't know lineage so no point?!? not much we can say?!?
    }elsif( &check_pair(\%pp,$pp1,$pp2,'Aunt-Uncle','Parent')       || &check_pair(\%pp,$pp1,$pp2,'Grandparent','Parent')                 ||
              &check_pair(\%pp,$pp1,$pp2,'Grandparent','Aunt-Uncle')  || &check_pair(\%pp,$pp1,$pp2,'Sibling','Great great aunt-uncle')   ||
              &check_pair(\%pp,$pp1,$pp2,'Sibling','Niece-Nephew')    || &check_pair(\%pp,$pp1,$pp2,'Sibling','Great aunt-uncle')         ||
              &check_pair(\%pp,$pp1,$pp2,'Parent','Great aunt-uncle') || &check_pair(\%pp,$pp1,$pp2,'Half sibling','Parent')              ||
              &check_pair(\%pp,$pp1,$pp2,'Half sibling','Grandparent') || 
              ### whatever?!? some we should 'perhaps' check?!?
              &check_pair(\%pp,$pp1,$pp2,'Sibling','First cousin') || 
              &check_pair(\%pp,$pp1,$pp2,'Grandparent','Great aunt-uncle') || 
              &check_pair(\%pp,$pp1,$pp2,'Great grandparent','Great aunt-uncle') || 
              &check_pair(\%pp,$pp1,$pp2,'Great grandparent','Grandparent') || 
              &check_pair(\%pp,$pp1,$pp2,'Great grandparent','Aunt-Uncle') || 
              &check_pair(\%pp,$pp1,$pp2,'Great grandparent','Parent') || 
              &check_pair(\%pp,$pp1,$pp2,'Great grandparent','Sibling') || 
              &check_pair(\%pp,$pp1,$pp2,'Grandparent','Sibling') || 
              &check_pair(\%pp,$pp1,$pp2,'Aunt-Uncle','Sibling') || 
              &check_pair(\%pp,$pp1,$pp2,'Aunt-Uncle','Great aunt-uncle') || 
              &check_pair(\%pp,$pp1,$pp2,'Great niece-nephew','Grandchild') || 
              &check_pair(\%pp,$pp1,$pp2,'Child','First cousin') ||  
              &check_pair(\%pp,$pp1,$pp2,'N/A') ##### leave this for sample_id matches above?!?
    ) { 
        next; 

    ####### don't allow for Proband-Proband-other or Proband-Proband as they SHOULD have same sample_id/chgvid!?!
    }elsif( 
        &check_pair(\%pp,$pp1,$pp2,'Monozygotic twin','Proband') || 
        &check_pair(\%pp,$pp1,$pp2,'Proband-other tissue','Proband')
    ) {

        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin < 0.45 || $kin > 0.5 );                      next;

    }elsif( $pp{$pp1}{FamilyRelationProband} ne 'Proband' && $pp{$pp2}{FamilyRelationProband} ne 'Proband' ) { 

        print Dumper $pp{$pp1};
        print Dumper $pp{$pp2};
        die; 
        
    ##### all rest require a proband - so this is wasteful!?!?!? should have been handled by silly bit above?!?
    }elsif( &check_for_either(\%pp,$pp1,$pp2,['First cousin once removed','First cousin twice removed','Fourth cousin','Great-great grandchild',
              'Great-great grandparent','Half first cousin','Second cousin','Second cousin once removed','Second cousin twice removed',
              'Spouse','Step child','Step parent','Step sibiling','Third cousin','Third cousin once removed'])
    ) {
        ##### this seems absurd?!?
        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if( $kin < -1000|| $kin > 0.1 );                     next;

    }elsif( &check_for_either(\%pp,$pp1,$pp2,['First cousin','Great aunt-uncle','Great grandchild','Great grandparent','Great great aunt-uncle',
              'Great great niece-nephew','Great niece-nephew'])
    ) {
        # don't allow completely unrelated?!?
        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if( $kin < 0.0 || $kin > 0.1 );                      next;

    }elsif( &check_for_either(\%pp,$pp1,$pp2,['Aunt-Uncle','Grandchild','Grandparent','Half sibling','Niece-Nephew']) ) {

        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if($kin < 0.088388348 || $kin > 0.176776695 );       next;

    }elsif( &check_for_either(\%pp,$pp1,$pp2,['Child','Consanguineous parent','Dizygotic twin','Parent','Sibling-other tissue','Sibling']) ) {

        &rel_prob(\%PROBS,\%pp,$pp1,$pp2,$s1,$s2,$kin,__LINE__) if( $kin < 0.176776695 || $kin > 0.353553391 );      next;

    }

    print Dumper $pp{$pp1};
    print Dumper $pp{$pp2};
    die qq{what're these?!?\n};

}
}

# http://www.perlmonks.org/?node_id=2461
my@OKAY=grep { !defined $PROBS{$_} } keys%TOTAL;

print q{THERE ARE }.scalar(@OKAY).qq{ SAMPLES WITHOUT WARNINGS\n};
print q{THERE ARE }.scalar(keys%PROBS).qq{ SAMPLES WITH WARNINGS\n};

open(my$ffh,'>','/nfs/fastq_temp/dsth/relatedness/full.txt')||die;

for my$K (keys%PROBS){
    die if(scalar(keys%{$PROBS{$K}})!=1);
    my($K2)=keys%{$PROBS{$K}};
    my%warnings;
    my$V=$PROBS{$K}{$K2};
    my$silly=0;
    for my$w (@{$V}){ 
        ++$silly;
        $warnings{$w->[0]}++;
    }

    print qq{HAVE $K : $K2 (warnings=$silly)\n\n};

    my$full=q{};
    my$summary=q{};
    my$short=q{};

    # if($silly>2){ # short just gets counts?!? do we order them at all?!?
    for my$k (sort {$b cmp $a}keys%warnings) {
        $summary.=$k.':'.$warnings{$k}.';';
    }
    #### do we still put in non-cryptic ones?!?
    # for my$w (@{$V}){ 
        # if($silly>2){ 
            # next if(substr($w->[0],0,length('Cryptic')) eq 'Cryptic');
        # }
        # $short.=$w->[0].':'.$w->[1].':'.$w->[2].';';
        # if(++$added>=3) {
            # $short.='...;';
            # last;
        # }
    # }
    # bored and lazy?!?
    my$added=0;
    my$lim=3;
    ### should sort and give ranges...
    for my$w (@{$V}){ 
        $full.=$w->[0].':'.$w->[1].':'.$w->[2].';';
        if(++$added<$lim) {
            $short.=$w->[0].':'.$w->[1].':'.$w->[2].';';
        }
        $short.='...;' if($added==$lim);
            # last;
        # }
    }
    # print qq{short=$short\nfull=$full\n};

    $short=~s{Cryptic_Relatedness}{CR}g;
    substr($summary,-1,1)=q{};;
    substr($full,-1,1)=q{};;
    substr($short,-1,1)=q{};;

    print qq{summary=$summary\nshort=$short\nfull=$full\n};

    if(
        $pp{$K}->{relatedness_summary} ne $summary ||
        $pp{$K}->{relatedness_report} ne $short
    ) {
        &ARGH::mu("update dragen_qc_metrics set relatedness_summary = '$summary', relatedness_report = '$short' where pseudo_prepid = $K");
    }


    print $ffh qq{$K\t$K2\t$silly\t}.$pp{$K}{PercentContamination}.qq{\t$summary\t$short\t$full\n};

    print Dumper $pp{$K};
    # print Dumper \%warnings;
    # print Dumper $V;

}

for my$K (@OKAY){
    print qq{No Warnings $K\n};
    print Dumper $pp{$K};
    if(
        $pp{$K}->{relatedness_summary} ne 'Pass' ||
        $pp{$K}->{relatedness_report} ne 'N/A' 
    ) {
        &ARGH::mu("update dragen_qc_metrics set relatedness_summary = 'Pass', relatedness_report = 'N/A' where pseudo_prepid = $K");
    }
}

print qq{bye\n};
# print Dumper \%PROBS;
