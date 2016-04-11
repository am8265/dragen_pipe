#!/nfs/goldstein/software/perlbrew/perls/perl-5.20.2/bin/perl -w

use strict;

use lib '/nfs/goldstein/goldsteinlab/software/bioperl-1.2.3';
use lib '/nfs/goldstein/software/ensembl_74/modules';
use lib '/nfs/goldstein/software/ensembl-variation_74/modules';
use lib '/nfs/goldstein/software/perlbrew/perls/perl-5.20.2/lib/site_perl/5.20.2';

use Bio::EnsEMBL::Registry;
use Digest::MD5 qw(md5_hex);
use Config::Simple;

my $cfg = new Config::Simple("/nfs/goldstein/software/dragen/.my.cnf");

my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(
	-host => $cfg->param("clienthomo_sapiens_variation_74_37.host"),
	-user => $cfg->param("clienthomo_sapiens_variation_74_37.user"),
	-pass => $cfg->param("clienthomo_sapiens_variation_74_37.password"),
	-db_version => '74'
);

my $transcript_adaptor = $reg->get_adaptor('Human', 'Core', 'Transcript');
my $dbh = DBI->connect("DBI:mysql:homo_sapiens_core_74_37;mysql_read_default_file=/nfs/goldstein/software/annodb/.my.cnf;mysql_read_default_group=clienthomo_sapiens_variation_74_37", undef, undef);
my $query = "SELECT stable_id FROM transcript";
my $sth = $dbh->prepare($query);
$sth->execute();
my $results = $sth->fetchall_arrayref;
$sth->finish();
my $dbh2 = DBI->connect("DBI:mysql:homo_sapiens_variation_74_37;mysql_read_default_file=/nfs/goldstein/software/annodb/.my.cnf;mysql_read_default_group=clienthomo_sapiens_variation_74_37", undef, undef);
my $md5_query = "SELECT translation_md5_id FROM translation_md5 WHERE translation_md5 = ?";
my $md5_sth = $dbh2->prepare($md5_query);
my $stable_id = "";
my $transcript;
my $tl;
my $seq;
my $md5;

foreach my $row (@{$results}) {
	$stable_id = $row->[0];
	$transcript = $transcript_adaptor->fetch_by_stable_id($stable_id);
	if(!$transcript) {
		print "${stable_id}\n";
		next;
	}
	$tl = $transcript->translation;
	if(!defined($tl)) {
		print STDERR "${stable_id}\n";
		next;
	}
	$seq = $tl->seq;
	$md5 = md5_hex($seq);
	$md5_sth->execute($md5);
	my $row_ref = $md5_sth->fetchrow_arrayref;
	if ($row_ref) {
		my $translation_md5_id = $row_ref->[0];
		print STDOUT "${stable_id}\t${translation_md5_id}\n";
	}
	else {
		print STDERR "${stable_id}\n";
	}
}
$dbh->disconnect();
$dbh2->disconnect();
