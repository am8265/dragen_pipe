use strict;
use lib '/home/dh2880/.bin/';
use ARGH;
use Data::Dumper;
    my$cmd=q{cd ~/pipeline/dragen_genotyping ; export PYTHONPATH=/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness ; /nfs/goldstein/software/python2.7.7/bin/luigi --module create_ped_luigi_dragen CreatePed --local-scheduler --markers /nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/data/filtered.37MB.master.training.map };
while(1) {
# my@a=&ARGH::mq("select sample_name, d.pseudo_prepid, sample_type, upper(sample_type) STD, p.*, q.alignseqfileloc from dragen_sample_metadata d join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid join dragen_ped_status p on d.pseudo_prepid=p.pseudo_prepid where create_ped = 1 order by d.pseudo_prepid desc");
# print Dumper \@a;
# for my$x (@a) {
    my@a=&ARGH::mq("select sample_name, d.pseudo_prepid d_pp, p.pseudo_prepid p_pp, genotyping_rate, sample_type, upper(sample_type) STD, q.alignseqfileloc from dragen_sample_metadata d join dragen_qc_metrics q on d.pseudo_prepid=q.pseudo_prepid left join dragen_ped_status p on d.pseudo_prepid=p.pseudo_prepid where p.pseudo_prepid is null order by d_pp  desc limit 1");
    my$x=$a[0];
    $x||die;
    print Dumper $x;
    die if($x->{p_pp} ne 'NULL');
    # die if($x->{d_pp} eq 'NULL');
    # create_ped' => '1',
    #           'STD' => 'GENOME_AS_FAKE_EXOME',
    #                     'pseudo_prepid' => '400',
    #                               'predict' => '1',
    #                                         'sample_name' => 'b57mod159',
    #                                                   'sample_type' =>
    #                                                   'Genome_As_Fake_Exome',
    #                                                             'append_ped' => '0'
    #                                                                     };
    #
    my$pp=$x->{d_pp};
    # my$pp=$x->{pseudo_prepid};
    my$cmd=$cmd
      .' --pseudo-prepid '.$pp
      .' --sample-name '.$x->{sample_name}
      .' --sample-type '.$x->{sample_type}
      .' --output-directory /nfs/fastq_temp2/dsth/genotyping/'.$x->{sample_name}.'.'.$pp.'.'.$x->{sample_type}
      .' --align-loc '.$x->{alignseqfileloc}.'/'.$x->{sample_name}.'.'.$pp;
    print qq{$cmd\n};
    system($cmd) && die;
}
# --align-loc /nfs/archive/p2018/ALIGNMENT/BUILD37/DRAGEN/EXOME/aa10295.165375/ 


