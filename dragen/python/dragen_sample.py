#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Define a class describing samples going through the dragen pipeline
"""

import argparse
import MySQLdb
import os
import pdb
import subprocess
import sys
import time
from glob import glob

def get_prepid(curs, sample):
    """Retrieve qualifying prepids"""
    query = ("SELECT prepid FROM pseudo_prepid WHERE pseudo_prepid={0}"
            ).format(sample["pseudo_prepid"])
    #print query
    curs.execute(query)
    prepids = curs.fetchall()
    prepids = [x[0] for x in prepids]
    return prepids

def get_priority(curs,sample):
    query = ("SELECT priority "
            "FROM SampleT s "
            "JOIN prepT p ON s.CHGVID=p.CHGVID "
            "JOIN SeqType st ON st.prepid=p.prepid "
            "WHERE p.CHGVID='{sample_name}' "
            "AND st.seqtype='{sample_type}' "
            "AND exomekit='{capture_kit}' "
            "ORDER BY priority ASC "
            "LIMIT 1"
            ).format(sample_name=sample['sample_name'],
                    sample_type=sample['sample_type'],
                    capture_kit=sample['capture_kit'])
    #print query
    curs.execute(query)
    priority = curs.fetchone()
    return priority[0]

def get_pseudo_prepid(curs,sample):
    query = ("SELECT DISTINCT pseudo_prepid "
            "FROM pseudo_prepid pp "
            "JOIN prepT p ON pp.prepid=p.prepid "
            "JOIN SeqType s ON s.prepid=p.prepid "
            "WHERE CHGVID='{sample_name}' "
            "AND seqtype='{sample_type}' "
            "AND exomekit='{capture_kit}' "
            "AND failedprep=0"
            ).format(sample_name=sample['sample_name'],
                    sample_type=sample['sample_type'],
                    capture_kit=sample['capture_kit'])
    curs.execute(query)
    pseudo_prepid = curs.fetchone()
    return pseudo_prepid[0]

def get_bed_file_loc(curs,capture_kit):
    """Retrieve the bed file location for a given capture_kit name"""
    query = (("SELECT region_file_lsrc FROM captureKit WHERE name='{0}' "
        "and chr = 'all'"
        ).format(capture_kit))
    curs.execute(query)
    bed_file_loc = curs.fetchone()
    return bed_file_loc[0]

def get_fastq_loc(curs, sample):
    locs = []
    #correcting the downstream consequences of "custom capture" as the sample_type
    corrected_sample_type = sample['sample_type'].upper().replace(" ","_")
    #pdb.set_trace()
    for prepid in sample["prepid"]:
        query = ("SELECT seqsataloc,FCIllumID FROM Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            "AND l.prepID={0} AND f.fail=0 "
            "GROUP BY f.fcid"
            ).format(prepid)
        curs.execute(query)
        seqsatalocs = curs.fetchall()
        """For cases where there is not flowell information the sequeuncing
        will have to found be manually.  There will be two types of samples
        that under this catergory: Old Casava 1.8 sample (pre-seqDB) and
        Casava 1.7 samples sequenced by the Illumina GAII."""
        if seqsatalocs:
            """For externally submitted sample they have a fake Flowcell entry
            in the database.  The Flowcell.FCIllumID begins with 'X' always.
            Typically when processing the external sample each read group is
            enumerated 1,2,3...etc.

            Secondly when external samples are archived sometimes the FCIllumID
            is preserved otherwise its enumerated."""
            #print sample
            #for externally submitted samples
            if seqsatalocs[0][1][0] == 'X':
                print "Looking for external sample..."
                if 'SRR' in sample['sample_name']: #specifically for SRR samples
                    fastq_loc = glob(('/nfs/seqscratch10/SRA/{}/1'
                                ).format(sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif 'CGNDHDA' in sample['sample_name'] or 'FA000000' in sample['sample_name'] or 'NEU' in sample['sample_name']:
                    fastq_loc = glob(('/nfs/seqscratch09/tx_temp/tx_2390/CGND_11418-fastq/Project_CGND_11418_B01_GRM_WGS.2016-03-30/{}/1'
                                ).format(sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif 'pgm' in sample['sample_name'][0:3] or 'PGM' in sample['sample_name'][0:3]:
                    fastq_loc = glob(('/nfs/fastq1*/PGM/{}/[0-9]'
                                ).format(sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/igmdata[0-9]/{}/{}/[0-9]'
                    ).format(corrected_sample_type,sample['sample_name'])) != []: #for external igmdata samples with enumerated folders
                    print 1
                    fastq_loc = glob(('/nfs/fastq1[0-9]/{}/{}/[0-9]'
                                ).format(corrected_sample_type,sample['sample_name']))
                    if fastq_loc:
                        for flowcell in fastq_loc:
                            locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/igmdata[0-9]/{}/{}/*XX'
                    ).format(corrected_sample_type,sample['sample_name'])) != []: #for external igmdata samples with actual flowcell name
                    fastq_loc = glob(('/nfs/fastq1[0-9]/{}/{}/*XX'
                            ).format(corrected_sample_type,sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/fastq1[0-9]/{}/{}/[0-9]'
                    ).format(corrected_sample_type,sample['sample_name'])) != []: #for external fastq16 samples with enumerated folders
                    fastq_loc = glob(('/nfs/fastq1[0-9]/{}/{}/[0-9]'
                                ).format(corrected_sample_type,sample['sample_name']))
                    if fastq_loc:
                        for flowcell in fastq_loc:
                            locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/fastq1[0-9]/{}/{}/*XX'
                    ).format(corrected_sample_type,sample['sample_name'])) != []: #for external fastq16 samples with actual flowcell name
                    fastq_loc = glob(('/nfs/fastq1[0-9]/{}/{}/*XX'
                            ).format(corrected_sample_type,sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/[0-9]'
                    ).format(sample['sample_name'])) != []: # for external seqsata samples with enumerated folders
                    fastq_loc = glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/[0-9]'
                        ).format(sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                elif glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/*XX'
                    ).format(sample['sample_name'])) != []: #for external seqsata samples with actual flowcell name
                    fastq_loc = glob(('/nfs/seqsata*/seqfinal/whole_genome/{}/*XX'
                        ).format(sample['sample_name']))
                    for flowcell in fastq_loc:
                        locs.append(os.path.realpath(flowcell))
                else:
                    raise Exception, 'fastq files not found!'
            else: #For regular samples
                for flowcell in seqsatalocs:
                    if 'igmdata' in flowcell[0] or 'fastq' in flowcell[0]: # for igmdata## or fastq_temp##
                        fastq_loc = ('/nfs/{0}/{1}/{2}/{3}'
                                ).format(flowcell[0],corrected_sample_type,
                                        sample['sample_name'],flowcell[1])
                    elif 'seqsata' in flowcell[0]: # for seqsata## drives
                        fastq_loc = ('/nfs/{0}/seqfinal/whole_genome/{1}/{2}'
                                ).format(flowcell[0],sample['sample_name'],flowcell[1])
                    else:
                        raise Exception, "fastqloc does not within seqsata or fastq drive!"
                    read = glob(os.path.realpath(fastq_loc))
                    if read != []:
                        locs.append(os.path.realpath(fastq_loc))
                    else:
                        # For samples in the database but stored on the quantum and 
                        # have not had their location properly restored
                        fastq_loc = glob('/nfs/stornext/seqfinal/casava1.8/whole_{0}/{1}/*XX'.format(
                            corrected_sample_type.lower(),sample['sample_name']))
                        if fastq_loc:
                            for flowcell in fastq_loc:
                                locs.append(os.path.realpath(flowcell))

        else: #For samples with flowcells not in sequenceDB
            if glob('/nfs/fastq16/{}/{}/*XX'.format(corrected_sample_type,sample['sample_name'])):
                fastq_loc = glob('/nfs/fastq16/{}/{}/*XX'.format(corrected_sample_type,sample['sample_name']))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/fastq16/{}/{}/[0-9]'.format(corrected_sample_type,sample['sample_name'])) != []:
                fastq_loc = glob('/nfs/fastq16/{}/{}/[0-9]'.format(corrected_sample_type,sample['sample_name']))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/stornext/seqfinal/casava1.8/whole_{0}/{1}/*XX'.format(
                corrected_sample_type.lower(),sample['sample_name'])):
                fastq_loc = glob('/nfs/stornext/seqfinal/casava1.8/whole_{0}/{1}/*XX'.format(corrected_sample_type.lower(),sample['sample_name']))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            elif glob('/nfs/seqsata*/seqfinal/whole_genome/{}/*XX'.format(sample['sample_name'])) != []:
                fastq_loc = glob('/nfs/seqsata*/seqfinal/whole_genome/{}/*XX'.format(sample['sample_name']))
                for flowcell in fastq_loc:
                    locs.append(os.path.realpath(flowcell))
            else:
                raise Exception, "Sample {0} Fastq files not found!".format(sample['sample_name'])

    """For samples in the database we can exclude any samples that only have
    R1 data however for sampels that predate the database we have to manually
    check for R2 existance"""
    locs = check_fastq_locs(list(set(locs)))
    return locs

def check_fastq_locs(locs):
    """Determine if loc has read 2 fastq files"""

    valid_locs = []
    for loc in locs:
        read2 = glob("{loc}/*R2_[0-9]*fastq*".format(loc=loc))
        if read2 != []:
            valid_locs.append(loc)
        else:
            print 'Did not find fastq mate pair for: {}!'.format(loc)
    return valid_locs

def get_output_dir(curs,sample):
    """Generate ouput directory for Dragen results.  Dependent on seqtype"""

    # Custom capture samples need to be partitioned by capture_kit or 
    # pseudo_prepid since they are often sequenced with multiple capture kits.
    # Example: EpiMIR and SchizoEpi samples
    query = ("SELECT seqscratch_drive "
             "FROM dragen_sample_metadata "
             "WHERE pseudo_prepid = {} "
            ).format(sample['pseudo_prepid'])
    curs.execute(query)
    seqscratch_drive = curs.fetchone()
    if seqscratch_drive is None:
        seqscratch_drive = 'seqscratch_ssd'
    else:
        seqscratch_drive = seqscratch_drive[0]
    output_dir = ('/nfs/{3}/ALIGNMENT/BUILD37/DRAGEN/{0}/{1}.{2}/'
        ).format(sample['sample_type'].upper(),
                 sample['sample_name'],
                 sample['pseudo_prepid'],
                 seqscratch_drive)
    return output_dir

def get_lanes(curs,sample):
    """retrieve all qualifying lanes for the prepids associated with the sample"""
    lanes = []
    """For cases where there is not flowell information the sequeuncing
    will have to be manually.  There will be two types of samples that
    under this catergory: Old Casava 1.8 sample (pre-seqDB) and Casava 1.7
    samples sequenced by the Illumina GAII."""

    for prepid in sample['prepid']:
        """ Queries Lane table for lanes that do not have a failing lane for read1
        or read two since the dragen cannot use a mix of single and paired end reads
        """
        query = ("SELECT lanenum,FCIllumID from Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL and FailR2 IS NULL) "
            "AND l.prepID={0} AND p.failedPrep=0 "
            ).format(prepid)
        curs.execute(query)
        lane = curs.fetchall()
        if lane:
            lanes.append(lane)
        else:
            for flowcell in sample['fastq_loc']:
                fastqs = glob(flowcell + '/*fastq.gz')
                lane = []
                for fastq in fastqs:
                    lane_num = fastq[fastq.find('L00')+3]
                    lane.append(lane_num)
                lane_nums = set(sorted(lane))
                for lane_num in lane_nums:
                    lanes.append((lane_num,flowcell.split('/')[-1]))
            lanes = (lanes,)
    return lanes

class dragen_sample:
    # store all relevant information about a sample in a dictionary
    def __init__(self, sample_name, sample_type, pseudo_prepid, capture_kit, curs):
        self.metadata = {}
        self.metadata['sample_name'] = sample_name
        self.metadata['sample_type'] = sample_type.lower()
        self.metadata['pseudo_prepid'] = pseudo_prepid
        if sample_type.lower() != 'genome':
            self.metadata['capture_kit'] = capture_kit
            self.metadata['bed_file_loc'] = get_bed_file_loc(curs,self.metadata['capture_kit'])
        else:
            self.metadata['capture_kit'] = ''
            #Genome samples are set using the most current capture kit for any case which requires a target region.
            self.metadata['bed_file_loc'] = '/nfs/goldsteindata/refDB/captured_regions/Build37/65MB_build37/SeqCap_EZ_Exome_v3_capture.bed'
        self.metadata['prepid'] = get_prepid(curs,self.metadata)
        self.metadata['priority'] = get_priority(curs,self.metadata)
        self.metadata['fastq_loc'] = get_fastq_loc(curs,self.metadata)
        self.metadata['lane'] = get_lanes(curs,self.metadata)
        self.metadata['output_dir'] = get_output_dir(curs,self.metadata)
        self.metadata['script_dir'] = self.metadata['output_dir']+'scripts'
        self.metadata['fastq_dir'] = self.metadata['output_dir']+'fastq'
        self.metadata['log_dir'] = self.metadata['output_dir']+'logs'
        self.metadata['conf_file'] = "{script_dir}/{sample_name}.{pseudo_prepid}.dragen.conf".format(**self.metadata)
        self.metadata['dragen_stdout'] = "{log_dir}/{sample_name}.{pseudo_prepid}.dragen.out".format(**self.metadata)
        self.metadata['dragen_stderr'] = "{log_dir}/{sample_name}.{pseudo_prepid}.dragen.err".format(**self.metadata)

    def get_attribute(self, attribute):
        """return the value requested if present, otherwise raise a TypeError"""
        if attribute in self.metadata:
            return self.metadata[attribute]
        else:
            raise TypeError("{attribute} is not defined for {CHGVID}".format(
                attribute=attribute, CHGVID=self.CHGVID))

    def get_dict(self):
        """return a dict of the values of this sample"""
        return self.metadata

    def set(self, attribute, value):
        """set the specified attribute to the given value"""
        self.metadata[attribute] = value
