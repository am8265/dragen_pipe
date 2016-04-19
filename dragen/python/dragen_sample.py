#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Define a class describing samples going through the dragen pipeline
"""

import argparse
import MySQLdb
import os
import sys


def get_prepid(curs, sample, pseudo_prepid):
    query = ("SELECT prepid FROM pseudo_prepid where pseudo_prepid={0}"
            ).format(sample["pseudo_prepid"])

    curs.execute(query)
    prepids = curs.fetchall()
    return prepids[0]

def get_bed_file_loc(curs,capture_kit):
    query = (("SELECT region_file_lsrc FROM captureKit WHERE name='{0}' "
        "and chr = 'all'"
        ).format(capture_kit))
    curs.execute(query)
    bed_file_loc = curs.fetchone()
    return bed_file_loc[0]

def get_fastq_loc(curs, sample):
    fastq_locs = []
    for prepid in sample["prepid"]:
        query = ("SELECT seqsataloc,FCIllumID FROM Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL or FailR2 IS NULL) "
            "AND l.prepID={0} AND f.fail=0 "
            "GROUP BY f.fcid"
            ).format(prepid)
        curs.execute(query)
        seqsatalocs = curs.fetchall()

        locs = []
        for flowcell in seqsatalocs:
            if 'seqsata' in flowcell[0]:
                fastq_loc = ('/nfs/{0}/seqfinal/whole_genome/{1}/{2}'
                        ).format(flowcell[0],sample['sample_name'],flowcell[1])
            elif 'fastq' in flowcell[0]:
                fastq_loc = ('/nfs/{0}/{1}/{2}/{3}'
                        ).format(flowcell[0],sample['sample_type'].upper(),
                                sample['sample_name'],flowcell[1])
            else:
                raise Exception, "fastqloc does not within seqsata or fastq drive!"
            locs.append(os.path.realpath(fastq_loc))

        return locs
def get_output_dir(sample):
    """Custom capture samples need to be partitioned by capture_kit since they are often
       sequenced with multiple capture kits.  Example: EpiMIR and SchizoEpi
    """
    if sample['sample_type'] == 'custom_capture':
        output_dir = ('/nfs/fastq_temp/ALIGNMENT/BUILD37/DRAGEN/{0}/{1}/{2}'
            ).format(sample['sample_type'].upper(),
                    sample['capture_kit'],sample['sample_name'])
    else:
        output_dir = ('/nfs/fastq_temp/ALIGNMENT/BUILD37/DRAGEN/{0}/{1}'
            ).format(sample['sample_type'],sample['sample_name'])

    return output_dir

def get_lanes(curs,sample):
    lanes = []

    for prepid in sample['prepid']:
        """ Queries Lane table for lanes that do not have a failing lane for read1
        or read two since the dragen cannot use a mix of single and paired end reads
        """
        query = ("SELECT lanenum,FCIllumID from Flowcell f "
            "JOIN Lane l ON f.FCID=l.FCID "
            "JOIN prepT p ON l.prepID=p.prepID "
            "WHERE (FailR1 IS NULL or FailR2 IS NULL) "
            "AND l.prepID={0} AND p.failedPrep=0 "
            ).format(prepid)
        curs.execute(query)
        lane = curs.fetchall()
        lanes.append(lane)
    return lanes

#class prep(self, sample_name, curs):

class dragen_sample:
    # store all relevant information about a sample in a dictionary
    def __init__(self, sample_name, sample_type, pseudo_prepid, capture_kit, curs):
        self.metadata = {}
        self.metadata['sample_name'] = sample_name
        self.metadata['sample_type'] = sample_type
        self.metadata['pseudo_prepid'] = pseudo_prepid
        self.metadata['capture_kit'] = capture_kit
        if self.metadata['sample_type'] != 'genome':
            self.metadata['bed_file_loc'] = get_bed_file_loc(curs,self.metadata['capture_kit'])
        else:
            self.metadata['bed_file_loc'] = '/nfs/goldsteindata/refDB/captured_regions/hs37d5.bed'
        self.metadata['prepid'] = get_prepid(curs, self.metadata, pseudo_prepid)
        self.metadata['fastq_loc'] = get_fastq_loc(curs, self.metadata)
        self.metadata['output_dir'] = get_output_dir(self.metadata)
        self.metadata['script_dir'] = self.metadata['output_dir']+'/scripts'
        self.metadata['conf_file'] = "{0}/{1}.dragen.conf".format(self.metadata['script_dir'],self.metadata['sample_name'])
        self.metadata['conf_file'] = "{script_dir}/{sample_name}.dragen.conf".format(**self.metadata)
        self.metadata['fastq_dir'] = self.metadata['output_dir']+'/fastq'
        self.metadata['log_dir'] = self.metadata['output_dir']+'/logs'
        self.metadata['dragen_stdout'] = "{log_dir}/{sample_name}.out".format(**self.metadata)
        self.metadata['dragen_stderr'] = "{log_dir}/{sample_name}.err".format(**self.metadata)

        self.metadata['lane'] = get_lanes(curs,self.metadata)


    def get_attribute(self, attribute):
        """return the value requested if present, otherwise raise a TypeError
        """
        if attribute in self.metadata:
            return self.metadata[attribute]
        else:
            raise TypeError("{attribute} is not defined for {CHGVID}".format(
                attribute=attribute, CHGVID=self.CHGVID))

    def get_dict(self):
        """return a dict of the values of this sample
        """
        return self.metadata

    def set(self, attribute, value):
        """set the specified attribute to the given value
        """
        self.metadata[attribute] = value
