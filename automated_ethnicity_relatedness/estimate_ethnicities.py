#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Calculate a given sample's ethnicity probabilities
"""

import luigi
from luigi.contrib.sge import SGEJobTask
import os
import sys
sys.path.append(os.path.join(os.path.dirname(
    os.path.dirname(os.path.realpath(__file__))), "import", "src"))
from waldb_globals import *
from ethnicity_db_statements import *

class EthnicityStatus(luigi.Target):
    def __init__(self, prep_id, sample_name, sample_type, field, fn=None):
        self.prep_id = prep_id
        self.sample_name = sample_name
        self.sample_type = sample_type
        self.field = field
        self.fn = fn

    def exists(self):
        if self.fn:
            if not os.path.isfile(self.fn):
                return False
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            cur.execute(ETHNICITY_STATUS.format(
                self.field, self.prep_id, self.sample_name, self.sample_type))
            row = cur.fetchone()
            if row:
                return row[0] == 1
            else:
                cur.execute(INSERT_STATUS.format(
                    self.prep_id, self.sample_name, self.sample_type))
                db.commit()
                return False
        finally:
            if db.open:
                db.close()

def get_family_id(sample_name):
    db = get_connection("seqdb")
    try:
        cur = db.cursor()
        cur.execute(GET_FAMILY_ID.format(sample_name))
        row = cur.fetchone()
        return sample_name if row[0] == "N/A" else row[0]
    finally:
        if db.open:
            db.close()

class CreatePed(SGEJobTask):
    sample_name = luigi.Parameter(description="the sample's name identifier")
    prep_id = luigi.IntParameter(description="the sample's prep ID")
    sample_type = luigi.Parameter(description="the sample's sequencing type")
    vcf = luigi.InputFileParameter(description="the sample's VCF")
    bam = luigi.InputFileParameter(description="the sample's BAM")
    output_directory = luigi.Parameter(description="the output directory")
    markers = luigi.InputFileParameter(
        default="/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/data/"
        "filtered.37MB.master.training.map",
        description="the list of SNPs tested as part of the model")
    ped_script = luigi.InputFileParameter(
        default="/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/"
        "create_ped_map.py", description="the script used to create the PED file")
    task_name_format = luigi.Parameter(
        default="{task_family}.{sample_name}.{prep_id}.{sample_type}",
        significant=False)

    def __init__(self, *args, **kwargs):
        super(CreatePed, self).__init__(*args, **kwargs)
        kwargs["output_directory"] = kwargs["output_directory"].format(**kwargs)
        super(CreatePed, self).__init__(*args, **kwargs)
        if not os.path.isdir(self.output_directory):
            os.makedirs(self.output_directory)
        self.log_file = os.path.join(
            self.output_directory, "{}.{}.createped.log".format(
                self.sample_name, self.prep_id))
        self.family_id = get_family_id(self.sample_name)
        self.ped = ("{output_directory}/{sample_name}.{prep_id}.ped".format(
            output_directory=self.output_directory,
            sample_name=self.sample_name, prep_id=self.prep_id))

    def work(self):
        run_command("/nfs/goldstein/software/python2.7.7/bin/python {} --vcf {} "
                    "--sample_id {} --markers {} --seqtype {} --pseudo_prepid {} "
                    "--bam {} --stem {}.{} -out {} --logfile {} --family_id {}".
                    format(
                        self.ped_script, self.vcf, self.sample_name,
                        self.markers, self.sample_type, self.prep_id, self.bam,
                        self.sample_name, self.prep_id, self.output_directory,
                        self.log_file, self.family_id), self.output_directory,
                    "CreatePed")
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            cur.execute(UPDATE_PED_STATUS.format(
                self.prep_id, self.sample_name, self.sample_type))
            db.commit()
        finally:
            if db.open:
                db.close()

    def output(self):
        return EthnicityStatus(
            self.prep_id, self.sample_name, self.sample_type, "create_ped",
            self.ped)

class CreatePedTest(CreatePed):
    pass

@CreatePed.event_handler(luigi.Event.FAILURE)
def do_some_stuff(task, exception):
    with open("/home/bc2675/blah.txt", "w") as o:
        o.write("abc\n123\nxyz\nijk\n")

class CalculateProbabilities(SGEJobTask):
    """Run the classifier on the given sample
    """
    sample_name = luigi.Parameter(description="the sample's name identifier")
    prep_id = luigi.IntParameter(description="the sample's prep ID")
    sample_type = luigi.Parameter(description="the sample's sequencing type")
    vcf = luigi.InputFileParameter(description="the sample's VCF")
    bam = luigi.InputFileParameter(description="the sample's BAM")
    output_directory = luigi.Parameter(description="the output directory")
    training_model = luigi.InputFileParameter(
        default="/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/data/"
        "37MB_markers_model.obj", description="the already trained model")
    markers = luigi.InputFileParameter(
        default="/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/data/"
        "filtered.37MB.master.training.map",
        description="the list of SNPs tested as part of the model")
    ethnicity_script = luigi.InputFileParameter(
        default="/nfs/goldstein/software/dragen_pipe/master/automated_ethnicity_relatedness/"
        "model_ancestry.py",
        description="the script used to calculate the probabilities")
    task_name_format = luigi.Parameter(
        default="{task_family}.{sample_name}.{prep_id}.{sample_type}",
        significant=False)

    def __init__(self, *args, **kwargs):
        super(CalculateProbabilities, self).__init__(*args, **kwargs)
        kwargs["output_directory"] = kwargs["output_directory"].format(**kwargs)
        super(CalculateProbabilities, self).__init__(*args, **kwargs)
        if not os.path.isdir(self.output_directory):
            os.makedirs(self.output_directory)
        self.ped = os.path.join(
            self.output_directory, "{sample_name}.{prep_id}.ped".format(
                sample_name=self.sample_name, prep_id=self.prep_id))
        self.output_probs = os.path.join(
            self.output_directory, "ethnicity.probs.txt")

    def work(self):
        run_command("/nfs/goldstein/software/python2.7.7/bin/python "
                    "{ethnicity_script} --output-prob-file {output_probs} "
                    "--testing-ped {ped} --mapfile {map} --input-model {model}".
                    format(
                        ethnicity_script=self.ethnicity_script,
                        output_probs=self.output_probs, ped=self.ped,
                        map=self.markers, model=self.training_model),
                    self.output_directory, "CalculateProbabilities")
        db = get_connection("seqdb")
        try:
            cur = db.cursor()
            with open(self.output_probs) as fh:
                fields = fh.next().split("\t")
            sample_metadata = [self.prep_id, self.sample_name, self.sample_type]
            fields.extend(sample_metadata)
            cur.execute(UPDATE_PROBS.format(*fields))
            cur.execute(UPDATE_PREDICT_STATUS.format(*sample_metadata))
            db.commit()
        finally:
            if db.open:
                db.close()

    def requires(self):
        return self.clone(CreatePedTest)

    def output(self):
        return EthnicityStatus(
            self.prep_id, self.sample_name, self.sample_type, "predict")

if __name__ == "__main__":
    sys.exit(luigi.run())
