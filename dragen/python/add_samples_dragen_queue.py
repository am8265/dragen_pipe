#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Takes a sample and imports them into dragen_queue
"""
import argparse
import MySQLdb
import traceback
from ConfigParser import SafeConfigParser
from dragen_sample import dragen_sample
from glob import glob

def main(samples, sample_type, capture_kit, debug, priority, curs):
    for sample_name in samples:
        sample = dragen_sample(sample_name,sample_type,0,capture_kit, curs)
        pseudo_prepid = insert_pseudo_prepid(sample,debug)
        #insert_dragen_queue(sample,debug,priority,pseudo_prepid)
        insert_seqdbClone(sample,debug,priority,pseudo_prepid)

def insert_pseudo_prepid(sample,debug):
    #INSERT first prep and get the new pseudo_prepid
    prepID = getPrepID(sample,debug)
    query = ("INSERT INTO pseudo_prepid "
                "(pseudo_prepid,prepid) "
                "VALUES (NULL,{}) "
            ).format(prepID)
    if debug:
        print query
        print sample.metadata['prepid']
    curs.execute(query)

    pseudo_prepid_query = "SELECT LAST_INSERT_ID()"
    curs.execute(pseudo_prepid_query)
    pseudo_prepid = curs.fetchone()

    if debug:
       print pseudo_prepid[0]

    #Iterate over remaining prepids, if any inserting with the new pseudo_prepid
    if len(sample.metadata['prepid']) > 1:
        for prepid in sample.metadata['prepid'][1:]:
            multiple_query = ("INSERT INTO pseudo_prepid "
                    "(pseudo_prepid,prepid) "
                    "VALUES ({pseudo_prepid},{prepid}) ; "
            ).format(prepid=prepid,pseudo_prepid=pseudo_prepid[0])
            if debug:
               print multiple_query
            curs.execute(multiple_query)

    return pseudo_prepid[0]

def getPrepID(sample,debug):
    sample_name=sample.metadata['sample_name']
    sample_type=sample.metadata['sample_type']
    capture_kit=sample.metadata['capture_kit']

    query = ("SELECT p.prepID FROM prepT p "
            "JOIN SeqType st ON p.prepid=st.prepid "
            "WHERE CHGVID='{}' "
            "AND seqtype = '{}' "
            "AND exomekit = '{}' "
            "AND failedprep = 0"
            ).format(sample_name,sample_type,capture_kit)

    if debug:
        print query
    curs.execute(query)
    prepid = curs.fetchone()[0]
    #print prepid
    return prepid

def insert_dragen_queue(sample,debug,priority,pseudo_prepid):
    dragen_queue_query = ("INSERT INTO dragen_queue "
                        "(pseudo_prepid,sample_name,sample_type,capture_kit,priority) "
                        "VALUES ( "
                            "{pseudo_prepid},"
                            "'{sample_name}',"
                            "'{sample_type}',"
                            "'{capture_kit}',"
                            "{priority})"
                            ).format(pseudo_prepid=pseudo_prepid,
                                sample_name=sample.metadata['sample_name'],
                                sample_type=sample.metadata['sample_type'],
                                capture_kit=sample.metadata['capture_kit'],
                                priority=priority)

    if debug:
        print dragen_queue_query
    curs.execute(dragen_queue_query)

def insert_seqdbClone(sample,debug,priority,pseudo_prepid):
    for prepID in sample.metadata['prepid']:

        seqdbClone_query = ("UPDATE seqdbClone "
                            "SET pseudo_prepid={pseudo_prepid} WHERE prepID={prepID}"
                           ).format(pseudo_prepid=pseudo_prepid,prepID=prepID)
        if debug:
            print seqdbClone_query
        curs.execute(seqdbClone_query)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", nargs="+",
                        help="Specify sample")
    parser.add_argument("-t", "--sample_type", required=True,
                        help="Type of samples. Ex. Genome, Exome...")
    parser.add_argument("-k", "--capture_kit", required=True,
                        help="Type of capture_kit used by samples.")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
    parser.add_argument("-p", "--priority", default=4)
    args=parser.parse_args()

    CNF = "/nfs/goldstein/software/dragen/dragen.cnf" # defaults file for pipeline
    config_parser = SafeConfigParser()
    config_parser.read(CNF)
    parameters = {}
    for parameter in ("DB_CONFIG_FILE",):
        parameters[parameter] = config_parser.get("dragen_pipe", parameter.upper())

    if args.debug:
        print "Parameters used: {0}".format(parameters)

    if args.test:
        args.database="testdb"
    else:
        args.database="sequenceDB"
    db = MySQLdb.connect(db='sequenceDB',
            read_default_group="client{database}".format(database=args.database),
            read_default_file=parameters["DB_CONFIG_FILE"],
            )
    curs = db.cursor()
    curs.execute('BEGIN;')
    if args.sample:
        sample = args.sample
    else:
        sample = args.samples
    if args.capture_kit:
        capture_kit = args.capture_kit
    else:
        #Genome samples do have an empty string as their capture_kit
        capture_kit = 'N/A'

    """For unknown reasons custom capture samples are expressed as 
    "custom capture" in the SeqType table but in seqdbClone and other tables
    its "custom_capture" with the underscore.  Since prepids are determined 
    using a join of the prepT and SeqType table we use seqtype with spaces"""
    if args.sample_type.lower() == "custom_capture":
        sample_type = 'custom capture'
    else:
        sample_type = args.sample_type
    try:
        main(args.sample, sample_type, capture_kit, args.debug, args.priority, curs)
        if args.debug:
            print "MySQL COMMIT"
        curs.execute("COMMIT")
    except Exception, e:
        traceback.print_exc()
        curs.execute('ROLLBACK')
        curs.close()
        print "Import FAILURE!"
