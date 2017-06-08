#!/nfs/goldstein/software/python2.7/bin/python2.7.7
"""
Takes a list of samples and imports them into dragen_queue
"""
import argparse
import MySQLdb
import traceback
from ConfigParser import SafeConfigParser
from glob import glob

def main(sample_list, debug, curs):
    with open(sample_list) as infile:
        for sample in infile:
            sample_name,sample_type,capture_kit = sample.strip().split('\t')
            if capture_kit == 'N/A':
                #Genome samples do have an empty string as their capture_kit
                capture_kit = ''

            """For unknown reasons custom capture samples are expressed as
            "custom capture" in the SeqType table but in seqdbClone and other tables
            its "custom_capture" with the underscore.  Since prepids are determined
            using a join of the prepT and SeqType table we use seqtype with spaces"""
            if sample_type.lower() == "custom_capture":
                sample_type = 'custom capture'

            prepids = getPrepIDs(sample_name,sample_type,capture_kit,debug,curs)
            pseudo_prepid = insert_pseudo_prepid(prepids,sample_name,sample_type,capture_kit,debug,curs)
            insert_seqdbClone(sample_name,sample_type,capture_kit,debug,pseudo_prepid,prepids,curs)

def getPrepIDs(sample_name,sample_type,capture_kit,debug,curs):
    if sample_type.lower() == 'genome':
        query = ("SELECT p.prepID from prepT p "
                "JOIN SeqType st ON p.prepid=st.prepid "
                "WHERE CHGVID='{}' AND seqtype='{}' and failedPrep=0"
                ).format(sample_name,sample_type,capture_kit)
    else:
        query = ("SELECT p.prepID FROM prepT p "
                "JOIN captureKit c ON c.prepT_name=p.exomekit "
                "JOIN SeqType st ON p.prepid=st.prepid "
                "WHERE CHGVID='{}' AND seqtype='{}' AND "
                "name='{}' and chr='all' AND failedPrep=0"
                ).format(sample_name,sample_type,capture_kit)

    if debug:
        print query
    curs.execute(query)
    prepids = curs.fetchall()
    if debug:
        print prepids
    return prepids

def insert_pseudo_prepid(prepids,sample_name,sample_type,capture_kit,debug,curs):
    #INSERT first prep and get the new pseudo_prepid
    query = ("INSERT INTO pseudo_prepid "
                "(pseudo_prepid,prepid) "
                "VALUES (NULL,{prepID}) "
            ).format(prepID=prepids[0][0])
    if debug:
        print query
    curs.execute(query)

    pseudo_prepid_query = "SELECT LAST_INSERT_ID()"
    curs.execute(pseudo_prepid_query)
    pseudo_prepid = curs.fetchone()

    if debug:
       print pseudo_prepid[0]

    #Iterate over remaining prepids, if any inserting with the new pseudo_prepid
    if len(prepids) > 1:
        for prepid in prepids[1:]:
            multiple_query = ("INSERT INTO pseudo_prepid "
                    "(pseudo_prepid,prepid) "
                    "VALUES ({pseudo_prepid},{prepid}) ; "
            ).format(prepid=prepid[0],pseudo_prepid=pseudo_prepid[0])
            if debug:
               print multiple_query
            curs.execute(multiple_query)

    return pseudo_prepid[0]

def insert_seqdbClone(sample_name,sample_type,capture_kit,debug,pseudo_prepid,prepids,curs):
    for prepid in prepids:
        query = ("UPDATE seqdbClone "
                "SET pseudo_prepid={pseudo_prepid} WHERE prepID={prepID}"
                ).format(pseudo_prepid=pseudo_prepid,prepID=prepid[0])
        if debug:
            print query
        curs.execute(query)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample-list", required=True,
                        help="Specify sample list")
    parser.add_argument("-d", "--debug", default=False, action="store_true",
                        help="Verbose output for debugging")
    parser.add_argument("--test", default=False, action="store_true",
                        help="Query and updates to the database occur on the "
                        "test server")
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

    print args.database
    db = MySQLdb.connect(db=args.database,
            read_default_group="client{database}".format(database=args.database),
            read_default_file=parameters["DB_CONFIG_FILE"])

    curs = db.cursor()
    curs.execute('BEGIN;')
    try:
        main(args.sample_list, args.debug, curs)
        if args.debug:
            print "MySQL COMMIT"
        curs.execute("COMMIT")
    except Exception, e:
        traceback.print_exc()
        curs.execute('ROLLBACK')
        curs.close()
        print "Import FAILURE!"
