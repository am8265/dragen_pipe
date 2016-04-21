#Yield_Sum.sh
#Takes list of sample IDs and seqtype, and outputs run summary information
IFS=$'\n'
INPUT_LIST=$1
SEQTYPE=$2
EXOMEKIT=$3
SAMPLES=`less $INPUT_LIST`

sdb() {
        mysql --defaults-group-suffix=sequencedb "$@"
}


#Deletes previous run summary information

if [ "$EXOMEKIT" = 'roche' ]; then
    EXOMEKIT="Roche SeqCap EZ V3"
fi

#Ensures SEQTYPE is viable input
if [ "$SEQTYPE" = 'Custom' ]; then
    SEQTYPE='Custom Capture'
    abbrevSeqType='cc'
    fastqCheckType='genome'
fi
if [ "$SEQTYPE" = 'custom' ]; then
    SEQTYPE='Custom Capture'
    abbrevSeqType='cc'
    fastqCheckType='genome'

fi
if [ "$SEQTYPE" = 'custom_capture' ]; then
    SEQTYPE='Custom Capture'
    abbrevSeqType='cc'
    fastqCheckType='genome'

fi

if [ "$SEQTYPE" = 'genome' ]; then
    SEQTYPE='Genome'
    abbrevSeqType='wg'
    fastqCheckType='genome'

fi
if [ "$SEQTYPE" = 'exome' ]; then
    SEQTYPE='Exome'
    abbrevSeqType='ex'
    fastqCheckType='exome'

fi
if [ "$SEQTYPE" = 'rnaseq' ]; then
    SEQTYPE='RNASeq'
    abbrevSeqType='rs'
    fastqCheckType='genome'

fi


if [ -z "$SEQTYPE" ]; then
    echo Please enter seqtype
    exit 1
fi

#Loops over all samples in Input list and finds most recent prepID for specified SeqType. Then uses that prepID and gathers run summary information for each sample.

for s in $SAMPLES ; do

prepIDs=`sdb -e "select s.prepID from SeqType s join prepT p on p.prepID=s.prepID where s.DBID=(select DBID from SampleT where CHGVID='$s') and s.SeqType='$SEQTYPE' and p.failedPrep='0' and exomekit='$EXOMEKIT' ORDER BY p.prepDate DESC" -NB`

	#echo  "select s.prepID from SeqType s join prepT p on p.prepID=s.prepID where s.DBID=(select DBID from SampleT where CHGVID='$s') and s.SeqType='$SEQTYPE' and p.failedPrep='0' and exomekit='$EXOMEKIT' ORDER BY p.prepDate DESC"
        #echo $s $prepIDs
        pseudo_prepid=$(sdb -e "INSERT INTO pseudo_prepid (pseudo_prepid,prepid) VALUES (NULL,'$prepIDs') ; SELECT LAST_INSERT_ID();" -NB)

        for prepID in $prepIDs ; do
                sdb -e "UPDATE seqdbClone set pseudo_prepid='$pseudo_prepid' WHERE prepid=$prepID"
                sdb -e "INSERT INTO dragen_queue (pseudo_prepid,sample_name,sample_type,capture_kit,priority) VALUES ($pseudo_prepid,'$s','$(echo $SEQTYPE | sed 's/ /_/g')','$EXOMEKIT',4)"
	done
done
