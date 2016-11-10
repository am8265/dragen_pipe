"""
Definition for an abstract class running a given pipeline step for a maximum
amount of time given
"""
import sys
from Queue import PriorityQueue
from time import sleep
from threading import Thread
from BreakHandler import BreakHandler
import signal
import Command
import re

def preexec_function():
    # ignore SIGINT by setting to SIG_IGN
    signal.signal(signal.SIGINT, signal.SIG_IGN)

class ProcessSamples():
    def __init__(self, max_samples_concurrently, qdel_jobs=True, **kwargs):
        self.max_samples_concurrently = max_samples_concurrently
        self.qdel_jobs = qdel_jobs
        self.kwargs = kwargs
        # keep a set of samples processed so we don't try to run it again and
        # again
        self.samples_already_queued = set()
        self.samples_queue = PriorityQueue()
        self.all_threads = []
        self.bh = BreakHandler()
        self.bh.enable()

    def _delete_finished_threads(self):
        to_delete = []
        for x, thread in enumerate(self.all_threads):
            if not thread.is_alive():
                to_delete.append(x)
        for x in to_delete[::-1]:
            del self.all_threads[x]

    def _get_samples(self):
        # overwrite to get samples from DB, check files, etc.
        # should return a list of tuples of:
        # (sample_name, sample_type, capture_kit, prep_id (pseudo_prep_id), priority)
        pass

    def _get_command(self, sample_name, sample_type, capture_kit, prep_id):
        # overwrite to format the actual command to run, and specify the timeout
        pass

    def get_sge_job_ids(self, sample_name, sample_type, capture_kit, prep_id):
        """Pass in a pattern= parameter upon initialization if this method is to
        be used"""
        p = re.compile(r"^(\d+) ")
        job_name = self.pattern.format(
            sample_name=sample_name, sample_type=sample_type,
            capture_kit=capture_kit, prep_id=prep_id, **self.kwargs)
        proc = subprocess.Popen(["qstat", "-r", "-ne"], stdout=subprocess.PIPE,
                                preexec_fn=preexec_function)
        qstat_output = proc.communicate()[0]
        job_ids = []
        for x, line in enumerate(qstat_output):
            m = p.match(line)
            if m:
                job_id = m.group(1)
                if job_name == qstat_output[x + 1].split()[-1]:
                    job_ids.append(job_id)
        return job_ids

    def run_command(self, sample_name, sample_type, capture_kit, prep_id):
        # overwrite to run one or more commands
        cmd, timeout = self._get_command(
            sample_name=sample_name, sample_type=sample_type,
            capture_kit=capture_kit, prep_id=prep_id)
        command = Command.Command(cmd, preexec_fn=preexec_function)
        try:
            exit_code = command.run(timeout=timeout)
        except Command.TimeoutException:
            if self.qdel_jobs:
                for sge_job in self.get_sge_job_ids(
                    sample_name, sample_type, capture_kit, prep_id):
                    p = subprocess.Popen(
                        ["qdel", sge_job], preexec_fn=preexec_function)
                    p.communicate()
            sys.stderr("timed out importing {sample_name}\n".format(
                sample_name=sample_name))
        if exit_code:
            sys.stderr("failed importing {sample_name}\n".format(
                sample_name=sample_name))

    def process_samples(self):
        done = False
        while True:
            while True:
                self._delete_finished_threads()
                if self.bh.trapped:
                    sys.stderr.write("CTRL+C pressed: will exit when all "
                                     "currently running samples are done\n")
                    # need to confirm all threads are done prior to exiting
                    all_threads_done = False
                    while self.all_threads:
                        sys.stderr.write(
                            "Waiting 10 seconds for current threads to close...\n")
                        sleep(10)
                        self._delete_finished_threads()
                    done = True
                    break
                if len(self.all_threads) == self.max_samples_concurrently:
                    # wait until a slot is freed up to submit next sample
                    sleep(1)
                else:
                    break
            if done:
                break
            for sample_name, sample_type, capture_kit, prep_id, priority in (
                self._get_samples()):
                sample_metadata = (
                    sample_name, sample_type, capture_kit, prep_id)
                if sample_metadata not in self.samples_already_queued:
                    self.samples_queue.put((priority, sample_metadata))
                    self.samples_already_queued.add(sample_metadata)
            if self.samples_queue.empty():
                # sleep for a while since we don't have any samples to run
                sleep(300)
            else:
                priority, sample_metadata = self.samples_queue.get()
                thread = Thread(
                    target=self.command, args=sample_metadata)
                self.all_threads.append(thread)
                thread.start()
