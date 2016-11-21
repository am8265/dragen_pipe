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
import os
import subprocess
import logging

#def preexec_function():
#    # ignore SIGINT by setting to SIG_IGN
#    signal.signal(signal.SIGINT, signal.SIG_IGN)
logger = logging.getLogger(__name__)

class ProcessSamples(object):
    def __init__(self, max_samples_concurrently, qdel_jobs=True,
                 stdout=sys.stdout, stdout_mode="w", stderr=sys.stderr,
                 stderr_mode="w", **kwargs):
        self.max_samples_concurrently = max_samples_concurrently
        self.qdel_jobs = qdel_jobs
        self.stdout = stdout
        self.stdout_mode = stdout_mode
        self.stderr = stderr
        self.stderr_mode = stderr_mode
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
        # (sample_name, priority, (...))
        # an example would be (sample_type, priority, (capture_kit, prep_id (pseudo_prep_id)))
        pass

    def _get_command(self, sample_name, *args):
        # overwrite to format the actual command to run, and specify the timeout
        pass

    def get_sge_job_ids(self, sample_name, *args):
        """Pass in a pattern= parameter upon initialization if this method is to
        be used"""
        p = re.compile(r"^(\d+) ")
        job_name = re.compile(self.pattern.format(
            sample_name=sample_name, *args, **self.kwargs))
        proc = subprocess.Popen(["qstat", "-r", "-ne"], stdout=subprocess.PIPE,
                                preexec_fn=os.setpgrp)
        qstat_output = proc.communicate()[0].splitlines()
        job_ids = []
        for x, line in enumerate(qstat_output):
            m = p.match(line)
            if m:
                job_id = m.group(1)
                if job_name.match(qstat_output[x + 1].split()[-1]):
                    job_ids.append(job_id)
        return job_ids

    def run_command(self, sample_name, *args):
        # overwrite to run one or more commands
        cmd, timeout = self._get_command(sample_name, *args)
        print(cmd)
        if type(self.stdout) is file:
            close_stdout = False
            stdout = self.stdout
        else:
            close_stdout = True
            stdout = open(self.stdout, self.stdout_mode)
        if type(self.stderr) is file:
            close_stderr = False
            stderr = self.stderr
        else:
            close_stderr = True
            stderr = open(self.stderr, self.stderr_mode)
        command = Command.Command(
            cmd, preexec_fn=os.setpgrp, stdout=stdout, stderr=stderr)
        try:
            exit_code = command.run(timeout=timeout)
        except Command.TimeoutException:
            if self.qdel_jobs:
                for sge_job in self.get_sge_job_ids(sample_name, *args):
                    logger.debug("Killing job ID: {}".format(sge_job))
                    p = subprocess.Popen(
                        ["qdel", sge_job], preexec_fn=os.setpgrp)
                    p.communicate()
            sys.stderr.write("timed out processing {sample_name}\n".format(
                sample_name=sample_name))
            if close_stdout:
                stdout.close()
            if close_stderr:
                stderr.close()
            return
        if close_stdout:
            stdout.close()
        if close_stderr:
            stderr.close()
        if exit_code:
            if self.qdel_jobs:
                for sge_job in self.get_sge_job_ids(sample_name, *args):
                    p = subprocess.Popen(
                        ["qdel", sge_job], preexec_fn=os.setpgrp)
                    p.communicate()
            sys.stderr.write("failed processing {sample_name}\n".format(
                sample_name=sample_name))
        else:
            sys.stderr.write("succeeded processing {sample_name}\n".format(
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
            for sample_name, priority, sample_metadata in self._get_samples():
                if ((sample_name, sample_metadata)
                    not in self.samples_already_queued):
                    self.samples_queue.put(
                        (priority, (sample_name, sample_metadata)))
                    self.samples_already_queued.add(
                        (sample_name, sample_metadata))
            if self.samples_queue.empty():
                # sleep for a while since we don't have any samples to run
                print("waiting for 300 seconds for new samples...")
                sleep(300)
            else:
                priority, sample_metadata = self.samples_queue.get()
                sample_metadata = tuple(
                    [sample_metadata[0]] + list(sample_metadata[1]))
                thread = Thread(
                    target=self.run_command, args=sample_metadata)
                self.all_threads.append(thread)
                thread.start()
        self.bh.disable()
