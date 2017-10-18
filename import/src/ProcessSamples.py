"""
Definition for an abstract class running a given pipeline step for a maximum
amount of time given
"""
import sys
from Queue import Queue, PriorityQueue
from time import sleep, time
from threading import Thread
from BreakHandler import BreakHandler
import signal
import Command
import re
import os
import subprocess
import logging
import signal

logger = logging.getLogger(__name__)

class ProcessSamples(object):
    def __init__(self, max_samples_concurrently, ignore_priority=False,
                 force_failed_samples=False, run_locally=False, qdel_jobs=True,
                 local_scheduler=False, stdout=sys.stdout, stdout_mode="w",
                 stderr=sys.stderr, stderr_mode="w", **kwargs):
        self.max_samples_concurrently = max_samples_concurrently
        self.force_failed_samples = force_failed_samples
        self.run_locally = run_locally
        self.qdel_jobs = qdel_jobs
        self.local_scheduler = local_scheduler
        self.stdout = stdout
        self.stdout_mode = stdout_mode
        self.stderr = stderr
        self.stderr_mode = stderr_mode
        self.kwargs = kwargs
        self.__dict__.update(kwargs)
        # keep a set of samples processed so we don't try to run it again and
        # again
        self.samples_already_queued = set()
        self.ignore_priority = ignore_priority
        self.samples_queue = Queue() if ignore_priority else PriorityQueue()
        self.all_threads = []
        # if we reach the interrupt threshold and we're deleting 
        self.commands = set()
        self._original_handler = signal.getsignal(signal.SIGINT)
        signal.signal(signal.SIGINT, self._cleanup_processes)
        self.bh = BreakHandler()
        self.bh.enable()

    def _cleanup_processes(self, signum, frame):
        signal.signal(signal.SIGINT, self._original_handler)
        self.bh.enable()
        commands = set(self.commands)
        ps_check = subprocess.check_output(
            "ps ux | tail -n +2 | awk '{printf $2\"\\t\"; printf $11; "
            "for(i=12; i<=NF; i++) {printf \" \"$i;} printf \"\\n\"}'", shell=True)
        for p, COMMAND in [line.split("\t") for line in ps_check.splitlines()]:
            for command in commands:
                if command in COMMAND:
                    subprocess.call(["kill", p])
        sys.exit(1)

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
        # (pseudo_prepid, sample_name, priority, (...))
        # an example would be (sample_type, priority, (capture_kit, prep_id (pseudo_prep_id)))
        pass

    def _get_command(self, pseudo_prepid, sample_name, *args):
        # overwrite to format the actual command to run, and specify the timeout
        pass

    def get_sge_job_ids(self, pseudo_prepid, sample_name, *args):
        """Pass in a pattern= parameter upon initialization if this method is to
        be used"""
        p = re.compile(r"^(\d+) ")
        job_name = re.compile(self.pattern.format(
            pseudo_prepid=pseudo_prepid, sample_name=sample_name,
            *args, **self.kwargs))
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

    def run_command(self, pseudo_prepid, sample_name, *args):
        # overwrite to run one or more commands
        cmd, timeout = self._get_command(pseudo_prepid, sample_name, *args)
        logger.info(cmd)
        start_time = time()
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
            command.start()
            self.commands.add(cmd)
            exit_code = command.join(timeout=timeout)
            self.commands.remove(cmd)
        except Command.TimeoutException:
            if self.qdel_jobs:
                for sge_job in self.get_sge_job_ids(
                    pseudo_prepid, sample_name, *args):
                    logger.debug("Killing job ID: {}".format(sge_job))
                    p = subprocess.Popen(
                        ["qdel", sge_job], preexec_fn=os.setpgrp)
                    p.communicate()
            logger.warning("timed out processing {sample_name}\n".format(
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
        elapsed_seconds = int(round(time() - start_time))
        minutes, seconds = divmod(elapsed_seconds, 60)
        hours, minutes = divmod(minutes, 60)
        if exit_code:
            if self.qdel_jobs:
                for sge_job in self.get_sge_job_ids(
                    pseudo_prepid, sample_name, *args):
                    p = subprocess.Popen(
                        ["qdel", sge_job], preexec_fn=os.setpgrp)
                    p.communicate()
            logger.error("failed processing "
                         "{sample_name}.{pseudo_prepid} after "
                         "{h}:{m:0>2}:{s:0>2}".format(
                             sample_name=sample_name, h=hours,
                             m=minutes, s=seconds, pseudo_prepid=pseudo_prepid))
        else:
            logger.info("succeeded processing {sample_name}.{pseudo_prepid} after "
                        "{h}:{m:0>2}:{s:0>2}".format(
                            sample_name=sample_name, h=hours, m=minutes,
                            s=seconds, pseudo_prepid=pseudo_prepid))

    def process_samples(self):
        done = False
        while True:
            while True:
                self._delete_finished_threads()
                if self.bh.trapped:
                    logger.info("CTRL+C pressed: will exit when all "
                                "currently running samples are done")
                    # need to confirm all threads are done prior to exiting
                    while self.all_threads:
                        logger.info(
                            "Waiting 10 seconds for current threads to close...")
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
            for pseudo_prepid, sample_name, priority, sample_metadata in self._get_samples():
                if ((pseudo_prepid, sample_name, sample_metadata)
                    not in self.samples_already_queued):
                    if self.ignore_priority:
                        record = (pseudo_prepid, sample_name, sample_metadata)
                    else:
                        record = (priority,
                                  (pseudo_prepid, sample_name, sample_metadata))
                    self.samples_queue.put(record)
                    self.samples_already_queued.add(
                        (pseudo_prepid, sample_name, sample_metadata))
            if self.samples_queue.empty():
                # sleep for a while since we don't have any samples to run
                logger.info("waiting for 300 seconds for new samples...")
                sleep(300)
            else:
                if self.ignore_priority:
                    sample_metadata = self.samples_queue.get()
                else:
                    _, sample_metadata = self.samples_queue.get()
                args = list(sample_metadata[:2])
                for arg in sample_metadata[2]:
                    args.append(arg)
                thread = Thread(target=self.run_command, args=args)
                self.all_threads.append(thread)
                thread.start()
        self.bh.disable()
