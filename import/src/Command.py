import subprocess
import threading
from shlex import split as sxsplit

class TimeoutException(Exception):
    """This exception is to be raised when the timeout is reached for a command
    """
    pass

class Command():
    """Define a class to run a subprocess command with arbitrary arguments to
    Popen that runs for the maximum specified time
    http://stackoverflow.com/a/4825933/5441617
    """
    def __init__(self, cmd, **kwargs):
        self.cmd = cmd
        self.process = None
        self.kwargs = kwargs

    def start(self):
        """run the previously defined command for at most timeout seconds
        """
        def target():
            if "shell" in self.kwargs and self.kwargs["shell"]:
                cmd = self.cmd
            else:
                cmd = sxsplit(self.cmd)
            self.process = subprocess.Popen(cmd, **self.kwargs)
            self.process.communicate()

        self.thread = threading.Thread(target=target)
        self.thread.start()

    def join(self, timeout):
        self.thread.join(timeout)

        if self.thread.is_alive():
            # timeout reached
            self.process.terminate()
            self.thread.join()
            raise TimeoutException
        return self.process.returncode
