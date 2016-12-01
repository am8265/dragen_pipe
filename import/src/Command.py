import subprocess
import threading

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

    def run(self, timeout):
        """run the previously defined command for at most timeout seconds
        """
        def target():
            self.process = subprocess.Popen(self.cmd, shell=True, **self.kwargs)
            self.process.communicate()

        thread = threading.Thread(target=target)
        thread.start()
        thread.join(timeout)

        if thread.is_alive():
            # timeout reached
            self.process.terminate()
            thread.join()
            raise TimeoutException
        return self.process.returncode
