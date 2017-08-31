#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Run Dan's script and parse its output to test for Avere load before proceeding
with I/O intensive tasks.
"""

import subprocess
import os
import logging
import sys
import re

SCRIPT = "/home/dh2880/avere/avere-load"
MAX_LOAD_ANY_NODE = 0.9 # don't proceed if any node has a load >= this value
MAX_LOAD_ALL_NODES = 0.7 # don't proceed if all nodes have a load >= this value
MAX_NODES_DOWN = 1 # don't proceed if > this many nodes are inaccessible
MAX_SECONDS_FOR_UPDATE = 150 # don't proceed if the last update was > this many seconds ago

LOGGING_LEVELS = {
    "CRITICAL":logging.CRITICAL, "ERROR":logging.ERROR,
    "WARNING":logging.WARNING, "INFO":logging.INFO, "DEBUG":logging.DEBUG}


def check_avere(script_fn=SCRIPT, max_load_any_node=MAX_LOAD_ANY_NODE,
                max_load_all_nodes=MAX_LOAD_ALL_NODES,
                max_nodes_down=MAX_NODES_DOWN, logging_level=logging.WARNING):
    """
    Return True if Avere's status as per the parameters specified suggest that
    it's ok to proceed with I/O, False otherwise

    @script_fn - the path to Dan's avere load script
    @max_load_any_node - return False if any node has load >= this value
    @max_load_all_nodes - return False if all nodes have load >= this value
    @max_nodes_down - return False if > this number of nodes are unreachable
    """
    if not os.path.isfile(script_fn):
        raise OSError("{script_fn} does not exist".format(script_fn=script_fn))
    if not (0.0 <= max_load_any_node <= 1.0):
        raise ValueError("max_load_any_node ({}) must be within [0.0, 1.0]".
                         format(max_load_any_node))
    if not (0.0 <= max_load_all_nodes <= 1.0):
        raise ValueError("max_load_all_nodes ({}) must be within [0.0, 1.0]".
                         format(max_load_all_nodes))
    if max_nodes_down not in (0, 1, 2):
        raise ValueError("max_nodes_down ({}) must be in [0, 1, 2]".format(
            max_nodes_down))

    logger = logging.getLogger(__name__)
    logger.setLevel(logging_level)
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging_level)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    logger.debug("Executing {}".format(script_fn))
    p = subprocess.Popen(
        [script_fn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.splitlines()
    err = err.splitlines()
    logger.debug("Done with {}".format(script_fn))

    if p.returncode:
        if out:
            logger.debug("Stdout: {}".format(out))
        if err:
            logger.debug("Stderr: {}".format(err))
        raise subprocess.CalledProcessError(p.returncode, script_fn)
    status_update = re.compile(r"^Avere status logs updated (\d+)s ago$")
    found = False
    for line in out:
        m = status_update.match(line)
        if m:
            found = True
            seconds = int(m.group(1))
            if seconds > MAX_SECONDS_FOR_UPDATE:
                logger.warning("Max update time exceeded: {}".format(seconds))
                return False
    if not found:
        raise ValueError("Didn't find update time")
    logger.debug("Last update was {} seconds ago".format(seconds))
    load_p = re.compile(r"^Load: node1=(\d+)%, node2=(\d+)%, node3=(\d+)%")
    load = []
    found = False
    for line in out:
        m = load_p.match(line)
        if m:
            for node in xrange(3):
                load.append(float(m.group(1 + node)) / 100.0)
            found = True
            break
    if not found:
        raise ValueError("Didn't find load data")
    logger.debug("Load is {}".format(load))
    clients_found = [False, False, False]
    node_client_p = re.compile(r"^node(\d)_clients=")
    nodes_down = set()
    for line in out:
        m = node_client_p.match(line)
        if m:
            client = int(m.group(1)) - 1
            clients_found[client] = True
            if len(line.split("=")[1].split(",")) < 2:
                nodes_down.add(client)
                load[client] = max_load_all_nodes
    if not all(clients_found):
        raise ValueError("Didn't find clients for all nodes")
    if nodes_down:
        for node in nodes_down:
            logger.warning("Node {} is down".format(node))
    if len(nodes_down) > max_nodes_down:
        logger.warning("Max # nodes down exceeded: {}".format(len(nodes_down)))
        return False
    if any([value >= max_load_any_node for value in load]):
        logger.info("Max load exceeded for one node")
        return False
    if all([value >= max_load_all_nodes for value in load]):
        logger.info("Max load exceeded across all nodes")
        return False
    return True

if __name__ == "__main__":
    import argparse
    class DereferenceKeyAction(argparse.Action):
        """Define a class for automatically converting the key specified from
        choices to its corresponding value
        """
        def __init__(self, option_strings, dest, nargs=None, default=None,
                     choices=None, **kwargs):
            if nargs is not None:
                raise ValueError("nargs should not be specified")
            if type(choices) is not dict:
                raise TypeError("choices must be a dict")
            super(DereferenceKeyAction, self).__init__(
                option_strings, dest, choices=choices, **kwargs)
            if default:
                self.default = choices[default]

    def __call__(self, parser, namespace, values, option_string):
        setattr(namespace, self.dest, self.choices[values])
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--script", default=SCRIPT,
                        help="the path to the Avere checking script")
    parser.add_argument("--max_load_any_node", default=MAX_LOAD_ANY_NODE,
                        type=float,
                        help="return False if any node has a load >= this")
    parser.add_argument("--max_load_all_nodes", default=MAX_LOAD_ALL_NODES,
                        type=float,
                        help="return False if all nodes have a load >= this")
    parser.add_argument("--max_nodes_down", default=MAX_NODES_DOWN,
                        type=int, choices=(0, 1, 2),
                        help="return False if > this many nodes are inaccessible")
    parser.add_argument("-l", "--level", default="DEBUG",
                        action=DereferenceKeyAction, choices=LOGGING_LEVELS,
                        help="the logging level to use")
    args = parser.parse_args()
    avere_ok = check_avere(
        args.script, args.max_load_any_node, args.max_load_all_nodes,
        args.max_nodes_down, args.level)
    sys.exit(0 if avere_ok else 1)
