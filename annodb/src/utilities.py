import argparse
import sys
import os
import gzip
from socket import gethostname
from time import ctime

LOG_FORMAT = "Running on {host} @{ctime}: {command}"

def get_script_path():
    """return the path of the script being executed
    """
    return os.path.dirname(os.path.realpath(sys.argv[0]))

def log_output(command):
    """return the string to put in the log file that pertains to the command
    being run
    """
    return LOG_FORMAT.format(
        host=gethostname(), ctime=ctime(), command=command) + "\n"

def is_gzipped(file_name):
    """is the specified file gzipped?
    """
    with open(file_name, "rb") as fh:
        magic_number = fh.read(2)
    return magic_number == "\x1f\x8b"

def fh(file_name):
    """return a file handle to the file whether it's gzipped or not
    """
    if is_gzipped(file_name):
        return gzip.open(file_name)
    else:
        return open(file_name)

class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    """multiple inheritance of two argparse formatters
    """
    pass

def file_exists(arg):
    """check if the given file exists
    """
    if os.path.isfile(arg):
        return os.path.realpath(arg)
    else:
        raise argparse.ArgumentTypeError(arg + " does not exist")

def create_directory_for_file_if_needed(arg):
    """create the directory to contain the specified output file if necessary
    """
    rp = os.path.realpath(arg)
    directory_name = os.path.dirname(rp)
    if not os.path.isdir(directory_name):
        try:
            os.makedirs(directory_name)
        except OSError:
            raise arpgarse.ArgumentTypeError(
                "could not create directory {}".format(directory_name))
    return rp

def make_directory_for_new_file(arg):
    """verify file doesn't exist and make directory for the file if needed
    """
    if os.path.isfile(arg):
        raise argparse.ArgumentTypeError(arg + " already exists")
    elif os.path.isdir(arg):
        raise argparse.ArgumentTypeError(arg + " is a directory")
    else:
        rp = os.path.realpath(arg)
        make_directory_if_not_exists(os.path.dirname(rp))
        return rp

def make_directory_if_not_exists(arg):
    """make any directories in the given path if possible
    """
    rp = os.path.realpath(arg)
    if not os.path.isdir(rp):
        try:
            os.makedirs(rp)
        except Exception:
            raise argparse.ArgumentTypeError("cannot create " + rp)
    return rp

def input_file_list(arg):
    """a list of input arguments, one per line
    comments begin with # and are stripped out
    """
    rp = file_exists(arg)
    arg_list = []
    with open(rp) as fh:
        for line in fh:
            line = line.strip()
            entry = line.split("#")[0]
            if entry:
                arg_list.append(entry)
    return arg_list

def simplify_REF_ALT_alleles(REF, ALT):
    """take a potentially complex representation of a pair of alleles introduced
    into multiallelic sites and simplify to the canonical form
    http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
    """
    # strip shared suffix
    strip = 0
    for x in xrange(1, min(len(REF), len(ALT))):
        if REF[-x] == ALT[-x]:
            strip -= 1
        else:
            break
    if strip:
        REF = REF[:strip]
        ALT = ALT[:strip]
    # strip shared prefix
    strip = 0
    for x in xrange(0, min(len(REF), len(ALT)) - 1):
        if REF[x] == ALT[x]:
            strip += 1
        else:
            break
    # return simplified REF, ALT, and position offset
    return REF[strip:], ALT[strip:], strip

def create_INFO_dict(INFO):
    """take the INFO field from a VCF, parse and return a dict of values
    """
    return dict(entry.split("=", 1) for entry in INFO.split(";") if "=" in entry)

def merge_dicts(*dict_list):
    """merge an arbitrary number of dictionaries into a single dictionary and
    return it
    """
    ndicts = len(dict_list)
    if ndicts == 0:
        raise ValueError("no dicts passed in!")
    else:
        new_dict = dict_list[0].copy()
        for d in dict_list[1:]:
            new_dict.update(d)
        return new_dict

def binary_search_intervals(intervals, value, extension=0, debug=False):
    """search a sorted list of intervals [x, y] and return the index and
    interval with the maximum x such that x <= value
    extension allows for extending arbitrary number of positions past the
    intervals
    """
    min_index = 0
    nvalues = len(intervals)
    if not nvalues:
        raise TypeError("The intervals list is empty.")
    max_index = nvalues
    while min_index < max_index:
        mid_index = (min_index + max_index) / 2
        if debug:
            print("min_index: {}, mid_index: {}, max_index: {}, interval: {}".
                  format(min_index, mid_index, max_index, intervals[mid_index]))
        compare_value = intervals[mid_index][0] - extension
        if compare_value <= value:
            min_index = mid_index + 1
        else:
            max_index = mid_index
    if min_index < nvalues and value == (intervals[min_index][0] - extension):
        return min_index, intervals[min_index]
    elif min_index == 0:
        raise ValueError(
            str(value) + " is less than any value in the intervals.")
    else:
        return min_index - 1, intervals[min_index - 1]

def strip_suffix(string, suffix):
    """strip the suffix if present
    """
    return string[:-len(suffix)] if string.endswith(suffix) else string

def natural_number(arg):
    """return the integer value of arg if it is a natural number, otherwise
    raise an exception
    """
    try:
        value = int(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "{arg} is not a valid integer".format(arg=arg))
    if arg < 0:
        raise argparse.ArgumentTypeError(
            "{arg} is negative".format(arg=arg))
    return value
