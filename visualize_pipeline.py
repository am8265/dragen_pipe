#!/nfs/goldstein/software/python2.7.7/bin/python
"""
Take a python luigi script and a command for how it would be run, then generate
a DAG visualizing the dependency tree
"""
import argparse
import os
import sys
import luigi
from graphviz import Digraph
from collections import defaultdict, Counter

def build_graph(parent_class, task, graph, tasks, class_counter):
    task_exists = False
    for task_id, node_name in tasks:
        if task.task_id == task_id:
            task_exists = True
            # N.B. if the task_id matches, we should use this node_name to add
            # the edge
            break
    if not task_exists:
        class_name = task.__class__.__name__
        node_name = class_name
        if class_name in class_counter:
            node_name += str(class_counter[class_name])
        class_counter[class_name] += 1
        tasks.add((task.task_id, node_name))
        graph.node(node_name)
    # add the edge if the node is encountered again from a different source
    # class, but don't process its dependencies again
    graph.edge(node_name, parent_class)
    if not task_exists:
        for dep in task.deps():
            build_graph(node_name, dep, graph, tasks, class_counter)

def visualize_pipeline(luigi_args, output_fn, file_format):
    with luigi.cmdline_parser.CmdlineParser.global_instance(luigi_args) as ctxt:
        core_module = ctxt.known_args.core_module
        base_task = ctxt.get_task_obj()
    if not output_fn:
        output_fn = core_module
    tasks = set()
    class_counter = Counter()
    class_name = base_task.__class__.__name__
    tasks.add((base_task.task_id, class_name))
    class_counter[class_name] += 1

    graph = Digraph(format=file_format)
    graph.node(class_name)
    for dep in base_task.deps():
        build_graph(class_name, dep, graph, tasks, class_counter)
    graph.render(filename=output_fn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument("-o", "--output", help="specify an output file name "
                        "(otherwise will append .{format} to the script name)")
    parser.add_argument("-f", "--format", default="svg",
                        help="the output file type")
    args, luigi_args = parser.parse_known_args()
    visualize_pipeline(luigi_args, args.output, args.format)
