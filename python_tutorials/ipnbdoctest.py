#!/usr/bin/env python
"""
simple example script for running and testing notebooks.

Usage: `ipnbdoctest.py foo.ipynb [bar.ipynb [...]]`

Each cell is submitted to the kernel, and the outputs are compared with those stored in the notebook.

This is based on a gist under https://gist.github.com/jhprinz/f2a6dfcb4adf903b3669, which is in turn based on
a gist under https://gist.github.com/minrk/2620735

I have just tried to make things python 2/3 compatible, and updated it to reflect changes to the
ipython messaging specifications 5.0

TODO: It would be nice to add some more output about failed diffs. So far I leave it this way
"""

import os,sys
import base64
import re
import argparse
import platform

import difflib

printdiffs = False # default is not to print diff of differences between test and reference cells

try:
    #python 2
    from Queue import Empty
except:
    #python 3
    from queue import Empty

try:
    from IPython.kernel import KernelManager
except ImportError:
    from IPython.zmq.blockingkernelmanager import BlockingKernelManager as KernelManager
 
from IPython.nbformat.current import reads, NotebookNode

class TravisConsole(object):
    """
    A wrapper class to allow easier output to the console especially for travis
    """
    def __init__(self):
        self.stream = sys.stdout
        self.linebreak = '\n'
        self.fold_count = dict()
        self.fold_stack = dict()

    def fold_open(self, name):
        if name not in self.fold_count:
            self.fold_count[name] = 0
            self.fold_stack[name] = []

        self.fold_count[name] += 1
        fold_name = name.lower() + '.' + str(self.fold_count[name])
        self.fold_stack[name].append(fold_name)

        self.writeln("travis_fold:start:" + fold_name)

    def fold_close(self, name):
        fold_name = self.fold_stack[name].pop()
        self.writeln("travis_fold:end:" + fold_name)

    def _indent(self, s, num = 4):
        lines = s.splitlines(True)
        lines = map(lambda s: ' ' * num + s, lines)
        return ''.join(lines)

    def writeln(self, s, indent = 0):
        self.write(s, indent)
        if s[-1] != '\n':
            # make a line break if there is none present
            self.br()

    def br(self):
        """
        Write a linebreak
        """
        self.stream.write(self.linebreak)

    def write(self, s, indent = 0):
        if indent > 0:
            self.stream.write(self._indent(s, indent))
        else:
            self.stream.write(s)


    def compare_png(a64, b64):
        """compare two b64 PNGs (incomplete)"""
        try:
            import Image
        except ImportError:
            pass
        adata = base64.decodestring(a64)
        bdata = base64.decodestring(b64)
        return True

    def red(self, s):
        RED = '\033[31m'
        DEFAULT = '\033[39m'
        return RED + s + DEFAULT

    def green(self, s):
        GREEN = '\033[32m'
        DEFAULT = '\033[39m'
        return GREEN + s + DEFAULT

    def blue(self, s):
        BLUE = '\033[36m'
        DEFAULT = '\033[39m'
        return BLUE + s + DEFAULT

class IPyTestConsole(TravisConsole):
    """
    Add support for different output results
    """
    def __init__(self):
        super(IPyTestConsole, self).__init__()

        self.default_results = {
            'success' : True,       # passed without differences
            'kernel' : False,       # kernel (IPYTHON) error occurred
            'error' : False,        # errors during execution
            'timeout' : True,       # kernel run timed out
            'diff' : True,          # passed, but with differences in the output
            'skip' : True,          # cell has been skipped, not even tried to execute
            'ignore' : True         # cell has been executed, but not compared
        }

        self.pass_count = 0
        self.fail_count = 0

        self.result_count = { key : 0 for key in self.default_results.keys() }

    def write_result(self, result, okay_list = None):

        my_list = self.default_results.copy()
        if okay_list is not None:
            my_list.update(okay_list)

        if my_list[result]:
            tv.write(self.green('ok'))
            self.pass_count += 1
        else:
            tv.write(self.red('fail'))
            self.fail_count += 1

        tv.writeln(' [' + result + ']')
        self.result_count[result] += 1


class IPyKernel(object):
    """
    A simple wrapper class to run cells in an IPython Notebook.

    Notes
    -----
    Use `with` construct to properly instantiate
    """

    def __init__(self, console = None):
        # default timeout time is 60 seconds
        self.default_timeout = 60
        self.extra_arguments = ['--pylab=inline']

    def __enter__(self):
        self.km = KernelManager()
        self.km.start_kernel(extra_arguments=self.extra_arguments, stderr=open(os.devnull, 'w'))

        if platform.system() == 'Dwarin':
            sleep(1)

        self.kc = self.km.client()
        self.kc.start_channels()
        try:
            self.kc.wait_for_ready()
        except AttributeError:
            # Ipython < 3
            self._wait_for_ready_backport()

        # run %pylab inline, because some notebooks assume this
        # even though they shouldn't
        self.kc.execute("pass")
        reply = self.kc.get_shell_msg()        
        while True:
            try:
                msg = self.kc.get_iopub_msg(timeout=1)
            except Empty:
                break

        return self

    def _wait_for_ready_backport(self):
        """Backport BlockingKernelClient.wait_for_ready from IPython 3"""
        # Wait for kernel info reply on shell channel
        self.kc.kernel_info()
        while True:
            msg = self.kc.get_shell_msg(block=True, timeout=30)
            if mesg['msg_type'] == 'kernel_info_reply':
                break

        # Flush IOPub channel
        while True:
            try:
                msg = self.kc.get_iopub_msg(block=True, timeout=0.2)
            except Empty:
                break

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.kc.stop_channels()
        self.km.shutdown_kernel()
        del self.km

    def run(self, cell, timeout = None):
        use_timeout = self.default_timeout
        if timeout is not None:
            use_timeout = timeout
        self.kc.execute(cell.input)
        reply = self.kc.get_shell_msg(timeout=use_timeout)
        outs = []

        while True:
            try:
                msg = self.kc.get_iopub_msg(timeout=0.5)
            except Empty:
                break
            msg_type = msg['msg_type']
            if msg_type in ('status', 'pyin', 'execute_input'): # pyin renamed to execute_input in ipython message specification 5.0
                continue
            elif msg_type == 'clear_output':
                outs = []
                continue

            content = msg['content']
            out = NotebookNode(output_type=msg_type)

            if msg_type == 'stream':
                out.stream = content['name']
                if 'text' in content:
                    out.text = content['text']
                else:
                    out.text = content['data']
            elif msg_type in ('display_data', 'pyout'):
                out['metadata'] = content['metadata']
                for mime, data in content['data'].items():
                    attr = mime.split('/')[-1].lower()
                    # this gets most right, but fix svg+html, plain
                    attr = attr.replace('+xml', '').replace('plain', 'text')
                    setattr(out, attr, data)
                if msg_type == 'pyout':
                    out.prompt_number = content['execution_count']
            elif msg_type in ('pyerr', 'error'): # pyerr changed to error in ipython msg spec 5.0
                out.ename = content['ename']
                out.evalue = content['evalue']
                out.traceback = content['traceback']
            else:
                print("unhandled iopub msg:", msg_type)

            outs.append(out)
        return outs

    def sanitize(self, s):
        """sanitize a string for comparison.

        fix universal newlines, strip trailing newlines, and normalize likely random values (memory addresses and UUIDs)
        """
        if not isinstance(s, (str, bytes)):
            return s
        # normalize newline:
        s = s.replace('\r\n', '\n')

        # ignore trailing newlines (but not space)
        s = s.rstrip('\n')

        # normalize hex addresses:
        s = re.sub(r'0x[a-f0-9]+', '0xFFFFFFFF', s)

        # normalize UUIDs:
        s = re.sub(r'[a-f0-9]{8}(\-[a-f0-9]{4}){3}\-[a-f0-9]{12}', 'U-U-I-D', s)

        return s

    def compare_outputs(self, test, ref, skip_compare=('png', 'traceback', 'latex', 'prompt_number', 'svg', 'html')):
        for key in ref:
            if key not in test:
                if printdiffs: 
                    print("missing key: %s != %s" % (test.keys(), ref.keys()))
                return False
            elif key not in skip_compare:
                s1 = self.sanitize(test[key])
                s2 = self.sanitize(ref[key])
                if s1 != s2:
                    if printdiffs:
                        print("mismatch %s:" % key)
                    expected=s1.splitlines(1)
                    actual=s2.splitlines(1)
                    if printdiffs:
                        diff=difflib.unified_diff(expected, actual)

                    if printdiffs:
                        print(''.join(diff))
                    return False
        return True

    def get_commands(self, cell):
        commands = {}
        if hasattr(cell, 'input'):
            lines = cell.input.splitlines()
            if len(lines) > 0:
                first_line = lines[0]
                if first_line.startswith('#!'):
                    txt = first_line[2:].strip()

                    parts = txt.split(',')
                    for part in parts:
                        subparts = part.split(':')
                        if len(subparts) == 1:
                            commands[subparts[0].strip().lower()] = True
                        elif len(subparts) == 2:
                            commands[subparts[0].strip().lower()] = subparts[1]

        return commands

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run all cells in an ipython notebook as a test and check whether these successfully execute and ' +
                    'optionally compare their output to the ones found inside the notebook.')

    parser.add_argument('file', metavar='file.ipynb', help='the notebook to be checked.', type=str)

    parser.add_argument('--timeout', dest='timeout',
                   type=int, default=300,
                   help='the default timeout time in seconds for a cell evaluation. Default is 300s.')

    parser.add_argument('--strict', dest='strict', action='store_true',
                   default=False,
                   help='if set to true then the default test is that cells have to have matching output' +
                        'otherwise only a fail in execution will be considered a failed test.')

    parser.add_argument('--fail-if-timeout', dest='no_timeout', action='store_true',
                   default=False,
                   help='if set to true then a timeout is considered a failed test.')

    parser.add_argument('--diff', dest='diff', action='store_true',
                    default=False, 
                    help='if set to true, will output the diff when the test\'s cells do not match the reference outputs.')

    args = parser.parse_args()
    ipynb = args.file

    tv = IPyTestConsole()

    if args.diff:
        printdiffs = True

    if args.strict:
        tv.default_results['diff'] = False

    if args.no_timeout:
        tv.default_results['timeout'] = False

    tv.writeln('testing ipython notebook : "%s"' % ipynb)
    tv.write("starting kernel ... ")

    with open(ipynb) as f:
        nb = reads(f.read(), 'json')

    with IPyKernel() as ipy:
        ipy.default_timeout = args.timeout
        tv.writeln("ok")

        nbs = ipynb.split('.')
        nb_class_name = nbs[1] + '.' + nbs[0].replace(" ", "_")

        tv.br()

        for ws in nb.worksheets:
            for cell in ws.cells:
                if cell.cell_type == 'markdown':
                    for line in cell.source.splitlines():
                        # only tv.writeln(headlines in markdown
                        if line.startswith('#'):
                            tv.writeln(line)

                if cell.cell_type == 'heading':
                    tv.writeln('#' * cell.level + ' ' + cell.source)

                if cell.cell_type != 'code':
                    continue

                # If code cell then continue with checking it

                if hasattr(cell, 'prompt_number'):
                    tv.write(nb_class_name + '.' + 'In [%3i]' % cell.prompt_number + ' ... ')
                else:
                    tv.write(nb_class_name + '.' + 'In [???]' + ' ... ')

                commands = ipy.get_commands(cell)

                result = 'success'

                timeout = ipy.default_timeout

                if 'skip' in commands:
                    tv.write_result('skip')
                    continue

                try:
                    if 'timeout' in commands:
                        outs = ipy.run(cell, timeout=int(commands['timeout']))
                    else:
                        outs = ipy.run(cell)

                except Exception as e:
                    # Internal IPython error occurred (might still be that the cell did not execute correctly)
                    if repr(e) == 'Empty()':
                        # Assume it has been timed out!
                        tv.write_result('timeout')
                        # tv.writeln('>>> TimeOut (%is)' % args.timeout)
                    else:
                        tv.write_result('kernel')
                        tv.fold_open('ipynb.kernel')
                        tv.writeln('>>> ' + out.ename + ' ("' + out.evalue + '")')
                        tv.writeln(repr(e), indent=4)
                        tv.fold_close('ipynb.kernel')
                    continue

                failed = False
                diff = False

                if outs != [] and cell.outputs == []: # check case where ref has no output
                    for out in outs:
                        if out.output_type in ('pyerr', 'error'): 
                            # An python error occurred. Cell is not completed correctly
                            tv.write_result('error')
                            tv.fold_open('ipynb.fail')
                            tv.write('>>> ' + out.ename + ' ("' + out.evalue + '")')
                            for idx, trace in enumerate(out.traceback):
                                tv.writeln(trace, indent=4)
                            tv.fold_close('ipynb.fail')
                            failed = True
                        else:
                            if not ipy.compare_outputs(out, ref):
                                # Output is different than the one in the notebook. This might be okay.
                                diff = True

                for out, ref in zip(outs, cell.outputs):
                    if out.output_type in ('pyerr', 'error'): 
                        # An python error occurred. Cell is not completed correctly
                        tv.write_result('error')
                        tv.fold_open('ipynb.fail')
                        tv.write('>>> ' + out.ename + ' ("' + out.evalue + '")')
                        for idx, trace in enumerate(out.traceback):
                            tv.writeln(trace, indent=4)
                        tv.fold_close('ipynb.fail')
                        failed = True
                    else:
                        if not ipy.compare_outputs(out, ref):
                            # Output is different than the one in the notebook. This might be okay.
                            diff = True
                
                if diff:
                    if 'strict' in commands:
                        # strict mode means a difference will fail the test
                        tv.write_result('diff', okay_list={ 'diff' : False })
                    elif 'ignore' in commands:
                        # ignore mode means a difference will pass the test
                        tv.write_result('diff', okay_list={ 'diff' : True })
                    else:
                        # use defaults
                        tv.write_result('diff')

                if not failed and not diff:
                    tv.write_result('success')

        tv.br()
        tv.writeln("testing results")
        tv.writeln("===============")
        if tv.pass_count > 0:
            tv.writeln("    %3i cells passed [" % tv.pass_count + tv.green('ok') + "]" )
        if tv.fail_count > 0:
            tv.writeln("    %3i cells failed [" % tv.fail_count + tv.red('fail') + "]" )

        tv.br()
        tv.writeln("    %3i cells have been successfully replicated [success]" % tv.result_count['success'])
        tv.writeln("    %3i cells had mismatched outputs [diff]" % tv.result_count['diff'])
        tv.writeln("    %3i cells timed out during execution [time]" % tv.result_count['timeout'])
        tv.writeln("    %3i cells ran with python errors [fail]" % tv.result_count['error'])
        tv.writeln("    %3i cells have been executed without comparison [ignore]" % tv.result_count['ignore'])
        tv.writeln("    %3i cells failed to even execute (IPython error) [kernel]" % tv.result_count['kernel'])
        tv.writeln("    %3i cells have been skipped [skip]" % tv.result_count['skip'])

        tv.br()
        tv.write("shutting down kernel ... ")

    tv.writeln('ok')

    if tv.fail_count != 0:
        tv.writeln(tv.red('some tests not passed.'))
        exit(1)
    else:
        tv.writeln(tv.green('all tests passed.'))
        exit(0)
