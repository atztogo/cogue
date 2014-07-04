import os
import time
import datetime
from cogue.task import TaskSet
from cogue.qsystem import EmptyQueue

def date():
    return datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")

class AutoCalc:
    def __init__(self, name=None, log_name=None, verbose=False):
        if name is None:
            self._name = "autocalc"
            self._taskset = TaskSet(name="autocalc")
        else:
            self._name = name
            self._taskset = TaskSet(directory=name)

        if log_name is None:
            self._log_name = self._name
        else:
            self._log_name = log_name

        self._queue = EmptyQueue()
        self._tid_count = 0
        self._f_log = None
        self._cwd = None
        self._verbose = verbose
        self._dot_count = 1

    def set_queue(self, queue):
        self._queue = queue

    def append(self, directory, task):
        taskset = TaskSet(directory)
        taskset.append(task)
        self._taskset.append(taskset)

    def run(self, check_period=10):  # in second
        self._begin()

        while True:
            time.sleep(check_period)
            self._write_log("-" * 40 + "> %s\n" % date())
            self._queue.qstat()
            self._deep_run(self._taskset)
            self._overwrite_settings()
            self._write_log("<" + "-" * 40 + "\n")
            self._write_dot()
            self._write_qstatus()
            if self._taskset.done():
                break

        self._end()

    def _begin(self):
        if self._verbose:
            self._f_log = open("%s.log" % self._log_name, 'w')
        self._cwd = os.getcwd()
        self._deep_begin(self._taskset)

    def _end(self):
        if self._verbose:
            self._f_log.close()
        os.chdir(self._cwd)

    def _deep_begin(self, task):
        directory = task.get_directory()
        if not directory == None:
            if not os.path.exists(directory):
                os.mkdir(directory)

        cwd = self._chdir_in(directory)

        task.set_tid(self._tid_count)
        self._tid_count += 1
        task.begin()
        subtasks = task.get_tasks()
        if subtasks: # Task-set
            for subtask in task.get_tasks():
                self._deep_begin(subtask)
        else: # Execution task
            self._queue.register(task)

        self._chdir_out(cwd, task.get_status())
        
    def _deep_run(self, task):
        orig_cwd = self._chdir_in(task.get_directory())

        subtasks = task.get_tasks()
        if subtasks: # Task-set
            for subtask in subtasks:
                if not subtask.done():
                    self._deep_run(subtask)

        else: # Execution task
            self._queue.submit(task)
            
        task.set_status()
        if task.done():
            try:
                next_subtasks = task.next()
            except StopIteration:
                task.end()
            else:
                for next_subtask in next_subtasks:
                    self._deep_begin(next_subtask)

        if task.get_log():
            self._write_log(task.get_log() + "\n")
        self._chdir_out(orig_cwd, task.get_status())

    def _chdir_in(self, directory_in):
        if directory_in == None:
            return None
        else:
            cwd = os.getcwd()
            os.chdir(directory_in)
            self._write_log("--> %s\n" %
                            os.getcwd().replace(self._cwd, '').lstrip('/'))
            return cwd

    def _chdir_out(self, cwd, status): 
       if not cwd == None:
           self._write_log("        [ %s ]\n" % status)
           directory = cwd.replace(self._cwd, '').lstrip('/')
           if directory == "":
               directory = '.'
           self._write_log("    %s <--\n" % directory)
           os.chdir(cwd)

    def _overwrite_settings(self):
        filename = "%s.cogue" % self._name
        if os.path.exists(filename):
            print "%s is found." % filename
            import yaml
            data = yaml.load(open(filename))
            if 'max_jobs' in data:
                self._queue.set_max_jobs(data['max_jobs'])
                print "Overwrite max number of jobs by %d." % data['max_jobs']

            print "File %s was renamed to %s.done." % (filename, filename)
            if os.path.exists("%s.done" % filename):
                os.remove("%s.done" % filename)
            os.rename("%s" % filename, "%s.done" % filename)

    def _write_log(self, log):
        if self._verbose > 1:
            self._f_log.write(log)

    def _write_dot(self):
        if self._verbose:
            f_dot = open("%s.dot" % self._log_name, 'w')
            self._dot_count += 1
            f_dot.write("digraph %s {\n" %
                        self._name.replace('-', '_').replace('.', '_'))
            f_dot.write("graph [ rankdir = \"LR\" ] ;\n")
            for task in self._taskset.get_tasks():
                self._write_dot_labels(task, f_dot)
                self._write_dot_tids(task, f_dot)
            f_dot.write("}\n")
            f_dot.close()

    def _write_dot_labels(self, task, f_dot):
        tid = task.get_tid()
        name = task.get_name()
        status = task.get_status()
        comment = task.get_comment()
        if status == "done":
            color = "lightblue2"
        elif status == "terminate":
            color = "pink"
        else:
            color = "khaki1"

        f_dot.write("n%d ;\n" % tid)
        if status == None:
            f_dot.write("n%d [label=\"%s\"] ;\n" % (tid, name))
        else:
            if comment:
                f_dot.write("n%d [color=%s, style=filled, "
                            "shape=box, label=\"[%d] %s\\n%s\\n%s\"] ;\n" %
                            (tid, color, tid, status, name, comment))
            else:
                f_dot.write("n%d [color=%s, style=filled, "
                            "shape=box, label=\"[%d] %s\\n%s\"] ;\n" %
                            (tid, color, tid, status, name))
        if task.get_tasks():
            for t in task.get_tasks():
                self._write_dot_labels(t, f_dot)

    def _write_dot_tids(self, task, f_dot):
        tid = task.get_tid()
        if task.get_tasks():
            for t in task.get_tasks():
                f_dot.write("n%d -> n%d ;\n" % (tid, t.get_tid()))
                self._write_dot_tids(t, f_dot)
    
    def _write_qstatus(self):
        if self._verbose:
            self._queue.write_qstatus(self._log_name)
