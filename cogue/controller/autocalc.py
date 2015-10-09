import os
import time
import yaml
import datetime
from cogue.task import TaskSet
from cogue.qsystem.queue import EmptyQueue

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
        self._cwd = None
        self._verbose = verbose
        self._dot_count = 1
        self._log = []

    def set_queue(self, queue):
        self._queue = queue

    def append(self, directory, task):
        taskset = TaskSet(directory)
        taskset.append(task)
        self._taskset.append(taskset)

    def get_tasks(self):
        return [x.get_tasks()[0] for x in self._taskset.get_tasks()]

    def run(self, check_period=10):  # in second
        self._begin()

        while True:
            self._queue.qstat()
            time.sleep(check_period)
            self._log.append("-" * 40 + "> %s" % date())
            self._deep_run(self._taskset)
            self._log.append("<" + "-" * 40)
            self._overwrite_settings()
            self._write_log()
            self._write_dot()
            self._write_qstatus()
            if self._taskset.done():
                break

        self._end()

    def _begin(self):
        self._cwd = os.getcwd()
        self._deep_begin(self._taskset)

    def _end(self):
        os.chdir(self._cwd)

    def _deep_begin(self, task):
        directory = task.get_directory()
        if directory is not None:
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
            for next_taskset in task:
                for next_task in next_taskset:
                    self._deep_begin(next_task)
                break

        log = task.get_log().rstrip()
        if log:
            self._log.append(log)
            task.set_log("")

        self._chdir_out(orig_cwd, task.get_status())

    def _chdir_in(self, directory_in):
        if directory_in is None:
            return None
        else:
            cwd = os.getcwd()
            os.chdir(directory_in)
            self._log.append("--> %s" %
                             os.getcwd().replace(self._cwd, '').lstrip('/'))
            return cwd

    def _chdir_out(self, cwd, status): 
       if cwd is not None:
           self._log.append("        [ %s ]" % status)
           directory = cwd.replace(self._cwd, '').lstrip('/')
           if directory == "":
               directory = '.'
           self._log.append("    %s <--" % directory)
           os.chdir(cwd)

    def _overwrite_settings(self):
        filename = "%s.cogue" % self._name
        if os.path.exists(filename):
            print "%s is found." % filename

            with open(filename) as f:
                data = yaml.load(f)
                if 'max_jobs' in data:
                    max_jobs = data['max_jobs']
                    self._queue.set_max_jobs(max_jobs)
                    print "Overwrite max number of jobs by %d." % max_jobs

            print "File %s was renamed to %s.done." % (filename, filename)
            if os.path.exists("%s.done" % filename):
                os.remove("%s.done" % filename)
            os.rename("%s" % filename, "%s.done" % filename)

    def _write_log(self):
        if self._verbose > 1:
            with open("%s.log" % self._log_name, 'a') as w:
                w.write("\n".join(self._log))
                self._log = []

    def _write_dot(self):
        if self._verbose:
            with open("%s.dot" % self._log_name, 'w') as f_dot:
                self._dot_count += 1
                f_dot.write("digraph %s {\n" %
                            self._name.replace('-', '_').replace('.', '_'))
                f_dot.write("graph [ rankdir = \"LR\" ] ;\n")
                for task in self._taskset.get_tasks():
                    self._write_dot_labels(task, f_dot)
                    self._write_dot_tids(task, f_dot)
                f_dot.write("}\n")

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
        if status is None:
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
