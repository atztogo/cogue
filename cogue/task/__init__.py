import os

class TaskBase:
    def __init__(self):
        self._status = ""
        self._comment = ""
        self._log = ""
        self._tid = None
    
    def __iter__(self):
        return self

    def done(self):
        return True

    def begin(self):
        pass

    def next(self):
        raise StopIteration

    def end(self):
        pass

    def get_name(self):
        return self._name

    def get_type(self):
        return self._task_type

    def get_directory(self):
        return self._directory

    def get_log(self):
        return self._log

    def get_status(self):
        return self._status

    def get_comment(self):
        return self._comment

    def get_tasks(self):
        return self._tasks

    def set_tid(self, tid):
        self._tid = tid
        
    def get_tid(self):
        return self._tid


class TaskElement(TaskBase):
    def __init__(self):

        TaskBase.__init__(self)
        self._job = None
        self._traverse = False # Do not submit job

    def set_job(self, job):
        if (isinstance(self._job, list) or
            isinstance(self._job, tuple)):
            self._job = list(job)
        else:
            self._job = job
        
    def get_job(self):
        return self._job

    def set_traverse(self, traverse):
        self._traverse = traverse

    def get_traverse(self):
        return self._traverse

    def _overwrite_settings(self):
        if os.path.exists(".coguerc"):
            import yaml
            data = yaml.load(open(".coguerc"))
            if 'max_iteration' in data:
                self._max_iteration = data['max_iteration']
            if 'min_iteration' in data:
                self._min_iteration = data['min_iteration']
            if 'traverse' in data:
                self._traverse = data['traverse']

class TaskSet(TaskBase):
    def __init__(self,
                 directory=None,
                 name=None):

        TaskBase.__init__(self)

        self._directory = directory
        if not name:
            self._name = directory
        else:
            self._name = name
        self._tasks = []
        self._task_type = "task_set"

    def append(self, task):
        self._tasks.append(task)

    def set_status(self):
        done = True
        for task in self._tasks:
            done &= task.done()
        if done:
            self._status = "done"
        self._write_yaml()

    def done(self):
        return (self._status == "done")

    def begin(self):
        self._status = "begin"

    def _write_yaml(self):
        if self._directory:
            w = open("%s.yaml" % self._directory, 'w')
            w.write("status: %s\n" % self._status)
            w.write("tasks:\n")
            for task in self._tasks:
                if task.get_status():
                    w.write("- name:   %s\n" % task.get_name())
                    w.write("  status: %s\n" % task.get_status())
            w.close()
