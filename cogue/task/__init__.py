import os

class TaskBase:
    def __init__(self):
        self._status = ""
        self._comment = ""
        self._log = ""
        self._tid = None
        self._name = None
        self._task_type = None
        self._directory = None
        self._tasks = None
    
    def __iter__(self):
        return self

    def done(self):
        return True

    def begin(self):
        pass

    def next(self):
        raise StopIteration

    def get_name(self):
        return self._name

    def get_type(self):
        return self._task_type

    def get_directory(self):
        return self._directory

    def set_log(self, log):
        self._log = log

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

    def get_yaml_lines(self):
        lines = []
        if self._name:
            lines.append("name:      %s" % self._name)
        if self._task_type:
            lines.append("type:      %s" % self._task_type)
        if self._tid:
            lines.append("id:        %s" % self._tid)
        if self._status:
            lines.append("status:    %s" % self._status)
        if self._directory:
            lines.append("directory: %s" % self._directory)
        return lines

    def __str__(self):
        return "\n".join(self.get_yaml_lines())

    def _write_yaml(self, filename=None):
        if filename:
            w = open(filename)
        else:
            w = open("%s.yaml" % self._task_type, 'w')
        w.write("\n".join(self.get_yaml_lines()))
        w.close()

class TaskElement(TaskBase):
    def __init__(self):
        TaskBase.__init__(self)
        self._job = None
        self._traverse = False # Do submit job

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

    def get_yaml_lines(self):
        lines = TaskBase.get_yaml_lines(self)
        if self._traverse is True:
            lines.append("traverse:  True")
        elif self._traverse is not False:
            lines.append("traverse:  %s" % self._traverse)
        return lines

    def _overwrite_settings(self):
        if os.path.exists(".coguerc"):
            import yaml
            data = yaml.load(open(".coguerc"))
            if 'max_iteration' in data:
                if '_max_iteration' in self.__dict__:
                    self._max_iteration = data['max_iteration']
            if 'min_iteration' in data:
                if '_min_iteration' in self.__dict__:
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
        terminate = False
        for task in self._tasks:
            done &= task.done()
            if task.get_status() == "terminate":
                terminate = True
        if done:
            if terminate:
                self._status = "terminate"
            else:
                self._status = "done"
        self._write_yaml()

    def done(self):
        return (self._status == "done" or
                self._status == "terminate")

    def begin(self):
        self._status = "begin"

    def get_yaml_lines(self):
        lines = TaskBase.get_yaml_lines(self)
        if self._tasks:
            lines.append("tasks:")
        for task in self._tasks:
            if not task:
                continue
            for i, line in enumerate(TaskBase.get_yaml_lines(task)):
                if i == 0:
                    lines.append("- " + line)
                else:
                    lines.append("  " + line)

        return lines

        
