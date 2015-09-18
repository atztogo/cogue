class JobBase:
    def __init__(self):
        self._jobname = None
        self._status = None

    def set_jobname(self, jobname):
        self._jobname = jobname

    def get_jobname(self):
        return self._jobname

    def set_status(self, status, jobid=None):
        if jobid:
            self._status = "%s (id:%d)" % (status, jobid)
        else:
            self._status = status
    
    def get_status(self):
        return self._status
    
    
