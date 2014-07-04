class EmptyQueue:
    def __init__(self):
        pass

    def register(self, task):
        pass

    def submit(self, task):
        pass

    def qstat(self):
        pass

    def set_max_jobs(self, max_jobs):
        pass

    def write_qstatus(self, name):
        pass

class QueueBase:
    def __init__(self, max_jobs=None):
        self._max_jobs = max_jobs
        self._qstatus = None
        self._tid_queue = []
        self._tid2jobid = {}
        self._shell = None
        self._shell_type = None

    def register(self, task):
        self._tid_queue.append(task.get_tid())
        job = task.get_job()
        job.set_status("preparing")

    def write_qstatus(self, name):
        f_qstat = open("%s.qstat" % name, 'w')
        f_qstat.write("%8s %8s %8s\n" % ("tid", "jobid", "status"))
        for tid in self._tid_queue:
            f_qstat.write("%8d %8s %8s\n" % (tid, 'None', 'Queued'))
            
        for tid, jobid in self._tid2jobid.iteritems():
            if jobid in self._qstatus:
                f_qstat.write("%8d %8d %8s\n" %
                              (tid, jobid, self._qstatus[jobid]))
        f_qstat.close()

    def submit(self):
        """To be implemented in specific queue"""
        pass
        
    def set_max_jobs(self, max_jobs):
        self._max_jobs = max_jobs

    def _set_job_status(self, job, tid):
        if "preparing" in job.get_status():
            if tid == self._tid_queue[0]:
                will_submit = True
                if self._max_jobs:
                    if len(self._tid2jobid) > self._max_jobs:
                        will_submit = False

                if will_submit:
                    job.set_status("ready")
        else:
            jobid = self._tid2jobid[tid]
            if jobid in self._qstatus:
                if self._qstatus[jobid] == 'Running':
                    job.set_status("running", jobid)
            else:
                del self._tid2jobid[tid]
                job.set_status("done")

