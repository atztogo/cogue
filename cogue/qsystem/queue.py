import os
import shlex
import shutil
import tarfile

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

class LocalQueueBase(QueueBase):
    def __init__(self,
                 max_jobs=None,
                 qsub_command="qsub"):
        try:
            import spur
        except ImportError:
            print "You need to install spur."
            exit(1)
        
        QueueBase.__init__(self, max_jobs=max_jobs)
        self._qsub_command = qsub_command
        self._shell = spur.LocalShell()

    def submit(self, task):
        job = task.get_job()
        tid = task.get_tid()
        self._set_job_status(job, tid)
        if "ready" in job.get_status():
            job.write_script()
            if task.get_traverse():
                jobid = None
            else:
                qsub_out = self._shell.run(
                    shlex.split(self._qsub_command + " " + "job.sh"),
                    cwd=os.getcwd()).output
                jobid = self._get_jobid(qsub_out)
            self._tid2jobid[tid] = jobid
            self._tid_queue.pop(0)
            job.set_status("submitted", jobid)


class RemoteQueueBase(QueueBase):
    def __init__(self,
                 ssh_shell,
                 temporary_dir,
                 max_jobs=None,
                 name=None,
                 qsub_command="qsub"):
        QueueBase.__init__(self, max_jobs=max_jobs)
        self._qsub_command = qsub_command
        self._shell = ssh_shell
        self._name = name

        if not os.path.exists(temporary_dir):
            self._shell.run(["mkdir", "-p", temporary_dir])
        
        if self._name is not None:
            self._working_dir = "%s/%s" % (temporary_dir, self._name)
            if not os.path.exists(self._working_dir):
                self._shell.run(["mkdir", "-p", self._working_dir])
        else:
            self._working_dir = "%s" % temporary_dir

    def submit(self, task):
        job = task.get_job()
        tid = task.get_tid()
        remote_dir = "%s/c%05d" % (self._working_dir, tid)
        self._set_job_status(job, tid)
        if "ready" in job.get_status():
            job.write_script()
            self._shell.run(["mkdir", "-p", remote_dir])
            tar = tarfile.open("cogue.tar", "w")
            for name in os.listdir("."):
                tar.add(name)
            tar.close()
            with open("cogue.tar", "rb") as local_file:
                with self._shell.open("%s/%s" % (remote_dir, "cogue.tar"),
                                      "wb") as remote_file:
                    shutil.copyfileobj(local_file, remote_file)
                    os.remove("cogue.tar")
                    self._shell.run(["tar", "xvf", "cogue.tar"], cwd=remote_dir)
                    self._shell.run(["rm", "cogue.tar"], cwd=remote_dir)
            if task.get_traverse():
                jobid = None
            else:
                qsub_out = self._shell.run(
                    shlex.split(self._qsub_command + " " + "job.sh"),
                    cwd=remote_dir).output
                jobid = self._get_jobid(qsub_out)
            self._tid2jobid[tid] = jobid
            self._tid_queue.pop(0)
            job.set_status("submitted", jobid)

        elif "done" in job.get_status():
            names = self._shell.run(["/bin/ls"], cwd=remote_dir).output.split()
            self._shell.run(["tar", "cvf", "cogue.tar"] + names, cwd=remote_dir)
            with self._shell.open("%s/%s" % (remote_dir, "cogue.tar"),
                                  "rb") as remote_file:
                with open("cogue.tar", "wb") as local_file:
                    shutil.copyfileobj(remote_file, local_file)
                    tar = tarfile.open("cogue.tar")
                    tar.extractall()
                    tar.close()
                    os.remove("cogue.tar")
                    self._shell.run(["rm", "cogue.tar"], cwd=remote_dir)

