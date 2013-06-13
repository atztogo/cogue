__all__ = ['queue', 'job']

import subprocess
import shlex
import spur
import os
import shutil
import tarfile

def queue(max_jobs=None,
          ssh_shell=None,
          temporary_dir=None):
    if ssh_shell is None:
        return LocalQueue(max_jobs=max_jobs)
    elif temporary_dir is not None:
        return RemoteQueue(ssh_shell,
                           temporary_dir,
                           max_jobs=max_jobs)
        

def job(script=None,
        shell=None,
        cwd=True,
        jobname=None,
        q=None,
        l=None,
        pe=None,
        stdout=None,
        stderr=None):
    return Job(script=script,
               shell=shell,
               cwd=cwd,
               jobname=jobname,
               q=q,
               l=l,
               pe=pe,
               stdout=stdout,
               stderr=stderr)

class Queue:
    """Store and control jobs

    1. A job is registered. Task-ID is used as the identifier.
       --> Task-ID is pushed to self._tid_queue. [preparing]

    2. The job is submitted to queueing system if number of submitted
       jobs are less then specified max number of jobs. [submitted]
       --> Task-ID is removed from self._tid_queue.
       --> Job-ID is mapped to the task-ID by self._tid2jobid.

    3. If the job-ID is not in the list of job-IDs in queueing system,
       the job is recognized as finished. [done]
       --> The key of task-ID in self._tid2jobid is removed.

    4. If the job-ID is in the list of job-IDs in queueing system and
       the status of the job in queueing system is 'R':
       the job is recognized as running. [running]

    """
    def __init__(self,
                 max_jobs = None,
                 qsub_command="qsub"):
        self._max_jobs = max_jobs
        self._qsub_command = qsub_command
        self._qstatus = None
        self._tid_queue = []
        self._tid2jobid = {}
        self._shell = None
        self._shell_type = None

    def register(self, task):
        self._tid_queue.append(task.get_tid())
        job = task.get_job()
        job.set_status("preparing")

    def get_qstatus(self):
        return self._qstatus

    def qstat(self):
        """qstatout

        Text of output of 'qstat'
        If this is None, qstat is exectued in python.

        """
        qstat_out = self._shell.run(["qstat"]).output.split('\n')
        self._qstatus = {}

        # subprocess.check_out is new in python 2.7
        # for line in subprocess.check_output(["qstat"]).split('\n'):
        for line in qstat_out:
            if len(line.split()) > 5:
                jobid = line.split()[0]
                if jobid.isdigit():
                    jobid = int(jobid)
                    s = line[40:44].strip()
                    self._qstatus[jobid] = s
                    if s == 'r':
                        self._qstatus[jobid] = 'Running'
                    elif s == 'qw':
                        self._qstatus[jobid] = 'Pending'

    
    def set_max_jobs(self, max_jobs):
        self._max_jobs = max_jobs

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

    def _upload_files(self):
        pass

    def _download_files(self):
        pass
        

class RemoteQueue(Queue):
    def __init__(self,
                 ssh_shell,
                 temporary_dir,
                 max_jobs=None,
                 qsub_command="qsub"):
        Queue.__init__(self,
                       max_jobs=max_jobs,
                       qsub_command=qsub_command)
        self._shell = ssh_shell
        self._temporary_dir = temporary_dir

    def submit(self, task):
        job = task.get_job()
        tid = task.get_tid()
        remote_dir = "%s/c%05d" % (self._temporary_dir, tid)
        self._set_job_status(job, tid)
        if "ready" in job.get_status():
            job.write_script()
            self._shell.run(["mkdir", "-p", remote_dir])
            tar = tarfile.open("cogue.tar", "w")
            for name in os.listdir("."):
                tar.add(name)
            tar.close()
            with open("cogue.tar", "rb") as local_file:
                with self._shell.open("%s/%s" % (remote_dir, "cogue.tar"), "wb") as remote_file:
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
                jobid = int(qsub_out.split()[2])
            self._tid2jobid[tid] = jobid
            self._tid_queue.pop(0)
            job.set_status("submitted", jobid)

        elif "done" in job.get_status():
            names = self._shell.run(["/bin/ls"], cwd=remote_dir).output.split()
            self._shell.run(["tar", "cvf", "cogue.tar"] + names, cwd=remote_dir)
            with self._shell.open("%s/%s" % (remote_dir, "cogue.tar"), "rb") as remote_file:
                with open("cogue.tar", "wb") as local_file:
                    shutil.copyfileobj(remote_file, local_file)
                    tar = tarfile.open("cogue.tar")
                    tar.extractall()
                    tar.close()
                    os.remove("cogue.tar")
                    self._shell.run(["rm", "cogue.tar"], cwd=remote_dir)

class LocalQueue(Queue):
    def __init__(self,
                 max_jobs=None,
                 qsub_command="qsub"):
        Queue.__init__(self,
                       max_jobs=max_jobs,
                       qsub_command=qsub_command)
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
                jobid = int(qsub_out.split()[2])
            self._tid2jobid[tid] = jobid
            self._tid_queue.pop(0)
            job.set_status("submitted", jobid)

class Job:
    def __init__(self,
                 script=None,
                 shell=None,
                 cwd=True,
                 jobname=None,
                 q=None,
                 l=None,
                 pe=None,
                 stdout=None,
                 stderr=None):

        if not script:
            print "Command not found"
            return False
        else:
            self._script = script

        if not shell:
            self._shell = "/bin/bash"
        else:
            self._shell = shell

        if not jobname:
            self._jobname = "cogue-job"
        else:
            self._jobname = jobname

        self._pe = pe
        self._q = q
        self._l = l

        self._cwd = cwd
        self._stdout = stdout
        self._stderr = stderr

        self._status = ""

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
    
    def copy(self, jobname):
        return Job(script=self._script,
                   shell=self._shell,
                   cwd=self._cwd,
                   jobname=jobname,
                   q=self._q,
                   l=self._l,
                   pe=self._pe,
                   stdout=self._stdout,
                   stderr=self._stderr)

    def write_script(self, filename="job.sh"):
        """
        #$ -S /bin/zsh
        #$ -cwd
        #$ -N NaCl
        #$ -pe mpi* 4
        #$ -q lowspeed
        #$ -l exclusive=false
        #$ -e err.log
        #$ -o std.log

        mpirun vasp5212mpi
        """
        
        w = open(filename, 'w')
        if not self._shell:
            w.write("#$ -S %s\n" % "/bin/bash")
        else:
            w.write("#$ -S %s\n" % self._shell)
        if self._cwd:
            w.write("#$ -cwd\n")
        if self._jobname:
            w.write("#$ -N %s\n" % self._jobname)
        if self._q:
            w.write("#$ -q %s\n" % self._q)
        if self._l:
            w.write("#$ -l %s\n" % self._l)
        if self._pe:
            w.write("#$ -pe %s\n" % self._pe)
        if self._stderr:
            w.write("#$ -e %s\n" % self._stderr)
        if self._stdout:
            w.write("#$ -o %s\n" % self._stdout)
        
        w.write("\n")

        w.write(self._script)

        w.close()
