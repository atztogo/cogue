__all__ = ['queue', 'job']

import subprocess
import shlex

def queue(max_jobs=None):
    return Queue(max_jobs=max_jobs)

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
        self._qsub_command = qsub_command

        self._qstatus = None
        self._tid_queue = []
        self._tid2jobid = {}
        self._max_jobs = max_jobs

    def qstat(self):
        """qstatout

        Text of output of 'qstat'
        If this is None, qstat is exectued in python.

        """
        self._qstatus = {}

        # subprocess.check_out is new in python 2.7
        # for line in subprocess.check_output(["qstat"]).split('\n'):
        for line in subprocess.Popen(["qstat"], stdout=subprocess.PIPE).communicate()[0].split('\n'):
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

    def register(self, task):
        self._tid_queue.append(task.get_tid())
        job = task.get_job()
        job.set_status("preparing")

    def get_qstatus(self):
        return self._qstatus

    def set_job_status(self, task):
        job = task.get_job()
        tid = task.get_tid()
        if "preparing" in job.get_status():
            if tid == self._tid_queue[0]:
                if self._max_jobs:
                    if len(self._tid2jobid) < self._max_jobs:
                        self._submit(tid, job, task.get_traverse())
                else:
                    self._submit(tid, job, task.get_traverse())
        else:
            jobid = self._tid2jobid[tid]
            if jobid in self._qstatus:
                if self._qstatus[jobid] == 'Running':
                    job.set_status("running", jobid)
            else:
                del self._tid2jobid[tid]
                job.set_status("done")

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

    def set_max_jobs(self, max_jobs):
        self._max_jobs = max_jobs
        
    def _submit(self, tid, job, traverse=False):
        if traverse:
            jobid = None
        else:
            jobid = self._qsub(job)
        self._tid2jobid[tid] = jobid
        self._tid_queue.pop(0)
        job.set_status("submitted", jobid)

    def _qsub(self, job, filename="job.sh"):
        # subprocess.check_out is new in python 2.7
        # stdout = subprocess.check_output(shlex.split(
        #         self._qsub_command + " " + filename))
        job.write_script()
        stdout = subprocess.Popen(shlex.split(self._qsub_command + " " + filename), stdout=subprocess.PIPE).communicate()[0]
        jobid = int(stdout.split()[2])
        return jobid

            

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
