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

__all__ = ['queue', 'job']

import sys
from cogue.qsystem.queue import QueueBase, LocalQueueBase, RemoteQueueBase
from cogue.qsystem.job import JobBase

def queue(max_jobs=None,
          ssh_shell=None,
          temporary_dir=None,
          name=None,
          sleep_time=None):
    if ssh_shell is None:
        return LocalQueue(max_jobs=max_jobs)
    elif temporary_dir is not None:
        return RemoteQueue(ssh_shell,
                           temporary_dir,
                           max_jobs=max_jobs,
                           name=name,
                           sleep_time=sleep_time)

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

class Qstat:
    def qstat(self):
        """qstatout

        Text of output of 'qstat'

        """
        qstat_out = self._shell.run(["qstat"]).output.split(b'\n')
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
        
def _parse_jobid(qsub_out):
    return int(qsub_out.split()[2])
                    
class LocalQueue(LocalQueueBase,Qstat):
    def __init__(self,
                 max_jobs=None,
                 qsub_command="qsub"):
        LocalQueueBase.__init__(self,
                                max_jobs=max_jobs,
                                qsub_command=qsub_command)

    def _get_jobid(self, qsub_out):
        return _parse_jobid(qsub_out)

class RemoteQueue(RemoteQueueBase,Qstat):
    def __init__(self,
                 ssh_shell,
                 temporary_dir,
                 max_jobs=None,
                 name=None,
                 sleep_time=None,
                 qsub_command="qsub"):
        RemoteQueueBase.__init__(self,
                                 ssh_shell,
                                 temporary_dir,
                                 max_jobs=max_jobs,
                                 name=name,
                                 sleep_time=sleep_time,
                                 qsub_command=qsub_command)

    def _get_jobid(self, qsub_out):
        return _parse_jobid(qsub_out)

class Job(JobBase):
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

        JobBase.__init__(self)

        if not script:
            print("Queue script not found.")
            sys.exit(1)
        else:
            self._script = script

        if shell is None:
            self._shell = "/bin/bash"
        else:
            self._shell = shell

        if jobname is None:
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

    def copy(self, jobname=None):
        if jobname is None:
            jobname_new = self._jobname 
        else:
            jobname_new = jobname
        return Job(script=self._script,
                   shell=self._shell,
                   cwd=self._cwd,
                   jobname=jobname_new,
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
        w.write("#$ -S %s\n" % self._shell)
        
        if self._cwd:
            w.write("#$ -cwd\n")
            
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
