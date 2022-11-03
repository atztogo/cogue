"""Queue classes."""
import datetime
import os
import shlex
import shutil
import sys
import tarfile
import time
import traceback
from subprocess import check_output


def get_time():
    """Return current time."""
    return datetime.datetime.today().strftime("%H:%M:%S")


class EmptyQueue:
    """EmptyQueue class."""

    def __init__(self):
        """Init method."""
        pass

    def register(self, task):
        """Register."""
        pass

    def submit(self, task):
        """Submit."""
        pass

    def qstat(self):
        """Qstat."""
        pass

    def set_max_jobs(self, max_jobs):
        """Set max jobs."""
        pass

    def write_qstatus(self, name):
        """Write qstatus."""
        pass


class QueueBase:
    """Queue base class."""

    def __init__(self, max_jobs=None):
        """Init method."""
        self._max_jobs = max_jobs
        self._qstatus = None
        self._tid_queue = []
        self._tid2jobid = {}
        self._shell = None
        self._shell_type = None

    def register(self, task):
        """Register."""
        if task.get_traverse() is False:
            self._tid_queue.append(task.get_tid())
            job = task.get_job()
            job.set_status("preparing")

    def write_qstatus(self, name):
        """Write qstatus."""
        with open("%s.qstat" % name, "w") as f_qstat:
            f_qstat.write("%8s %8s %8s\n" % ("tid", "jobid", "status"))
            for tid in self._tid_queue:
                f_qstat.write("%8d %8s %8s\n" % (tid, "None", "Queued"))

            for tid in self._tid2jobid:
                jobid = self._tid2jobid[tid]
                if jobid in self._qstatus:
                    f_qstat.write("%8d %8d %8s\n" % (tid, jobid, self._qstatus[jobid]))

    def submit(self):
        """Submit."""
        pass

    def set_max_jobs(self, max_jobs):
        """Set max jobs."""
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
                if self._qstatus[jobid] == "Running":
                    job.set_status("running", jobid)
            else:
                del self._tid2jobid[tid]
                job.set_status("done")


class LocalQueueBase(QueueBase):
    """LocalQueue base class."""

    def __init__(self, max_jobs=None, qsub_command="qsub"):
        """Init method."""
        try:
            import spur
        except ImportError:
            print("You need to install spur.")
            exit(1)

        QueueBase.__init__(self, max_jobs=max_jobs)
        self._qsub_command = qsub_command
        self._shell = spur.LocalShell()

    def submit(self, task):
        """Submit."""
        if task.get_traverse() is not False:
            return

        job = task.get_job()
        tid = task.get_tid()
        self._set_job_status(job, tid)
        if "ready" in job.get_status():
            job.write_script()
            if task.get_traverse():
                jobid = None
            else:
                qsub_out = self._shell.run(
                    shlex.split(self._qsub_command + " " + "job.sh"), cwd=os.getcwd()
                ).output
                jobid = self._get_jobid(qsub_out)
            self._tid2jobid[tid] = jobid
            self._tid_queue.pop(0)
            job.set_status("submitted", jobid)


class RemoteQueueBase(QueueBase):
    """RemoreQueue base class."""

    def __init__(
        self,
        ssh_shell,
        temporary_dir,
        max_jobs=None,
        name=None,
        sleep_time=None,
        qsub_command="qsub",
    ):
        """Init method."""
        QueueBase.__init__(self, max_jobs=max_jobs)
        self._qsub_command = qsub_command
        self._shell = ssh_shell
        self._name = name
        if sleep_time is None:
            self._sleep_time = 0.1
        else:
            self._sleep_time = sleep_time

        if not os.path.exists(temporary_dir):
            self._shell_run(["mkdir", "-p", temporary_dir])

        if self._name is not None:
            self._working_dir = "%s/%s" % (temporary_dir, self._name)
            if not os.path.exists(self._working_dir):
                self._shell_run(["mkdir", "-p", self._working_dir])
        else:
            self._working_dir = "%s" % temporary_dir

    def submit(self, task):
        """Submit."""
        if task.get_traverse() is not False:
            return

        job = task.get_job()
        tid = task.get_tid()
        self._set_job_status(job, tid)

        if "ready" in job.get_status():
            self._submit(task)
        elif "done" in job.get_status():
            self._collect(task)

    def _submit(self, task):
        job = task.get_job()
        tid = task.get_tid()
        remote_dir = "%s/c%05d" % (self._working_dir, tid)
        task_log = task.get_log()

        job.write_script()
        self._shell_run(["mkdir", "-p", remote_dir])

        with tarfile.open("cogue.tar", "w") as tar:
            for name in os.listdir("."):
                tar.add(name)

        for i in range(20):
            with open("cogue.tar", "rb") as local_file:
                with self._shell.open(
                    "%s/%s" % (remote_dir, "cogue.tar"), "wb"
                ) as remote_file:
                    shutil.copyfileobj(local_file, remote_file)
                    time.sleep(self._sleep_time)

            shasum_l = check_output(["shasum", "cogue.tar"]).split()[0]
            shasum_r = self._shell_run(
                ["shasum", "cogue.tar"], cwd=remote_dir
            ).output.split()[0]
            task_log += "    copy local %s -> remote %s (%s tid-%05d)\n" % (
                shasum_l[:8],
                shasum_r[:8],
                get_time(),
                tid,
            )

            if shasum_r == shasum_l:
                break
            else:
                err_log = (
                    "    copying to remote, waiting for 10s..."
                    " (%s tid-%05d)\n" % (get_time(), tid)
                )
                sys.stderr.write(err_log)
                task_log += err_log
                time.sleep(10)

        os.remove("cogue.tar")
        self._shell_run(["tar", "xvf", "cogue.tar"], cwd=remote_dir)
        self._shell_run(["rm", "cogue.tar"], cwd=remote_dir)

        if task.get_traverse():
            jobid = None
        else:
            qsub_out = self._shell_run(
                shlex.split(self._qsub_command + " " + "job.sh"), cwd=remote_dir
            ).output
            jobid = self._get_jobid(qsub_out)
        self._tid2jobid[tid] = jobid
        self._tid_queue.pop(0)
        job.set_status("submitted", jobid)

        task.set_log(task_log)

    def _collect(self, task):
        tid = task.get_tid()
        remote_dir = "%s/c%05d" % (self._working_dir, tid)
        task_log = task.get_log()

        names = [
            filename.decode("utf-8")
            for filename in self._shell_run(["/bin/ls"], cwd=remote_dir).output.split()
        ]
        self._shell_run(["tar", "cvf", "cogue.tar"] + names, cwd=remote_dir)

        for i in range(20):
            with open("cogue.tar", "wb") as local_file:
                with self._shell.open(
                    "%s/%s" % (remote_dir, "cogue.tar"), "rb"
                ) as remote_file:
                    shutil.copyfileobj(remote_file, local_file)
                    time.sleep(self._sleep_time)

            shasum_r = self._shell_run(
                ["shasum", "cogue.tar"], cwd=remote_dir
            ).output.split()[0]

            shasum_l = check_output(["shasum", "cogue.tar"]).split()[0]
            task_log += "    copy local %s <- remote %s (%s tid-%05d)\n" % (
                shasum_l[:8],
                shasum_r[:8],
                get_time(),
                tid,
            )

            if shasum_r == shasum_l:
                break
            else:
                err_log = (
                    "    copying from remote, waiting for 10s..."
                    " (%s tid-%05d)\n" % (get_time(), tid)
                )
                sys.stderr.write(err_log)
                task_log += err_log
                time.sleep(10)

        with tarfile.open("cogue.tar") as tar:

            def is_within_directory(directory, target):

                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)

                prefix = os.path.commonprefix([abs_directory, abs_target])

                return prefix == abs_directory

            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):

                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")

                tar.extractall(path, members, numeric_owner=numeric_owner)

            safe_extract(tar)

        os.remove("cogue.tar")
        self._shell_run(["rm", "cogue.tar"], cwd=remote_dir)

        task.set_log(task_log)

    def _shell_run(self, command, cwd=None):
        import spur

        shell_done = False
        for i in range(10):
            try:
                if cwd is None:
                    return self._shell.run(command)
                else:
                    return self._shell.run(command, cwd=cwd)
                time.sleep(self._sleep_time)
                shell_done = True
            except spur.results.RunProcessError:
                date = datetime.datetime.today().strftime("%Y-%m-%d %H:%M:%S")
                print("%s: SshShell.run() failed (%d)." % (date, i + 1))
                print("command: %s" % command)
                if cwd:
                    print("cwd: %s" % cwd)
                if i < 9:
                    err_log = traceback.format_exc()
                    sys.stderr.write(err_log)
                    print("    copying from remote, waiting for 10s...")

                    time.sleep(10)
                else:
                    raise RuntimeError

            if shell_done:
                break
