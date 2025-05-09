import os
from datetime import datetime
import subprocess
from rich import print as rprint
from .utils.add_log import add_log
import time
import re
from scwf import TASK_RECORD


def generate_slurm_script(work_name, mem, p, command, work_dir):
    slurm_template = f"""#!/bin/bash
#SBATCH -J {work_name}
#SBATCH -p debug
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={p}
#SBATCH --time 87600:00:00
#SBATCH --mem={mem}
#SBATCH -o {work_dir}/{work_name}_std.%j.out
#SBATCH -e {work_dir}/{work_name}_std.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zhliuccc@163.com

{command}
"""
    return slurm_template


def parse_submit_info_get_jobIDs(text):
    job_id = re.match(r'Submitted batch job (\d+)', text).group(1) # get job id')
    return job_id

class CommandExecutor:
    def __init__(self, 
                 command, 
                 module_name,
                 verbose = True, 
                 cpu = None, memory_gb = None):
        self.command = command
        self.verbose = verbose
        self.module_name = module_name
        self.cpu = cpu
        self.memory_gb = memory_gb
        self.now = datetime.now()
        self.sufx = f"{self.now.year:04d}{self.now.month:02d}{self.now.day:02d}_{self.now.hour:02d}{self.now.minute:02d}{self.now.second:02d}"
        self.sufx2 = f"{self.now.year:04d}-{self.now.month:02d}-{self.now.day:02d}-{self.now.hour:02d}:{self.now.minute:02d}:{self.now.second:02d}"
    
    @add_log
    def mk_workdir(self):
        os.system(f'mkdir -p ./_workdir_/{self.module_name}_{self.sufx}')
        self.wd = f"./_workdir_/{self.module_name}_{self.sufx}"
        self.logfile = f"./_workdir_/{self.module_name}_{self.sufx}/_run.log"
        self.mk_workdir.logger.info(f'workdir: {self.wd}')
        self.mk_workdir.logger.info(f'logfile: {self.logfile}')

    @add_log
    def exec_bash(self):
        if self.verbose in [False, 'False']:
            cmd = f"{self.command} &> {self.logfile}"
        if self.verbose in [True, 'True']:
            cmd = f"{self.command} &> {self.logfile}"
        
        with open(f"{self.wd}/cmd.sh", 'w') as fd:
            fd.write(cmd)

        self.exec_bash.logger.info(f'cmd: {cmd}')
        self.exec_bash.logger.info(f"{self.wd}/cmd.sh")

        try:
            subprocess.check_call(cmd, shell=True)
        except subprocess.CalledProcessError:
            self.exec_bash.logger.error(f"exec ERROR!")
            self.exec_bash.logger.error(f"WD: {self.wd}")
            with open(self.logfile) as fd:
                rprint(fd.read())

    @add_log
    def exec_submit(self, p, m, work_name = None):
        
        if work_name in [None, 'None']:
            work_name = f"{self.module_name}_{self.sufx}"

        script = generate_slurm_script(
            work_name=work_name,
            mem=m,
            p=p,
            command=self.command,
            work_dir=os.path.abspath(self.wd)
        )
        with open(f"{os.path.abspath(self.wd)}/submit.slurm", 'w') as fd:
            fd.write(script)

        # 获取提交信息
        sbatch_command = f'sbatch {os.path.abspath(self.wd)}/submit.slurm'
        result = subprocess.run(sbatch_command, shell=True, capture_output=True, text=True)

        jobID = parse_submit_info_get_jobIDs(result.stdout)
        with open(f"{TASK_RECORD}", 'a') as fd:
            fd.write(f"{self.sufx2}\t{jobID}\t{work_name}\t{os.path.abspath(self.wd)}\n")
            

    def exec(self, method = 'bash', p = 4, mem = '60G', work_name = None):
        self.mk_workdir()
        if method == 'bash':
            self.exec_bash()
        if method == 'slurm':
            self.exec_submit(p, mem, work_name)
            