import glob
import os
import click
import subprocess
from .utils.add_log import add_log
from .exec import CommandExecutor
import json
from rich.console import Console
from rich.markdown import Markdown
from rich.table import Table
from rich import print as rprint
from rich.text import Text
from collections import deque
import re
from scwf import TASK_RECORD
import scwf
from scwf import ROOT_DIR

#ROOT_DIR = os.path.dirname(scwf.__file__)

def find_analysis_method(name):
    sc_dir = ROOT_DIR
    if not os.path.exists(f"{sc_dir}/ana_module/{name}"):
        raise FileNotFoundError()
    else:
        return f"{sc_dir}/ana_module/{name}"

def load_env(env_p, soft = 'singularity', sub_cmd = 'exec -B /public/home:/public/home --nv'):
    if not os.path.exists(env_p):
        raise FileExistsError(f"Not Found {env_p}")
    return ' '.join([soft, sub_cmd, os.path.abspath(env_p)])


def json_to_dict(js):
    with open(js, 'r', encoding='utf-8') as file:
            parsed_dict = json.load(file)
    return parsed_dict

def parse_job_STATE(text):
    text = text.strip()
    out_lines = text.split("\n")

    # 移除 ANSI 转义码（如果有）
    ansi_escape = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
    out_lines = [ansi_escape.sub('', line).strip() for line in out_lines if line.strip()]

    if len(out_lines) > 2:
        state = out_lines[2]
        return state
    else:
        return f"None" 
    
def get_file_size_gb(file_path):
    size_bytes = os.path.getsize(file_path)
    size_gb = size_bytes / (1024 ** 3)
    return round(size_gb, 2)

class ModuleRunner:
    def __init__(
            self,
            module_name,
            args_dict = None,
            env = 'auto',
            exec = 'auto'
            ):
        self.module_name = module_name
        self.env = env
        self.args_dict = args_dict
        self.exec = exec

    def set_env(self, p):
        self.env = p
    
    def set_sc(self, p):
        self.start_sc = p
    
    def set_exec(self, e):
        self.exec = e

    def get_module_file(self):
        self.module_dir = find_analysis_method(self.module_name)
        
        self.exec = json_to_dict(f"{self.module_dir}/info.json")['exec']

        if (self.exec == 'python'):
            self.start_sc = f"{self.module_dir}/start.py"
            self.exec = 'python'

        elif (self.exec == 'R'):
            self.start_sc = f"{self.module_dir}/start.r"
            self.exec = 'Rscript'
        
        self.env = f"{ROOT_DIR}/libs/{self.module_name}.sif"
        if not os.path.exists(self.env):
            self.env = f"{ROOT_DIR}/libs/_common.sif"
        self.readme = f"{self.module_dir}/README.md"

    def show_readme(self):
        console = Console()
        with open(self.readme) as readme:
            markdown = Markdown(readme.read())
        console.print(markdown)

    def get_runner_env(self):
        self.runner_env = load_env(self.env) + ' ' + self.exec
    
    def get_run_args(self):
        args_lst = []
        for k, v in self.args_dict.items():
            args_lst.append(f"--{k} {v}")
        self.args_str = " ".join(args_lst)
        

    def mk_run_cmd(self, no_env= False):
        if no_env:
            self.cmd = ' '.join([self.exec, self.start_sc, self.args_str])
        else:
            self.cmd = ' '.join([self.runner_env, self.start_sc, self.args_str])
    
    def get_cmd(self):
        return self.cmd

            


@click.group()
def cli():
    pass


@cli.command(context_settings=dict(ignore_unknown_options=True))
def ls():
    console = Console()
    table = Table(show_header=True, header_style="bold magenta")
    table.title = "单细胞转录组分析流程"
    table.add_column("序号", style="dim", width=5)
    table.add_column("流程名称", style="dim", width=20)
    table.add_column("路径")

    inx = 0
    for fname in os.listdir(f"{ROOT_DIR}/ana_module/"):
        inx += 1
        table.add_row(str(inx), fname, f"{ROOT_DIR}/ana_module/{fname}/")
    
    console.print(table)


@cli.command(context_settings=dict(ignore_unknown_options=True))
def lsenv():
    console = Console()
    table = Table(show_header=True, header_style="bold magenta")
    table.title = "Analysis envirnment"
    table.add_column("序号", style="dim", width=5)
    table.add_column("env", style="dim", width=20)
    table.add_column("size(GB)")

    inx = 0
    total_size = 0
    env_all = glob.glob(f"{ROOT_DIR}/libs/*.sif")
    for file in env_all:
        inx += 1
        size = get_file_size_gb(file)
        table.add_row(str(inx), file.split("/")[-1], str(size))
        total_size += size
    
    console.print(table)
    rprint(f"Total size: {total_size}G")
    rprint(f"Env dir: {ROOT_DIR}/libs/")



@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.argument('module')
def help(module):
    r = ModuleRunner(module)
    r.get_module_file()
    r.show_readme()


@cli.command()
@click.option('--ntask','-n', default=30, required=False)
def task_stat(ntask):

    if not os.path.exists(TASK_RECORD):
        rprint('no task stat!')
    else:
        with open(TASK_RECORD) as fd:
            console = Console()
            table = Table(show_header=True, header_style="#39C5BB")
            table.add_column("date")
            table.add_column("task_name")
            table.add_column("jobid")
            table.add_column("stat")
            table.add_column("work_dir")

            last_lines = deque(fd, maxlen=ntask)

            for line in last_lines:
                line = line.strip('\n')
                time, jobid, jobname, wd = line.split('\t')
                command = f"sacct -o State -j {jobid} "
                result = subprocess.run(command, shell=True, capture_output=True, text=True)
                job_status = parse_job_STATE(result.stdout)

                if job_status == 'COMPLETED':
                    job_status = Text(job_status, style="#FF69B4")
                elif job_status == 'RUNNING':
                    job_status = Text(job_status, style="#1E90FF")
                elif job_status == 'FAILED':
                    job_status = Text(job_status, style="#FF4500")
                elif job_status == 'CANCELLED+':
                    job_status = Text(job_status, style="#FFD700")
                else:
                    job_status = Text(job_status, style="white")

                table.add_row(time, jobname, jobid, job_status,  wd)
        console.print(table)


@add_log
@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option('--module', help='')
@click.option('--sc_type', help='')
@click.option('--add_colordb', help='', default = 'False')
def mk_pipline(module, sc_type, add_colordb):
    os.system(f"mkdir -p ./{module}  && touch ./{module}/README.md && mkdir ./{module}/demo/")
    if sc_type == 'R':
        os.system(f"touch ./{module}/start.r")
    if sc_type == 'python':
        os.system(f"touch ./{module}/start.py")

    with open(f'./{module}/info.json', 'w') as file:
            json.dump({"exec":sc_type}, file, indent=4)

    if add_colordb in [True, 'True']:
        pass



@add_log
@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option('--module_path', help='')
@click.option('--env_path', help='', default = 'None')
def install_pipline(module_path, env_path):
    module_name = module_path.strip('/').split('/')[-1]
    env_name = env_path.split('/')[-1]

    module_list = [fname for fname in os.listdir(f"{ROOT_DIR}/ana_module/")]
    env_list = [fname for fname in os.listdir(f"{ROOT_DIR}/libs/")]
    if module_name in module_list:
        install_pipline.logger.error(f'duplicate  module_name: {module_name}')
        raise KeyError('')
    
    if env_name in env_list:
        install_pipline.logger.error(f'duplicate  env_name: {env_name}')
        raise KeyError('')
    
    if env_path in [None, 'None']:
        install_pipline.logger.warning('Not Found env path!')
    
    cmd = f"mv {module_path} {ROOT_DIR}/ana_module/ "
    subprocess.check_call(cmd, shell=True)

    if env_path not in [None, 'None']:
        cmd2 = f"cp {env_path} {ROOT_DIR}/libs/ "
        subprocess.check_call(cmd2, shell=True)
    
    install_pipline.logger.info(f'install finish!')
    

    



@add_log
@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option('--module', help='')
@click.option('--env', default='auto')
@click.option('--exec', default='auto')
@click.option('--sc', default='auto')
@click.option('--no_env', default='False')
@click.option('--verbose', default='True')
@click.option('--exec_method', default='bash')
@click.option('--cpu', default='4')
@click.option('--mem', default='60G')
@click.option('--wn', default=None)
@click.argument('args', nargs=-1, type=click.UNPROCESSED)
def run(module, env, exec,sc, args, no_env, verbose, exec_method, cpu, mem, wn):
    run.logger.info(f"RUN Module Name: {module}")
    run.logger.info(f"RUN env: {env}")
    run.logger.info(f"RUN exec: {exec}")
    run.logger.info(f"RUN exec_method: {exec_method}")
    run.logger.info(f"RUN script: {sc}")
    
    
    arg_dict = {}
    i = 0
    while i < len(args):
        if args[i].startswith('--'):
            key = args[i].lstrip('-')
            if '=' in key:
                k, v = key.split('=', 1)
                arg_dict[k] = v
            elif i + 1 < len(args) and not args[i + 1].startswith('--'):
                arg_dict[key] = args[i + 1]
                i += 1  # skip value
            else:
                arg_dict[key] = True  # flag-style, e.g., --debug
        i += 1

    run.logger.info(f"RUN args: {arg_dict}")

    r = ModuleRunner(module, args_dict = arg_dict)
    r.get_module_file()

    if env != 'auto':
        r.set_env(env)
    
    if exec != 'auto':
        r.set_exec(exec)
    
    if sc != 'auto':
        r.set_sc(sc)

    r.get_runner_env()
    r.get_run_args()
    
    if no_env in ['True', True]:
        r.mk_run_cmd(no_env=True)
    else:
        r.mk_run_cmd()
    
    ce = CommandExecutor(
        r.get_cmd(), 
        module,
        verbose = verbose)
    ce.exec(method = exec_method, p = cpu, mem = mem, work_name = wn)



if __name__ == '__main__':
    cli()
