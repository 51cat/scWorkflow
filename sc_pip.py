#!/public/home/zhliu/PIPLINE/singlecell/conda_env/sc_start/bin/python

import glob
import os
import click
import subprocess
from utils.add_log import add_log
import json
from rich.console import Console
from rich.markdown import Markdown
from rich.console import Console
from rich.table import Table
from rich import print as rprint


ROOT_DIR = os.path.dirname(__file__)

def find_analysis_method(name):
    sc_dir = ROOT_DIR
    if not os.path.exists(f"{sc_dir}/ana_module/{name}"):
        raise FileNotFoundError()
    else:
        return f"{sc_dir}/ana_module/{name}"

def load_env(env_p, soft = 'singularity', sub_cmd = 'exec --nv'):
    if not os.path.exists(env_p):
        raise FileExistsError(f"Not Found {env_p}")
    return ' '.join([soft, sub_cmd, os.path.abspath(env_p)])


def json_to_dict(js):
    with open(js, 'r', encoding='utf-8') as file:
            parsed_dict = json.load(file)
    return parsed_dict


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
    
    @add_log
    def run_cmd(self):
        self.cmd += " &> ./__run.log__ "
        self.run_cmd.logger.info(f'RUN CMD: {self.cmd}')
        try:
            subprocess.check_call(self.cmd, shell=True)
        except subprocess.CalledProcessError:
            with open("./__run.log__") as fd:
                self.run_cmd.logger.error(f"RUN: Fail!")
                rprint(fd.read())
                exit()

            


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
@click.argument('module')
def help(module):
    r = ModuleRunner(module)
    r.get_module_file()
    r.show_readme()



@add_log
@cli.command(context_settings=dict(ignore_unknown_options=True))
@click.option('--module', help='')
@click.option('--sc_type', help='')
def mk_pipline(module, sc_type):
    os.system(f"mkdir -p ./{module}  && touch ./{module}/README.md && mkdir ./{module}/demo/")
    if sc_type == 'R':
        os.system(f"touch ./{module}/start.r")
    if sc_type == 'python':
        os.system(f"touch ./{module}/start.py")

    with open(f'./{module}/info.json', 'w') as file:
            json.dump({"exec":sc_type}, file, indent=4)



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
@click.argument('args', nargs=-1, type=click.UNPROCESSED)
def run(module, env, exec,sc, args, no_env):
    run.logger.info(f"RUN Module Name: {module}")
    run.logger.info(f"RUN env: {env}")
    run.logger.info(f"RUN exec: {exec}")
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
    
    r.run_cmd()


if __name__ == '__main__':
    cli()
