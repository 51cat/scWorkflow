import click
import scwf
import os
import subprocess

ROOT_DIR = os.path.dirname(scwf.__file__)

@click.group()
def main():
    pass

@main.command()
def install_module():
    cmd = f'cp -a ./data/ana_module/ {ROOT_DIR}/'
    subprocess.check_call(cmd, shell=True)

@main.command()
def mk_config():
    cmd = f'cp -a ./data/configs/ {ROOT_DIR}/'
    subprocess.check_call(cmd, shell=True)


@main.command()
@click.option('--env_p', default='./libs/', help='')
def install_env(env_p):
    cmd = f'cp -a {env_p} {ROOT_DIR}/'
    subprocess.check_call(cmd, shell=True)


@main.command()
@click.option('--method', type=click.Choice(['oss', 'bypy']), default='bypy', help='')
def download_env(method):
    if method == 'oss':
        print('remove!')

    elif method == 'bypy':
        config_dir = f"{ROOT_DIR}/configs/.bypy"
        cmd1 = f"bypy --config-dir {config_dir} downfile libs.tar.gz"
        subprocess.check_call(cmd1, shell=True)

        cmd2 = f'tar -zxvf libs.tar.gz'
        subprocess.check_call(cmd2, shell=True)


@main.command()
def export_env():
    lib_dir = os.path.join(ROOT_DIR, 'libs')
    if not os.path.exists(lib_dir):
        click.echo(f"libs 目录不存在: {lib_dir}")
        return

    tar_cmd = f'tar -zcvf ./libs.tar.gz -C {ROOT_DIR} libs'
    subprocess.check_call(tar_cmd, shell=True)
    click.echo("环境打包完成：libs.tar.gz")


@main.command()
def export_module():
    module_dir = os.path.join(ROOT_DIR, 'ana_module')
    if not os.path.exists(module_dir):
        click.echo(f"module_dir 目录不存在: {module_dir}")
        return

    tar_cmd = f'tar -zcvf ./modules.tar.gz -C {ROOT_DIR} ana_module'
    subprocess.check_call(tar_cmd, shell=True)
    click.echo("分析模块打包完成：modules.tar.gz")


@main.command()
@click.option('--method', type=click.Choice(['oss', 'bypy']), default='bypy', help='')
def upload_env(method):
    if method == 'oss':
        print('remove!')

    elif method == 'bypy':
        config_dir = f"{ROOT_DIR}/configs/.bypy"
        cmd1 = f"bypy --config-dir {config_dir} upload libs.tar.gz"
        subprocess.check_call(cmd1, shell=True)


@main.command()
@click.option('--method', type=click.Choice(['oss', 'bypy']), default='bypy', help='')
def upload_module(method):
    if method == 'oss':
        print('remove!')

    elif method == 'bypy':
        config_dir = f"{ROOT_DIR}/configs/.bypy"
        cmd1 = f"bypy --config-dir {config_dir} upload modules.tar.gz"
        subprocess.check_call(cmd1, shell=True)