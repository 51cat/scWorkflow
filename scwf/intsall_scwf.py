import click
import scwf
import os
import subprocess

ROOT_DIR = os.path.dirname(scwf.__file__)
Bucket_name = 'yikebucket'

@click.group()
def main():
    pass

@main.command()
def module():
    cmd = f'cp -a ./data/ana_module/ {ROOT_DIR}/'
    subprocess.check_call(cmd, shell=True)


@main.command()
@click.option('--method', type=click.Choice(['oss', 'bypy']), default='oss', help='')
def download_env(method):
    if method == 'oss':
        #cmd1 = f'obsutil cp obs://{Bucket_name}/lib.tar.gz ./'
        cmd2 = f'tar -zxvf ./libs.tar.gz'
        cmd3 = f'cp -a ./libs {ROOT_DIR}/'
        cmd4 = f'rm -rf ./lib/ ./libs.tar.gz'

        #subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)
        subprocess.check_call(cmd4, shell=True)

    elif method == 'bypy':
        click.echo("bypy 下载方式尚未实现")


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
@click.option('--method', type=click.Choice(['oss', 'bypy']), default='oss', help='')
def upload_env(method):
    if not os.path.exists('./libs.tar.gz'):
        click.echo("libs.tar.gz 文件不存在，请先运行 export_env 命令")
        return

    if method == 'obss':
        cmd = f'obsutil cp ./libs.tar.gz obs://{Bucket_name}/'
        subprocess.check_call(cmd, shell=True)
        click.echo("环境上传完成：libs.tar.gz 上传到 OBS")

    elif method == 'bypy':
        click.echo("bypy 上传方式尚未实现")