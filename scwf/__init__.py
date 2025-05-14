import os
import scwf

TASK_RECORD = f"{os.environ["HOME"]}/__job_ID__.log"
ROOT_DIR = os.path.dirname(__file__)


def get_env(name, load = False, other_arg = '-B /public/home:/public/home'):
    if load:
        return f"singularity exec {other_arg} --nv {ROOT_DIR}/libs/{name}.sif"
    else:
        return f"{ROOT_DIR}/libs/{name}.sif"