import os

from ebi_eva_common_pyutils.config import cfg


def load_config(*args):
    """
    Load a config file from any path provided.
    If none are provided then read from a file path provided in the environment variable ASSEMBLYCONFIG.
    If not provided then default to .assembly_config.yml place in the current users' home
    """
    cfg.load_config_file(
        *args,
        os.getenv('ASSEMBLYCONFIG'),
        os.path.expanduser('~/.assembly_config.yml'),
    )


def get_nextflow_config_flag():
    """
    Return the commandline flag for Nextflow to use the config provided in environment variable ASSEMBLY_NEXTFLOW_CONFIG.
    If not provided, return an empty string, which allows Nextflow to use the default precedence as described here:
    https://www.nextflow.io/docs/latest/config.html
    """
    env_val = os.getenv('ASSEMBLY_NEXTFLOW_CONFIG')
    if env_val:
        return f'-c {env_val}'
    return ''
