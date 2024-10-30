import argparse
from pathlib import Path
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description='set a seed and run the pipeline')
    parser.add_argument('--seed', type=int, help='seed')
    parser.add_argument(
        '--config_template', type=str, help='Provide path to template config file'
    )
    return parser.parse_args()


def main(args):
    with open(args.config_template, 'r') as file:
        config_content = file.read()

    config_content = config_content.replace('__seed__', str(args.seed))

    with open(
        Path(args.config_template).parent / f'config_{args.seed}.yaml', 'w'
    ) as file:
        file.write(config_content)

    command = (
        'source $(conda info --base)/etc/profile.d/conda.sh && conda activate snakemake; '
        f'sbatch --time=24:00:00 --wrap="snakemake --cores 10 --configfile config/config_{args.seed}.yaml '
        '--software-deployment-method conda --rerun-incomplete  -p  --keep-going --profile profiles/slurm/"'
    )

    subprocess.run(command, shell=True, check=True, executable='/bin/bash')


if __name__ == '__main__':
    args = parse_args()
    main(args)
