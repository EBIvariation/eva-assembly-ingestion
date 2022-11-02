import glob
from distutils.core import setup
from os.path import join, abspath, dirname

base_dir = abspath(dirname(__file__))
requirements_txt = join(base_dir, 'requirements.txt')
requirements = [l.strip() for l in open(requirements_txt) if l and not l.startswith('#')]

version = open(join(base_dir, 'VERSION')).read().strip()

setup(
    name='eva_assembly_ingestion',
    packages=['eva_assembly_ingestion'],
    package_data={'eva_assembly_ingestion': ['VERSION', 'nextflow/*']},
    scripts=glob.glob(join(dirname(__file__), 'bin', '*.py')),
    version=version,
    license='Apache',
    description='Tools for ingesting a new assembly for a species',
    url='https://github.com/EBIvariation/eva-assembly-ingestion',
    keywords=['ebi', 'eva', 'python'],
    install_requires=requirements,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3'
    ]
)
