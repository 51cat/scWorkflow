from setuptools import setup, find_packages
from glob import glob

setup(
    name="scwf",
    version='0.1',
    description='',
    url='',
    packages=find_packages(),
    keywords=['keyword'],
    install_requires=[],
    python_requires=">=3",
    entry_points="""
        [console_scripts]
        scwf=scwf.sc_pip:cli
        scwf_tk=scwf.scwf_tk:main
    """,
    include_package_data=True,
    zip_safe=False
)
