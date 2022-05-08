import subprocess
import sys
import setuptools
from setuptools.command.install import install

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

class Install(install):
    """Customized setuptools install command - starts Makefile ensuring Jellyfish presence."""
    def run(self):
        jf_command_unix = ["make", "jellyfish"]
        jf_command_mac = ["make", "jellyfish-mac"]
        if subprocess.call(jf_command_unix) != 0:
            sys.exit(-1)
        install.run(self)


setuptools.setup(
    name="sbat-alex_skyslakova",
    version="0.0.3",
    author="Alexandra Skyslakova",
    author_email="alexandra.skyslakova@gmail.com",
    description="A tool for strand bias analysis of NGS data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/alex-skyslakova/strand-bias-analysis-tool",
    #project_urls={
    #    "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    #},
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Operating System :: MacOS",
        "Operating System :: Unix",
    ],
    package_dir={"": "sbat"},
    packages=setuptools.find_packages(where="sbat"),
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            'sbat=sbat.main:main'
        ]
    },
)
