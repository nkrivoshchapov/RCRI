import subprocess,os
from setuptools import setup

def runcmd(cmd):
    subprocess.call(cmd.split(" "))

runcmd("mkdir build")
os.chdir("./build")
runcmd("cmake ..")
runcmd("make")
runcmd("make install")
os.chdir("..")
setup(
    name='rcri',
    version='0.0.1',
    author='Nikolai Krivoshchapov',
    packages=['rcrilib'],
    package_data={'rcrilib': ['CoreParsers/*py', 'Helpers/*py', 'Solvers/*py', 'Solvers/libtlc.so']},
    install_requires=[
            'networkx',
            'numpy',
            'scipy'
        ]
)
