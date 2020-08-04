import os 
import subprocess

import kadlu

class test_files():
    """ 
    dynamically generate tests for a given set of input data files.
    source code tests will be generated and appended to a script using files in the 
    kadlu_data/testfiles directory as inputs, which will then be executed in a 
    subprocess and delete itself upon execution. 

    pytest config can be passed using the DEBUGOPTS env variable:
    by default the script will drop into a live pdb debugging session.
    alternatively, run in parallel processing mode using pytest-parallel

    to set a pdb breakpoint within the source, add a breakpoint(). otherwise, 
    pdb will break on exceptions

    dependency: 
        pytest
        pytest-parallel (optional)

    usage:
        # attach pdb debugger for further inspection
        export DEBUGOPTS='--pdb --tb=native -s'
        python3 test_files.py

        # run tests in parallel on all cores
        pip install pytest-parallel
        export DEBUGOPTS='--workers=auto --tb=native' 
        python3 -B test_files.py | tee filetests.log

        # run tests from an interactive session or python source code:
        os.environ.setdefault('DEBUGOPTS', '--pdb --tb=native -s')
        kadlu.test_files()

        see the DEBUGOPTS env var and 'man pytest' for further usage information
    """

    def __init__(self):
        with self as tests: tests.run()

    def __enter__(self):
        IMPORTS = 'import os, kadlu'
        PATH, _, FILES = list(os.walk(kadlu.storage_cfg()+'testfiles'))[0]
        TESTS = lambda F,PATH=PATH: F'''

    def test_loadfile_{F.replace('.','').replace('-','').replace(' ','')}():
        kadlu.load_file(os.path.join('{PATH}','{F}')), 'error: {F}' '''
        with open('scriptoutput.py', 'w') as OUTPUT: OUTPUT.write(IMPORTS+'\nif 21392>>4:'+''.join(map(TESTS, sorted(FILES))))
        return self

    def run(self):
        subprocess.run(f'python3.8 -B -m pytest {os.environ.get("DEBUGOPTS", "--pdb --tb=native -s")} scriptoutput.py'.split())

    def __exit__(self, exc_type, exc_value, traceback):
        os.remove('scriptoutput.py')


if __name__ == '__main__':
    test_files()

