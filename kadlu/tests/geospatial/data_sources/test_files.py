""" 
dynamically generate tests for a given set of input data files.
source code tests will be generated and appended to a script using files in the 
kadlu_data/testfiles directory as inputs, which will then be executed in a 
subprocess and delete itself upon execution. 

pytest config can be passed using the DEBUGOPTS env variable:
by default the script will drop into a live pdb debugging session.
alternatively, run in parallel processing mode using pytest-parallel

dependency: 
    pytest
    pytest-parallel (optional)

usage:
    # attach pdb debugger for further inspection
    export DEBUGOPTS='--pdb --tb=native -s'
    python3 test_files.py

    # run tests in parallel on all cores
    pip install pytest-parallel
    export DEBUGOPTS='--workers=auto --tb=line' 
    python3 -B test_files.py | tee filetests.log

    # run tests from an interactive session or python source code:
    os.environ.setdefault('DEBUGOPTS', '--pdb --tb=native -s')
    kadlu.test_files()

    see the DEBUGOPTS env var and 'man pytest' for further usage information
"""
import os, kadlu, subprocess

class test_files():
    def __init__(self):
        IMPORTS = 'import os, kadlu'
        CORPUS, _, CFILES = list(os.walk(kadlu.storage_cfg()+'testfiles'))[0]
        EXPORTS = lambda C,PATH=CORPUS: f'def test_loadfile_{C.replace(".","").replace("-","").replace(" ","")}():\n\tres = kadlu.load_file(os.path.join("{PATH}","{C}")); print(f"\\n{C}:\\n{{res}}\\n"); assert res, "error: {C}"\n\n'
        with open('scriptoutput.py', 'w') as OUTPUT: OUTPUT.write(IMPORTS+'\n\n'+''.join(map(EXPORTS, CFILES))+'os.remove("scriptoutput.py")')
        subprocess.run(f'python3 -B -m pytest scriptoutput.py {os.environ.get("DEBUGOPTS","--pdb --tb=native -s")}'.split())

