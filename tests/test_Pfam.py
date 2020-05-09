# import
## batteries
import os
import sys
import pytest

# test/data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

def test_help(script_runner):
    ret = script_runner.run('traitar', 'pfam', '-h')
    assert ret.success

@pytest.mark.dependency()
def test_pfam(script_runner, tests_tmpdir):
    out_dir = os.path.join(tests_tmpdir, 'pfam')
    ret = script_runner.run('traitar', 'pfam', out_dir)
    assert ret.success
    
