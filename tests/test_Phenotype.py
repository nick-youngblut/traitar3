# import
## batteries
import os
import sys
import pytest

# test/data dir
test_dir = os.path.join(os.path.dirname(__file__))
data_dir = os.path.join(test_dir, 'data')

def test_help(script_runner):
    ret = script_runner.run('traitar', 'phenotype', '-h')
    assert ret.success
    
@pytest.mark.dependency(depends=['tests/test_Pfam.py::test_pfam[inprocess]'], scope='session')    
def test_phenotype(script_runner, tests_tmpdir):
    pfam_dir = os.path.join(tests_tmpdir, 'pfam')
    sample_file = os.path.join(data_dir, 'samples.txt')
    out_dir = os.path.join(tests_tmpdir, 'phenotype_out')
    ret = script_runner.run('traitar', 'phenotype', '--overwrite',
                            pfam_dir, data_dir, sample_file, 'from_genes', out_dir)
    assert ret.success

