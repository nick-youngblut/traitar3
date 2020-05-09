import pytest
import tempfile

@pytest.fixture(scope='session')
def tests_tmpdir():
    return tempfile.mkdtemp()
