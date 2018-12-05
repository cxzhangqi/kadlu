
import pytest
import os

path_to_assets = os.path.join(os.path.dirname(__file__),"assets")

@pytest.fixture
def one():
    return 1

