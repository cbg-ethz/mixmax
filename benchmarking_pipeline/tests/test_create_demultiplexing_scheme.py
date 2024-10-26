import unittest
from pathlib import Path
import importlib.util

script_directory = Path(__file__).resolve().parent

path = (script_directory.parent / "workflow" / "scripts" / "create_demultiplexing_scheme.py").as_posix()
spec = importlib.util.spec_from_file_location("create_demultiplexing_scheme", path)
target = importlib.util.module_from_spec(spec)
spec.loader.exec_module(target)
multiplexing_scheme_format2pool_format = target.multiplexing_scheme_format2pool_format
select_samples_for_pooling = target.select_samples_for_pooling
define_demultiplexing_scheme_optimal_case = target.define_demultiplexing_scheme_optimal_case

class TestCreateDemultiplexingScheme(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCreateDemultiplexingScheme, self).__init__(*args, **kwargs)
        self.multiplexing_scheme = {1: (0, 1, 0), 2: (0, 2, 0), 3: (1, 2, 0), 4: (0, 0, 0), 5: (1, 1, 0), 6: (2, 2, 0), 7: (0, 1, 1), 8: (0, 2, 1), 9: (1, 2, 1), 10: (0, 0, 1), 11: (1, 1, 1), 12: (2, 2, 1)}
        self.multiplexing_scheme_pool_format = {(0, 0): [1,2, 4,-4], (1, 0): [-1,3,5,-5], (2, 0): [-2, -3, 6, -6], (0, 1): [7, 8, 10, -10], (1,1): [-7, 9, 11, -11], (2, 1): [-8, -9, 12, -12]}
        self.mock_samples = ['/test/folder/loom1.loom', '/test/folder/loom2.loom', '/test/folder/loom3.loom']
        self.mock_path = '/another/test/folder'
        self.another_pool_scheme = {(0,0): [1,-1,3], (1,0): [2,-2,-3]}

    def test_multiplexing_scheme_format2pool_format(self):
        self.assertEqual(multiplexing_scheme_format2pool_format(self.multiplexing_scheme), self.multiplexing_scheme_pool_format)
    
    def test_select_samples_for_pooling(self):
        self.assertEqual(select_samples_for_pooling(self.another_pool_scheme, self.mock_path, self.mock_samples), {'(0.0)': ['/another/test/folder/loom1_split1.loom', '/another/test/folder/loom1_split2.loom', '/another/test/folder/loom3_split1.loom'],
                                                                                                             '(1.0)': ['/another/test/folder/loom2_split1.loom', '/another/test/folder/loom2_split2.loom', '/another/test/folder/loom3_split2.loom']})

    def test_define_demultiplexing_scheme_optimal_case(self):
        self.assertEqual(define_demultiplexing_scheme_optimal_case(maximal_number_of_samples = 3, maximal_pool_size = 2, n_samples = 3), {1: (0, 1, 0), 2: (0, 0, 0), 3: (1, 1, 0)})

if __name__ == '__main__':
    unittest.main()