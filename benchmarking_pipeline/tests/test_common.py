import unittest
from pathlib import Path
import importlib.util

script_directory = Path(__file__).resolve().parent

path = (script_directory.parent / "workflow" / "rules" / "common.smk").as_posix()
spec = importlib.util.spec_from_file_location("common", path)
target = importlib.util.module_from_spec(spec)
spec.loader.exec_module(target)
get_split_files = target.get_split_files
determine_number_of_different_donors = target.determine_number_of_different_donors


class TestCommon(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestCommon, self).__init__(*args, **kwargs)
        self.sample_list = script_directory / "sample_list_test.txt"
        
    def test_get_split_files(self):
        split_files = get_split_files(self.sample_list)
        self.assertEqual(split_files, [self.sample_list.parent / "loom1_split1.loom", self.sample_list.parent / "loom1_split2.loom", self.sample_list.parent / "loom2_split1.loom",self.sample_list.parent / "loom2_split2.loom"])

    def test_determine_number_of_different_donors(self):
        self.assertEqual(target.determine_number_of_different_donors(script_directory / 'pools_test.txt'), "(0.0)", 4)

if __name__ == '__main__':
    unittest.main()