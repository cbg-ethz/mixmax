import unittest
from pathlib import Path
import importlib.util
import logging

import numpy as np
import loompy

logging.basicConfig(format="{asctime} - {levelname} - {message}", style="{", datefmt="%Y-%m-%d %H:%M",level=logging.DEBUG)

script_directory = Path(__file__).resolve().parent

split_loom_file_path = (script_directory.parent / "workflow" / "scripts" / "split_loom_files.py").as_posix()
spec = importlib.util.spec_from_file_location("split_loom_file", split_loom_file_path)
target = importlib.util.module_from_spec(spec)
spec.loader.exec_module(target)
split_loom_file = target.split_loom_file


class TestSplitLoomFiles(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TestSplitLoomFiles, self).__init__(*args, **kwargs)
        self.filename = "split_loom_data_test.loom"
        #matrix = np.arange(4).reshape(2,2)
        #row_attrs = { "SomeRowAttr": np.arange(2), "OtherRowAttr": ['A','B'] }
        #col_attrs = { "SomeColAttr": np.arange(2), "OtherColAttr": ['C','D'] }
        #other_layer = np.arange(4, 8).reshape(2,2)
        #loompy.create((script_directory / self.filename).as_posix(), layers = {"":matrix, "other": other_layer}, row_attrs = row_attrs, col_attrs = col_attrs)

    def test_split_loom_files(self):
        loom_file = script_directory / self.filename
        split_loom_file(0.5, script_directory / self.filename, loom_file.parent, seed = 42)
        with loompy.connect(script_directory / 'temp' / f'{Path(self.filename).stem}_split1.loom') as ds:
            self.assertEqual(ds.shape[1], 1)
            self.assertEqual(ds.shape[0], 2)
            self.assertTrue(np.all(ds[:,:] == np.array([[0], [2]])))
            self.assertTrue(np.all(ds.layers["other"][:,:] == np.array([[4], [6]])))
            self.assertTrue(np.all(ds.ca["SomeColAttr"] == np.array([0])))
            self.assertTrue(np.all(ds.ca["OtherColAttr"] == np.array(['C'])))
            self.assertTrue(np.all(ds.ra["SomeRowAttr"] == np.array([0, 1])))
            self.assertTrue(np.all(ds.ra["OtherRowAttr"] == np.array(['A', 'B'])))

        with loompy.connect(script_directory / 'temp' / f'{Path(self.filename).stem}_split2.loom') as ds:
            self.assertEqual(ds.shape[1], 1)
            self.assertEqual(ds.shape[0], 2)
            self.assertTrue(np.all(ds[:,:] == np.array([[1], [3]])))
            self.assertTrue(np.all(ds.layers["other"][:,:] == np.array([[5], [7]])))
            self.assertTrue(np.all(ds.ca["SomeColAttr"] == np.array([1])))
            self.assertTrue(np.all(ds.ca["OtherColAttr"] == np.array(['D'])))
            self.assertTrue(np.all(ds.ra["SomeRowAttr"] == np.array([0, 1])))
            self.assertTrue(np.all(ds.ra["OtherRowAttr"] == np.array(['A', 'B'])))

if __name__ == '__main__':
    unittest.main()

