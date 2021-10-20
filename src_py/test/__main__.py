# PYTHON Wrapper for Clothoids
# 
# License MIT - See LICENSE file
# 
# 2019 Matteo Ragni, Claudio Kerov Ghiglianovich,
#      Enrico Bertolazzi, Marco Frego

import os
import unittest


if __name__ == "__main__":
    test_path = os.path.normpath(os.path.join(__file__, ".."))
    suite = unittest.defaultTestLoader.discover(test_path)    
    unittest.TextTestRunner(verbosity=2).run(suite)
