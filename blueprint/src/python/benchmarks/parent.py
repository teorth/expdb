
# Hack for import python files from parent directory
import os
import sys
sys.path.insert(1, os.path.join(sys.path[0], '..'))