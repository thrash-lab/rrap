import subprocess
import os

base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sum_arg.py")
print(base_path)
subprocess.run("python3 {0} 1 2 --sum".format(base_path), shell=True)