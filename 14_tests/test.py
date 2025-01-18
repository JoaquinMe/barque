import pandas as pd
import os
import subprocess
import re

result=subprocess.run(["shasum",f"14_tests/test.py"],capture_output=True)
checksum=re.match(r"^(\w*) ",result.stdout.decode("utf-8"))
print([checksum.group(1)])
