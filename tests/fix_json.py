import os
import json
from qcparse import parse

for f in os.listdir():
    if f.endswith(".json"):
        with open(f, 'r') as js:
            d = json.load(js)
        name = f.split(".")[0]
        tcout = name + ".out"
        if not os.path.isfile(tcout):
            stdout = d['stdout']
        else:
            with open(tcout, 'r') as out:
                stdout = ''.join(out.readlines())
            os.remove(tcout)
        with open(name + ".out", 'w') as out:
            out.write("TeraChem vTestSuite\n")
            out.write("Git Version: TestSuite\n")
            out.write(stdout)
        print(tcout)
        results = parse(tcout, "terachem")
        os.remove(f)
        #os.rename(f, name + "_old.json")
        #results.save(f)
