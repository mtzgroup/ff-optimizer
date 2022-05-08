import pytest

def test_findRepeatIndex():
    pass

def test_readCharges_out():
    os.chdir(os.path.join(os.path.dirname(__file__),"resp"))
    options = {}
    options['sampleDir'] = "test"
    options['mol2'] = "dasa.mol2"
    options['mode'] = 1
    respPriors = RespPriors(options)
    with open("resp.out",'r') as f:
        lines = f.readlines()
    esp, resp = respPriors(lines)

    

