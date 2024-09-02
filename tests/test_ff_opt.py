import os

from ff_optimizer import ff_opt, model, optengine

from .test_inputs import getDefaults


def monkeyDoMMSampling(self, i):
    pass


def monkeyDoQMCalculations(self, i):
    pass


def monkeyMakeFBTargets(self, i):
    return []


def monkeyCopyLeapFiles(self, f):
    pass


def monkeyGetSampleFolders(self, i):
    pass


def monkeyCopyQMResults(self, f1, f2):
    pass


def monkeyCopyResults(self, i):
    pass


def monkeySetupInputFiles(self, i):
    pass


def monkeyRunValidPrevious(self, i):
    self.validPrevious.append(1.0)


def monkeyRunTraining(self, i):
    self.train.append(0.5)
    with open("cycle.txt", "w") as f:
        f.write(str(i))


def monkeyCopyResults(self, i):
    pass


def monkeyRunValid(self, i):
    self.valid.append(1.0)


def monkeyInitModel(self, inp):
    self.converged = False
    self.optEngine = self.initializeOptEngine(inp)
    self.inp = inp
    self.optdir = inp.optdir
    self.home = self.optdir


def monkeyInitOpt(self, inp):
    self.converged = False
    self.inp = inp
    self.determineRestart()


def monkeyInitialCycle(self):
    pass


def monkeyGetFinalValidations(self, j):
    return 1


def monkeyCopyFinalResults(self, best):
    pass


def test_optimize_converge(monkeypatch):
    os.chdir(os.path.dirname(__file__))
    os.chdir("model")
    os.chdir("test1")
    inp = getDefaults()
    inp.validinitial = False
    inp.initialtraining = False
    inp.skipchecks = True
    inp.optdir = "."
    inp.toYaml("temp.yaml")
    monkeypatch.setattr(model.Model, "doMMSampling", monkeyDoMMSampling)
    monkeypatch.setattr(model.Model, "doQMCalculations", monkeyDoQMCalculations)
    monkeypatch.setattr(model.Model, "makeFBTargets", monkeyMakeFBTargets)
    monkeypatch.setattr(model.Model, "copyLeapFiles", monkeyCopyLeapFiles)
    monkeypatch.setattr(model.Model, "getSampleFolders", monkeyGetSampleFolders)
    monkeypatch.setattr(model.Model, "copyQMResults", monkeyCopyQMResults)
    monkeypatch.setattr(model.Model, "__init__", monkeyInitModel)
    monkeypatch.setattr(model.Model, "initialCycle", monkeyInitialCycle)
    monkeypatch.setattr(optengine.OptEngine, "__init__", monkeyInitOpt)
    monkeypatch.setattr(optengine.OptEngine, "copyResults", monkeyCopyResults)
    monkeypatch.setattr(optengine.OptEngine, "setupInputFiles", monkeySetupInputFiles)
    monkeypatch.setattr(optengine.OptEngine, "runValidPrevious", monkeyRunValidPrevious)
    monkeypatch.setattr(optengine.OptEngine, "runTraining", monkeyRunTraining)
    monkeypatch.setattr(optengine.OptEngine, "copyResults", monkeyCopyResults)
    monkeypatch.setattr(optengine.OptEngine, "runValid", monkeyRunValid)
    monkeypatch.setattr(
        optengine.OptEngine, "getFinalValidations", monkeyGetFinalValidations
    )
    monkeypatch.setattr(optengine.OptEngine, "copyFinalResults", monkeyCopyFinalResults)
    ff_opt.optimize("temp.yaml")
    with open("cycle.txt", "r") as f:
        cycle = int(f.readline().split()[0])
    os.remove("temp.yaml")
    os.remove("cycle.txt")
    assert cycle == 5
