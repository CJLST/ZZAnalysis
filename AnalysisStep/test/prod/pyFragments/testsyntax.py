#!/usr/bin/env python
import os

class DummyThingThatHasAnyAttribute(object):
    def __getattr__(self, attr):
        return self
    def __call__(self, *args, **kwargs):
        return self

#this means that something like
#process.ZZTree.lheProbabilities.extend(theLHEProbabilities)
#will succeed
process = DummyThingThatHasAnyAttribute()
#for cms.untracked.VLuminosityBlockRange
cms = DummyThingThatHasAnyAttribute()

for filename in os.listdir("."):
    if not filename.endswith(".py"):
        continue
    if filename == os.path.basename(__file__):
        continue #otherwise infinite loop!
    execfile(filename)

print "All files are valid, as far as I can tell."
