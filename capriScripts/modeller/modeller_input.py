from modeller import *
from modeller.automodel import *
env = environ()
a = automodel(env,alnfile  = 'mdl_vs_targt.ali',knowns='targt',sequence = 'target')
a.starting_model= 1
a.ending_model  = 1
a.make()
