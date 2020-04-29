#!python
#!/usr/bin/env python

import cubit
import cubit2dg2d
import os
import sys
from math import *
from numpy import *


reload(cubit2dg2d)

cubit.init([""])
cubit.cmd('set duplicate block elements on')
command = 'open "testtri.cub"'
cubit.cmd(command)

cubit2dg2d.mesh()
