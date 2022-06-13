# NOT WORKING!!!

# Initialize HEAsoft from within a Python script like so:
# >> import heasoft
# >> heasoft.init()
#
import os
import subprocess
thisdir     = os.path.dirname( os.path.abspath(__file__) )
start_xspec = os.path.join(thisdir, 'start_xspec.sh')
system      = "localhost"
cmd         = [start_xspec, system]
popen       = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

#from xspec import *
