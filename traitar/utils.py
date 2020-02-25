#!/usr/bin/env python
from __future__ import print_function
# import
## batteries
import os
import sys
import logging
import subprocess
import multiprocessing as mp


def execute_commands(commands, parallel = 1, joblog = None):
    """ (parallel) execution of >=1 subcommand """
    if len(commands) == 0:            
        return
    if parallel > 1:
        pool = mp.Pool(parallel)
        pool.map(execute_command, commands)
    else:
        map(execute_command, commands)
        
def execute_command(cmd):
    """ subcommand execution """
    logging.info('CMD: {}'.format(cmd))
    p = subprocess.Popen(cmd,
                         stderr = subprocess.PIPE,
                         stdout = subprocess.PIPE,
                         shell = True) 
    output, err = p.communicate()
    rc = p.returncode
    if rc != 0:
        print('--- Subprocess ERROR ---')
        print(output.decode())
        print(err.decode())
        print('------------------------')
        sys.exit(1)
            
