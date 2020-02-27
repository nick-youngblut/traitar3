#!/usr/bin/env python
from __future__ import print_function
# import
## batteries
import os
import sys
import time
import signal
import logging
import subprocess
import multiprocessing as mp
from functools import wraps

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def execute_commands(commands, parallel = 1, joblog = None):
    """ (parallel) execution of >=1 subcommand """
    if len(commands) == 0:            
        return
    if parallel > 1:
        pool = mp.Pool(parallel, init_worker)
        pool.map(execute_command, commands)
    else:
        for cmd in commands:
            execute_command(cmd)
        
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
            

def retry(ExceptionToCheck, tries=4, delay=3, backoff=2, logger=None):
    """Retry calling the decorated function using an exponential backoff.

    http://www.saltycrane.com/blog/2009/11/trying-out-retry-decorator-python/
    original from: http://wiki.python.org/moin/PythonDecoratorLibrary#Retry

    :param ExceptionToCheck: the exception to check. may be a tuple of
        exceptions to check
    :type ExceptionToCheck: Exception or tuple
    :param tries: number of times to try (not retry) before giving up
    :type tries: int
    :param delay: initial delay between retries in seconds
    :type delay: int
    :param backoff: backoff multiplier e.g. value of 2 will double the delay
        each retry
    :type backoff: int
    :param logger: logger to use. If None, print
    :type logger: logging.Logger instance
    """
    def deco_retry(f):

        @wraps(f)
        def f_retry(*args, **kwargs):
            mtries, mdelay = tries, delay
            while mtries > 1:
                try:
                    return f(*args, **kwargs)
                except ExceptionToCheck as e:
                    msg = '{}; Retrying in {} seconds...'.format(str(e), mdelay)
                    logging.warning(msg)
                    time.sleep(mdelay)
                    mtries -= 1
                    mdelay *= backoff
            return f(*args, **kwargs)

        return f_retry  # true decorator

    return deco_retry
