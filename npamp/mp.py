
# Copyright (C) 2012 Victor Semionov
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import sys
import os

import signal
import threading
import multiprocessing

import numpy as np

import params
import cfg


# The MP IO subsystem has some flaws and needs to be redesigned.
# It contains some concurrency issues (race conditions after the occurrence of errors) and is not completely reliable.
# As a workaround for any eventual issues, use the command-line interface (as opposed to the GUI).


class MPOutput(object):
    
    def __init__(self, send):
        self.send = send
        self.buff = str()
    
    def write(self, s):
        buff = self.buff + s
        lines = buff.split('\n')
        self.buff = lines[-1]
        for line in lines[:-1]:
            self.send(line)
    
    def flush(self):
        if self.buff:
            self.send(self.buff)
            self.buff = str()


class _sentinel_type(object): pass
def _input_thread(io_queue, outfile, oldfiles):
    try:
        while True:
            s = io_queue.get()
            if type(s) is _sentinel_type:
                break # workaround (for python bugs?)
            print >>outfile, s
    finally:
        sys.stdout, sys.stderr = oldfiles

def monitor_thread(in_conn):
    try:
        in_conn.recv()
    except EOFError:
        os.kill(os.getpid(), signal.SIGTERM)

def _task_init(io_queue, (in_conn, out_conn), conf):
    mpout = MPOutput(io_queue.put)
    sys.stdout = sys.stderr = mpout
    
    out_conn.close()
    thr = threading.Thread(target=monitor_thread, args=(in_conn,))
    thr.daemon = True
    thr.start()
    
    params.__dict__.update(conf)

def _dispatch_task((task, args)):
    return task(*args)

class TaskPool(object):
    def __init__(self, num_tasks, parallel, task_params):
        num_tasks = num_tasks or None
        if not parallel:
            num_tasks = 1
        self._pool = None
        self.num_tasks = num_tasks or multiprocessing.cpu_count()
        if num_tasks != 1:
            conf = cfg.copy_conf(task_params)
            self._io_queue = multiprocessing.Queue()
            self._oldfiles = sys.stdout, sys.stderr
            self._input_thread = threading.Thread(target=_input_thread, args=(self._io_queue, sys.stdout, self._oldfiles))
            self._input_thread.daemon = True
            self._input_thread.start()
            mpout = MPOutput(self._io_queue.put)
            sys.stdout = sys.stderr = mpout
            try:
                in_conn, out_conn = multiprocessing.Pipe(False)
                self._out_conn = out_conn
                self._pool = multiprocessing.Pool(processes=num_tasks, initializer=_task_init, initargs=(self._io_queue, (in_conn, out_conn), conf))
            except:
                sys.stdout, sys.stderr = self._oldfiles
                raise
    
    def close(self):
        if self._pool is not None:
            self._pool.close()
            self._pool.terminate()
            self._pool.join()
            
            sys.stdout, sys.stderr = self._oldfiles # avoid a race condition
            self._io_queue.put(_sentinel_type())
            self._input_thread.join()
        self._pool = None
    
    def __del__(self):
        self.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, t, v, tb):
        self.close()
    
    def parallel_task(self, task, coords, var_args, const_args):
        ndims = len(coords)
        to_array = [
            lambda res: np.array(res),
            lambda res: np.array(res).reshape(count_x, count_y),
            ]
        results = []
        if ndims == 2:
            X, Y = coords
            count_x, count_y = len(X), len(Y)
            count_tot = count_x * count_y
            indices = np.indices((count_x, count_y))
            i_args = indices[0].reshape(count_tot)
            j_args = indices[1].reshape(count_tot)
            MG = np.meshgrid(X, Y)
            x_args = MG[0].T.reshape(count_tot)
            y_args = MG[1].T.reshape(count_tot)
            va_list = [va.reshape(count_tot) for va in var_args]
            ca_list = zip(*([const_args] * count_tot))
            task_args = zip(zip(i_args, j_args), zip(x_args, y_args), *(va_list + ca_list))
        elif ndims == 1:
            X, = coords
            count_x = len(X)
            count_tot = count_x
            i_args = range(count_x)
            x_args = X
            va_list = list(var_args)
            ca_list = zip(*([const_args] * count_tot))
            task_args = zip(i_args, x_args, *(va_list + ca_list))
        else:
            assert False, "invalid number of coordinates"
        dispatch_args = zip((task,) * count_tot, task_args)
        mapper = self._pool.map if self._pool is not None else map
        task_results = mapper(_dispatch_task, dispatch_args)
        if task_results:
            if type(task_results[0]) in [tuple, list]:
                results = map(to_array[ndims-1], zip(*task_results))
            else:
                results = to_array[ndims-1](task_results)
        return results
