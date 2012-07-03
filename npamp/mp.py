
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
#  * Neither the name of the copyright holder nor the names of the contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
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
import core


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

def monitor_thread(in_conn):
    try:
        in_conn.recv()
    except EOFError:
        os.kill(os.getpid(), signal.SIGTERM)

def input_thread(io_queue):
    while True:
        s = io_queue.get()
        print s

def task_init(io_queue, (in_conn, out_conn), conf):
    mpout = MPOutput(io_queue.put)
    sys.stdout = mpout
    sys.stderr = mpout
    out_conn.close()
    thr = threading.Thread(target=monitor_thread, args=(in_conn,))
    thr.daemon = True
    thr.start()
    params.__dict__.update(conf)

def dispatch_task((task, args)):
    return task(*args)

class TaskPool(object):
    def __init__(self, num_tasks, parallel, task_params):
        num_tasks = num_tasks or None
        if not parallel:
            num_tasks = 1
        self.pool = None
        self.num_tasks = num_tasks or multiprocessing.cpu_count()
        if num_tasks != 1:
            conf = core.copy_conf(task_params)
            io_queue = multiprocessing.Queue()
            thr = threading.Thread(target=input_thread, args=(io_queue,))
            thr.daemon = True
            thr.start()
            in_conn, out_conn = multiprocessing.Pipe(False)
            self._out_conn = out_conn
            self.pool = multiprocessing.Pool(processes=num_tasks, initializer=task_init, initargs=(io_queue, (in_conn, out_conn), conf))
    
    def close(self):
        if self.pool is not None:
            self.pool.close()
            self.pool.terminate()
            self.pool.join()
        self.pool = None
    
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
        mapper = self.pool.map if self.pool is not None else map
        task_results = mapper(dispatch_task, dispatch_args)
        if task_results:
            if type(task_results[0]) in [tuple, list]:
                results = map(to_array[ndims-1], zip(*task_results))
            else:
                results = to_array[ndims-1](task_results)
        return results
