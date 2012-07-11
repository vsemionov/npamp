
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


import re

from decimal import Decimal


quantities = (
    ("diameter", 2.0),
)

units = {
    "cm": 1.0e-2,
    "mm": 1.0e-3,
    "um": 1.0e-6,
    "nm": 1.0e-9,
    "pm": 1.0e-12,
    "fm": 1.0e-15,
    "ms": 1.0e-3,
    "us": 1.0e-6,
    "ns": 1.0e-9,
    "ps": 1.0e-12,
    "fs": 1.0e-15,
    "mJ": 1.0e-3,
    "uJ": 1.0e-6,
    "nJ": 1.0e-9,
    "cm^2": 1.0e-4,
    "cm^-3": 1.0e6,
    "J/cm^2": 1.0e4,
    "W/cm^3": 1.0e6,
    "s^-1": 1.0,
    "cm^-3 s^-1": 1.0e6,
    "%": 1.0e-2,
    
    "cm^{-3}": 1.0e6,
    "{\mu}s": 1.0e-6,
    "cm^{-3}s^{-1}": 1.0e6,
    "s^{-1}": 1.0,
    r"\%": 1.0e-2,
}


def convert_from_input(label, value):
    value = Decimal(value)
    for external, coef in quantities:
        if re.search(r"\b%s\b.*\[.+\]\s*:\s*$" % external, label):
            value /= Decimal(repr(coef))
            break
    m = re.search(r"\[\s*([^\[\]]+)\s*\]\s*:\s*$", label)
    if m:
        unit = m.group(1)
        if unit in units:
            coef = units[unit]
            value *= Decimal(repr(coef))
    value = float(value)
    return value

def convert_to_input(label, value):
    value = Decimal(repr(value))
    for external, coef in quantities:
        if re.search(r"\b%s\b.*\[.+\]\s*:\s*$" % external, label):
            value *= Decimal(repr(coef))
            break
    m = re.search(r"\[\s*([^\[\]]+)\s*\]\s*:\s*$", label)
    if m:
        unit = m.group(1)
        if unit in units:
            coef = units[unit]
            value /= Decimal(repr(coef))
    value = float(value)
    return value

def output_scale(label):
    if label is None:
        return 1.0
    scale = 1.0
    for external, coef in quantities:
        if re.search(r"\b%s\b.*\[([^\[]+)\]\s*\$?\s*$" % external, label):
            scale *= coef
            break
    m = re.search(r"\[([^\[]+)\]\s*\$?\s*$", label)
    if m:
        unit = m.group(1)
        if unit in units:
            coef = units[unit]
            scale /= coef
    return scale

def print_result(fmt, runits, rvalues):
    assert len(rvalues) >= len(runits), "wrong result length"
    punits = runits
    if len(runits) > 0:
        if len(runits) == 1:
            runits = runits * len(rvalues)
        for external, coef in quantities:
            if re.search(r"\b%s\b.*\[.+\]" % external, fmt):
                rvalues = tuple([float(Decimal(repr(value)) * Decimal(repr(coef))) if value is not None else None for value in rvalues])
                break
        rvalues = tuple([(v / (units[u] if u in units else 1.0)) if v is not None else None for (u, v) in zip(runits, rvalues)])
    print fmt.format(*(punits + rvalues))
