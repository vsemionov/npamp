
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


import numpy as np

import matplotlib as mpl
mpl.rcParams["backend"] = "Agg" # workaround: the worker process's interpreter crashes on windows after executing the program when graphs are enabled

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D; Axes3D = Axes3D # suppress warning

import params
import unitconv


dat_file_format = "dat"

line_styles = ["-", "o-", "s-", "^-", "--", "o--", "s--", "^--", ":", "o:", "s:", "^:"]


def output_data_2d(basename, X, Y):
    if type(Y) not in [tuple, list]:
        X = [X]
        Y = [Y]
    xlist = []
    for Xel in X:
        xlist += list(Xel)
    Xu = sorted(list(set(xlist)))
    Yu = [np.interp(Xu, Xel, Yel) for (Xel, Yel) in zip(X, Y)]
    with open(basename + '.' + dat_file_format, "wt") as f:
        for i, x in enumerate(Xu):
            line = "%e" % x
            for k in range(len(Y)):
                y = Yu[k][i]
                line += " %e" % y
            f.write(line + "\n")

def output_data_3d(basename, X, Y, Z):
    Z = Z.T
    assert Z.shape == (X.shape + Y.shape)
    with open(basename + '.' + dat_file_format, "wt") as f:
        for i, x in enumerate(X):
            for j, y in enumerate(Y):
                z = Z[i, j]
                f.write("%e %e %e\n" % (x, y, z))
            f.write("\n")

def output_data_3d_ext(basename, XY, YX, Z):
    assert XY.shape == YX.shape == Z.shape
    with open(basename + '.' + dat_file_format, "wt") as f:
        for i in range(Z.shape[0]):
            for j in range(Z.shape[1]):
                x = XY[i, j]
                y = YX[i, j]
                z = Z[i, j]
                f.write("%e %e %e\n" % (x, y, z))
            f.write("\n")

def plot_data(basename, title, (X, xscale, xlim, xlabel), (Y, yscale, ylim, ylabel), legend=None, style="-", vlines=False, grid="both", xvals=None, yvals=None, callback=None):
    if params.output_data:
        output_data_2d(basename, X, Y)
    uxscale, uyscale = unitconv.output_scale(xlabel), unitconv.output_scale(ylabel)
    plt.clf()
    if xscale:
        plt.xscale(xscale)
    if yscale:
        plt.yscale(yscale)
    if xlim is not None:
        xlim = [t * uxscale for t in xlim]
        plt.xlim(xlim)
    if ylim is not None:
        ylim = [t * uyscale for t in ylim]
        plt.ylim(ylim)
    if params.output_plot_titles:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if grid:
        plt.grid(axis=grid)
    kwargs = dict(color = "black") if not params.output_color else {}
    if type(Y) in [tuple, list]:
        for i, (x, y, label) in enumerate(zip(X, Y, legend)):
            x = x * uxscale
            y = y * uyscale
            style = line_styles[i] if params.output_styled_lines else "-"
            marker_stride = max(len(x) // params.out_markers_step_divisor, 1)
            plt.plot(x, y, style, markevery=marker_stride, label=label, **kwargs)
        plt.legend(loc=0)
    else:
        X = X * uxscale
        Y = Y * uyscale
        if vlines:
            plt.vlines(X, 0.0, Y, color="black")
        plt.plot(X, Y, style, **kwargs)
    if xvals is not None:
        for x, label in xvals:
            x *= uxscale
            plt.axvline(x, color="black", linestyle="--")
            if label:
                plt.text(x, 0.0, label)
    if yvals is not None:
        for y, label in yvals:
            y *= uyscale
            plt.axhline(y, color="black", linestyle="--")
            if label:
                plt.text(0.0, y, label)
    if callback is not None:
        callback()
    plt.savefig(basename + '.' + params.output_file_format)

def plot_projection(basename, title, (XY, xlim, xlabel), (YX, ylim, ylabel), (Z, zlim, zlabel), (elev, azim), (stride_x, stride_y)):
    if params.output_data:
        output_data_3d_ext(basename, XY, YX, Z)
    uxscale, uyscale, uzscale = unitconv.output_scale(xlabel), unitconv.output_scale(ylabel), unitconv.output_scale(zlabel)
    plt.clf()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(elev, azim)
    if params.output_plot_titles:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    XY = XY * uxscale
    YX = YX * uyscale
    Z = Z * uzscale
    color = "blue" if params.output_color else "black"
    ax.plot_surface(XY, YX, Z, rstride=stride_x, cstride=stride_y, color="white", edgecolor=color, shade=False)
    if xlim is not None:
        xlim = [t * uxscale for t in xlim]
        ax.set_xlim(xlim)
    if ylim is not None:
        ylim = [t * uyscale for t in ylim]
        ax.set_ylim(ylim)
    if zlim is not None:
        zlim = [t * uzscale for t in zlim]
        ax.set_zlim(zlim)
    plt.savefig(basename + '.' + params.output_file_format)

def plot_error(basename, title, (X, xscale, xlim, xlabel), ((exact, approx), yscale, ylim, ylabel), **kwargs):
    error = np.fabs((exact - approx) / exact)
    plot_data(basename, title, (X, xscale, xlim, xlabel), (error, yscale, ylim, ylabel), **kwargs)

def plot_color(basename, title, (X, xscale, xlim, xlabel), (Y, yscale, ylim, ylabel), (Z, zlim, zlabel), auto_contours=None, contours=None, extra_contours=None, grid="both", xvals=None, yvals=None, colorbar=True, ccolor="black", callback=None):
    # begin patch
    if ccolor == "black":
        output_data = params.output_data
        params.output_data = False
        plot_color(basename, title, (X, xscale, xlim, xlabel), (Y, yscale, ylim, ylabel), (Z, zlim, zlabel), auto_contours=auto_contours, contours=contours, extra_contours=extra_contours, grid=grid, xvals=xvals, yvals=yvals, colorbar=colorbar, ccolor="white", callback=callback)
        params.output_data = output_data
    if ccolor == "white":
        import os
        dirname, filename = os.path.split(basename)
        dirname = os.path.join(dirname, "cc")
        basename = os.path.join(dirname, filename)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
    # end patch
    if params.output_data:
        output_data_3d(basename, X, Y, Z)
    cextend = lambda Arr: np.append(Arr, Arr[-1] + (Arr[-1] - Arr[0]) / (len(Arr) - 1))
    uxscale, uyscale, uzscale = unitconv.output_scale(xlabel), unitconv.output_scale(ylabel), unitconv.output_scale(zlabel)
    plt.clf()
    if xscale:
        plt.xscale(xscale)
    if yscale:
        plt.yscale(yscale)
    if xlim is not None:
        xlim = [t * uxscale for t in xlim]
        plt.xlim(xlim)
    if ylim is not None:
        ylim = [t * uyscale for t in ylim]
        plt.ylim(ylim)
    if params.output_plot_titles:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    X = X * uxscale
    Y = Y * uyscale
    Z = Z * uzscale
    cmap = None if params.output_color else getattr(plt.cm, "gray")
    plt.pcolor(cextend(X), cextend(Y), Z, cmap=cmap)
    if zlim is not None:
        zlim = [t * uzscale for t in zlim]
        plt.clim(zlim)
    if colorbar:
        cb = plt.colorbar()
        if zlabel:
            cb.set_label(zlabel)
    if auto_contours:
        plt.contour(X, Y, Z, auto_contours, colors="gray")
    if contours is not None:
        contours = [c * uzscale for c in contours]
        cs = plt.contour(X, Y, Z, contours, colors=ccolor)
        plt.clabel(cs, fmt=str)
    if extra_contours:
        for Z_extra, val, label in extra_contours:
            cs = plt.contour(X, Y, Z_extra, [val], colors=ccolor)
            if label:
                plt.clabel(cs, fmt = lambda x: label)
    if xvals is not None:
        for x, label in xvals:
            x *= uxscale
            plt.axvline(x, color=ccolor, linestyle="--")
            if label:
                plt.text(x, 0.0, label, color=ccolor)
    if yvals is not None:
        for y, label in yvals:
            y *= uyscale
            plt.axhline(y, color=ccolor, linestyle="--")
            if label:
                plt.text(0.0, y, label, color=ccolor)
    if grid:
        plt.grid(axis=grid)
    if callback is not None:
        callback()
    plt.savefig(basename + '.' + params.output_file_format)

def plot_contour(basename, title, (X, xscale, xlim, xlabel), (Y, yscale, ylim, ylabel), contours, grid="both", xvals=None, yvals=None, callback=None):
    uxscale, uyscale = unitconv.output_scale(xlabel), unitconv.output_scale(ylabel)
    plt.clf()
    if xscale:
        plt.xscale(xscale)
    if yscale:
        plt.yscale(yscale)
    if xlim is not None:
        xlim = [t * uxscale for t in xlim]
        plt.xlim(xlim)
    if ylim is not None:
        ylim = [t * uyscale for t in ylim]
        plt.ylim(ylim)
    if params.output_plot_titles:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    X = X * uxscale
    Y = Y * uyscale
    for Z, val, label in contours:
        cs = plt.contour(X, Y, Z, [val], colors="black")
        if label:
            plt.clabel(cs, fmt = lambda x: label)
    if xvals is not None:
        for x, label in xvals:
            x *= uxscale
            plt.axvline(x, color="black", linestyle="--")
            if label:
                plt.text(x, 0.0, label, rotation="vertical")
    if yvals is not None:
        for y, label in yvals:
            y *= uyscale
            plt.axhline(y, color="black", linestyle="--")
            if label:
                plt.text(0.0, y, label)
    if grid:
        plt.grid(axis=grid)
    if callback is not None:
        callback()
    plt.savefig(basename + '.' + params.output_file_format)
