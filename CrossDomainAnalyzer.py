#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class CrossDomainAnalyzer(object):
    def __init__(self, fname):
        self._fname = fname
        self._data  = []
        for l in open(fname, 'r'): self._data.append(l[:-1].split(','))
        self._center = float(self._data[365][1])
        self._span   = float(self._data[366][1])
        self._ref    = float(self._data[397][1])
        self._att    = float(self._data[616][1])
        self._unit   = self._data[791][3]
        self._point  = [int(x) for x in np.array(self._data[792:1793])[:,1]]
        self._level  = [float(x) for x in np.array(self._data[792:1793])[:,3]]
    def center(self):
        return self._center
    def span(self):
        return self._span
    def ref(self):
        return self._ref
    def att(self):
        return self._att
    def unit(self):
        return self._unit
    def point(self):
        return self._point
    def level(self):
        return self._level
    def freq(self):
        start = self._center - self._span/2
        width = self._span/(len(self._point)-1)
        return [start + width*p for p in self._point]
    def baseline(self, N=100):
        """
        Mask +- N(=100) point of signal peak
        """
        lvl = list(self.level())
        ind = lvl.index(max(lvl))
        if ind <= N:
            b = self.level()[ind+N:]
        elif ind <= len(self._point)-1 - N:
            b = self.level()[:ind-N]
        else:
            b = self.level()[:ind-N] + self.level()[ind+N:]
        return np.average(b)
    def fwhm(self):
        lvl = list(self.level())
        ind_max = lvl.index(max(lvl))
        ind = []
        for i, lv in enumerate(self.level()):
            if lv >= lvl[ind_max]-3:
                ind.append(i)
            if len(ind) < 2:
                return self.freq()[ind_max+1] - self.freq()[ind_max-1]
        return self.freq()[ind[-1]] - self.freq()[ind[0]]
    def average(self, denominator=10):
        ret = []
        f_sum = 0
        l_sum = 0
        i = 0
        for f, l in zip(self.freq(), self.level()):
            f_sum += f
            l_sum += l
            i += 1
            if i == denominator:
                ret.append([f_sum/denominator, l_sum/denominator])
                f_sum = 0
                l_sum = 0
                i = 0
        return np.array(ret)
    def signal(self):
        lvl = list(self.level())
        ind = lvl.index(max(lvl))
        return self.freq()[ind], max(lvl)

def plot_csv(fname, save=False):
    cda = CrossDomainAnalyzer(fname)
    fig, ax = plt.subplots()
    ax.plot(cda.freq(), cda.level())
    ax.plot(cda.average()[:,0], cda.average()[:,1], c='red', lw=1.5)
    ax.set_xlim(min(cda.freq()), max(cda.freq()))
    ax.set_ylim(min(cda.level())-3, max(cda.level())+3)
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power [%s]' % cda.unit())
    ax.grid()
    if save:
        plt.savefig('%s.png' % fname[:-4])
    print """Parameters:
    Center frequency : %e Hz
    Span             : %e Hz
    Reference        : %.2f dBm
    Attenuate        : %.1f dB
    Signal frequency : %e Hz
    Signal power     : %.2f dBm
    Baseline         : %.2f dBm
    FWHM             : %e Hz
    """ % (cda.center(), cda.span(), cda.ref(), cda.att(), cda.signal()[0], cda.signal()[1], cda.baseline(), cda.fwhm())
