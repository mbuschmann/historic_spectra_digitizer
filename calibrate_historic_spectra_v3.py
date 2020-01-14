#!/usr/bin/python3
from __future__ import print_function, division
import os, pickle, sys
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
from scipy.optimize import curve_fit
import pdb


class Calibrator():
    #def read_params(self, lbl):
    #    self.params[lbl] = {}
    #    fname = os.path.join(self.parampath, 'digitizing_params_'+lbl+'.dat')
    #    with open(fname, 'r') as f:
    #        ll = f.readlines()
    #    for l in ll:
    #        lsp = l.split(' ')
    #        if len(lsp)==2:
    #            self.params[lbl][lsp[0]] = float(lsp[1])
    #        elif len(lsp)==3:
    #            self.params[lbl][lsp[0]] = (float(lsp[1]), float(lsp[2]))
    #        else:
    #            pass
    #    for k in ['endpoints', 'subset', 'yoffset']:
    #        if k in self.params[lbl].keys():
    #            pass
    #        else:
    #            if k=='endpoints':
    #                self.params[lbl]['endpoints'] = (2859, 2927)
    #            elif k=='subset':
    #                self.params[lbl]['subset'] = (0, len(self.data[lbl]['data']))
    #            elif k=='yoffset':
    #                self.params[lbl]['yoffset'] = 0
    #            else: pass
    
    #def write_params(self, event):
    #    lbl = self.currlabel
    #    fname = os.path.join(self.parampath, 'digitizing_params_'+lbl+'.dat')
    #    with open(fname, 'a') as f:
    #        p = 'subset'
    #        tmp = p+' '+str(self.params[lbl][p][0])+' '+str(self.params[lbl][p][1])+'\n'
    #        f.write(tmp)
    #        p = 'yoffset'
    #        tmp = p+' '+str(self.params[lbl][p])+'\n'
    #        f.write(tmp)
    #        p = 'endpoints'
    #        tmp = p+' '+str(self.params[lbl][p][0])+' '+str(self.params[lbl][p][1])+'\n'
    #        f.write(tmp)
    #    print('changed: ',fname)
    #    #print(type(self.params['d']['endpoints']))
    #    #print(type(self.params['d']['endpoints'])==tuple)
    #    #print(type(self.params['d']['yoffset']))
    #    #print(type(self.params['d']['yoffset'])==float)
    
    def print_params(self, event):
        for lbl in self.labellist:
            for p in self.params[lbl].keys():
                print(lbl, p, self.params[lbl][p])
    
    def xshift(self, val):
        a = self.aslider.val
        b = self.bslider.val
        self.s, self.e = (a, b)
        y = self.y
        if self.radio.value_selected=='linear lamda':
            xn = 1e4/np.linspace(1e4/a, 1e4/b, len(y))
        else: 
            xn = np.linspace(a, b, len(y))
        o = self.yoffset
        self.digiline.set_xdata(xn)
        self.digiline.set_ydata((y-o)/np.max(y-o))
        self.fig.canvas.draw_idle()
        
    def zshift(self, val):
        self.yoffset = self.zslider.val
        self.y2 = (self.y-self.yoffset)/np.max(self.y-self.yoffset)
        self.digiline.set_ydata(self.y2)
        self.fig.canvas.draw_idle()

#    def draw_sliders(self):
#        #self.axas.clear()
#        #self.axbs.clear()
#        #self.axzs.clear()
#        #self.aslider = Slider(self.axas, 'L', self.s-(self.e-self.s), self.s+(self.e-self.s), valinit=self.s, facecolor='gray')
#        #self.bslider = Slider(self.axbs, 'R', self.e-(self.e-self.s), self.e+(self.e-self.s), valinit=self.e, facecolor='gray')
#        #self.zslider = Slider(self.axzs, 'z', 0, 3000, valinit=self.yoffset, facecolor='gray')
#        #self.aslider.on_changed(self.xshift)
#        #self.bslider.on_changed(self.xshift)
#        #self.zslider.on_changed(self.zshift)
#        pass

#    def changesrc(self, lbl):
#        self.currlabel = self.radio.value_selected
#        e0 = self.params[self.currlabel]['endpoints'][0]
#        e1 = self.params[self.currlabel]['endpoints'][1]
#        self.draw_sliders()
#        s0 = int(self.params[self.currlabel]['subset'][0])
#        s1 = int(self.params[self.currlabel]['subset'][1])
#        self.x2 = np.linspace(e0, e1, len(self.data[self.currlabel]['data'][s0:s1]))
#        self.y2 = data[self.currlabel]['data'][s0:s1]-self.params[self.currlabel]['yoffset']
#        self.digiline.set_xdata(self.x2)
#        self.digiline.set_ydata((self.y2)/np.max(self.y2))
#        self.fig.canvas.draw_idle()

    def reset(self, event):
        self.aslider.reset()
        self.bslider.reset()
    
    def get_new_xy(self):
        self.x2 = np.linspace(self.s, self.e, len(self.data.y))
        self.y2 = (self.y-self.yoffset)/np.max(self.y-self.yoffset)
        
    def print_sfit_readable_spectrum(self,event, sza=63.0, latlon=(46.55, 7.98), d=dt.datetime(1951, 4, 15, 7, 30, 0), res=0.25, apo='TRI', sn=100.0, rearth=6377.9857):
        spc = self.yreduced
        wvn_bounds = [np.min(self.xcal), np.max(self.xcal)]
        fname = os.path.join(self.path, self.fname[:-4]+'_cal_lines.dat')
        with open(fname, 'w') as f:
            for i, j in zip(self.ax1_lines, self.ax2_lines):
                f.write('%4.4f %4.4f\n'%(i,j))
        print('Wrote file', fname)
        fname = os.path.join(self.path, self.fname[:-4]+'_calibrated.dat')
        #pdb.set_trace()
        s = ' %4.2f  %8.4f  %4.2f  %5.2f  %5i\n'%(sza, rearth, latlon[0], latlon[1] , sn)
        s = s + d.strftime(' %Y %m %d %H %M %S\n')
        s = s + d.strftime(' %d/%m/%Y, %H:%M:%S')+', RES=%5.4f  APOD FN = %3s\n'%(res, apo)
        s = s + ' %7.3f %7.3f %11.10f %7i'%(wvn_bounds[0], wvn_bounds[1], (wvn_bounds[1]-wvn_bounds[0])/float(len(spc)), len(spc))
        for i in spc:
            s+='\n %8.5f'%(i)
        with open(fname, 'w') as f:
            f.write(s)
        print('Wrote file', fname)
    
    def read_cali_lines(self, event):
        fname = os.path.join(self.path, self.fname[:-4]+'_cal_lines.dat')
        with open(fname, 'r') as f:
            ll = f.readlines()
        self.ax1_lines = [float(i.split()[0]) for i in ll]
        self.ax2_lines = [float(i.split()[1]) for i in ll]
        print(self.ax1_lines)
        print(self.ax2_lines)
        print('Read calibration lines from', fname)
        self.plotlines()
        self.calibrate(None)
    
    def reduce_points(self, event, res=0.05):
        xo = self.xcal
        yo = self.y2
        yn = []
        xmin, xmax = np.min(self.xcal), np.max(self.xcal)
        xn = np.linspace(xmin, xmax, (xmax-xmin)/res)
        for x in xn:
            yn.append(np.median(yo[(xo>x-res/2) & (xo<x+res/2)]))
        self.yreduced = np.array(yn)
        self.xreduced = xn
        self.ax1.plot(self.xreduced, self.yreduced, '-')
        self.fig.canvas.draw_idle()

        
    def plotlines(self):
        self.ax1l = self.ax1.vlines(self.ax1_lines, 0,1)
        self.ax2l = self.ax2.vlines(self.ax2_lines, 0,1)
        self.fig.canvas.draw_idle()
        
    def onclick(self, event):
        #print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % ('double' if event.dblclick else 'single', event.button,event.x, event.y, event.xdata, event.ydata))
        #pdb.set_trace()
        #print('%s'%event.inaxes)
        if event.dblclick:
            if event.inaxes==self.ax1:
                self.ax1_lines.append(event.xdata)
            elif event.inaxes==self.ax2:
                self.ax2_lines.append(event.xdata)
            else:
                pass
            #print(self.ax1_lines, self.ax2_lines)
            self.plotlines()
    
    def calibrate(self, event):
        l1 = np.unique(np.array(self.ax1_lines))
        l2 = np.unique(np.array(self.ax2_lines))
        l1.sort()
        l2.sort()
        f = lambda x, a,b: a*x+b
        p, pcov = curve_fit(f, l2, l1)
        #print(p)
        #pdb.set_trace()
        self.xcal = f(self.x2, *p)
        self.ax1.plot(self.xcal, self.y2)
        self.calpax.plot(l2, l1, 'o')
        self.calpax.plot(l2, f(l2, *p), '--')
        self.fig.canvas.draw_idle()

    def __init__(self, path, fname):
        self.path = path
        self.fname = fname
        d = np.recfromtxt('data/S16428AC_025_trian_zf2.dpt', names=['w', 'i'], encoding='utf8')
        self.x1 = d.w
        self.y1 = d.i
        print('Reading ', os.path.join(self.path, self.fname), '...')
        print(os.path.join(self.path, self.fname))
        self.data = np.recfromtxt(os.path.join(self.path, self.fname), names=['y'], skip_header=0, encoding='utf8')
        self.y = self.data.y*(-1)+np.max(self.data.y)
        print(np.max(self.y), np.min(self.y))
        self.s = int(self.fname[-15:-11])
        self.e = int(self.fname[-10:-6])
        self.yoffset = 0
        print('Setting initial spectral range to: ', self.s, '--', self.e, '...')
        #
        #self.get_new_xy()
        self.y2 = (self.y-self.yoffset)/np.max(self.y-self.yoffset)
        self.x2 = np.arange(len(self.y2))
        #
        self.ax1_lines, self.ax2_lines = [], []#[2916.313432835821, 2916.9211087420044, 2917.6460554371006, 2918.712153518124, 2919.1385927505335, 2920.6417910447763, 2921.3347547974417, 2922.944562899787, 2905.7270788912588, 2906.6759061833695, 2907.336886993604, 2908.339019189766, 2909.2771855010665, 2910.567164179105, 2900.108742004265, 2896.985074626866, 2895.087420042645], [6973.705578609254, 7126.3022367882595, 7276.891044201751, 7529.880240656417, 7632.280629697592, 7979.638812131381, 8134.243321075899, 8465.538697385582, 4640.5829890828845, 4859.438722523826, 5004.003977640778, 5220.8518603162065, 5441.715444522662, 5740.885208584133, 3333.4721407337734, 2648.795029693763, 2241.2013242945777]
        self.fig = plt.figure('Spectrum calibration', figsize=(12,7))
        self.ax1, self.ax2 = self.fig.add_subplot(211), self.fig.add_subplot(212)
        self.fig.subplots_adjust(left=0.3, bottom=0.2, right=0.97)
        self.ax1.set_xlim(self.s, self.e)
        self.ax1.set_ylim(0,1)
        self.ax1.grid()
        self.ax2.grid()
        self.calpax = plt.axes([0.04, 0.5, 0.15, 0.2])
        self.plotlines()
        self.ax1l = self.ax1.vlines(self.ax1_lines, 0,1)
        self.ax2l = self.ax2.vlines(self.ax2_lines, 0,1)
        self.ax1.plot(self.x1[(self.x1>=2800) & (self.x1<3000)], (self.y1[(self.x1>=2800) & (self.x1<3000)]/np.max(self.y1[(self.x1>=2800) & (self.x1<3000)])))
        self.digiline, = self.ax2.plot(self.x2, (self.y2)/np.max(self.y2))
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        #self.axas = plt.axes([0.25, 0.1, 0.65, 0.03])
        #self.axbs = plt.axes([0.25, 0.15, 0.65, 0.03])
        #self.axzs = plt.axes([0.25, 0.90, 0.65, 0.03])
        #
        #self.draw_sliders()
        #self.resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        #self.button = Button(self.resetax, 'Reset', hovercolor='0.975', color='gray')
        #self.button.on_clicked(self.reset)
        #self.printax = plt.axes([0.6, 0.025, 0.1, 0.04])
        #self.printbutton = Button(self.printax, 'Print', hovercolor='0.975', color='gray')
        #self.printbutton.on_clicked(self.print_params)
        #self.saveax = plt.axes([0.4, 0.025, 0.1, 0.04])
        #self.savebutton = Button(self.saveax, 'Save', hovercolor='0.975', color='gray')
        #self.savebutton.on_clicked(self.write_params)
        #self.rax = plt.axes([0.025, 0.7, 0.15, 0.15])
        #self.radio = RadioButtons(self.rax, ['linear nu', 'linear lambda'], active=0)
        #self.radio.on_clicked(self.xshift)
        self.caliax = plt.axes([0.6, 0.025, 0.1, 0.04])
        self.buttonc = Button(self.caliax, 'Calibrate', hovercolor='0.975', color='gray')
        self.buttonc.on_clicked(self.calibrate)
        self.rcaliax = plt.axes([0.5, 0.025, 0.1, 0.04])
        self.buttonrc = Button(self.rcaliax, 'Read Cal.file', hovercolor='0.975', color='gray')
        self.buttonrc.on_clicked(self.read_cali_lines)
        self.reduceax = plt.axes([0.3, 0.025, 0.1, 0.04])
        self.buttonr = Button(self.reduceax, 'Reduce SPC', hovercolor='0.975', color='gray')
        self.buttonr.on_clicked(self.reduce_points)
        self.saveax = plt.axes([0.4, 0.025, 0.1, 0.04])
        self.buttons = Button(self.saveax, 'Save sfit4 SPC', hovercolor='0.975', color='gray')
        self.buttons.on_clicked(self.print_sfit_readable_spectrum)
        plt.show()

#print_sfit_readable_spectrum('fname', testy2r, (2858.7, 2926.7))


if __name__ == '__main__':
    if len(sys.argv)==2:
        fname = sys.argv[1]
        path = 'data/'
        Calibrator(path, fname)
    else:
        print('Add image filename as cmdline arg ....')




