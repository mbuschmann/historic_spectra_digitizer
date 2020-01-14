#!/usr/bin/python3
from __future__ import print_function, division
import matplotlib
matplotlib.use('Qt5Agg')
import os, pickle, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import imageio

class Digitizer():
    def digitize_jfj_spectrum(self, pic, rg_val=2.25, black_val=60.0, rg_val_min=1.5, bg_val=0.9, bg_val_min=0.8, bright_val=100):
        l, w = self.l, self.w
        W = np.arange(w)
        data = []
        datamin=[]
        datamax=[]
        for i in range(l):
            r = pic[i,:,0].astype(float)+0.01
            g = pic[i,:,1].astype(float)+0.01
            b = pic[i,:,2].astype(float)+0.01
            line_rg = r/g
            line_bg = b/g
            cond_rg = line_rg > rg_val
            cond_bg = line_bg > bg_val
            cond_rg_min = line_rg > rg_val_min
            cond_bg_min = line_bg > bg_val_min
            cond_bright = np.mean(pic[i,:,:], axis=1) > bright_val
            if len(W[(cond_bg & cond_rg) & cond_bright])>5:
                data.append(np.mean(W[line_rg > rg_val]))
                datamin.append(np.min(W[line_rg > rg_val]))
                datamax.append(np.max(W[line_rg > rg_val]))
            elif len(W[(cond_bg_min & cond_rg_min) & cond_bright])>5:
                data.append(np.mean(W[line_rg > rg_val_min]))
                datamin.append(np.min(W[line_rg > rg_val_min]))
                datamax.append(np.max(W[line_rg > rg_val_min]))
            else:
                data.append(-1.0)
                datamin.append(-1.0)
                datamax.append(-0.5)
        self.data = {'data':np.array(data), 'datamin': np.array(datamin), 'datamax': np.array(datamax)}

    def digitize_spectrum(self, event):
        print('.. digitizing ...')
        self.digitized_spectrum = []
        for i in range(self.img.shape[1]):
            self.digitized_spectrum.append(self.digitize_point(i))
        self.ax2.plot(self.digitized_spectrum, 'b-')
        self.fig.canvas.draw_idle()

    def digitize_point(self, s):
        x = np.arange(self.img.shape[0])
        #if s<9000:
        #    cond1 = (x<4000) & (x>1400)
        #else:
        #    cond1 = (x<4700) & (x>1400)
        y = self.img[:,s, 0]/(self.img[:,s, 1]+0.01)
        #cond = y[cond1]>(np.max(y[cond1])*0.80)
        cond = y>(np.max(y)*0.80)
        #return np.mean(x[cond1][cond])
        return np.mean(x[cond])
        
    def save_spectrum(self, event):
        if self.digitized_spectrum:
            with open(os.path.join(self.savepath, 'digitized_'+self.fname[:-3]+'dat'), 'w') as f:
                for i in self.digitized_spectrum:
                    f.write('%4.2f\n'%i)
        print('digitized spectrum saved as ', os.path.join(self.savepath, 'digitized_'+self.fname[:-3]+'dat'))

    def plot_figure(self, init):
        if init:
            self.fig, (self.ax, self.ax2) = plt.subplots(ncols=2, nrows=1, figsize=(15,6), gridspec_kw={'width_ratios':(1,3)})
            self.randomizeax = plt.axes([0.8, 0.025, 0.1, 0.04])
            self.buttonr = Button(self.randomizeax, 'Randomize', hovercolor='0.975', color='gray')
            self.buttonr.on_clicked(self.randomize)
            self.digitizeax = plt.axes([0.7, 0.025, 0.1, 0.04])
            self.buttond = Button(self.digitizeax, 'Digitize', hovercolor='0.975', color='gray')
            self.buttond.on_clicked(self.digitize_spectrum)
            self.saveax = plt.axes([0.6, 0.025, 0.1, 0.04])
            self.buttons = Button(self.saveax, 'Save Spectrum', hovercolor='0.975', color='gray')
            self.buttons.on_clicked(self.save_spectrum)
        else:
            pass
        self.ax.grid()
        self.ax2.imshow(self.img)
        for s in range(self.nsel):
            self.ax2.plot(2*[self.sel[s]], [0,self.img.shape[0]], self.clrs[s]+'-', linewidth=0.5, label=str(self.sel[s]))
            y = self.img[:,self.sel[s], 0]/self.img[:,self.sel[s], 1]
            y2 = self.img[:,self.sel[s], 2]/self.img[:,self.sel[s], 1] 
            #y3 = (self.img[:,self.sel[s], 2]+self.img[:,self.sel[s], 1]+self.img[:,self.sel[s], 0])/3/256
            fwhm = np.mean(y[y>np.max(y)*0.80])
            self.ax.plot(y, self.clrs[s]+'-', label=str(s), linewidth=0.5)
            self.ax.plot(y2, self.clrs[s]+'--', label=str(s), linewidth=0.5, alpha=0.5)
            self.ax.plot([0,self.img.shape[0]], [fwhm, fwhm], self.clrs[s]+'--', label=str(self.sel[s]), linewidth=0.5, alpha=0.7)
        self.ax2.legend(loc='upper left', ncol=5)

    def randomize(self, event):
        print('.. randomizing ...')
        self.sel = [int(i) for i in np.random.random(self.nsel)*self.img.shape[1]]#+[4680]
        self.ax2.clear()
        self.ax.clear()
        self.plot_figure(init=False)
        self.fig.canvas.draw_idle()
    
    def __init__(self, path, fname, savepath, nsel=5):
        self.path = path
        self.fname = fname
        self.savepath = savepath
        self.img = imageio.imread(os.path.join(self.path, self.fname))
        self.nsel = nsel
        self.clrs = ['r', 'g', 'b', 'k', 'm']
        self.l, self.w = self.img.shape[0], self.img.shape[1]
        self.sel = [int(i) for i in np.random.random(nsel)*self.img.shape[1]]
        self.s = 0
        #
        self.plot_figure(init=True)
        self.digitize_spectrum(None)
        #
        plt.show()

if __name__ == '__main__':
    if len(sys.argv)==2:
        fname = sys.argv[1]
        path = 'images/'
        savepath = 'data/'
        Digitizer(path, fname, savepath)
    else:
        print('Add image filename as cmdline arg ....')





