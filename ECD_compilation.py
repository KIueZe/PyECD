import math
import time
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

def gaussian_func(sigma, E, Ek, Rk):
    cons = float(2.4572717053474 * np.power(10, 2)) # 10, 38
    # cons = 1/22.96 * np.power(np.pi, -1/2) * np.power(10, 40)
    sig_re = float(1/sigma)
    # reciprocal of sigma (half of the width of the CD band at 1/e peak height)
    expower = float(-np.power((E-Ek) * sig_re, 2))

    epsilon = float(cons * sig_re * Ek * Rk * np.power(np.e, expower))
    return epsilon

def gauss(sigma, E, Ek, Rk):
    cons = float(2.4572717053474 * math.pow(10, -3)) # 10, 38
    # cons = 1/22.96 * np.power(np.pi, -1/2) * np.power(10, 40)
    sig_re = float(1/sigma)
    # reciprocal of sigma (half of the width of the CD band at 1/e peak height)
    expower = float(-math.pow((E-Ek) * sig_re, 2))

    epsilon = float(cons * sig_re * Ek * Rk * math.pow(math.e, expower))
    return epsilon


def nm_ev(nm):
    # convert ev to nm, or vice versa
    ev = float(1239.85/float(nm))
    return ev

def read_bil(cd_bil):
    # return a dict, key for Ek, value for Rk
    ER_dict = {}
    with open(cd_bil, 'r') as f:
        # f.readline()
        line = f.readlines()
        for i in line:
            readline = i.split()
            if readline:
                Ek = nm_ev(readline[0]) # convert nm to eV
                Rk = float(readline[1])
                ER_dict[Ek] = Rk
    return ER_dict

def data_of_draw():
    ER_dict = read_bil('test.cd.bil')
    data_range = np.linspace(190, 400, num=2101)
    ret = {}
    for nm in data_range:
        E = nm_ev(nm)
        g = 0
        for Ek in ER_dict:
            g += gauss(0.16, E, Ek, ER_dict[Ek])
        ret[nm] = g
    return ret  

def draw_fig():
    # import data
    g = data_of_draw()
    listx = []
    listy = []
    for x in g:
        listx.append(x)
        listy.append(g[x])

    # set fig size
    fig = plt.figure(figsize=(8,4.9), dpi=300)
    ax = fig.add_subplot(111)

    ################################################################
    # set curve
    plt.plot(listx, listy,':r', linewidth=2, label='Calcd ECD of ' + '1')
    # plt.plot(listx, [i * -1 for i in listy], '-b', linewidth=2, label='Exptl ECD')
    
    ################################################################

    # set zeroline
    plt.plot(range(190,400), [0 for j in range(190,400)],':k',linewidth=1)
    
    # set tick labelsize
    plt.tick_params(labelsize=10)
    ticklabels = ax.get_xticklabels() + ax.get_yticklabels()
    [ticklabel.set_fontname('Times New Roman') for ticklabel in ticklabels]

    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    # ax.xaxis.set_major_locator(MultipleLocator(40))

    # set x, y label title
    plt.xlabel(r'Wavelength (nm)', family='Times new roman', weight='bold', fontsize=16)
    plt.ylabel(r'$\Delta\epsilon$', family='Times new roman',  weight='bold',fontsize=16)
    # plt.title(r'ECD', fontsize=20)

    # set legend
    legendfont = {'family':'Times new roman',
                'weight': 'bold',   # 'ultralight', 'light', 'normal', 'bold', 'heavy', 'extra bold', 'black'
                'style': 'normal',  # supported values are 'normal', 'italic', 'oblique'
                'size':'12'
    }
    plt.legend(loc=0, prop=legendfont,frameon=False)

    # set frame (spine) visiblity
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # save figure
    plt.savefig('ECD.png')
    plt.show()



if __name__ == '__main__':
    draw_fig()
