from matplotlib import colors
import numpy as np
import pylab as pl

cdict = {'red':   ((0.0, 1.0,1.0,),
                   (0.4, 0.5,0.5,),
                   (0.5, 0.2,0.2,),
                   (0.6, 0.0,0.0,),
                   (1.0, 0.0,0.0,)),
         'green': ((0.0, 0.0,0.0,),
                   (0.3, 0.2,0.2,),
                   (0.5, 1.0,1.0,),
                   (0.7, 0.2,0.2,),
                   (1.0, 0.0,0.0,)),
         'blue':  ((0.0, 0.0,0.0,),
                   (0.4, 0.0,0.0,),
                   (0.5, 0.2,0.2,),
                   (0.6, 0.5,0.5,),
                   (1.0, 1.0,1.0,)),
         'alpha': ((0.0, 1.0,1.0,),
                   (0.4, 0.7,0.7,),
                   (0.5, 0.5,0.5,),
                   (0.6, 0.7,0.7,),
                   (1.0, 1.0,1.0,)),
        }

VeloCmap = colors.LinearSegmentedColormap('VeloCmap',cdict,256)
VeloCmap.set_bad('w')
VeloCmap.set_under('w')
VeloCmap.set_over('w')

if __name__ == "__main__":
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    pl.clf()
    ax = pl.gca()
    ax.imshow(gradient, aspect='auto', cmap=VeloCmap)
    pl.show()
