import pyspeckit
import paths

cube1 = pyspeckit.Cube(paths.dpath('W51_H2CO11_taucube_supersampled.fits'))
cube1.xarr.refX_units='GHz'
cube1.xarr.refX = 4.829659400
E = cube1.cube[cube1.xarr.as_unit('km/s') < 0].std(axis=0)
cube1.errorcube = np.repeat(np.reshape(E,(1,)+E.shape),cube1.shape[0],axis=0)

parcubefilename = paths.dpath('H2CO11_h2cofit_parameters.fits')
cube1.fiteach(guesses=[1,60,1,1,70,1],
        absorption=False,
        integral=False,
        fittype='formaldehyde',
        multicore=1,
        signal_cut=3,
        parlimited=[(True,False), (True,True), (True,False)]*2, 
        parlimits=[(0,0), (40,65), (0,0)]+[(0,0), (65,75), (0,0)],
        start_from_point=[95,85])
cube1.write_fit(parcubefilename,clobber=True)
