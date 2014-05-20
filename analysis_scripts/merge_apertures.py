import pyregion

apregfiles = ['dense_filament_spectral_apertures.reg',
              'maus_spectral_apertures.reg',
              'filament_leftside_spectral_apertures.reg',
              'middlechunk_spectral_apertures.reg',
              'w51main_spectral_apertures.reg', ]

regions = pyregion.ShapeList(reduce(pyregion.ShapeList.__add__, [pyregion.open(r) for r in apregfiles]))

regions.write('merged_spectral_apertures.reg')
