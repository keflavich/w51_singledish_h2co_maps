@aoinit
.r reduce_session
.r reduce_map

flist = get_scanlist(0,136,obsdate='20120910',machine='eta',bsg='b2s1g1')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213coW_spectra_0910.fits',obsdate='20120910',machine='eta',line='h213coW',percentile=10,/doplot
flist = get_scanlist(0,136,obsdate='20120910',machine='eta',bsg='b1s1g0')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213co_spectra_0910.fits',obsdate='20120910',machine='eta',line='h213co',percentile=10,/doplot

flist = get_scanlist(0,134,obsdate='20120911',machine='eta',bsg='b2s1g1')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213coW_spectra_0911.fits',obsdate='20120911',machine='eta',line='h213coW',percentile=10,/doplot
flist = get_scanlist(0,134,obsdate='20120911',machine='eta',bsg='b1s1g0')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213co_spectra_0911.fits',obsdate='20120911',machine='eta',line='h213co',percentile=10,/doplot

flist = get_scanlist(0,134,obsdate='20120912',machine='eta',bsg='b2s1g1',projid='a2705')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213coW_spectra_0912.fits',obsdate='20120912',machine='eta',line='h213coW',percentile=10,/doplot
flist = get_scanlist(0,134,obsdate='20120912',machine='eta',bsg='b1s1g0',projid='a2705')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213co_spectra_0912.fits',obsdate='20120912',machine='eta',line='h213co',percentile=10,/doplot

flist = get_scanlist(3,145,obsdate='20120915',machine='eta',bsg='b2s1g1',projid='a2705')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213coW_spectra_0915.fits',obsdate='20120915',machine='eta',line='h213coW',percentile=10,/doplot
flist = get_scanlist(3,145,obsdate='20120915',machine='eta',bsg='b1s1g0',projid='a2705')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h213co_spectra_0915.fits',obsdate='20120915',machine='eta',line='h213co',percentile=10,/doplot
