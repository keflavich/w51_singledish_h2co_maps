.r reduce_session
.r reduce_map
;82-126
; SKIP - first round was corrupted, so this is pointless flist = get_scanlist(0,134,obsdate='20120915',machine='eta',projid='a2705')
; SKIP - first round was corrupted, so this is pointless accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_spectra_0915.fits',velocities=[-10,100],obsdate='20120915',machine='eta',projid='a2705'
;127-171
;flist = get_scanlist(127,44,obsdate='20120910',machine='eta')
;accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_spectra_0910.fits',velocities=[0,20]

; file read errors?
flist = get_scanlist(3,145,obsdate='20120915',machine='eta',projid='a2705',bsg='b5s1g1')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h108a_spectra_0915.fits',velocities=[-30,-10,100,130],obsdate='20120915',machine='eta',projid='a2705',line='h108a',percentile=10,/doplot,/do_mask_line

flist = get_scanlist(0,145,obsdate='20120915',machine='eta',bsg='b3s1g0',projid='a2705')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0915.fits',velocities=[-30,-5,10,37,72,115],obsdate='20120915',machine='eta',line='h2coW',percentile=10,/doplot,/do_mask_line

flist = get_scanlist(0,145,obsdate='20120915',machine='eta',bsg='b0s1g0',projid='a2705')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h2co_spectra_0915.fits',velocities=[-30,-5,10,37,72,115],obsdate='20120915',machine='eta',line='h2co',percentile=10,/doplot,/do_mask_line