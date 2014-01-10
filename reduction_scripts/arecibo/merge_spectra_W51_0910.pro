.r reduce_session
.r reduce_map
;82-126
; CORRUPTED flist = get_scanlist(0,136,obsdate='20120910',machine='eta')
; CORRUPTED accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_spectra_0910.fits',velocities=[-10,100],obsdate='20120910',machine='eta'
;accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_spectra_0910.fits',velocities=[0,0],obsdate='20120910',machine='eta'
;127-171
;flist = get_scanlist(127,44,obsdate='20120910',machine='eta')
;accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_spectra_0910.fits',velocities=[0,20]

; there is RFI in this band
flist = get_scanlist(0,136,obsdate='20120910',machine='eta',bsg='b5s1g1')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h108a_spectra_0910.fits',velocities=[-30,-10,100,130],obsdate='20120910',machine='eta',line='h108a',percentile=10,/doplot,/do_mask_line

flist = get_scanlist(0,136,obsdate='20120910',machine='eta',bsg='b3s1g0')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0910.fits',velocities=[-10,-2,10,45,72,115],obsdate='20120910',machine='eta',line='h2coW',percentile=10,/doplot,/do_mask_line

; this one was tuned incorrectly
flist = get_scanlist(0,136,obsdate='20120910',machine='eta',bsg='b0s1g0')
accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h2co_spectra_0910.fits',velocities=[-10,-2,10,45,72,115],obsdate='20120910',machine='eta',line='h2co',percentile=10,/doplot,/do_mask_line


; flist = get_scanlist(99,4,obsdate='20120910',machine='eta',bsg='b3s1g0')
; accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0910.fits',velocities=[-10,100],obsdate='20120910',machine='eta',line='h2coW',/debug,percentile=10
