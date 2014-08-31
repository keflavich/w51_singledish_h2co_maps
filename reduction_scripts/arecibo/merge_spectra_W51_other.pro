; created Thurs, Sep 26, 2013
; see arecibo_bsg_freqref.txt

readcol,'arecibo_bsg_freqref.txt',name,minfr,maxfr,bw,line,restfrq,format='(A,F,F,F,A,F)'

obsdates = ['20120910','20120911','20120912','20120915']
; redo these with fixed projids
;obsdates = ['20120912','20120915']
scanranges = [[0,136],[0,134],[0,134],[3,145]]
;scanranges = [[0,134],[0,145]]


for ii=0,n_elements(name)-1 do begin 
    bsg = name[ii] 
    linename = line[ii]
    restfreq = restfrq[ii]

    inrange = (restfreq-(50/3e5*restfreq) gt minfr[ii]) and (restfreq+(150/3e5*restfreq) lt maxfr[ii])

    if restfreq ne 0 and inrange then begin
        for jj=0,n_elements(obsdates)-1 do begin
            firstscan = scanranges[0,jj]
            lastscan = scanranges[1,jj]
            obsdate = obsdates[jj]
            projid = 'a2584'
            if obsdate eq '20120912' or obsdate eq '20120915' then projid = 'a2705'

            ; don't let failures break the loop...
            ;catch,error_status
            ;IF Error_status NE 0 THEN BEGIN	; This statement begins the error handler.
            ;    PRINT, 'Error index: ', Error_status
            ;    PRINT, 'Error message:', !ERR_STRING
            ;    continue
            ;ENDIF

            print,obsdate,firstscan,lastscan,bsg,linename,restfreq
            flist = get_scanlist(firstscan,lastscan,obsdate=obsdate,machine='eta',bsg=bsg,projid=projid)
            prefix = '/Users/adam/observations/arecibo/'+strtrim(obsdate)
            accum_map,flist,savefile=prefix+'/W51_'+strtrim(linename)+'_spectra_'+strmid(obsdate,4,4)+'.fits',$
                velocities=[-30,30,100,130],obsdate=obsdate,machine='eta',line=linename,percentile=10,/do_mask_line,$
                restfreq=restfreq,domedsmooth=domedsmooth
        endfor
    endif
endfor

end



; flist = get_scanlist(99,4,obsdate='20120910',machine='eta',bsg='b3s1g0')
; accum_map,flist,savefile='/Users/adam/observations/arecibo/20120910/W51_h2coW_spectra_0910.fits',velocities=[-10,100],obsdate='20120910',machine='eta',line='h2coW',/debug,percentile=10

