

function masflist,scannum,fbase=fbase,bsg=bsg,obsdate=obsdate,machine=machine,projid=projid
    if ~keyword_set(obsdate) then obsdate='20101111'
    if machine eq 'macbook' then begin
        if ~keyword_set(fbase) then fbase = getenv('AODATAROOT')+'/'+obsdate+'/'
    endif else if machine eq 'eta' then begin
        if ~keyword_set(fbase) then fbase = getenv('AODATAROOT')+'/'+obsdate+'/'
    endif else if machine eq 'ao' then begin
        pdev = (long(strmid(bsg,1,1))+1)+(long(strmid(bsg,5,1))*7)    
        if ~keyword_set(fbase) then fbase = '/share/pdata'+string(pdev,format="(I0)")+'/pdev/'
    endif else message,'wrong machine '+machine
    if ~keyword_set(bsg) then bsg = 'b0s1g0'
    if ~keyword_set(projid) then projid = 'a2584'
    fprefix = projid+'.'+obsdate+'.'+bsg+'.'
    suffix = '00.fits'
    flist = [fbase+fprefix+string(scannum,format='(I03)')+suffix,$
        fbase+fprefix+string(scannum+1,format='(I03)')+suffix,$
        fbase+fprefix+string(scannum+2,format='(I03)')+suffix,$
        fbase+fprefix+string(scannum+3,format='(I03)')+suffix]
    return,flist
end

function tsysmodel,za,p=p,tsys_scale=tsys_scale
    ; use fit_tsys fit, which fit (scaled) tsys for all cal observations from 2006 to present
    ; ASSUMES TSYS(ZA<10) = 27.2385
    if n_elements(p) eq 0 then p=[26.9417,0.0508974,-0.0269482,0.00190561]
    if n_elements(tsys_scale) eq 0 then tsys_scale=1

    model = p[0] + za*p[1] + za^2*(za gt 14)*p[2] + za^3*(za gt 14)*p[3]
    return,model*tsys_scale
end

pro wcs_line,header,line=line,crpix=crpix,cdelt=cdelt,crvalL=crvalL,vel_lsr=vel_lsr,crvalT=crvalT,$
    restfreq=restfreq
    if keyword_set(restfreq) then begin
        if ((restfreq lt 4000) or (restfreq gt 6000)) then begin
            message,'Rest frequency is not within C-band.  Must be specified in MHz.  Was: '+strtrim(restfreq)
        endif
    endif

    if ~keyword_set(line) then begin
        if ~keyword_set(restfreq) then restfreq = 4829.6594
    endif else if line eq 'h2co' or line eq 'h2coW' then begin
        restfreq = 4829.6594
    endif else if line eq 'h213co' or line eq 'h213coW' then begin
        restfreq = 4593.09
    endif else if line eq 'hdco' then begin
        restfreq = 4.48908*1e3
    endif else if line eq 'ccs' then begin
        restfreq = 5.40262*1e3
    endif else if line eq 'nh2cho' then begin
        restfreq = 4.61856*1e3
    endif else if line eq 'ch2nh' then begin
        restfreq = 5.29012*1e3
    endif else if line eq 'h108a' then begin
        restfreq = 5.14870*1e3
    endif else if line eq 'h135b' then begin
        restfreq = 5.22912*1e3
    endif else if line eq 'h2cn' or line eq 'h2cnW' then begin
        restfreq = 4.77777*1e3
    endif else if line eq 'OH' then begin
        restfreq = 4.66024*1e3
    endif else if line eq 'ch3oh' then begin
        restfreq = 5.00533*1e3
    endif else if line eq 'hcooh' then begin
        restfreq = 4.91632*1e3
    endif

    speedoflight=2.99792458e8

    ;x = masfreq(header)
    velo_topo = masfreq(header,/retvel,restFreq=restfreq,velcrdsys='T') ;(x-restfreq)/restfreq * speedoflight
    velo_bary = masfreq(header,/retvel,restfreq=restfreq,velcrdsys='B') ;velo_topo - header.vel_bary
    ra  = header.crval2
    dec = header.crval3
    vec = anglestovec3(ra*!dtor,dec*!dtor)
    offset_lsr=velLsrProj(vec,header.vel_bary/speedoflight)*speedoflight 
    ;velo_lsr1 = (x-restfreq)/restfreq * speedoflight+offset_lsr
    velo_lsr = masfreq(header,/retvel,restfreq=restfreq,velcrdsys='L') ;velo_topo - header.vel_bary
    ;plot,velo_lsr,velo_lsr1

    ;plot,velo_lsr/1e3,spec.d,psym=10
    vel_lsr = offset_lsr

    m = min(abs(velo_lsr),wherezero)
    
    cdelt = mean(velo_lsr[1:8191] - velo_lsr[0:8190])
    crvalL = velo_lsr[wherezero]
    crvalT = velo_topo[wherezero]
    crpix = wherezero+1

end

function add_lines_to_header,header,freqrange,spechead,linedict=linedict,machine=machine,projid=projid
    if machine eq 'macbook' or machine eq 'eta' then begin
        readcol,getenv('AODATAROOT')+'/arecibo_lines_named.txt',names,freq,tupp,format='(A,F,F)',/silent
    endif else if machine eq 'ao' then begin
        if ~keyword_set(projid) then projid='a2584'
        readcol,'/share/obs4/usr/'+projid+'/arecibo_lines_named.txt',names,freq,tupp,format='AFF',/silent
    endif else message,'wrong machine '+machine

    csuffix = ['A','B','C','D','E','G','H','I','J','K','L','M','N','O','P','Q','R','S','U','W','X','Y','Z']

    jj = 0
    ;linedict = [{''}]
    for ii=0,n_elements(names)-1 do begin 
        ;linedict = [linedict,{name:names[ii],freq:freq[ii]}]
        if (freqrange[0] lt freq[ii]*1e3) and (freqrange[1] gt freq[ii]*1e3) then begin
            wcs_line,spechead,crpix=crpix,cdelt=cdelt,crvalL=crvalL,vel_lsr=vel_lsr,crvalT=crvalT,restfreq=freq[ii]*1e3
            header = [header,'CRPIX1'+csuffix[jj]+' =  ' + string(strcompress(crpix),format='(A65)')]
            header = [header,'CRVAL1'+csuffix[jj]+' =  ' + string(strcompress(CRVALL),format='(A65)')]
            header = [header,'CDELT1'+csuffix[jj]+' =  ' + string(strcompress(CDELT),format='(A65)')]
            header = [header,'CTYPE1'+csuffix[jj]+' =  ' + string("'VRAD-LSR'",format='(A65)')]
            header = [header,'CUNIT1'+csuffix[jj]+' =  ' + string("'km/s'",format='(A65)')]
            header = [header,'RESTFRQ'+csuffix[jj]+'=  ' + string(strcompress(freq[ii]*1e9),format='(A65)')]
            header = [header,'LINENM'+csuffix[jj]+' =  ' + string("'"+names[ii]+"'",format='(A65)')]
            header = [header,string(names[ii],format='(A-8)')+"=  '" + string(csuffix[jj],format='(A)')+"'"]
            jj += 1
            if jj eq n_elements(csuffix) then begin
                print,"Can't add any more WCS to header!"
                ii=n_elements(names)
                return,header
            endif
            ;print,"Added line ",names[ii]," at frequency",freq[ii]," GHz to header"
        endif
    endfor

    return,header
end

pro reduce_session,line=line,output_prefix=output_prefix,offsmooth=offsmooth,obsdate=obsdate,$
    dostop=dostop,scanstart=scanstart,nscans=nscans,machine=machine,projid=projid

    if ~keyword_set(obsdate) then obsdate='20101111'
    if ~keyword_set(line) then line = 'h2co'
    if ~keyword_set(machine) then machine = 'macbook'
    if ~keyword_set(projid) then projid='a2584'
    if machine eq 'macbook' then begin
        if ~keyword_set(output_prefix) then output_prefix = getenv('AODATAROOT')+'/'+obsdate+'/'
    endif else if machine eq 'eta' then begin
        if ~keyword_set(output_prefix) then output_prefix = getenv('AODATAROOT')+'/'+obsdate+'/'
    endif else if machine eq 'ao' then begin
        ;if ~keyword_set(output_prefix) then output_prefix = '/home/aginsbur/nov'+strmid(obsdate,6,2)+'/'
        if ~keyword_set(output_prefix) then output_prefix = '/share/'+projid+'dat/'+obsdate+'/'
    endif else message,'wrong machine '+machine
    if ~keyword_set(offsmooth) then offsmooth=-64

    if line eq 'h2co' then begin
        bsg = 'b0s1g0'
    endif else if line eq 'h213co' then begin
        bsg = 'b1s1g0'
    endif else if line eq 'nh2cho' then begin
        bsg = 'b1s1g1'
    endif else if line eq 'ccs' then begin
        bsg = 'b6s1g0'
    endif else if line eq 'ch2nh' then begin
        bsg = 'b6s1g1'
    endif else if line eq 'h108a' then begin
        bsg = 'b5s1g1'
    endif else if line eq 'h135b' then begin
        bsg = 'b5s1g0'
    endif else if line eq 'h2coW' then begin
        bsg = 'b3s1g0'
    endif else if line eq 'h2cn' then begin
        bsg = 'b0s1g1'
    endif else if line eq 'h2cnW' then begin
        bsg = 'b3s1g1'
    endif else if line eq 'h213coW' then begin
        bsg = 'b2s1g1'
    endif else if line eq 'OH' then begin
        bsg = 'b2s1g0'
    endif else if line eq 'b4s1g1' then begin
        bsg = 'b4s1g1'
    endif else if line eq 'b4s1g0' then begin
        bsg = 'b4s1g0'
    endif else if line eq 'hcooh' then begin
        bsg = 'b4s1g1'
    endif else if line eq 'ch3oh' then begin
       bsg = 'b4s1g0'
    endif

    if n_elements(nscans) eq 0 then nscans = 24
    if n_elements(scanstart) eq 0 then scanstart = 0
    scans = indgen(nscans)*4+scanstart

    print,"Beginning reduction.  Obsdate=",obsdate," machine=",machine," line=",line," bsg=",bsg
    for ii=0,n_elements(scans)-1 do begin
        t0 = systime(/sec)
        scannum=scans[ii]

        flist = masflist(scannum,bsg=bsg,obsdate=obsdate,machine=machine,projid=projid)

        istat = masposonoff(flist,spec,specon,specoff,/sclcal,calI=calI,sclmask=1)
        istat = masposonoff(flist,specsmm,speconsmm,specoffsmm,/sclcal,calI=calIsm,smooff=offsmooth,sclmask=1)
        stop
        specOnDOff = masmath(specon,specoff,/div)
        speconrms=masrms(specon,/rembase)
        specoffrms=masrms(specoff,/rembase)
        specOnDOffrms=masrms(specondoff,/rembase)
        speconrmssmm=masrms(speconsmm,/rembase)
        specOnDOffsmm = masmath(speconsmm,specoffsmm,/div)
        specOnDOffsmmrms=masrms(specondoffsmm,/rembase)
        ;specoffrmssmm=masrms(specoffsmm,/rembase)
        speconmed = median(specon.d,dim=3)
        specoffmed = median(specoff.d,dim=3)
        ;speconsmmmed = median(speconsmm.d,dim=3)
        ;specoffsmmmed = median(specoffsmm.d,dim=3)

        speclength = n_elements(spec.d[*,0])
        nints = n_elements(specon)
        sublength = speclength / 16
        noisearr = dblarr(16,2)
        for pol=0,1 do begin
            for jj=0,15 do begin
                noisearr[jj,pol]=stddev(spec.d[sublength*jj:sublength*(jj+1)-1,pol]) 
            endfor 
        endfor
        medianrms = median(noisearr,dim=1)
        tsys = medianrms/(total(specondoffrms.d,1)/sqrt(nints)/double(speclength))
        spec.h.tsys = mean(tsys)

        tsysvector = replicate(1d,speclength) # tsys

        ; error propagation:
        ; T_A^* = (On - Off) / Off * C
        ; C = (CalOn - CalOff) / CalOff * Tsys * const
        ; C is assumed to have no error
        ;
        ; Assuming ZERO correlation between on & off
        ; sigma_T / T = (sigma_on^2-sigma_off^2)/(on-off)^2 + sigma_off^2/off^2
        ; (note that masposonoff has already scaled all values by C)
        ;
        ; alternately, T_A^* = on/off - 1
        ; sigma_T_A* = sqrt( (sigma_on/on)^2 * T_A^2 + (sigma_off/off)^2 * T_A^2 )
        errspec = specondoffrms.d/sqrt(nints)*tsysvector ;sqrt(speconrms.d^2/speconmed^2 + specoffrms.d^2/specoffmed^2)*abs(spec.d);*(replicate(1,speclength)#cali.calval)
        errspecavg = total(errspec,2)/2.0/sqrt(2d)
        errspecsm = specondoffsmmrms.d/sqrt(nints)*tsysvector ;sqrt(speconrmssmm.d^2/speconsmmmed^2 + specoffrmssmm.d^2/specoffsmmmed^2)*abs(specsmm.d);*(replicate(1,speclength)#calism.calval)
        errspecsmavg = total(errspecsm,2)/2.0/sqrt(2d)

        medianonrms = median(speconrms.d,dim=1)
        rmsonrms = mad(speconrms.d,dim=1) ;[stddev(speconrms.d[*,0]),stddev(speconrms.d[*,1])]
        if line eq 'ch2nh' then $
            badvals = where(speconrms.d gt ((medianonrms+3*rmsonrms)##replicate(1,speclength)),nbad) $
        else $
            badvals = where(speconrms.d gt ((medianonrms+5*rmsonrms)##replicate(1,speclength)),nbad)
        badarr = fltarr(speclength)
        if nbad gt 0 then begin
            print,"File ",bsg," line ",line,": Flagged ",nbad," bad points on rms"
            spec.d[badvals] = !values.f_nan
            badarr[badvals] = 1
        endif
        avgspec = total(spec.d,2)/2.0
        avgsmspec = total(specsmm.d,2)/2.0

        hor;,4828.5,4830.5
        ver
        masplot,spec,title=spec.h.object+line+" "+bsg
        print,"Scan number ",scannum," object ",spec.h.object," line ",line," bsg ",bsg
        ; there was a change to the Arecibo version of some component that required spec.h.restfrq to be set forcibly
        ; the change occured between 9/5/2011 and 10/10/2011
        wcs_line,spec.h,crpix=crpix,cdelt=cdelt,crvalL=crvalL,vel_lsr=vel_lsr,crvalT=crvalT,line=line,restfreq=restfreq
        spec.h.restfrq = restfreq * 1e6
        tagnames = tag_names(spec.h)
        fitsheader = strarr(n_elements(tagnames)) 
        for jj=0,n_elements(tagnames)-1 do begin
            szval = size(spec.h.(jj))
            if (szval[1] eq 1 or szval[1] eq 2) then begin
                printval = string(spec.h.(jj),format='(I66)')
            endif else if (szval[1] eq 4 or szval[1] eq 5) then begin
                printval = string(spec.h.(jj),format='(F66)')
            endif else begin
                printval = string("'"+strcompress(spec.h.(jj))+"'",format='(A66)')
            endelse
            fitsheader[jj] = string(tagnames[jj],format='(A8)') + "= " + printval
        endfor
        daycnv, spec.h.MJDXXOBS+2400000.5D, yr, mn, day, hr
        ; not needed daystart_str = long(string([yr-2000,mn-0,day],format='(I02,I02,I02)'))
        ; not needed dayend_str   = long(string([yr-2000,mn+1,day],format='(I02,I02,I02)'))
        ; not needed n=mmgetarchive(daystart_str,dayend_str,mm,rcvnum=9)
        ; best fit for all data from 1/1/2006 - 8/25/2011
        etamb_polyfit = poly((90-spec.h.elevatio),[0.491544,0.00580397,-0.000341992])
        fitsheader = [fitsheader,'DATEOBS =  ' + string([yr,mn,day,hr],format='(I04,"/",I02,"/",I02," ",F6.4)')]
        fitsheader = [fitsheader,'ETAMB   =  ' + string(etamb_polyfit,format='(F)')]
        exposuretime = total(specon.h.exposure)
        fitsheader = [fitsheader,'EXPOSURE=  ' + string(exposuretime)]
        fitsheader = [fitsheader,'TIMETOT =  ' + string(exposuretime)]
        radiff_seconds = (min(specoff.h.MJDXXOBS)-min(specon.h.MJDXXOBS))*86400
        waittime = (radiff_seconds-exposuretime)
        fitsheader = [fitsheader,'WAITTIME=  ' + string(waittime)]
        istatgain = masgainget(spec.h,gainval)
        fitsheader = [fitsheader,'GAINVAL =  ' + string(gainval)]
        fitsheader = [fitsheader,'CNTSTOKA=  ' + string(cali.cntstok[0])]
        fitsheader = [fitsheader,'CNTSTOKB=  ' + string(cali.cntstok[1])]
        fitsheader = [fitsheader,'RAOFFPOS=  ' + string(spec.h.crval2+radiff_seconds/3600.*15/cos(spec.h.crval3*!dpi/180.0))]
        tsys_fromfit = tsysmodel(90-spec.h.elevatio)
        fitsheader = [fitsheader,'TSYSFIT =  ' + string(strcompress(tsys_fromfit),format='(F)')]
        fitsheader = [fitsheader,'NINTS   =  ' + string(strcompress(nints),format='(A65)')]
        fitsheader = [fitsheader,'TSYSA   =  ' + string(strcompress(tsys[0]),format='(A65)')]
        fitsheader = [fitsheader,'TSYSB   =  ' + string(strcompress(tsys[1]),format='(A65)')]
        fitsheader = [fitsheader,'CRPIX1V =  ' + string(strcompress(crpix),format='(A65)')]
        fitsheader = [fitsheader,'CRVAL1V =  ' + string(strcompress(CRVALL),format='(A65)')]
        fitsheader = [fitsheader,'CDELT1V =  ' + string(strcompress(CDELT),format='(A65)')]
        fitsheader = [fitsheader,'CTYPE1V =  ' + string("'VRAD-LSR'",format='(A65)')]
        fitsheader = [fitsheader,'CUNIT1V =  ' + string("'km/s'",format='(A65)')]
        fitsheader = [fitsheader,'RESTFRQV=  ' + string(strcompress(spec.h.restfrq),format='(A65)')]
        fitsheader = [fitsheader,'VEL_LSR =  ' + string(vel_lsr,format='(F65)')]
        fitsheader = [fitsheader,'CRPIX1T =  ' + string(strcompress(crpix),format='(A65)')]
        fitsheader = [fitsheader,'CRVAL1T =  ' + string(strcompress(CRVALT),format='(A65)')]
        fitsheader = [fitsheader,'CDELT1T =  ' + string(strcompress(CDELT),format='(A65)')]
        fitsheader = [fitsheader,'CTYPE1T =  ' + string("'VRAD-TOP'",format='(A65)')]
        fitsheader = [fitsheader,'CUNIT1T =  ' + string("'km/s'",format='(A65)')]
        fitsheader = [fitsheader,'CTYPE1  =  ' + string("'FREQ'",format='(A65)')]
        fitsheader = [fitsheader,'CUNIT1  =  ' + string("'Hz'",format='(A65)')]
        fitsheader = [fitsheader,'BUNIT   =  ' + string("'K'",format='(A65)')]
        fitsheader = [fitsheader,'LINE1   =  ' + string("'Average Spectrum'",format='(A65)')]
        fitsheader = [fitsheader,'LINE2   =  ' + string("'Error on Avg Spectrum'",format='(A65)')]
        fitsheader = [fitsheader,'LINE3   =  ' + string("'Pol A average'",format='(A65)')]
        fitsheader = [fitsheader,'LINE4   =  ' + string("'Pol B average'",format='(A65)')]
        fitsheader = [fitsheader,'LINE5   =  ' + string("'Pol A error'",format='(A65)')]
        fitsheader = [fitsheader,'LINE6   =  ' + string("'Pol B error'",format='(A65)')]
        fitsheader = [fitsheader,'LINE7   =  ' + string("'Pol A on median'",format='(A65)')]
        fitsheader = [fitsheader,'LINE8   =  ' + string("'Pol B on median'",format='(A65)')]
        fitsheader = [fitsheader,'LINE9   =  ' + string("'Pol A on rms'",format='(A65)')]
        fitsheader = [fitsheader,'LINE10  =  ' + string("'Pol B on rms'",format='(A65)')]
        fitsheader = [fitsheader,'LINE11  =  ' + string("'Pol A off median'",format='(A65)')]
        fitsheader = [fitsheader,'LINE12  =  ' + string("'Pol B off median'",format='(A65)')]
        fitsheader = [fitsheader,'LINE13  =  ' + string("'Pol A off rms'",format='(A65)')]
        fitsheader = [fitsheader,'LINE14  =  ' + string("'Pol B off rms'",format='(A65)')]
        fitsheader = [fitsheader,'LINE15  =  ' + string("'BADFLAGS'",format='(A65)')]
        fitsheader = [fitsheader,'LINE16  =  ' + string("'Average spec, off smoothed by "+strcompress(offsmooth)+"'",format='(A65)')]
        fitsheader = [fitsheader,'LINE17  =  ' + string("'Error spec smooth off'",format='(A65)')]
        fitsheader = [fitsheader,'NBAD    =  ' + string(nbad,format='(A65)')]
        freq = masfreq(spec.h)
        fitsheader = add_lines_to_header(fitsheader,[min(freq),max(freq)],spec.h,machine=machine,projid=projid)
        fitsheader = [fitsheader,'END     ']
        outdata = [[avgspec],[errspecavg],[spec.d],[errspec],[speconmed],$
            [speconrms.d],[specoffmed],[specoffrms.d],[badarr],[avgsmspec],[errspecsmavg]]
        fits_write,output_prefix+spec.h.object+"_"+string(line,format='(A0)')+".fits",outdata,fitsheader
        t1 = systime(/sec)
        print,spec.h.object+line+" "+bsg+" Reduction took ",t1-t0," seconds"
        if keyword_set(dostop) then stop
    endfor

    end


pro reduce_session_wrapper,line=line,obsdate=obsdate,scanstart=scanstart,nscans=nscans,dostop=dostop,machine=machine,projid=projid
    if ~keyword_set(machine) then machine='eta'

    ;reduce_session
    ;read,obsdate
    print,line,obsdate,scanstart,nscans
    if total(obsdate eq ['20101110','20101111','20101112','20101113','20101114','20110222','20110223','20110224','20110830','20110901','20110902','20110903','20110904','20110905','20111010','20111011','20111012','20111013','20111014','20111122','20130912']) eq 1 then begin
        if line eq 'all' then begin
            reduce_session,line='h2co',obsdate=obsdate  ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h2coW',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h213co',obsdate=obsdate,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h213coW',obsdate=obsdate,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h2cn',obsdate=obsdate  ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h2cnW',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='hcooh',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='ch3oh',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='ch2nh',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='nh2cho',obsdate=obsdate,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='ccs',obsdate=obsdate   ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='OH',obsdate=obsdate    ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h108a',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            reduce_session,line='h135b',obsdate=obsdate ,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
            print,"Done"
        endif else reduce_session,line=line,obsdate=obsdate,dostop=dostop,machine=machine,scanstart=scanstart,nscans=nscans,projid=projid
    endif else print,"Compiled, but obsdate was invalid: ",obsdate
end
