;function map_flist,scanstart,scanend,bsg=bsg,obsdate=obsdate,machine=machine

function get_scanlist,scanstart,nscans,fbase=fbase,bsg=bsg,obsdate=obsdate,machine=machine,projid=projid
    if ~keyword_set(obsdate) then obsdate='20101111'
    if machine eq 'macbook' or machine eq 'eta' then begin
        if ~keyword_set(fbase) then fbase = getenv('AODATAROOT')+'/'+obsdate+'/'
    endif else if machine eq 'ao' then begin
        pdev = (long(strmid(bsg,1,1))+1)+(long(strmid(bsg,5,1))*7)    
        if ~keyword_set(fbase) then fbase = '/share/pdata'+string(pdev,format="(I0)")+'/pdev/'
    endif else message,'wrong machine '+machine
    if ~keyword_set(bsg) then bsg = 'b0s1g0'
    if ~keyword_set(projid) then projid='a2584'
    fprefix = projid+'.'+obsdate+'.'+bsg+'.'
    suffix = '00.fits'
    flist = [""] ; initialize as string array (could use list for IDL 8...)
    for jj=scanstart,nscans+scanstart-1 do begin
        ;if (jj-scanstart) mod 3 eq 2 then begin
            flist = [flist,fbase+fprefix+string(jj,format='(I03)')+suffix]
        ;endif
    endfor
    return,flist[1:*]
end

function mask_line,velo,spec,velocities,fitorder,velocity_buffer=velocity_buffer,$
         pp=pp
    ; velocities should be specified as a set of ranges to fit
    ; in between will be interpolated over
    ; e.g.:
    ; [0,10,20,30] will fit 0-10,20-30 and will replace 10-20 with the fit

    ; velocity_buffer gives the range around the masked region to include in the polyfit
    ;if n_elements(velocity_buffer) eq 0 then velocity_buffer = 20 

    fit_mask = replicate(0,n_elements(velo))
    replace_mask = replicate(0,n_elements(velo))
    ; n_e(v)-2 because elements 0,-1 are the endpoints
    for ii=0,n_elements(velocities)-2 do begin
        vmin = velocities[ii]
        vmax = velocities[ii+1]
        wh = where((velo gt vmin)*(velo lt vmax),nw)
        if nw gt 0 then begin
            if ii mod 2 eq 0 then fit_mask[wh] = 1 $
                else replace_mask[wh] = 1
        endif else message,'No matches found to velocity range '+string(vmin,vmax,format="(F,'-',F)")
    endfor
    fit_region = where(fit_mask,nfit)
    replace_region = where(replace_mask,nreplace)
    if nfit eq 0 or nreplace eq 0 then message,'Did not find enough points to fit / replace'

    pp = poly_fit(velo[fit_region],spec[fit_region],fitorder)
    spec[replace_region] = poly(velo[replace_region],pp)

    ; mask out selected velocities (e.g., where line is)
    ; then, replace with best-fit across that region
    ;mask_region = where((velo lt max(velocities))*(velo gt min(velocities)),nmask,complement=not_mask)
    ;if nmask gt 0 then begin
    ;    fitvelocities = [min(velocities)-velocity_buffer,max(velocities)+velocity_buffer]
    ;    fit_region = where((velo lt max(fitvelocities))*(velo gt max(velocities))+$
    ;        (velo gt min(fitvelocities))*(velo lt min(velocities)),nfit,complement=not_fit)
    ;    pp = poly_fit(velo[fit_region],spec[fit_region],fitorder)
    ;    spec[mask_region] = poly(velo[mask_region],pp)
    ;endif

    return,spec
end

; plot things:
; for ii=0,5 do begin & i1 = masopen(flist[ii],desc) & i2 = masavg(desc, 2, avspec) & masclose,desc & masplot,avspec,pollist=[0],/over,/norm & endfor

function make_off,flist,exclude_middle=exclude_middle,percentile=percentile,doplot=doplot,$
    continuum_fitpars=continuum_fitpars, fitorder_off=fitorder_off,scanmeans=scanmeans,dostop=dostop,$
    speccube=speccube, scaled_speccube=speccube_scaled, savefile=savefile, do_mask_line=do_mask_line, $
    fitorder_line=fitorder_line, velocities=velocities, velocity_buffer=velocity_buffer, $
    restfreq=restfreq,hack_obsmode=hack_obsmode,domedsmooth=domedsmooth
    ; continuum_fitpars : output array with polyfit parameters for each pol
    ; scanmeans : output array of accumulated spectrum means

    if keyword_set(doplot) then begin
        window,0
        window,1
    endif

    if n_elements(exclude_middle) eq 0 then exclude_middle = 1 ; try to exclude the central "spike" from the Mock
    if n_elements(percentile) eq 0 then percentile=10
    if n_elements(fitorder_off) eq 0 then fitorder_off = 8
    if n_elements(fitorder_line) eq 0 then fitorder_line = 6
    speclist = []
    speccube = []
    accumtotals = []
    inds_by_scan1 = []
    inds_by_scan2 = []
    za = []
    clearline = string(13b) ;"
    scannum = 0
    print,format='("Creating OFF using ",I," files")',n_elements(flist)
    for jj=0,n_elements(flist)-1 do begin
        ;print,"Loading file ",jj,scannum,clearline,format="($,A,I4,' scan ',I4,A)"
        file = flist[jj]
        istat1 = masopen(file,desc)
        istat2 = masgetfile(desc,spec)
        istat3 = masaccum(spec,accumspec,/new,/avg)
        istat4 = mascalval(accumspec.h,calI)
        print,"Loaded file ",jj,scannum,accumspec.h.obsmode,clearline,format="($,A,I4,' scan ',I4,' obsmode: ',A,A)"
        if accumspec.h.obsmode ne 'CAL' or $
                (keyword_set(hack_obsmode) and (jj mod 4 eq 2 or jj mod 4 eq 3)) then begin

            ; if both mock boxes are used, maybe?
            if size(spec.d,/n_dim) gt 3 then begin
                ; average over 3rd dimension of 4
                ; total is 1-indexed, IDL is otherwise 0-indexed
                specdata = total(spec.d,3)/(size(spec.d,/dim))[2]
            endif else specdata = spec.d

            ; check for bad data
            if size(accumtotals,/n_dim) eq 3 then begin
                ; skip empty speccube
                if (size(accumtotals,/dim))[2] ne (size(specdata,/dim))[2] then begin
                    print,"Failed to load scan",jj,"with shape",size(specdata)
                    ;stop
                    continue
                endif
            endif
            
            speccube = [[[[speccube]]],[[specdata]]]
            speclist = [speclist,spec]
            velo = masfreq(accumspec.h,velcrdsys='L',/retvel,restfreq=restfreq)
            za = [za,90-accumspec.h.elevatio]
            ; number of spectral elements
            specsize = (size(specdata,/dim))[0]
            ; totals are summed across spectral dimension, leaving npols x nscans
            if exclude_middle then begin
                totals = total(specdata[0.2*specsize:4090,*,*],1) + $
                         total(specdata[4102:0.8*specsize,*,*],1)
                npix = 8192*0.6-12
            endif else begin
                totals = total(specdata[0.1*specsize:0.9*specsize,*,*],1)
                npix = 8192*0.8
            endelse
            ;print,"Totals size: ",size(totals,/dim)
            accumtotals = [[accumtotals],[totals]]
            ; compute best-fit line with time
            npts = (size(totals,/dim))[1]
            xinds = findgen(npts)
            ;linpars1 = linfit(xinds,totals[0,*],yfit=linfit1)
            ;linpars2 = linfit(xinds,totals[1,*],yfit=linfit2)
            ; sort each pol
            sorted1 = sort(totals[0,*]);-linfit1)
            sorted2 = sort(totals[1,*]);-linfit2)
            ; nsorted = nscans
            nsorted = npts
            npercentile = floor(nsorted*percentile/100.) ; number above percentile
            firstN1 = sorted1[0:npercentile]
            firstN2 = sorted2[0:npercentile]
            inds_by_scan1 = [[inds_by_scan1],[firstN1 + npts * scannum]]
            inds_by_scan2 = [[inds_by_scan2],[firstN2 + npts * scannum]]
            scannum += 1
        endif
        masclose,desc
    endfor
    scanmeans = accumtotals / float(npix)
    
    if size(scanmeans,/n_dim) gt 2 then message,'Error: scanmeans have the wrong size and shape.'

    nscans_thr = nsorted*scannum
    nscans = (size(speccube,/dim))[2]
    if nscans ne nscans_thr then message,"Warning: nscans="+strtrim(nscans)+" but should be "+strtrim(nscans_thr),/info
    xinds = findgen(nscans)

    ; polyfit to get "background level"
    ; (didn't work so well)
    pfit1 = poly_fit(inds_by_scan1,scanmeans[0,inds_by_scan1],fitorder_off,yfit=yfit1) 
    pfit2 = poly_fit(inds_by_scan2,scanmeans[1,inds_by_scan2],fitorder_off,yfit=yfit2) 
    pol1 = poly(xinds,pfit1)
    pol2 = poly(xinds,pfit2)

    ; 0'th order polynomial - flat baseline instead of weird curvy junk
    pfit1 = [min(scanmeans[0,*])] 
    pfit2 = [min(scanmeans[1,*])] 
    continuum_fitpars = [pfit1,pfit2] ; output parameter

    ; any integrations with strong continuum are scaled down to the fitted lower
    ; surface for later averaging
    ; scale_factors = transpose([[pol1/accumtotals[0,*]],[pol2/accumtotals[1,*]]])
    ; That's not right... they should just be scaled to all have the same mean
    scale_factors = ([1./scanmeans[0,*],1./scanmeans[1,*]])

    speccube_scaled = speccube*0.
    ; to save memory, do this in a for loop rather than broadcasting scale_factors
    for ii=0,specsize-1 do begin
        speccube_scaled[ii,0:1,0:nscans-1] = speccube[ii,0:1,0:nscans-1] * scale_factors
    endfor
    avg_off_A = total(speccube_scaled,3) / nscans

    if keyword_set(do_mask_line) then begin
        pol1 = avg_off_A[*,0]
        pol2 = avg_off_A[*,1]
        if keyword_set(domedsmooth) then begin
            pol1 = medsmooth(pol1, 5)
            pol2 = medsmooth(pol2, 5)
        endif
        off0 = mask_line(velo, pol1, velocities, fitorder_line, velocity_buffer=velocity_buffer,pp=pp1)
        off1 = mask_line(velo, pol2, velocities, fitorder_line, velocity_buffer=velocity_buffer,pp=pp2)
        avg_off = [[off0],[off1]]
    endif else begin
        avg_off = avg_off_A
    endelse

    if keyword_set(doplot) then begin
        wset,1
        hor,0,nsorted*scannum
        ver,min([pol1,pol2]),max([pol1,pol2]);6e8,1.4e9 ;min(accumtotals),max(accumtotals)
        plot,[0,nsorted],[min(scanmeans),max(scanmeans)],/nodata,xtitle="Integration #",ytitle="Total inner 80%"
        oplot,scanmeans[0,*],color='00FF00'x
        oplot,scanmeans[1,*],color='0000FF'x
        ;oplot,xinds,linpars1,color='000F0F'x
        ;oplot,xinds,linpars2,color='0F000F'x
        oplot,inds_by_scan1,scanmeans[0,inds_by_scan1],color='FFFF00'x,psym=3
        oplot,inds_by_scan2,scanmeans[1,inds_by_scan2],color='FF00FF'x,psym=3
        oplot,xinds,pol1
        oplot,xinds,pol2,color='FFFF00'x
        oplot,xinds,replicate(pfit1,nscans),color='00FFFF'x
        oplot,xinds,replicate(pfit2,nscans),color='CCCCFF'x
        outimg = strjoin(strsplit(flist[0],".fits",/regex,/extract,/preserve_null),"_continuum.png")
        write_png,outimg,tvrd(/true)

        wset,0
        hor
        ;ver,min(avg_off_A),max(avg_off_A)
        ver,0.95,1.05 ; normalized...
        plot,velo,avg_off[*,0],/nodata
        oplot,velo,avg_off_A[*,0],color='AA22AA'x
        oplot,velo,avg_off_A[*,1],color='FFFF00'x
        oplot,velo,avg_off[*,0],color='00FF00'x
        oplot,velo,avg_off[*,1],color='0000FF'x
        outimg = strjoin(strsplit(flist[0],".fits",/regex,/extract,/preserve_null),"_offspectra.png")
        write_png,outimg,tvrd(/true)
    endif

    if keyword_set(dostop) then stop

    if keyword_set(savefile) then begin
        outfile = strjoin(strsplit(savefile,".fits",/regex,/extract,/preserve_null),"_velo.fits")
        writefits,outfile,velo
        outfile = strjoin(strsplit(savefile,".fits",/regex,/extract,/preserve_null),"_offspectra.fits")
        writefits,outfile,avg_off
        if keyword_set(do_mask_line) then begin
            outfile = strjoin(strsplit(savefile,".fits",/regex,/extract,/preserve_null),"_offspectra_nomask.fits")
            writefits,outfile,avg_off_A
        endif
        outfile = strjoin(strsplit(savefile,".fits",/regex,/extract,/preserve_null),"_scanmeans.fits")
        writefits,outfile,scanmeans
    endif

    return,avg_off
end


; should change: flist should be scanstart, nscans and should call get_scanlist;
; looks like this was awkwardly hacked together...
pro accum_map,flist,savefile=savefile,line=line,output_prefix=output_prefix,offsmooth=offsmooth,obsdate=obsdate,$
    dostop=dostop,scanstart=scanstart,nscans=nscans,machine=machine,velocities=velocities,restfreq=restfreq,$
    fitorder=fitorder,projid=projid,debug=debug,flatten=flatten,percentile=percentile,exclude_middle=exclude_middle,$
    velocity_buffer=velocity_buffer,zoomregion=zoomregion,doplot=doplot,do_mask_line=do_mask_line,hack_obsmode=hack_obsmode,$
    domedsmooth=domedsmooth

    if ~keyword_set(obsdate) then obsdate='20101111'
    if ~keyword_set(line) then line = 'h2co'
    if ~keyword_set(machine) then machine = 'macbook'
    if machine eq 'macbook' or machine eq 'eta' then begin
        if ~keyword_set(output_prefix) then output_prefix = getenv('AODATAROOT')+'/'+obsdate+'/'
    endif else if machine eq 'ao' then begin
        ;if ~keyword_set(output_prefix) then output_prefix = '/home/aginsbur/nov'+strmid(obsdate,6,2)+'/'
        if ~keyword_set(output_prefix) then output_prefix = '/share/'+projid+'dat/'+obsdate+'/'
    endif else message,'wrong machine '+machine
    ; not used if ~keyword_set(offsmooth) then offsmooth=-64
    if ~keyword_set(fitorder) then fitorder=2
    if n_elements(exclude_middle) eq 0 then exclude_middle = 1 ; try to exclude the central "spike" from the Mock
    ; number of km/s above/below the specified velocities to include when fitting
    if n_elements(velocity_buffer) eq 0 then velocity_buffer = 20 
    ; frequency in MHz to zoom in on when plotting "final" reduced spec
    if n_elements(zoomregion) eq 0 then zoomregion = [4828,4831]
    if n_elements(debug) eq 0 then debug=0

    if line eq 'h2co' then begin
        bsg = 'b0s1g0'
        if ~keyword_set(restfreq) then restfreq = 4829.6594
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
        if ~keyword_set(restfreq) then restfreq = 5.14870*1e3
    endif else if line eq 'h135b' then begin
        bsg = 'b5s1g0'
    endif else if line eq 'h2coW' then begin
        bsg = 'b3s1g0'
        if ~keyword_set(restfreq) then restfreq = 4829.6594
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

    print,"Beginning reduction.  Obsdate=",obsdate," machine=",machine," line=",line

    ; why flist[2] ? is that when I am not cal?
    ; I think so: you do calON, calOFF, then the rest is scans
    file = flist[2]
    istat1 = masopen(file,desc)
    istat2 = masgetfile(desc,spec)
    istat3 = masgetfile(desc,avspec,/avg)
    ;istat3 = masavg(desc,1,avspec)
    masclose,desc

    wcs_line,avspec.h,crpix=crpix,cdelt=cdelt,crvalL=crvalL,vel_lsr=vel_lsr,crvalT=crvalT,line=line,restfreq=restfreq
    tagnames = tag_names(avspec.h)
    ;fxbhmake,fitsheader,n_elements(flist)*n_elements(spec),'data',/initialize,/date
    fxhmake,fitsheader,/extend,/initialize,/date
    fitsheader = fitsheader[0:4]
    for jj=0,n_elements(tagnames)-1 do begin
        if stregex(tagnames[jj],"TDIM",/bool) then continue
        szval = size(avspec.h.(jj))
        if (szval[1] eq 1 or szval[1] eq 2) then begin
            printval = string(avspec.h.(jj),format='(I66)')
        endif else if (szval[1] eq 4 or szval[1] eq 5) then begin
            printval = string(avspec.h.(jj),format='(F66)')
        endif else begin
            printval = string("'"+strcompress(avspec.h.(jj))+"'",format='(A66)')
        endelse
        fitsheader = [fitsheader,string(tagnames[jj],format='(A8)') + "= " + printval]
    endfor
    daycnv, spec.h.MJDXXOBS+2400000.5D, yr, mn, day, hr
    dateobs = []
    for qq=0,n_elements(yr)-1 do begin
        dateobs = [dateobs,string([yr[qq],mn[qq],day[qq],hr[qq]],format='(I04,"/",I02,"/",I02," ",F8.4)')]
    endfor
    etamb_polyfit = poly((90-spec.h.elevatio),[0.491544,0.00580397,-0.000341992])
    fitsheader = [fitsheader,'DATEOBS =  ' + dateobs]
    fitsheader = [fitsheader,'ETAMB   =  ' + string(etamb_polyfit,format='(F)')]
    exposuretime = total(spec.h.exposure)
    fitsheader = [fitsheader,'EXPOSURE=  ' + string(exposuretime)]
    fitsheader = [fitsheader,'TIMETOT =  ' + string(exposuretime)]
    istatgain = masgainget(spec.h,gainval)
    fitsheader = [fitsheader,'GAINVAL =  ' + string(gainval)]
    tsys_fromfit = tsysmodel(90-spec.h.elevatio)
    fitsheader = [fitsheader,'TSYSFIT =  ' + string(strcompress(tsys_fromfit),format='(F)')]
    fitsheader = [fitsheader,'CRPIX1V =  ' + string(strcompress(crpix),format='(A65)')]
    fitsheader = [fitsheader,'CRVAL1V =  ' + string(strcompress(CRVALL),format='(A65)')]
    fitsheader = [fitsheader,'CDELT1V =  ' + string(strcompress(CDELT),format='(A65)')]
    fitsheader = [fitsheader,'CTYPE1V =  ' + string("'VRAD-LSR'",format='(A65)')]
    fitsheader = [fitsheader,'CUNIT1V =  ' + string("'km/s'",format='(A65)')]
    fitsheader = [fitsheader,'VEL_LSR =  ' + string(vel_lsr,format='(F65)')]
    fitsheader = [fitsheader,'RESTFRQV=  ' + string(strcompress(spec.h.restfrq),format='(A65)')]
    fitsheader = [fitsheader,'CRPIX1T =  ' + string(strcompress(crpix),format='(A65)')]
    fitsheader = [fitsheader,'CRVAL1T =  ' + string(strcompress(CRVALT),format='(A65)')]
    fitsheader = [fitsheader,'CDELT1T =  ' + string(strcompress(CDELT),format='(A65)')]
    fitsheader = [fitsheader,'CTYPE1T =  ' + string("'VRAD-TOP'",format='(A65)')]
    fitsheader = [fitsheader,'CUNIT1T =  ' + string("'km/s'",format='(A65)')]
    fitsheader = [fitsheader,'CTYPE1  =  ' + string("'FREQ'",format='(A65)')]
    fitsheader = [fitsheader,'CUNIT1  =  ' + string("'Hz'",format='(A65)')]
    fitsheader = [fitsheader,'RESTFRQ =  ' + string(restfreq,format='(F65)')]
    fitsheader = [fitsheader,'BUNIT   =  ' + string("'K'",format='(A65)')]
    fitsheader = [fitsheader,'BTYPE   =  ' + string("'T_MB'",format='(A65)')]
    fitsheader = [fitsheader,'ZA      =  ' + string(90-spec.h.elevatio)]
    freq = masfreq(avspec.h,restfreq=restfreq)
    fitsheader = add_lines_to_header(fitsheader,[min(freq),max(freq)],avspec.h,machine=machine)
    fitsheader = [fitsheader,'END     ']
    fxwrite,savefile,fitsheader
    ;print,fitsheader[2:3]
    ;help,fitsheader

    nrows = n_elements(flist)/3*n_elements(spec)
    specarr = reform(fltarr(avspec.nchan[0])) ;fltarr(avspec.nchan[0],nrows)
    ra      = dblarr(1) ;fltarr(nrows)
    dec     = dblarr(1) ;fltarr(nrows)
    glon    = dblarr(1) ;fltarr(nrows)
    glat    = dblarr(1) ;fltarr(nrows)
    
    ; generate the observation-averaged 'off' position
    offspec_avg = make_off( flist, exclude_middle=exclude_middle, percentile=percentile, $
        doplot=doplot, continuum_fitpars=continuum_fitpars,  fitorder_off=fitorder_off, $
        scanmeans=scanmeans, dostop=dostop, savefile=savefile, do_mask_line=do_mask_line, $
        velocities=velocities, velocity_buffer=velocity_buffer, $
        restfreq=restfreq,hack_obsmode=hack_obsmode,domedsmooth=domedsmooth)

    help,offspec_avg
    help,specarr

    fxbhmake,fitsheader,nrows,'AreciboMap','Individual spectra for mapping'
    fxbaddcol,1,fitsheader,specarr,'SPECTRA'
    fxbaddcol,2,fitsheader,ra,'RA'
    fxbaddcol,3,fitsheader,dec,'DEC'
    fxbaddcol,4,fitsheader,glon,'GLON'
    fxbaddcol,5,fitsheader,glat,'GLAT'
    fxbaddcol,6,fitsheader,ra,'CRVAL2'
    fxbaddcol,7,fitsheader,dec,'CRVAL3'
    fxbaddcol,8,fitsheader,dec,'CALA'
    fxbaddcol,9,fitsheader,dec,'CALB'
    fxbaddcol,10,fitsheader,dblarr(1),'CRPIX1'
    fxbaddcol,11,fitsheader,dblarr(1),'CRVAL1'
    fxbaddcol,12,fitsheader,dblarr(1),'CDELT1'
    fxbaddcol,13,fitsheader,dblarr(1),'ETAMB'
    fxbaddcol,14,fitsheader,[1],'SCANNUM'
    fxbaddcol,15,fitsheader,dblarr(1),'SCANMIN'
    fxbaddcol,16,fitsheader,dblarr(1),'INTMEAN'
    ;fxbaddcol,14,fitsheader,dblarr(1),'SCANMEAN'
    ;            fxbwrite,fitsnumber,offspec_avg,13,1
    ;            fxbwrite,fitsnumber,continuum_fitpars,14,1
    ;            fxbwrite,fitsnumber,scanmeans,15,1

    fxbcreate,fitsnumber,savefile,fitsheader,extension
    print,"fitsnumber: ",fitsnumber," extension: ",extension
    ;print,flist
    fxbhelp,fitsnumber

    ;istat3 = masaccum(spec,accumspec,/new)

    ; debug stuff?  keep track of cal numbers
    calonarr = []
    caloffarr = []
    cntstokarr = []

    ; for status bar things
    clearline = string(13b) ;"

    ; normalized offs
    ; from Phil:
    ;   CntsToK trys to convert from spectrometer counts to kelvins. You can multiply it by any value  that is spectrometer counts.
    ; on*cntstoK, off*cntsToK, or (on-off)*cntsToK).
    ; Problems arise if you try to remove a baseline by say dividing an on / off.. This output is no longer in spectrometer counts...so you can't use
    ; cntsToK (Usually in this case you divide by a normalized off so on/Norm(off) still has spectrometer units)
    normoff0 = offspec_avg[*,0]
    normoff1 = offspec_avg[*,1]

    nspec = 0 ; debug parameter - count total # of spectra to make sure it matches expectations
    npts = 0 ; initialize npts to be zero for comparison later

    scannum = 0
    for ii=0,n_elements(flist)-1 do begin
        t0 = systime(/sec)

        file = flist[ii]
        istat1 = masopen(file,desc)
        istat2 = masgetfile(desc,spec)
        istat3 = masaccum(spec,accumspec,/new,/avg)
        istat4 = mascalval(accumspec.h,calI)
        masclose,desc

        ;print,file,accumspec.h.obsmode,accumspec.h.scantype
        ; cals are done at the end of each scan, and they should not cover the source,
        ; unlike at GBT where cals are done on-the-fly
        if ((accumspec.h.obsmode eq 'CAL' and accumspec.h.scantype eq 'ON' and ~keyword_set(hack_obsmode)) $
                or (keyword_set(hack_obsmode) and ((ii mod 4 eq 0)))) then begin 
            if debug then print,ii,"CAL ON"
            calon=accumspec.d 
            calonarr = [calonarr,calon]
            hascalon=1
        endif else if (accumspec.h.obsmode eq 'CAL' and accumspec.h.scantype eq 'OFF' and ~keyword_set(hack_obsmode) $
                or (keyword_set(hack_obsmode) and ((ii mod 4 eq 1)))) then begin 
            if debug then print,ii,"CAL OFF"
            caloff=accumspec.d
            caldiff = calon - caloff
            caloffarr = [caloffarr,caloff]
            hascaloff=1
        endif else begin
            if debug then print,ii,"OBSMODE: ",accumspec.h.obsmode," SCANTYPE: ",accumspec.h.scantype

            if not (hascalon and hascaloff) and ~keyword_set(hack_obsmode) then begin
                message,"No calibrations found for scan "+string(scannum)+" iteration "+string(ii)+" calon: "+$
                    string(hascalon)+" caloff: "+string(hascaloff)
            endif else begin
                ; reset to zero for next iter
                hascalon=0
                hascaloff=0
            endelse
            calI.cntsToK[0]=calI.calval[0]/(total(caldiff[calI.indused,0],1)/calI.npnts)
            calI.cntsToK[1]=calI.calval[1]/(total(caldiff[calI.indused,1],1)/calI.npnts)
            v = masfreq(accumspec.h,velcrdsys='L',/retvel,restfreq=restfreq)
            cntstokarr = [cntstokarr,[[calI.cntsToK]]]

            if n_elements(spec) ne npts and (ii gt 3) then message,'mismatch in spectra size'
            npts = n_elements(spec)

            ; debug tool: count total # of spectra
            nspec += npts

            if debug then begin
                print,format='("Off scaling at scannum ",I4," is ",E15,E15," for pol 1/2")',scannum,$
                    poly(scannum,continuum_fitpars[0,*]),poly(scannum,continuum_fitpars[1,*])
            endif else begin
                print,format='($,"Off scaling at scannum ",I4," is ",E15,E15," for pol 1/2",A)',scannum,$
                    poly(scannum,continuum_fitpars[0,*]),poly(scannum,continuum_fitpars[1,*]),$
                    clearline
            endelse

            etamb_polyfit = poly((90-spec.h.elevatio),[0.491544,0.00580397,-0.000341992])

            if debug then print,"Does ii/3 = scannum?",ii/3,scannum,ii/3 eq scannum

            ; subtract off the lowest value in each scan to get an approximate background subtraction
            ; (better than the unfortunate alternative of having higher power at high ZA)
            ; why does the total power increase with ZA?  Don't know!  it's as if the gain increases...
            scanmin1 = min(scanmeans[0,scannum*npts:(scannum+1)*npts-1])
            scanmin2 = min(scanmeans[1,scannum*npts:(scannum+1)*npts-1])
            scanmin = (scanmin1+scanmin2)/2.

            ; scale to the fitted location...
            ; off0 = normoff0 * poly(scannum,continuum_fitpars[0,*])
            ; off1 = normoff1 * poly(scannum,continuum_fitpars[1,*])
            ; scale to the minimum measured scanmean
            off0 = normoff0 * scanmin1
            off1 = normoff1 * scanmin2

            offspec = accumspec
            offspec.d = [[off0],[off1]]


            for jj=0,npts-1 do begin
                ; create an index for the FITS file (starts at 1)
                ; ii/3 = scannum, I hope
                index = jj+scannum*npts+1
                ; based on reduce_session, should be (on-off)/off * c
                ;specA = (spec[jj].d[*,0] - off0) * calI.cntsToK[0]
                ;specB = (spec[jj].d[*,1] - off1) * calI.cntsToK[1]
                specA = (spec[jj].d[*,0] - off0)/normoff0 * calI.cntsToK[0] / etamb_polyfit[jj]
                specB = (spec[jj].d[*,1] - off1)/normoff1 * calI.cntsToK[1] / etamb_polyfit[jj]
                ;plot,v,specA,/ys,xrange=[-20,50]
                fxbwrite,fitsnumber,float(total([[specA],[specB]],2)/2.0),1,index
                fxbwrite,fitsnumber,spec[jj].h.crval2,2,index
                fxbwrite,fitsnumber,spec[jj].h.crval3,3,index
                fxbwrite,fitsnumber,spec[jj].h.crval2,6,index
                fxbwrite,fitsnumber,spec[jj].h.crval3,7,index
                fxbwrite,fitsnumber,spec[jj].h.glon,4,index
                fxbwrite,fitsnumber,spec[jj].h.glat,5,index
                fxbwrite,fitsnumber,double(cali.cntstok[0]),8,index
                fxbwrite,fitsnumber,double(cali.cntstok[1]),9,index
                fxbwrite,fitsnumber,double(CRPIX),10,index
                fxbwrite,fitsnumber,double(CRVALL),11,index
                fxbwrite,fitsnumber,double(CDELT),12,index
                fxbwrite,fitsnumber,double(etamb_polyfit[jj]),13,index
                fxbwrite,fitsnumber,scannum,14,index
                fxbwrite,fitsnumber,double(scanmin),15,index
                fxbwrite,fitsnumber,double(mean(scanmeans[*,index-1])),16,index
            endfor

            ; only count scans that are not end-of-scan calibration
            scannum += 1

            if keyword_set(dostop) then stop

        endelse

    endfor

    fxbhelp,fitsnumber
    fxbfinish,fitsnumber

    print,"Done.  stop statement included just so that you can check out some local variables: "
    help,scanmeans, calonarr, caloffarr , cntstokarr 
    print,"Standard deviation/mean of counts-to-kelvin conversion: ",$
        stddev(cntstokarr[0:*:2])/mean(cntstokarr[0:*:2]),stddev(cntstokarr[1:*:2])/mean(cntstokarr[1:*:2])

    if keyword_set(dostop) then stop

end
