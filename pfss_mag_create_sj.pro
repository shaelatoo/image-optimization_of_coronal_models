pro pfss_mag_create_sj,magout,magtype,nlat0,gridtype,file=file, $
   quiet=quiet,adapt_ind=adapt_ind,no_zero=no_zero

; this is a revised version of the routine pfss_mag_create 
;   written by Mark DeRosa.  the modification is to correct a problem
;    under the case magtype=2.  the original assumes the GONG
;   synoptic magnetogram starts at longitude=1, which is not
;   correct for the files i'm using.  i'm altering the code 
;   to parse the initial longitude from the ".mrbqs" header file - 2015?
;
; update: 6/27/2017 - modified to adjust the image to have zero net 
;   flux, before resampling (sij)
; update: 8/29/2017 - modified to add magtype = 5 option: ADAPT flux 
;   map
; update:  9/6/2017 - added adapt_ind keyword to select which member
;   of the adapt ensemble should be chosen as the "magnetogram", 
;   default=0
; update:  1/19/2018 - added no_zero keyword to allow to create maps
;   where the net flux has not been removed; this may be the best
;   approach for the ADAPT maps
; update: ??/2018 - altered to adjust net flux *after* resampling 
;   instead of before
;
;+
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;
;  pfss_mag_create.pro - This procedure resamples a magnetic map of the full
;                        solar surface (such as a synoptic map) onto either a
;                        Legendre grid or an equally spaced lat-lon grid.  The
;                        Legendre grid is compatible with the PFSS field
;                        extrapolation software included in this package.
;
;  usage: pfss_mag_create,magout,magtype,nlat,gridtype,file=file,/quiet
;         where magout = output magnetogram image in (lon,lat) format
;               magtype = 0=IDL save file or structure containing flux
;                           concentrations
;                         1=Wilcox text file (see below)
;                         2=GONG 360x180 magnetogram in FITS format
;                         3=MDI 3600x1080 synoptic charts in FITS format
;                         4=HMI 3600x1440 diachronic charts in FITS format
;                         5= ? (others can be added in the future...)
;               nlat = number of gridpoints in latitude (default=48)
;               gridtype = 1=regularly spaced grid in lat and lon
;                          2=Legendre grid (default)
;               file = filename of input magnetic data (or in the case of
;                      magtype=0 it can be a structure containing the data)
;               quiet = set if minimal screen output is desired
;
;  notes: -This procedure merely resamples the input data.  It does *NOT*
;          convert from line-of-sight magnetic fields typically measured in
;          synoptic maps to Br (as would be needed to do a PFSS
;          extrapolation).
;
;         -The total (unsigned) flux of the input magnetogram, proportional to
;              mean(abs(data))*4*!dpi*R^2,
;          should about equal the total flux of the resampled magnetogram,
;              total(total(abs(magout),1)/nlon*weights)*2*!dpi*R^2
;          or, equivalently (using the mean_dtheta function),
;              mean_dtheta(total(abs(magout),1)/nlon,cos(theta))*4*!dpi*R^2.
;          The percentage difference between the two (sometimes as large as
;          5%) arises due to the fact that bilinear interpolation is used to
;          resample from a coarse grid, and from the fact that the gridpoints
;          near the poles are populated with the flux density of the closest
;          nearest gridpoint in the input synoptic map.
;
;         -One can display the resulting magnetogram in orthographic
;          projection using
;            @pfss_data_block  ;  gets lat and lon arrays
;            scim,mag,m=4,ortho=[270,0],olon=lon,olat=lat
;
;         -magtype=0: IDL save file (or structure) contains the following four
;                     arrays:
;                       nflux = number of flux concentrations
;                       fluxs = array of fluxes in units of Mx
;                       phis = array of longitudes in radians
;                       thetas = array of colatitudes in radians
;                     To get Mx/cm^2 at each point, multiply by
;                       nlat*nlon/(4*!dpi*rsun^2) where rsun=6.959d10 (cm)
;         -magtype=1: Wilcox text files are available at
;                     http://wso.stanford.edu/synopticl.html, with the
;                     input data in units of microTesla.
;         -magtype=2: GONG synoptic charts available at http://gong.nso.edu
;         -magtype=3: MDI synoptic charts are available at
;                     http://sun.stanford.edu/synop, and by clicking on the
;                     "Magnetic Field and Synoptic Maps from SOHO/MDI" link.
;                     This routine assumes one uses the Mr (radial field) FITS
;                     files, and not the Ml (line of sight) versions.  Note,
;                     doing a PFSS extrapolation at full resolution is
;                     lengthy, and not recommended due to the inaccuracies
;                     associated with higher-degree spherical harmonic
;                     transforms
;         -magtype=4: HMI synoptic charts are available from the JSOC, for
;                     dataseries of hmi.Synoptic_Mr* or hmi_test.Synoptic_Mr.
;                     This routine assumes one uses the Mr (radial field) FITS
;                     files, and not the Ml (line of sight) versions.  Note,
;                     doing a PFSS extrapolation at full resolution is
;                     lengthy, and not recommended due to the inaccuracies
;                     associated with higher-degree spherical harmonic
;                     transforms
;

;  M.DeRosa - 30 Jan 2002 - converted from earlier script
;              5 Mar 2002 - now successfully reads in Wilcox magnetogram
;                           tables without proper spacing between columns
;             27 Jun 2002 - added quiet keyword
;              8 Nov 2002 - added capability to resample onto rectangular
;                           grid, added gridtype parameter
;             12 May 2003 - converted common block to PFSS package format
;             21 May 2003 - added magtype = 0, based on output from Karel's
;                           flux-transport model
;              9 Mar 2005 - now parses WSO text files that do not start at 0
;                           degrees longitude
;              9 Mar 2005 - fixed registration error when parsing WSO synoptic
;                           maps (longitudes formerly ranged from 0 to 355
;                           when they should have ranged from 5 to 360)
;  G.Petrie - 14 Jul 2006 - added magtype=2 (GONG synoptic maps)
;  M.DeRosa - 25 Sep 2007 - for magtype=0, enabled a structure contiaining the
;                           data to supplied in lieu of the filename
;             15 Apr 2009 - added magtype=3 (MDI synoptic charts)
;              4 Mar 2010 - for magtype 1 (Wilcox) now converts from
;                           line-of-sight to radial magnetic field
;             10 Mar 2010 - for magtype 1 (Wilcox) now uses the same
;                           longitudinal grid as is used to compute the
;                           harmonic coefficients on the WSO webpage
;             19 Mar 2012 - added magtype=4 (HMI diachronic charts)
;              1 Jun 2012 - added better treatment of missing data for MDI
;
;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
;-

;pro pfss_mag_create,magout,magtype,nlat0,gridtype,file=file,quiet=quiet

  ;  print usage message
  if n_params() eq 0 then begin
    print,'  pfss_mag_create_sj,magout,magtype,nlat,gridtype,file=file,/quiet'
    return
  endif
  
  ;  include common block
  @pfss_data_block
  
  ;  set gridsize
  if n_elements(nlat0) eq 0 then nlat=48 else nlat=round(nlat0(0))
  nlon=nlat*2
  
  ;  set up grid
  if n_elements(gridtype) gt 0 then gridtype=gridtype(0) else gridtype=2
  case gridtype of
    1: begin  ;  uniformly spaced grid
      lat=linrange(nlat+1,-90,90)
      lat=(lat(0:nlat-1)+lat(1:nlat))/2
      theta=(90-lat)*!dpi/180
      weights=sin(theta)
      weights(0)=.5
      weights(nlat-1)=0.5
    end
    else: begin  ;  Legendre grid
      gaussquad_legendre,nlat,cth,weights
      theta=acos(cth)  ;  radians
      lat=90-theta*180/!dpi
    end
  endcase
  lon=(linrange(nlon+1,0,360))(0:nlon-1)
  phi=lon*!dtor
  
  case magtype of
  
    0: begin  ;  IDL save file filled with flux concentrations
    
      ;  determine whether it's a structure or a file
      if n_elements(file) eq 0 then file=-1
      case size(file(0),/type) of
        7: restore,file(0)  ;  string, so it's a filename
        8: begin  ;  structure
          phis=file(0).phis
          thetas=file(0).thetas
          fluxs=file(0).fluxs
          nflux=file(0).nflux
        end
        else: begin
          print,'  ERROR in pfss_mag_create_sj:  file keyword not valid'
          return
        end
      endcase
      phis=(phis+2*!pi) mod (2*!pi)  ;  puts phis between 0 and 2pi
      
      ;  create magnetogram by adding in each flux
      mag=fltarr(nlon,nlat)
      if not keyword_set(quiet) then print, $
        '  pfss_mag_create_sj:  adding '+strcompress(nflux,/r)+' sources'
      for i=0l,nflux-1 do begin
      
        if not keyword_set(quiet) then $
          pfss_print_time,'  adding sources: ',i+1,nflux,tst,slen
          
        ;  find out where this source is on our grid
        thc=get_interpolation_index(lat,90-!radeg*thetas(i))
        thc1=fix(thc)  &  thc2=(thc1+1)<(nlat-1)
        phc=get_interpolation_index([phi,!dpi*2],phis(i))
        phc1=fix(phc)  &  phc2=(phc1+1) mod nlon
        
        ;  add flux to grid, and divide by areal factor
        bco=(thc-thc1)
        aco=(phc-phc1)
        mag(phc1,thc1)=mag(phc1,thc1)+(1-aco)*(1-bco)*fluxs(i)/weights(thc1)
        mag(phc2,thc1)=mag(phc2,thc1)+aco*(1-bco)*fluxs(i)/weights(thc1)
        mag(phc1,thc2)=mag(phc1,thc2)+(1-aco)*bco*fluxs(i)/weights(thc2)
        mag(phc2,thc2)=mag(phc2,thc2)+aco*bco*fluxs(i)/weights(thc2)
        
      endfor
      magout=mag*mean(weights)  ;  normalization
      
    end
    
    1: begin  ;  Wilcox line-of-sight synoptic map
    
      ;  idiot check
      if (n_elements(file) eq 0) then begin
        print,'  pfss_mag_create_sj: filename must be specified'
        return
      endif
      
      ;  parse file
      openr,lun,file,/g
      fs=fstat(lun)
      tablebyte=bytarr(fs.size)
      readu,lun,tablebyte
      table=string(tablebyte)
      ix=strpos(table,'CT')
      posix=strpos(table,'CT',ix+1)
      repeat begin
        ix=[ix,posix]
        posix=strpos(table,'CT',posix+1)
      endrep until posix eq -1
      longitudes=round(float(strmid(table,ix+7,3)))
      ix=ix+18
      point_lun,lun,0
      data=fltarr(72,30)  ;  assumes 72 lon bins, 30 slat bins
      for i=0,71 do begin
        point_lun,lun,ix(i)
        buff1=fltarr(6)
        readf,lun,buff1,f='(6f9.3)'
        buff2=fltarr(8)
        readf,lun,buff2,f='(8f9.3)'
        buff3=fltarr(8)
        readf,lun,buff3,f='(8f9.3)'
        buff4=fltarr(8)
        readf,lun,buff4,f='(8f9.3)'
        data(i,*)=[buff1,buff2,buff3,buff4]
      endfor
      free_lun,lun
      data=shift(data,-(where(longitudes eq 360))(0),0)  ;  align longitudes to 0
      data=reverse(data,1)  ;  make longitude increasing
      data=reverse(data,2)  ;  make latitude (instead of colatitude) increasing
      
      ;  convert from line-of-sight to radial field
      dlonix=linrange(72,0,355)+2.5  ;  72 longitude bins
      dslatix=linrange(30,14.5,-14.5)/15  ;  30 slat bins
      dlatix=asin(dslatix)*180/!dpi
      slatgrid=replicate(1,72)#dslatix       ; MISTAKE?  WHY dslatix INSTEAD OF dlatix?
print,'Possible mistake in pfss_mag_create_sj.  Investigate.'
stop
      data=data/sqrt(1-slatgrid*slatgrid)

      ; adjust net flux
      ; note: I think Wilcox does their own net flux subtraction before 
      ;   publishing the synoptic charts, so this is probably 
      ;   redundant
      ; note: zero_array will probably not be tested with Wilcox data
      ;   for a while, if ever - caveat emptor
;      if NOT(KEYWORD_SET(no_zero)) then begin
;        zerofluxdata=ZERO_ARRAY(data)
;      endif else zerofluxdata=data
      
      ;  remap onto our grid
      dlatinterp=get_interpolation_index(reverse(dlatix),lat)
      dloninterp=get_interpolation_index(dlonix,lon)
      interpmap=interpolate(data,dloninterp,dlatinterp,/grid)
      if KEYWORD_SET(no_zero) then magout=interpmap else $
           magout=ZERO_ARRAY(interpmap)
      
      ;  no need to take into account unequal-sized pixels of Legendre grid,
      ;  since both input and output data are in terms of flux densities
      ;  (e.g. Tesla or Gauss) instead of flux (e.g. Weber or Maxwell)
      
    end
    
    2: begin  ;  (GONG) 360x180 magnetogram in FITS format
    
      ;  idiot check
      if (n_elements(file) eq 0) then begin
        print,'  pfss_mag_create_sj: filename must be specified'
        return
      endif
      
      ;  read fits file
      mreadfits,file,hdr,data
      
;      ; adjust net flux
;      if NOT(KEYWORD_SET(no_zero)) then begin
;        zerofluxdata=ZERO_ARRAY(data)
;      endif else zerofluxdata=data

      
      ;  remap onto our grid
      dlatix=asin(linrange(180,89.5,-89.5)/90.0)*180/!dpi  ;  180 slat bins
      dlonix=linrange(360,1,360) ;  360 longitude bins
      initlong=DOUBLE(hdr.mapedge)
      if initlong gt 0 then begin
        temp=data
        longitudes=dlonix+initlong
        loc=WHERE(longitudes gt 360)
        data=[temp[loc[0]:*,*],temp[0:loc[0]-1,*]]
        longitudes=[longitudes[loc[0]:*]-360.,longitudes[0:loc[0]-1]]
        dlonix=longitudes
      endif
      dlatinterp=get_interpolation_index(reverse(dlatix),lat)
      dloninterp=get_interpolation_index(dlonix,lon)
      interpmap=interpolate(data,dloninterp,dlatinterp,/grid)
      if KEYWORD_SET(no_zero) then magout=interpmap else $
           magout=ZERO_ARRAY(interpmap)
      
    end
    
    3: begin  ;  MDI synoptic charts
    
      ;  idiot check
      if (n_elements(file) eq 0) then begin
        print,'  pfss_mag_create_sj: filename must be specified'
        return
      endif
      
      ;  read fits file
      data=readfits(file,hdr)
      
      ;  set out of bounds points to zero
      blank=sxpar(hdr,'BLANK')
      wh=where(abs(data-blank) lt (1e-4),nwh)
      if nwh gt 0 then data(wh)=0.0
      
;      ; adjust net flux
;      if NOT(KEYWORD_SET(no_zero)) then begin
;        zerofluxdata=ZERO_ARRAY(data)
;      endif else zerofluxdata=data

      ;  remap onto our grid
      dlatix=asin(linrange(1080,539.5,-539.5)/540)*180/!dpi  ;  1080 slat bins
      dlonix=linrange(3600,0.1,360)  ;  3600 longitude bins
      dlatinterp=get_interpolation_index(reverse(dlatix),lat)
      dloninterp=get_interpolation_index(dlonix,lon)
      interpmap=interpolate(data,dloninterp,dlatinterp,/grid)
      if KEYWORD_SET(no_zero) then magout=interpmap else $
        magout=ZERO_ARRAY(interpmap)
        
        
      ;  no need to take into account unequal-sized pixels of Legendre grid,
      ;  since both input and output data are in terms of flux densities
      ;  (e.g. Tesla or Gauss) instead of flux (e.g. Weber or Maxwell)
      
    end
    
    4: begin  ;  HMI synoptic charts
    
      ;  idiot check
      if (n_elements(file) eq 0) then begin
        print,'  pfss_mag_create_sj: filename must be specified'
        return
      endif
      
      ;  read fits file
      data=readfits(file,hdr)
      
;      ; adjust net flux
;      if NOT(KEYWORD_SET(no_zero)) then begin
;        zerofluxdata=ZERO_ARRAY(data)
;      endif else zerofluxdata=data

      ;  set out of bounds points to zero
      wh=where(data lt (-3e4),nwh)
      if nwh gt 0 then data(wh)=0.0
      wh=WHERE(FINITE(data,/nan),nwh)
      if nwh gt 0 then data(wh)=0.0
      
      ;  remap onto our grid
      dlatix=asin(linrange(1440,719.5,-719.5)/720)*180/!dpi  ;  1440 slat bins
      dlonix=linrange(3600,0.1,360)  ;  3600 longitude bins
      dlatinterp=get_interpolation_index(reverse(dlatix),lat)
      dloninterp=get_interpolation_index(dlonix,lon)
      interpmap=interpolate(data,dloninterp,dlatinterp,/grid)
      if KEYWORD_SET(no_zero) then magout=interpmap else $
        magout=ZERO_ARRAY(interpmap)

      
      ;  no need to take into account unequal-sized pixels of Legendre grid,
      ;  since both input and output data are in terms of flux densities
      ;  (e.g. Tesla or Gauss) instead of flux (e.g. Weber or Maxwell)
      
    end
    
    5: begin  ;  ADAPT flux map
    
      ;  idiot check
      if (n_elements(file) eq 0) then begin
        print,'  pfss_mag_create_sj: filename must be specified'
        return
      endif
      
      ;  read fits file
      fxread,file,data,hdr
      if NOT(KEYWORD_SET(adapt_ind)) then adapt_ind=0
      data=REFORM(data[*,*,adapt_ind])
      
      ; adjust net flux
;      if NOT(KEYWORD_SET(no_zero)) then begin
;        zerofluxdata=ZERO_ARRAY(data)
;      endif else zerofluxdata=data
          
      
      ;  remap onto our grid  -  assumes map has GONG resolution 
      ;     and first longitude bin corresponds to 0 degrees
      if nlat ne 180 or nlon ne 360 then begin
        dlatix=linrange(180,89.5,-89.5)  ;  180 latitude bins
        dlonix=linrange(360,0,359) ;  360 longitude bins
        dlatinterp=get_interpolation_index(reverse(dlatix),lat)
        dloninterp=get_interpolation_index(dlonix,lon)
        interpmap=interpolate(data,dloninterp,dlatinterp,/grid)
        if KEYWORD_SET(no_zero) then magout=interpmap else $
          magout=ZERO_ARRAY(interpmap)
      endif else magout=ZERO_ARRAY(data)
        
    end
    

    else: begin
      print,'  pfss_mag_create_sj: invalid magtype'
      return
    end
    
  endcase
  
  if not keyword_set(quiet) then print,'  pfss_mag_create_sj:  magnetogram created'
  
  lat0=lat
  lon0=lon
  
end

