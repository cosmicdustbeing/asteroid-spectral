pro stardiv_atran_prism,dir,asteroid,star,outdir,faint=faint
;
; DESCRIPTION:
;     performs telluric feature removal for an asteroid and star
;     PRISM spectrum by diving by a standard star
;
; INPUT:
;     astfits - asteroid spectrum (output from SpeXtool extraction steps)
;
; KEYWORD:
;     faint - set to bin data when plotting (useful for low S/N)
;
; OUTPUT:
;     this procedure automatically write the divided spectrum for the
;     asteroid and writes to a fits file with "_div_" inserted in the
;     filename in the 'outdir' directory
;

resolve_all
wholeband = 0

; search for asteroid and star files in directory
ast_files = file_search(dir+'/'+asteroid+'*prism*.fits')
star_files = file_search(dir+'/'+star+'*prism*.fits')

; initialize arrays
astbeam=strarr(n_elements(ast_files))
starbeam=strarr(n_elements(star_files))

astairmass=strarr(n_elements(ast_files))
starairmass=strarr(n_elements(star_files))

astzen=strarr(n_elements(ast_files))
starzen=strarr(n_elements(star_files))

pixshift=strarr(n_elements(star_files))

; loop through all fits files for asteroid and star
for i=0,n_elements(ast_files)+n_elements(star_files)-1 do begin

   ; read in fits files
   if i lt n_elements(ast_files) then file = ast_files[i] $
                                 else file = star_files[i-n_elements(ast_files)]
   spec = readfits(file,header,/silent)

   ; extract airmass and beam position
   airmass = fxpar(header,'TCS_AM')
   elevation = fxpar(header,'TCS_EL')
   zenith = nint(floor(91.-elevation))
   beam = fxpar(header,'BEAM')
   lam=reform(spec[*,0])

   ; perform atran telluric removal
   fitatran_brent,lam,zenith,reform(spec[*,1]),reform(spec[*,2]),reform(spec[*,3]),$
                  atrnout,h2o,atrancorr

   ; make structure to hold corrected spectra for asteroid and star frames
   if i eq 0 then allspec={asteroid:fltarr(n_elements(ast_files),$
                                           (size(spec))[1],(size(spec))[2]),$
                           astatran:fltarr(n_elements(ast_files),$
                                           (size(spec))[1],(size(spec))[2]),$
                           star:fltarr(n_elements(star_files),$
                                       (size(spec))[1],(size(spec))[2]),$
                           stratran:fltarr(n_elements(star_files),$
                                           (size(spec))[1],(size(spec))[2])}

   ; asteroid frames
   if i lt n_elements(ast_files) then begin
      allspec.asteroid[i,*,*] = reform(spec[*,0:3])
      allspec.astatran[i,*,0:2] = atrancorr & allspec.astatran[i,*,3] = reform(spec[*,3])
      astbeam[i] = beam & astairmass[i] = airmass
      if i eq 0 then asthead = strarr(n_elements(ast_files),n_elements(header))
      asthead[i,*] = header
   endif else begin
   ; star frames
      allspec.star[i-n_elements(ast_files),*,*] = reform(spec[*,0:3])
      allspec.stratran[i-n_elements(ast_files),*,0:2] = atrancorr
      allspec.stratran[i-n_elements(ast_files),*,3] = reform(spec[*,3])
      starbeam[i-n_elements(ast_files)] = beam
      starairmass[i-n_elements(ast_files)] = airmass
   endelse
endfor

; loop through all asteroid frames
for j=0,n_elements(ast_files)-1 do begin

   ; look at star files in the same position along the slit (A or B)
   ; also ignore if the asteroid/star airmass difference is more than 0.1
   index=where(starbeam eq astbeam[j]); and $
;              abs(float(starairmass)-float(astairmass[j])) lt 0.35)

   ; array to hold divided spectra and error
   div1 = make_array(n_elements(index),n_elements(spec[*,0]))
   err1 = make_array(n_elements(index),n_elements(spec[*,0]))
;   div2 = div1 & err2 = err1

   astframe = (strsplit(ast_files[j],'-.',/extract))$
              [n_elements((strsplit(ast_files[j],'-.',/extract)))-2]

   redo:
   for k=0,n_elements(index)-1 do begin

      ; divide asteroid by star
      div = allspec.astatran[j,*,1]/allspec.stratran[k,*,1]
      err = div*sqrt((allspec.astatran[0,*,2]/allspec.astatran[0,*,1])^2. + $
                     (allspec.star[k,*,2]/allspec.star[k,*,1])^2.)

      div1[k,*] = div/div[mean(where(allspec.astatran[j,*,0] gt 1.18 and $
                                     allspec.astatran[j,*,0] lt 1.22))]

      err1[k,*] = err*div1[k,*]/div

      if k eq 0 then window,0,ysize=150+100*nint(n_elements(index)-1),xsize=730
      if k eq 0 then begin
         cgplot,findgen(1),/nodata,yrange=[.7,1.5+((n_elements(index)-1)*.35)],$
                position=[0.08,0.1,0.98,0.95],xrange=[0.6,2.6],xstyle=1,ystyle=1, $
                xtitle='Wavelength (microns)',ytitle='Scaled Reflectance'

         cgcolorfill,[1.09,1.09,1.15,1.15],[0.7,1.5+((n_elements(index)-1)*.35),$
                                         1.5+((n_elements(index)-1)*.35),0.7],$
                     color='Light Gray'
         cgcolorfill,[1.32,1.32,1.46,1.46],[0.7,1.5+((n_elements(index)-1)*.35),$
                                         1.5+((n_elements(index)-1)*.35),0.7],$
                     color='Light Gray'
         cgcolorfill,[1.77,1.77,1.95,1.95],[0.7,1.5+((n_elements(index)-1)*.35),$
                                         1.5+((n_elements(index)-1)*.35),0.7],$
                     color='Light Gray'

         cgaxis,/save,xstyle=1
         cgaxis,/save,xstyle=1,/xaxis

         ; load color table
         cgloadct,13
      endif
      if keyword_set(faint) then $
         cgplot,allspec.asteroid[j,*,0],smooth(div1[k,*],5)+.37*k,yrange=[0.7,1.5],$
                symsize=0.7,/overplot,thick=2,psym=16,$
                color=60*(k mod 2 eq 1)+250*(k mod 2 eq 0) $
      else $
         cgplot,allspec.asteroid[j,*,0],div1[k,*]+.37*k,yrange=[0.7,1.5],psym=16,$
                symsize=0.7,/overplot,$
                color=60*(k mod 2 eq 1)+250*(k mod 2 eq 0)

   endfor

usespec=-99
print,'Working on asteroid frame #'+astframe
;print,'enter 0 to skip or ...'
read, 'enter the # of the spectrum to write [1-'+strtrim(n_elements(index),2)+']: ',$
     usespec

if usespec eq -1 then begin
   if wholeband eq 0 then wholeband=1 else wholeband=0
endif else wholeband=0
if usespec eq -1 then begin
   print,' '
   if wholeband eq 0 then print,'fitting to telluric edge...'
   if wholeband eq 1 then print,'fitting to whole band...'
   print,' '
   goto,redo
endif
if usespec eq 0 then goto,skip
if usespec le -2 then goto,quit

; build array for output
ast_spec = [[reform(allspec.asteroid[j,*,0])],$
            [reform(div1[usespec-1,*])],[reform(err1[usespec-1,*])],$
            [reform(allspec.asteroid[j,*,3])]]

; add history to header
hdr = mc_addhistlabel(reform(asthead[j,*]),'My Stardiv ATRAN Asteroid History',CANCEL=cancel)
addhist = 'The spectra were divided by an ATRAN model spectrum of the atmosphere and then by the spectra of the standard star '+ $
          (strsplit(star_files[index[usespec-1]],'/',/extract))$
          [n_elements((strsplit(star_files[index[usespec-1]],'/',/extract)))-1]+'. '+$
          'Therefore the spectrum is unitless. The asteroid spectrum was shifted by '+$
          strtrim(pixshift[usespec-1],2)+' pixels.'
history = mc_splittext(addhist,70,CANCEL=cancel)
sxaddhist,history,hdr

; write output asteroid and star spectra and final divided spectrum to fits
outstr = outdir+'/'+asteroid+'-'+astframe+'_tellcordiv_'+star+'-'+ $
         (strsplit(star_files[index[usespec-1]],'-.',/extract))$
         [n_elements((strsplit(star_files[index[usespec-1]],'-.',/extract)))-2]+'.fits'

writefits,outstr,ast_spec,hdr
print,'wrote: '+outstr
print,' '

skip:
if usespec eq 0 then begin
   print,'skipping frame...' & print,' ' & endif

endfor
quit:
if usespec eq -99 then print,'Exiting... bye!' else $
   print,'Done!'
end
