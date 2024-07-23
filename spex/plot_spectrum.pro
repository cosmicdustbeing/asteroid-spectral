pro plot_spectrum,infile,title,overplot=overplot,nonormal=nonormal,noerr=noerr,offset=offset,errcol=errcol
;
;
;
;
;

if not keyword_set(noerr) then begin

readcol,infile,wave,refl,/silent
;specerr,infile,'tempspec'
;readcol,'tempspec',wave,refl,err,format='F,F,F',/silent
;spawn,'rm -f tempspec'
err=fltarr(n_elements(wave))
endif else $
   readcol,infile,wave,refl,err,format='F,F,F',/silent
   wave = wave[where(refl gt 0.)]
   err = err[where(refl gt 0.)]
   refl = refl[where(refl gt 0.)]
if min(wave) le 0.451 then rng='VisNIR' else rng='NIR'
if n_elements(offset) lt 1 then offset=0

if rng eq 'NIR' then $
normrefl = refl/refl[mean(where(wave gt 1.18 and wave lt 1.22))] else $
normrefl = refl/mean(refl[where(wave gt 1.18 and wave lt 1.22)])

newerr = err*normrefl/refl

if keyword_set(nonormal) then begin

   normrefl = refl & newerr = err
endif

miny = (min(normrefl[where(wave lt 1.2 and finite(normrefl))])< 0.8); < min(refl)) > 0.2; - 0.03
maxy = (max(normrefl[where((wave lt 2.3 and wave gt 2.05) or wave lt 1.65 and $
                          finite(normrefl))]) +0.05); > (max(normrefl)) + 0.05)
minx = min(wave) & maxx = min(wave)

if not keyword_set(overplot) then begin
cgplot,wave,normrefl+offset,xtitle='Wavelength ($\mu$m)',ytitle='Normalized Reflectance',$
       yrange=[miny*0.96,maxy*1.04],xrange=[minx,maxx],xthick=3,ythick=3,charthick=3,$
       charsize=1.2,/nodata,title=title,xstyle=1

cgcolorfill,[1.09,1.09,1.15,1.15],[miny*.96,maxy*1.04,maxy*1.04,miny*.96],$
            color='Light Gray'
cgcolorfill,[1.32,1.32,1.46,1.46],[miny*.96,maxy*1.04,maxy*1.04,miny*.96],$
            color='Light Gray'
cgcolorfill,[1.77,1.77,1.95,1.95],[miny*.96,maxy*1.04,maxy*1.04,miny*.96],$
            color='Light Gray'
endif

oploterror,wave,normrefl+offset,newerr,psym=3,errthick=3.7,color=errcol,/nohat
cgplot,wave,normrefl+offset,psym=16,/overplot,symsize=0.3

if not keyword_set(overplot) then begin
cgaxis,xaxis=0,xrange=[minx,maxx],/save,xstyle=1,xthick=4,xtickformat='(A1)'
cgaxis,xaxis=1,xrange=[minx,maxx],/save,xstyle=1,xthick=4,xtickformat='(A1)'
cgaxis,yaxis=0,yrange=[miny*.96,maxy*1.04],/save,ystyle=1,ythick=4,ytickformat='(A1)'
cgaxis,yaxis=1,yrange=[miny*.96,maxy*1.04],/save,ystyle=1,ythick=4,charthick=3,charsize=1.2
endif
end
