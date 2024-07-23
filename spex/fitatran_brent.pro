pro shft_atrn,a,b,c,d
;
;shifts values around between these four variables
;
	a=b & b=c
	if (n_elements(d) eq 0) then c=a else c=d

end
;
;-------------------------------------------------------------------------
;
function atrnfunc,h2o,in,fin
;
;Shifts the object spectrum by h2o, divides the star spectrum,
; and calculates a chisq in 1.4 and 1.8 um bands.
;

        if h2o mod 0.1 ne 0 then begin

           geth2o = [[min(where(h2o lt in.h2o)-1)],$
                     [min(where(h2o lt in.h2o))]]
           useatran = reform(in.atran[geth2o[0],*]*10.*(.1-(h2o mod 0.1))) + $
                      reform(in.atran[geth2o[1],*]*10.*((h2o mod 0.1)))
        endif else begin

           geth2o = where(h2o eq in.h2o)
           useatran = reform(in.atran[geth2o,*])
        endelse

        ; find best shift
        objshift_brent,in.lam,in.spec,in.err,in.flags,useatran,useatran*10000.,out,$
                       dxmin,shiftspec ; AVOID wholeband keyword
        diverr=out.err[out.minsub,*]

        ; select telluric band to fit
        specsub = [shiftspec[in.sub1],shiftspec[in.sub2]]
        c1 = poly_fit(in.lam[in.injh],shiftspec[in.injh],3,/double)
        c2 = poly_fit(in.lam[in.inhk],shiftspec[in.inhk],3,/double)
        smthspecsub = [poly(in.lam(in.sub1),c1),poly(in.lam(in.sub2),c2)]

        chisq=total(((specsub-smthspecsub)/diverr)^2)

        if (n_elements(specsub) gt 1) then $
           chisq=chisq/((n_elements(specsub)-3.) > 1.) ;reduced chisq

	fin.h2o[fin.i]=h2o
	fin.chi[fin.i]=chisq
        fin.shift[fin.i]=dxmin
	fin.newspec[fin.i,*]=out.newspec[out.minsub,*]
	fin.div[fin.i,*]=shiftspec
        fin.err[fin.i,*]=diverr
	fin.i=fin.i+1
;	ins=where(fin.chi ne 0)
;	if (chisq le min(fin.chi[ins])) then fin.modflux=modflux
	return,chisq

end
;
;-------------------------------------------------------------------------
;
pro bracket_atrn,h2o1,h2o2,in,h2oarr,chiarr,fin
;
;*************************************************************************
;Provides initial bracketing of points before beginning actual 1-D
; minimization.  From routine in Numerical Recipes, 10.1.
;INPUTS:
;   h2o1   - one initial guess of h2o
;   h2o2   - a second initial guess of h2o
;   in   - a structure containing necessary information for the thermal
;          model, see main program for included parameters.
;OUTPUTS:
;   h2oarr   - 3 element array of h2o points that bracket the minimum
;   chiarr - 3 element array of associated chisquare values.
;*************************************************************************
;
	grat=1.618034	;ratio by which successive intervals are magnified
	glim=100.		;maximumum magnification allowed
	tiny=1e-20
;
;Set up two radius limits, one for A=0.99 and one for A=0.999
;
	h2olim1=0.01
	h2olim2=10.
;
;Call shifting function for initial h2o estimates to get chisquare values
;
	chisq1=atrnfunc(h2o1,in,fin)
	chisq2=atrnfunc(h2o2,in,fin)
;
;Make sure chisq1 > chisq2
;
	if (chisq2 gt chisq1) then begin
		shft_sh,tmp,h2o1,h2o2 & shft_sh,tmp,chisq1,chisq2
	endif
;
;Estimate h2o3 from golden ratio, then run thermal model again for it
;
	h2o3=(h2o2+grat*(h2o2-h2o1)) > h2olim1
        h2o3=(h2o2+grat*(h2o2-h2o1)) < h2olim2
	chisq3=atrnfunc(h2o3,in,fin)
;
;Now find 4th point (parabolic extrap), estimated at minimum.  Will check
; this point and do several bookkeeping things depending on relative value.
; Goal is to find 3 points such that chi1>chi2 and chi3>chi2.
;
	i=0
	while (chisq2 gt chisq3) do begin
		p=(h2o2-h2o1)*(chisq2-chisq3)	;these 3 lines are parabolic extrap
		q=(h2o2-h2o3)*(chisq2-chisq1)	;ru line also avoids divide by zero
		ru=h2o2-((h2o2-h2o3)*q-(h2o2-h2o1)*p)/(2.*(abs(q-p)>tiny)*$
			(1.-2.*((q-p) lt 0.)))
		ru=ru > h2olim1
                ru=ru < h2olim2
		rulim=(h2o2+glim*(h2o3-h2o2)); > h2olim2	;largest step allowed
;stop
		case 1 of
			((h2o2-ru)*(ru-h2o3)) gt 0.: begin	;ru between h2o2 and h2o3
					chisqu=atrnfunc(ru,in,fin)
					case 1 of
						chisqu lt chisq3: begin	;minimum between h2o2 and h2o3
							h2o1=h2o2 & h2o2=ru 		;put 2 into 1 and u into 2
							chisq1=chisq2 & chisq2=chisqu
							goto,done
						   end
						chisqu gt chisq2: begin ;min between h2o1 and h2o2
						    h2o3=ru & chisq3=chisqu
						    goto,done
						   end
						else:
					endcase
					ru=h2o3+grat*(h2o3-h2o2)	;parabolic no good, use grat and
;					while (ru gt h2olim1) do begin
;						grat=grat*.8 & ru=h2o3+grat*(h2o3-h2o2)
;					endwhile
					chisqu=atrnfunc(ru,in,fin)	; reset h2o3 and begin loop again
				  end
			((h2o3-ru)*(ru-rulim)) gt 0.: begin	;ru between h2o3 and rulim
					chisqu=atrnfunc(ru,in,fin)
;
;If this step doesn't produce min yet, set up so that the points propogated
; to the next loop are h2o3, ru, and a golden extrap from these two.
;(if chisqu gt chisq3, normal progression of program takes care of all)
;
					if (chisqu lt chisq3) then begin
						tmp=ru+grat*(ru-h2o3)
						while (tmp gt h2olim1) do begin
							grat=grat*.8 & tmp=ru+grat*(ru-h2o3)
						endwhile
						chitmp=atrnfunc(tmp,in,fin)
						shft_sh,h2o2,h2o3,ru,tmp
						shft_sh,chisq2,chisq3,chisqu,chitmp
					endif
				  end
			((ru-rulim)*(rulim-h2o3)) gt 0: begin	;ru beyond rulim, reset
					ru=rulim & chisqu=atrnfunc(ru,in,fin)
				  end
			else: begin		;last resort, golden extrapolation
					ru=h2o3+grat*(h2o3-h2o2)
;					while (ru gt h2olim1) do begin
;						grat=grat*.8 & ru=h2o3+grat*(h2o3-h2o2)
;					endwhile
					chisqu=atrnfunc(ru,in,fin)
				  end
		endcase
;
;Eliminate oldest point (h2o1) and continue.
;
;stop
		shft_sh,h2o1,h2o2,h2o3,ru
		shft_sh,chisq1,chisq2,chisq3,chisqu
		i=i+1
;		print,i
	endwhile

done:

	h2oarr=[h2o1,h2o2,h2o3]
	chiarr=[chisq1,chisq2,chisq3]

end
;
;-------------------------------------------------------------------------
;
pro brent_atrn,h2oarr,chiarr,in,h2omin,chimin,fin
;
;*************************************************************************
;Routine for Brent's method in 1-D.  Zeros in on minimum of function
; through successive calculation of function, and parabolic extrapolation,
; no derivatives. See Numerical Recipes, 10.2.
;INPUTS:
;   rarr   - 3 element array of radius points that bracket minimum
;   chiarr - 3 element array of corresponding chisq values
;   in     - structure of input parameters for thermal model
;OUTPUTS:
;   rmin   - minimum radius (km)
;   chimin - corresponding chisquare goodness of fit value
;   fin    - output structure
;*************************************************************************
;
	R=0.61803399	;golden ratio
	C=(1.0-R)
	tol=1e-4		;tolerance, i.e. level to which min is known
	tiny=1e-10
	e=0.
;
;Keeps track of 6 points.  Initally, points a and b are h2oarr(0) and h2oarr(2)
; (a & b bracket minimum, x is lowest function value so far, w is 2nd
;  lowest, v is previous value of w, u is point of most recent evaluation)
; x,v,w are initiallized as h2oarr(1).
;
	a=h2oarr(0) < h2oarr(2)
	b=h2oarr(0) > h2oarr(2)
	x=(w=(v=h2oarr(1)))
	fw=(fv=(fx=chiarr(1)))

	xm=0.5*(a+b)	;midpoint between a and b
	tol1=tol*abs(x)+tiny
	tol2=2.*tol1
;
;Begin loop which selects additional points, narrowing the gap between
; a and b until it is smaller than tolerance.
;
	while (abs(x-xm) gt tol2-0.5*(b-a)) do begin
;stop
		if (abs(e) gt tol1) then begin	;construct trial parabolic fit
			r=(x-w)*(fx-fv)
			q=(x-v)*(fx-fw)
			p=(x-v)*q-(x-w)*r
			q=2.0*(q-r)
			if (q gt 0.0) then p=-p
			q=abs(q)
			etmp=e & e=d
;
;First condition below rejects the parabolic fit and takes the golden
; section step into the larger of the two segments, secon accepts
; parabolic fit.
;
			if ((abs(p) ge abs(0.5*q*etmp)) or (p le q*(a-x)) or $
				(p ge q*(b-x))) then begin
					if (x ge xm) then e=a-x else e=b-x
					d=c*e
			endif else begin
				d=p/q & u=x+d
				if (u-a lt tol2 or b-u lt tol2) then $
					d=tol1*(xm-x)/abs(xm-x)
			endelse
		endif else begin	;golden section step
			if (x ge xm) then e=a-x else e=b-x
			d=c*e
		endelse

		if (abs(d) ge tol1) then u=x+d else u=x+tol1*d/abs(d)
		fu=atrnfunc(u,in,fin)
;		print,u,fu

		if (fu le fx) then begin	;fu new minimum, shift values over
			if (u ge x) then a=x else b=x
			shft_sh,v,w,x,u
			shft_sh,fv,fw,fx,fu
		endif else begin	;fu not new min, several possibilites
			if (u lt x) then a=u else b=u
			case 1 of
				(fu le fw or w eq x): begin
						v=w & w=u & fv=fw & fw=fu
					end
				(fu le fv or v eq x or v eq w): begin
						v=u & fv=fu
					end
				else:
			endcase
		endelse

		xm=0.5*(a+b)	;midpoint between a and b
		tol1=tol*abs(x)+tiny
		tol2=2.*tol1
	endwhile

	h2omin=x & chimin=fx

end
;
;-------------------------------------------------------------------------
;
pro fitatran_brent,lam,zenith,objspec,objerr,flags,fin2,h2omin,atrancorr,h2o1=h2o1
;
;*************************************************************************
;Wraps around sub-pixel shift to zero in on the best relative shift
; between the object and star spectrum (due to e.g. instrument flexure).
; Uses Golden Section Search as described in Numerical Recipes.
;
; INPUTS:
;   lam      - wavelength (same size as objspec)
;   zenith - zenith angle of the asteroid observation (integer)
;   objspec  - object spectrum (single array of flux values)
;   objerr   - object uncertainty (same size as objspec)
;
; OPTIONAL INPUTS:
;   h2o1   - first guess for water amount (in mm)
;
; OUTPUTS:
;   fin2   - structure containing all h2o attempts, chisquares, and
;             best fit divided spectra.
;   h2omin  - water amount at minimum chisquare
;   atrancorr - atran corrected spectrum at minimum chisquare
;*************************************************************************
;
;
  h2o1=6.6 & h2o2=0.1
  h2o = findgen(101)/10.

; set atran files path
atranpath = '/Volumes/mica/research/tools/irtf_idl_library/atran_prism/'

; obtain atran table
atran0 = read_atran_prism('A'+strtrim(zenith,2)+'.ALL2',atranpath)

; ATRAN files do not cover the full wavelength range of SpeX Prism mode
; so atmospheric transmission of 1 is assumed outside of the range
; provided in the tables
atran1 = fltarr(101,n_elements(lam))

; interpolate to asteroid wavelength values
for i=1,100 do begin

   atran1[i-1,*] = (interpol(atran0[i,*],atran0[0,*],lam,/spline)) < 1.
   atran1[i-1,where(lam gt max(atran0[0,*]))] = 1.
   atran1[i-1,where(lam lt min(atran0[0,*]))] = 1.
endfor
nlam = n_elements(lam)

sub1 = where(lam ge 1.35 and lam le 1.45 and flags lt 1)
sub2 = where(lam gt 1.77 and lam lt 1.95 and flags lt 1)
inj = where(lam gt 1.18 and lam lt 1.31 and flags lt 1)
inh = where(lam gt 1.51 and lam lt 1.74 and flags lt 1)
inks = where(lam gt 2.05 and lam lt 2.25 and flags lt 1)

;
;Call routine to find bracketing points before beginning minimization.
; First set up structure to easily pass all these variables.
;
	in={spec:objspec,err:objerr,flags:flags,atran:atran1,h2o:h2o,lam:lam,$
            sub1:sub1,sub2:sub2,injh:[inj,inh],inhk:[inh,inks]}
	fin={i:0,shift:fltarr(100),h2o:fltarr(100),chi:fltarr(100),$
             newspec:fltarr(100,nlam), div:fltarr(100,nlam),err:fltarr(100,nlam)}
	bracket_atrn,h2o1,h2o2,in,h2oarr,chiarr,fin
;
;Golden Section search.  Use the 3 radius points found by bracket,
; and only will need one of the chisq.  (They use none, but do a function
; call, I'll avoid that call).
;
	brent_atrn,h2oarr,chiarr,in,h2omin,chimin,fin

;	print,' '
	ins=where(fin.chi ne 0)
	mn=min(fin.chi(ins),minsub)
	sub=where(fin.chi eq chimin)
	h2omin=fin.h2o(minsub)
	newdiv=reform(fin.div(minsub,*))
        newlam=interpol(reform(lam),findgen(nlam),findgen(nlam)+fin.shift(minsub))

        atran = fin.newspec[minsub,*]/fin.div[minsub,*]
        atranerr = objerr/atran

	print,'h2o = '+strtrim(h2omin,2)+' mm';,chimin
        atrancorr = [[newlam],[newdiv],[atranerr]]
       fin2={h2o:fin.h2o(ins),chi:fin.chi(ins),newspec:fin.newspec(ins,*),$
             div:fin.div(ins,*),err:fin.err(ins,*),minsub:minsub}

done:

end
