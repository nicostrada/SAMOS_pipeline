pro pipeline_V2
!PATH=!PATH+';C:\Documents and Settings\robberto\My Documents\_SCIENCE\Projects\IRMOS\Commissioning\mychecks'
;PATH='C:\Documents and Settings\robberto\My Documents\_SCIENCE\Projects\IRMOS\Run6/night4/'
PATH='C:\Documents and Settings\robberto\My Documents\_SCIENCE\Projects\IRMOS\Run6\Night4\'
SOURCE_NAME='Orion_N4_'

; 1) LOAD FILES
; 1.1) Source
first_frame=22
last_frame=31
name_cube='source'
IRSPEC_load,PATH,SOURCE_NAME,first_frame,last_Frame,SOURCE_DCube
Nsource_frames=last_frame-first_frame+1

; 1.2) Source darks
first_frame=32
last_frame=34
IRSPEC_load,PATH,SOURCE_NAME,first_frame,last_Frame,DKSource_DCube

; 1.3) flats (lamp)
first_frame=36
last_frame=37
IRSPEC_load,PATH,SOURCE_NAME,first_frame,last_Frame,Flat_DCube

; 1.4) dark flats (just a single image in this case, so go with read_fits)
frame=35
fits_read,PATH+SOURCE_NAME+STRTRIM(STRING(frame),2)+'.fit',DKFlat_DCube

; 2.) First we do some work on the flats
;2.1) combine the flats with a median and subtract their dark
;medcomb,Flat_DCube,flat0
flat0 = MEDIAN(Flat_DCube, DIMENSION=3)
flat0=flat0-DKFlat_DCube
;2.2) remove the hotpixels. Cleaning comes always before rotation!
sigma=7
r_in=2
r_out=5
IRMOS_removeHOTPIX,flat0,flat0_c,sigma,r_in,r_out
;2.3) Rotate the flats to have vertical slits
im=reverse(rotate(flat0_c,3),2)
flat0_cr=rot(im,1.3,/interp)
;2.4) We use the flats for two reasons:
;2.4.1) Location of the individual spectra on the rotated flats
;2.4.2) Find the normalization factor of each flat
;For the first point you must visually inspect the image to define
n_targets=15  ;nr. of slits, if you don't know...
threshold=70  ;input to TARGET_FINDER. Start with this value, run
;TARGET_FInDER, inspect the plot with the spikes and make sure that the
;threshold is not too high/low. May need to iterate
TARGET_FINDER,flat0_cr,n_targets,threshold,x0_target,X1_target
;For the second point, we just go with a median value to build an
;array of flat field normalization factors NormFF
;.. NormFF=fltarr(n_targets)
;.. for ix = 0,N_targets-1 do NormFF[ix] = MEDIAN(flat0_cr[X0_Target[ix]:X1_Target[ix],*])

;3) BACK TO THE STARS
;we have some flacky guiding. Cannot coadd these source data directly, but we need to shift
;each frame to realign.
;
;3.1) Build the source dark
;medcomb,DKSource_DCube,dark0
dark0 = MEDIAN(DKSource_DCube, DIMENSION=3)
;3.2) Apply the flat field to the cube
SOURCE_DCube_FF = fltarr(1024,1024,NSource_frames)
for ix = 0,NSource_FRames-1 do begin & $
    Source_DCube_FF[*,*,ix] = (Source_DCube[*,*,ix]-dark0)/(flat0>0.5) & $
endfor
;3.3) normally here we would coadd, filter the result and rotate but in this case
;we need to realign before combining
;...not this:
;medcomb,Source_DCube_FF,Source_FF
;sigma=7
;r_in=2
;r_out=5
;IRMOS_removeHOTPIX,Source_FF,Source_FFc,sigma,r_in,r_out
;im=reverse(rotate(Source_FFc,3),2) & Source_FFcR=rot(im,1.3,/interp)
;... but this....
;-> We will realign and combine the rotate spectra. Before rotate we always
;have to clean them so the badpix do not spread out to the neighbouring ones.
Source_DCube_FFcr = fltarr(1024,1024,NSource_frames)
sigma=50 & r_in=2 & r_out=5
for ix = 0,NSource_FRames-1 do begin & $
    IRMOS_removeHOTPIX,Source_DCube_FF[*,*,ix],result,sigma,r_in,r_out & $
    im = reverse(rotate(result,3),2) & $
    Source_DCube_FFcr[*,*,ix] = rot(im,1.3,/interp) & $
endfor
;-> now that everything is clean and vertical,
;find on which columns are the stellar spectra. We refere to the first image
;so the loop starts from 1
v=median(Source_DCube_FFcr,DIMENSION=2)
shifts=[-3,-2,-1,0,1,2,3]
shiftOK=fltarr(NSource_Frames)
for ix=1,NSource_Frames-1 do begin  &    $
    result=c_correlate(v[250:*,0],v[250:*,ix],shifts)  &    $
    shiftOK=where(result EQ MAX(result))-3 &   $
    print,ix,shiftOK & $
    Source_DCube_FFcrs=shift(Source_DCube_FFcr,-shiftOK,0,0)    &  $
endfor
;3.4) finally we can coadd
medcomb,Source_DCube_FFcrs,SOURCE_FFcr

stop

end
pro aa


;4) Flat field the image
;FFimage=image/((lamp-lamp_dark)>1E-3)  ;not yet normalized!
;FFimage(WHERE(FFIMAGE GT 1 OR FFIMAGE LT -1))=!VALUES.F_NAN
;im=reverse(rotate(FFimage,3),2)
;spectra=rot(im,1.3,/interp)
SPECTRA=SOURCE_FFcr
;FIND OUT THE STARS ON THE MAP
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/Orion_N4_12.fit',field
im=reverse(rotate(field,3),2)
field=rot(im,1.3,/interp)
;
;NOW THE SPECTRA CAN BE EXTRACTED
Xsize=X1_target-X0_Target+1
is=0 ;slit6, y0=415
SPECTRA0=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=1 ;slit5, y0=387
SPECTRA1=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=2 ;slit4, y0=353
SPECTRA2=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=3 ;slit14
SPECTRA3=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=4 ;slit3 ;SKY
SPECTRA4=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=5 ;slit2
SPECTRA5=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=6 ;slit13
SPECTRA6=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=7 ;slit 1
SPECTRA7=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=8 ;slit12
SPECTRA8=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=9 ;slit11 ;SKY
SPECTRA9=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=10 ;slit 10
SPECTRA10=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=11 ;slit0
SPECTRA11=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=12 ;slit9
SPECTRA12=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=13 ;slit 8
SPECTRA13=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
is=14 ;slit7
SPECTRA14=SPECTRA[x0_target[is]:x1_target[is],*]*NormFF[is] ;ff norm is here!
;
s0=fltarr(1024) & for iy=0,1023 do begin s0[iy]=median(Spectra0[1:Xsize[0]-2,iy])
plot,s0
p1=850 & p2=950 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         902
p1=800 & p2=850 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         820
p1=700 & p2=750 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         728
p1=700 & p2=720 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         704
p1=650 & p2=700 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         680
p1=640 & p2=670 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         650
p1=600 & p2=650 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         626
p1=550 & p2=600 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         570
p1=270 & p2=320 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         299
p1=180 & p2=220 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         197
p1=130 & p2=180 & print,where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
;         145

s0_OHpix  = [ 902, 820, 728, 704, 680, 650, 626, 570, 299, 197, 145]
OHlines= [1.97665,2.00025,2.02709,2.0337,2.0406,2.04961,2.05561,2.07153,2.14988,2.1793,2.19492]
plot,OHlines,s0_OHpix,psym=2
h=linfit(ohlines,s0_ohpix)
print,h
oplot,ohlines,h[1]*ohlines+h[0]
wl=1.95+indgen(300)/1000.
plot,wl,h[1]*wl+h[0]
;conversion pix=>wl
pix=indgen(1024)
ll=(pix-h[0])/h[1]
plot,pix,ll
;check the Brgamma:
p1=220 & p2=270 & pixBg=where(s0[p1:p2] EQ MAX(S0[p1:p2]))+p1
ll=(pixBg-h[0])/h[1]
print,ll,s0

;s0 is at y=415. It is therefore
DIFF=415-s0_OHpix
;this goes in the calibration of the wl (calibrate_OH.pro)

;now on the second spectrum ay y=387
s1=fltarr(1024) & for iy=0,1023 do begin s1[iy]=median(Spectra1[1:4,iy])
plot,s1
calibrate_OH,s1,387,ll_1

s2=fltarr(1024) & for iy=0,1023 do begin s2[iy]=median(Spectra2[10:14,iy])
calibrate_OH,s2,353,ll_2
plot,ll_2,s2

s3=fltarr(1024) & for iy=0,1023 do begin s3[iy]=median(Spectra3[10:14,iy])
calibrate_OH,s3,412,ll_3
plot,ll_3,s3

s4=fltarr(1024) & for iy=0,1023 do begin s4[iy]=median(Spectra4[10:14,iy])
calibrate_OH,s4,405,ll_4
plot,ll_4,s4

s5=fltarr(1024) & for iy=0,1023 do begin s5[iy]=median(Spectra5[10:14,iy])
calibrate_OH,s5,395,ll_5
plot,ll_5,s5

s6=fltarr(1024) & for iy=0,1023 do begin s6[iy]=median(Spectra6[10:14,iy])
calibrate_OH,s6,88,ll_6
plot,ll_6,s6

s7=fltarr(1024) & for iy=0,1023 do begin s7[iy]=median(Spectra7[10:14,iy])
calibrate_OH,s7,455,ll_7
plot,ll_7,s7

s8=fltarr(1024) & for iy=0,1023 do begin s8[iy]=median(Spectra8[1:Xsize[8]-2,iy])
calibrate_OH,s8,77,ll_8
plot,ll_8,s8

s9=fltarr(1024) & for iy=0,1023 do begin s9[iy]=median(Spectra9[1:Xsize[9]-2,iy])
calibrate_OH,s9,70,ll_9
plot,ll_9,s9

s10=fltarr(1024) & for iy=0,1023 do begin s10[iy]=median(Spectra10[1:Xsize[10]-2,iy])
calibrate_OH,s10,77,ll_10
plot,ll_10,s10

s11=fltarr(1024) & for iy=0,1023 do begin s11[iy]=median(Spectra11[1:Xsize[11]-2,iy])
calibrate_OH,s11,447,ll_11
plot,ll_11,s11

s12=fltarr(1024) & for iy=0,1023 do begin s12[iy]=median(Spectra12[1:Xsize[12]-2,iy])
calibrate_OH,s12,208,ll_12
plot,ll_12,s12

s13=fltarr(1024) & for iy=0,1023 do begin s13[iy]=median(Spectra13[1:Xsize[13]-2,iy])
calibrate_OH,s13,290,ll_13
plot,ll_13,s13

s14=fltarr(1024) & for iy=0,1023 do begin s14[iy]=median(Spectra14[1:Xsize[14]-2,iy])
calibrate_OH,s14,357,ll_14
plot,ll_14,s14
;this completes the wavelength calibration for the sources

;source extraction
;s0: not clear wherer is the source
s0med=fltarr(Xsize[0],1024) & for i=0,1023 do begin s0med[*,i]=s0[i]
atv,spectra0-s0med<10>(-10)
;something at s0=9
final0=fltarr(1024)
for i=0,1023 do begin final0[i]=total(spectra0[8:10,i]-spectra0[11:13,i])/3.
plot,ll,final0

;s1
s1med=fltarr(Xsize[1],1024) & for i=0,1023 do begin s1med[*,i]=s1[i]
atv,spectra1-s1med<10>(-10)
;source at s0=11
;final1=fltarr(1024)
;for i=0,1023 do begin final1[i]=total(spectra1[10:12,i]-spectra1[2:4,i])/3.
S=fltarr(1024) & for i=0,1023 do begin S[i]=total(spectra1[10:12,i])/3.
B=fltarr(1024) & for i=0,1023 do begin B[i]=total(spectra1[2:4,i])/3.
iplot,ll_1,S-B
SP1=S-B

;s2
;s2med=fltarr(Xsize[1],1024) & for i=0,1023 do begin s2med[*,i]=s1[i]
;atv,spectra2-s2med<10>(-10)
V2=fltarr(Xsize[2]) & for i=0,Xsize[2]-1 do begin V2[i]=total(spectra2[i,*])
;source at s0=5; sky at 12
;final2=fltarr(1024)
S=fltarr(1024) & for i=0,1023 do begin S[i]=total(spectra2[4:6,i])/3.
B=fltarr(1024) & for i=0,1023 do begin B[i]=total(spectra2[11:13,i])/3.
iplot,ll_2,S-B


;5) WORK ON CALIBRATION OF THE Ne LAMP
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/Orion_N4_39.fit',a39,ha39
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/Orion_N4_40.fit',a40,ha40
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/Orion_N4_41.fit',a41,ha41
Neon=(a39+a40+a41)/3.
;... and its corresponding Dark
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/Orion_N4_42.fit',a42,ha42
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/Orion_N4_43.fit',a43,ha43
Neon_dark=(a42+a43)/2.

lines=reverse(rotate(Neon-Neon_dark,3),2)
lines=rot(lines,1.3)

ALLlines=fltarr(N_targets,max(xsize),1024)
i=0 ;for i=0,n_targets-1 do begin
    i_spec = lines[x0_target[i]:x1_Target[i],*]
    i_flat = flat[x0_target[i]:x1_Target[i],*]
    hist=histogram(i_flat)
    bins=min(i_flat)+FINDGEN(n_elements(hist))
    foo = MAX(hist,ind)
    mode = bins[ind]
  ; Normalize the flat field
    nflat = FLOAT(i_flat/mode)
    ALLlines[i,0:xsize[i],*] = i_spec/nflat
;endfor

;5) AND FOR THE STANDARD, Helias40335
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_7.fit',e7,he7
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_8.fit',e8,he8
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_9.fit',e9,he9
Elias=(e7+e8+e9)/3.
;... and its corresponding Dark
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_10.fit',e10,he10
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_11.fit',e11,he12
Elias_dark=(e10+e11)/2.
;... and its corresponding FLAT DARK
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_12.fit',e12,he12
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_13.fit',e13,he13
Elias_darklamp=(e12+e13)/2.
;... and its corresponding FLAT LAMP
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_14.fit',e14,he14
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_15.fit',e15,he15
Elias_lamp=(e14+e15)/2.


;find the flat normalization factor
flat=reverse(rotate(Elias_lamp-Elias_darklamp,3),2)
rot_flat=rot(flat,1.3,/interp)
mflat=median(rot_flat[520:537,*])

;standard
Elias=(Elias-Elias_dark)/(Elias_lamp-Elias_darklamp)*mflat
im=reverse(rotate(Elias,3),2)
rot_Elias=rot(im,1.3,/interp)

;;Check the rotation angle
;profile=fltarr(1024)
;for i=0,1023 do begin profile[i]=total(rot_elias[i,*])/1024. & endfor

OH=fltarr(1024)
for i=0,1023 do begin oh[i]=mean(rot_elias[535:537,I]) & endfor

PLOT,OH
;cursor,x,y,1,/data & print,x,y
OHlines= [1.97665,2.00025,2.02709,2.0337,2.0406,2.04961,2.05561,2.07153,2.14988,2.1793,2.19492]
OHpix  = 1024-[ 132.75, 216.15, 310.17,333.43,357.67, 388.68, 410.97, 469.14, 740.54, 844.26, 896.60]
plot,OHlines,OHpix,psym=2
h=linfit(ohlines,ohpix)
print,h
;      7814.04     -3502.80
oplot,ohlines,h[1]*ohlines+h[0]
wl=1.95+indgen(300)/1000.
plot,wl,h[1]*wl+h[0]
;conversion pix=>wl
pix=indgen(1024)
ll_std=(pix-h[0])/h[1]
plot,pix,ll_std
openps,'c:/Documents and Settings/robberto/my documents/projects/MacKenty-IRMOS/run6/OH_K3000_blue.ps'
plot,ll_std,OH,xtitle='wavelength (micron)',title='IRMOS, OH airglow on Kblue R=3000'
for i=0,10 do begin oplot,[OHlines[i],OHlines[i]],[-50,-10],linestyle=2
for i=0,10 do begin xyouts,2.22,80-i*8,STRING(OHlines[i])
closeps


;
;here is the spectrum of the standard star
Standard=fltarr(1024)
for i=0,1023 do begin standard[i]=mean(rot_elias[524:526,I]) & endfor
;
;plot the standard spectrum subtracting the oh airglow
plot,ll_std,standard
ill1=where(ABS(ll_std - 1.95) EQ MIN(ABS(ll_std - 1.95)))
ill2=where(ABS(ll_std - 2.20) EQ MIN(ABS(ll_std - 2.20)))
xl_std=ll_std(ill2:ill1)
Sp_Std=Standard(ill2:ill1)
plot,xl_std,Sp_std
BBFLUX=PLANCK_MIC(xl_std,10000.)
;oplot,xl_std,BBFLUX*Sp_STD(100)/BBFLUX(100)
plot,xl_std,SP_Std/(BBFLUX*Sp_STD(100)/BBFLUX(100))
;
;display one spectrum:
ill1_1=where(ABS(ll_1 - 1.95) EQ MIN(ABS(ll_1 - 1.95)))
ill1_2=where(ABS(ll_1 - 2.20) EQ MIN(ABS(ll_1 - 2.20)))
SP_Source1=SP1[ill1_2:ill1_1]
R=congrid(SP_Source1,NORM(ill2-ill1-1))
plot,R/Sp_Std


;CHECK THE NEON LAMP
;... Neon lamp
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_16.fit',e16,he16
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_17.fit',e17,he17
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_18.fit',e18,he18
Elias_Neonlamp=(e16+e17+e18)/3.
;... and its corresponding DARK LAMP
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_19.fit',e19,he19
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night4/HD_40335_20.fit',e20,he20
Elias_Neondark=(e19+e20)/2.
;Neon
Neon=(Elias_Neonlamp-Elias_Neondark)/(Elias_lamp-Elias_darklamp)*mflat
im=reverse(rotate(Neon,3),2)
rot_Neon=rot(im,1.3,/interp)
NeSpectrum=fltarr(1024)
for i=0,1023 do begin NeSpectrum[i]=mean(rot_Neon[520:537,I]) & endfor
openps,'c:/Documents and Settings/robberto/my documents/projects/MacKenty-IRMOS/run6/Neon_K3000_blue.ps'
plot,ll_std,Nespectrum,xtitle='wavelength (micron)',title='IRMOS, Ne lamp on Kblue R=3000'
oplot,[2.17081,2.17081],[-1000,-100],linestyle=2
oplot,[2.10413,2.10413],[-1000,-100],linestyle=2
oplot,[2.03502,2.03502],[-1000,-100],linestyle=2
oplot,[1.95771,1.95771],[-1000,-100],linestyle=2
xyouts,2.17081,-200,'2.17081'
xyouts,2.10413,-200,'2.10413'
xyouts,2.03502,-200,'2.03502'
xyouts,1.95771,-200,'1.95771'
closeps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



OPLOT,allSPECTRA[0,2:20,*]

end


pro cc
;sky0 : 230
sky0=fltarr(1024)
for i=0,1023 do begin sky0[i]=median(spectra[223:238,i])
;for i=0,1023 do begin resistant_mean,spectra[223:238,i],3,sky0[i]
for i=0,1023 do begin sky0[i]=mean(spectra[223:238,i])

;source 1: 283
source1=fltarr(1024) & sky1=fltarr(1024)
for i=0,1023 do begin source1[i]=mean(spectra[282:284,i])
for i=0,1023 do begin sky1[i]=mean(spectra[274:276,i])
spectrum1=source1-shift(sky1,0)
;ibad=where(spectrum1 LE -2000)
;spectrum1(ibad)=(spectrum1(ibad-1)+spectrum1(ibad+1))/2.

;source 2: 297
source2=fltarr(1024) & sky2=fltarr(1024)
for i=0,1023 do begin source2[i]=total(spectra[296:298,i])
for i=0,1023 do begin sky2[i]=total(spectra[304:306,i])
spectrum2=source2-sky2
ibad=where(spectrum1 LE -2000)
spectrum1(ibad)=(spectrum1(ibad-1)+spectrum1(ibad+1))/2.


end

pro test1


a12=reverse(a12,2)

a1r=reverse(ABS(rotate(a1,3)),2)
a2r=reverse(ABS(rotate(a2,3)),2)
cntrd,a1r,725,646,xref0,yref0
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/gridB1.fit',b1,hb1
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/gridB2.fit',b2,hb2
b12=ABS(rotate(b1-b2,3))
b12=reverse(b12,2)
b1r=reverse(ABS(rotate(b1,3)),2)
b2r=reverse(ABS(rotate(b2,3)),2)
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/gridB3.fit',b3,hb3
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/gridB4.fit',b4,hb4
b34=ABS(rotate(b3-b4,3))
b34=reverse(b34,2)

fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_50.fit',a1,ha1
a1=ABS(rotate(a1,3))
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_60.fit',a10,ha10
a5=ABS(rotate(a5,3))
a15=a1-a5
a15r=rot(a15,1.3)
;
spec=fltarr(1024)
sky=fltarr(1024)
for i=0,1023 do begin spec[i]=mean(r[443:445,i])
for i=0,1023 do begin sky[i]=mean(r[438:440,i])


fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_52.fit',a2,ha1
a2=ABS(rotate(a2,3))
a25=a2-a5
a25r=rot(a25,1.3)
;
spec=fltarr(1024)
sky=fltarr(1024)
for i=0,1023 do begin spec[i]=mean(r[443:445,i])
for i=0,1023 do begin sky[i]=mean(r[438:440,i])

fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_56.fit',a6,ha6
a6=ABS(rotate(a6,3))
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_57.fit',a7,ha7
a7=ABS(rotate(a7,3))
fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_58.fit',a8,ha8
a8=ABS(rotate(a8,3))

a68=a6-a8
a68r=rot(a68,1.3)
;
spec=fltarr(1024)
sky=fltarr(1024)
for i=0,1023 do begin spec[i]=mean(r[365:366,i])
for i=0,1023 do begin sky[i]=mean(r[356:358,i])

for i=0,1023 do begin spec[i]=mean(r[578:580,i])
for i=0,1023 do begin sky[i]=mean(r[569:571,i])

fits_read,'C:\Documents and Settings\robberto\My Documents\Projects\Mackenty-IRMOS\Run6/night2/M45_52.fit',a2,ha1
a2=ABS(rotate(a2,3))
a25=a2-a5
a25r=rot(a25,1.3)
;
spec=fltarr(1024)
sky=fltarr(1024)
for i=0,1023 do begin spec[i]=mean(r[443:445,i])
for i=0,1023 do begin sky[i]=mean(r[438:440,i])
end