pro cmp_ne_synop, ntomo, Ri, dmap,nmap,quiet=quiet,iw=iw,prf=prf,rhi=rhi,$
    den=den, mpg=mpg,Rin=Rin,Rout=Rout, crname=crname
;+purpose: create a synopic map of Ne at solar radius, Rin,
; from the Ne cube in Cartesian coordinate (x,y,z) with z
; along solar rotation axis, and xy for solar equatoral plane,
; Carrington longitude is related to positive x-axis increasing
; counter-clockwise.
; Input:
; dens(128,128,128) - Ne density cube in unit of cm^-3
; Ri - scalar, given radial distance for the synopic map
;     in units of Rsun
; Output:
; nemap(360,181) - synopic map of Ne at Ri
;   I-index for Carrington longtitude (0-359) deg
;   J-index for Carrington latitude (-90 to 90) deg
;-
 
if not keyword_set(ri) then Ri=2.5
if not keyword_set(crname) then crname='CR 2091'
 
sz=size(ntomo)  
 ntm=ntomo*1e7  ;cm^-3
 dir='./'
 if not keyword_set(den) then begin
    mdl=ntomo*1e7 ;verted into cm^-3 
    cnam='Same '
 endif else begin
    mdl=den*1e7  ;cm^-3
    cnam='Model:'
 endelse
ss=size(mdl)
n=ss[3]
if ss[3] gt sz[3] then begin
  mdl=rebin(mdl,sz[1],sz[2],sz[3])
  n=sz[3]
endif

 ;compare radial profiles averaged globally
 if keyword_set(prf) or keyword_set(mpg) then begin
  loadct,0
;  r0=1.5
   dr=0.1
;  nr=26
 if not keyword_set(Rin) then Rin=1.5
 if not keyword_set(Rout) then Rout=4.0
  nr=fix((Rout-Rin)/dr)+1
   rhi=findgen(nr)*dr+Rin
   prf=fltarr(nr,2)
   mov=bytarr(600,600,nr)
  ;img=bytarr(600,300)
   for i=0,nr-1 do begin
        model_ne_synop, mdl,   Rhi[i], nmap1, iw=0,name=cnam,/sm,crname=crname
        mov[0:599,300:*,i]=tvrd() 
        model_ne_synop, ntm,   Rhi[i], nmap2, iw=1,name='Tomography: ',/sm, crname=crname
        mov[0:599,0:299,i]=tvrd()
        prf[i,0]=average(nmap1>0)
        prf[i,1]=average(nmap2>0)
   endfor

 window,3,retain=2
 !p.multi=0
 linecolors
 m2=max(prf,min=m1)
 plot,rhi,prf[*,1],yr=[m1,m2],/yst,/ylog,thick=1,$
   charsize=1.4,xtit='Distance (Rsun)',ytit='Density (cm^-3)'
 oplot,rhi,prf[*,0],color=2
 legend2,['My Tomography',cnam],line=[0,0],colors=[255,2],/right,charsize=1.5
 mne1=string(max(prf[*,1]),format='(e10.3)')+' (cm!e-3!n)'
 mne0=string(max(prf[*,0]),format='(e10.3)')+' (cm!e-3!n)'
 xyouts,0.5,0.8,/norm,'Tomo N_max='+mne1,charsize=1.5
 xyouts,0.5,0.75,/norm,'Model N_max='+mne0,color=2,charsize=1.5

 loadct,0
 endif

 model_ne_synop, mdl, Ri, dmap, iw=0,name=cnam,/sm,crname=crname
 linecolors
if keyword_set(prf) then oplot,[0,360],[0,0],color=2
 model_ne_synop, ntm, Ri, nmap, iw=1,name='Tomography: ',/sm,crname=crname
 linecolors
if keyword_set(prf) then oplot,[0,360],[0,0],color=2

 window,4,retain=2
 device,decomposed=0
 m2=max([nmap[*,90],dmap[*,90]],min=m1)
 plot,nmap[*,90],yr=[m1,m2]
 linecolors
 oplot,dmap[*,90],color=2
 legend2,['My Tomo',cnam],line=[0,0],colors=[255,2],/right,charsize=1.5

 if keyword_set(mpg) then begin
    ss=size(mov)
    n=ss[3]
    loadct,0
    XINTERANIMATE, SET=[ss[1],ss[2],ss[3]], /SHOWLOAD
    FOR I=0,n-1 DO XINTERANIMATE, FRAME = I, IMAGE = mov[*,*,I]
    XINTERANIMATE

  oVid = IDLffVideoWrite('mov_Ne_synop.mp4', FORMAT='mp4')
  nx=ss[1]
  ny=ss[2]
  nfrm=ss[3]
  fps=12
  vidStream = oVid.AddVideoStream(nx, ny, fps)
  ;loadct,colortable
  image24 = BytArr(3, nx,ny)
  TVLCT, r, g, b, /Get
  xi=indgen(256)
   FOR indx=0, nfrm-1 DO BEGIN
         image24[0, 0:nx-1,0:ny-1] = r[mov[*,*,indx]]
         image24[1, 0:nx-1,0:ny-1] = g[mov[*,*,indx]]
         image24[2, 0:nx-1,0:ny-1] = b[mov[*,*,indx]]
      !NULL = oVid.Put(vidStream, image24)
   ENDFOR
 endif
end
