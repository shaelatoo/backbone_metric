pro show_density2D,file, ri,nmap,crname=crname
;+ read and show 3D density map at radial distabce, ri, from tomography reconstruction 
; using Cartesian grid, 128^3  
; ri in [2.2, 4.0] Rsun for 3-view reconstruction
 if n_elements(file) le 0 then begin 
   print,"please give parameters, e.g."
   print,"show_density2d,'Np-3view-d2f2-wt_cr2123_sz128_mu100t0.01.dat',2.5,nmap,crname='CR 2123'"
   return
 endif
 ntomo=read_nsolu_cart(file)
 if ri lt 2.2 or ri gt 4.0 then begin
   print,'Warning---: ri must be 2.2<ri<4.0 for 3-view reconstructions!'
 endif
 if n_elements(ri) le 0 then ri=3.0
 model_ne_synop, ntomo, ri, nmap, quiet=quiet,crname=crname
end
