function read_nsolu_cart,filef90
 nsz=0
 openr,lun,filef90,/get_lun,/f77_unformatted
 readu,lun,nsz
 nsolu=fltarr(nsz,nsz,nsz)
 readu,lun,nsolu
 free_lun,lun
 n3d=nsolu*1e7 ;cm^-3
 return,n3d
end

