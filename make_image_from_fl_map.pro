;USAGE
; .r /usr/local/cmod/codes/spectroscopy/ir/FLIR/make_image_from_fl_map.pro

restore,file='Xpt_view_fieldline_map_1150724011.sav',/verb
n_images=n_elements(fl_map(*).fl_shot)
fl_image=intarr(64,64,n_images)
smooth_param=7
;plot,indgen(10),xra=[0,64],yra=[0,64],/nodata,col=1,xst=1,ysty=1
for i=0,n_images-1 do begin
   loadct,45,/sil
   window,1,xsize=500,ysize=500
   plot,indgen(10),xra=[0,64],yra=[0,64],/nodata,col=1,xst=1,ysty=1
   plots,(fl_map(i).xpixfl)[0:fl_map(i).n_fl_seg-1],(fl_map(i).ypixfl)[0:fl_map(i).n_fl_seg-1],col=2
   spline_p,(fl_map(i).xpixfl)[0:fl_map(i).n_fl_seg-1],(fl_map(i).ypixfl)[0:fl_map(i).n_fl_seg-1],xr,yr,interval=1.
   plots,xr,yr,psym=2
;   stop
;endfor
loadct,3,/sil
;for i=0,n_images-1 do begin
;   n_pts=fl_map(i).n_fl_seg
n_pts=n_elements(xr)
   for j=0,n_pts-1 do begin
      for k=0,63 do begin
         for l=0,63 do begin
;            a=abs(k-(fl_map(i).xpixfl)[j]) & b=abs(l-(fl_map(i).ypixfl)[j])
            a=abs(k-xr(j)) & b=abs(l-yr(j))
;            if (a le 2. and  b le 2.) then fl_image(k,l,i)=1000./(sqrt(a^2+b^2) > 1.)
            if sqrt(a^2+b^2) le 2 then fl_image(k,l,i)=(1000./(sqrt(a^2+b^2) > 1.) < 1000.)
         endfor
      endfor
   endfor
   window,0,xsize=320,ysize=320
;   tvscl,rebin(fl_image(*,*,i),5*64,5*64) & xyouts,4,300,/dev,sval(fl_map(i).fl_R*100.,l=5)+','+sval(fl_map(i).fl_Z*100,l=6),charsiz=1.5,col=255 
   fl_image(*,*,i)=smooth(fl_image(*,*,i),smooth_param,/edge_trunc)
;   stop
   tvscl,rebin(fl_image(*,*,i),5*64,5*64) & xyouts,4,300,/dev,sval(fl_map(i).fl_R*100.,l=5)+','+sval(fl_map(i).fl_Z*100,l=6),charsiz=1.5,col=255 
;   stop
endfor
for i=0,n_images-1 do begin
   tvscl,rebin(fl_image(*,*,i),5*64,5*64) 
   xyouts,4,300,/dev,sval(fl_map(i).fl_R*100.,l=5)+','+sval(fl_map(i).fl_Z*100,l=6),charsiz=1.5,col=255 
   wait,0.3
endfor
shot=fl_map(0).fl_shot
time=fl_map(0).fl_time
fieldline_R=fl_map(*).fl_R
fieldline_Z=fl_map(*).fl_Z
fieldline_phi_start=fl_map(0).phipin
save,file='fl_images_'+strtrim(fl_map(0).fl_shot,2)+'.sav',shot,time,smooth_param,fieldline_R,fieldline_Z,fieldline_phi_start,fl_image
stop
end
