@/usr/local/cmod/codes/spectroscopy/ir/FLIR/get_ref_frame.pro
@/home/labombard/edge/modelling/geometry/xyz2rzphi.pro
@/usr/local/cmod/codes/spectroscopy/ir/FLIR/oplot_periscope_view.pro
@/usr/local/cmod/codes/spectroscopy/ir/FLIR/get_cmod_continuous_div_surface.pro
@/home/terry/idl_lib/ps.pro
@/home/terry/kn1d/sol_lc.pro
pro create_periscope_view,shot=shot,xbox,ybox,nseg,lseg,Xseg,Yseg,ColorSeg,LabelSeg,IndexSeg,Rpixel,Zpixel,Phipixel,$
                          cmod=cmod,plot=plot,showrays=showrays,zrot=zrot,xshift=xshift,yshift=yshift,$
                          Rpin=Rpin,Zpin=Zpin,phipin=phipin,show_roi=show_roi,mag=mag,rzphiplot=rzphiplot,$
                          windowmag=windowmag,camera_view=camera_view,roi=roi,set_roi_slot=set_roi_slot,$
                          view_ramped_tiles=view_ramped_tiles,show_field_line=show_field_line,proj_onto_pixels=proj_onto_pixels,$
                          ps=ps,locate_pt=locate_pt,continuous_div=continuous_div,shadow1=shadow1,shadow2=shadow2

common tiletile_outlinetile_outline_outline,_xbox,_ybox,_nseg,_lseg,_Xseg,_Yseg,_ColorSeg,_LabelSeg,_IndexSeg,_Rpixel,_Zpixel,_Phipixel,_roi,_cmod,_ramp_version
; USAGE
; .r /usr/local/cmod/codes/spectroscopy/ir/FLIR/create_periscope_view_Xpt.pro

; for the DIV camera view use Rpin=1.05,Zpin=-0.165, alpha=25.,beta=-15.,f_o_v_ang=38., zrot=-159
; create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-159.,yshift=.01,xshift=0.48,mag=0.795
; create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-159.,yshift=.01,xshift=0.48,mag=0.795,/locate_pt
; create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-159.,yshift=-.01,xshift=0.48,mag=0.795,camera_view='DIV'

; use Rpin=1.05,Zpin=0.164, alpha=25.,beta=-35.,f_o_v_ang=35., zrot=-155
; create_periscope_view,Rpin=1.05,Zpin=0.165,zrot=-135.,yshift=.28,xshift=0.013,mag=1.08

; use Rpin=1.05,Zpin=0.0, alpha=25.,beta=-25.,f_o_v_ang=35., zrot=-145
; create_periscope_view,Rpin=1.05,Zpin=0.0,zrot=-145.,yshift=.28,xshift=0.013,mag=1.08

; for a restricted f-o-v of 26 degrees use Rpin=1.05,Zpin=-0.165, alpha=23.,beta=-19.,f_o_v_ang=26., zrot=-148
; create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-148.,yshift=.01,xshift=0.48,mag=0.795,/plot,/continuous_div

; for the 14 degree f-o-v use Rpin=1.05,Zpin=-0.165, alpha=32.,beta=-19.,f_o_v_ang=14., zrot=-162
; create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-162.,yshift=.01,xshift=0.48,mag=0.795,/plot,/continuous_div
; create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-162.,yshift=.01,xshift=0.48,mag=0.795,/plot,/continuous_div,/locate_pt

; for the 50 degree f-o-v use Rpin=1.05,Zpin=-0.165, alpha=22.3,beta=-10.31.,f_o_v_ang=50., zrot=-162
;create_periscope_view,Rpin=1.05,Zpin=-0.165,zrot=-162.,yshift=.01,xshift=0.48,mag=0.795,/continuous_div,shadow1=70.,shadow2=45.


; to reproduce the 2012 IR camera view from the top of A-port, use Rpin=0.76,Zpin=0.44, alpha=-40.,beta=-45.,f_o_v_ang=20., zrot=15.
; create_periscope_view,Rpin=0.76,Zpin=0.44,zrot=15.0,yshift=-.007,xshift=0.132,mag=1.08,/plot,/continuous_div
; create_periscope_view,Rpin=0.76,Zpin=0.44,zrot=15.0,yshift=-.007,xshift=0.132,mag=1.08,/plot,/continuous_div,/locate_pt
; create_periscope_view,Rpin=0.76,Zpin=0.44,zrot=15.0,yshift=-.007,xshift=0.132,mag=1.08,camera_view='IR-A port',/view_ramped_tiles

; to create the view from H-port horizontal thru the 8" McPherson
; flange use  alpha=0.,beta=-7.75,f_o_v_ang=15.
; create_periscope_view,Rpin=3.65,Zpin=0.3045,zrot=-90.0,shadow1=45,shadow2=45.

; to reproduce the view of the actual lower X-pt view (2015)
; flange use  alpha=-33.,beta=-10.5,f_o_v_ang=13.
; create_periscope_view,Rpin=1.00,Zpin=-0.265,phipin=10.,zrot=-9.5,/plot,camera_view='X-pt'
; create_periscope_view,Rpin=1.00,Zpin=-0.265,phipin=10.,zrot=-9.5,/locate_pt
; create_periscope_view,Rpin=1.00,Zpin=-0.265,phipin=10.,zrot=-9.5,/show_field_line,camera_view='X-pt',/plot,/proj_onto_pixels

;
; KEYWORDS
; shot - specifying this sets the geometry of the ramped tiles to the
;        date of shot
; cmod - show the stick figure of the C-Mod PFCs
; show_background_tiles - not used
; view_ramped_tiles - load the J divertor ramped tiles if they
;                     will be in the periscope view
; camera_view - 'DIV' for the 2008-2012 A-port camera view or 'IR-A
;               port' for the IR view of J port from A port or
;               'X-PT' for the 2015 telescope view from the 
; showrays - if set, then rays from each ZYZseg point to the XYZ_eye
;            point are drawn (used for debugging only) 
; rzphiplot - plot the rzphi grids in the image plane
; zrot - in degrees
; xshift - in pixels assuming a xpix x ypix sized pixel image (used
;          when aligning with an actual camera image)
; yshift - in pixels assuming a xpix x ypix sized pixel image (used
;          when aligning with an actual camera image)
; locate_pt - if set, then you can use the cursor function to choose
;             up to 10 points in the image and determine the point-locations,
;             as well as angles and lengths of rays to those points
; continuous_div - if set, then the geometry of the 2014 continuous
;                  divertor (as designed 11/2012) is used in the view determination  
; mag - magnification of the stick plot 
; windowmag - keep at the default value (1) since otherwise it screws up the ROI
;             determinations
; show_field_line - shows the projections of specified fieldlines
;                   (set in lines 318-321)within the field-of-view
; roi - [not implemented]
; show_roi - [not implemented]
; set_roi_slot - [not implemented]



key_default,shot,1120904020L

if shot gt 1100601001 then ramp_version=2 else ramp_version=1
if shot eq -1 then ramp_version=2 

if ramp_version eq 1 then begin
  key_default,zrot,15.7
  key_default,xshift,0.013
  key_default,yshift,0.28
endif
if ramp_version eq 2 then begin
; use the 2010 reference frame 
   key_default,zrot,15.7
   key_default,xshift,-0.114
   key_default,yshift,0.092
endif
if ramp_version eq 1 then begin
;use the 2011 reference frame 
   key_default,zrot,15.7
   key_default,xshift,0.10
   key_default,yshift,0.166
endif
key_default,Rpin,1.05
key_default,Zpin,-0.165
key_default,phipin,0.
key_default,zrot,0
key_default,xshift,0.0
key_default,yshift,0.
key_default,plot,0
key_default,rzphiplot,0
key_default,showrays,0
key_default,cmod,1
key_default,mag,1.08
key_default,windowmag,1.0
key_default,camera_view,' '
key_default,show_roi,0
key_default,set_roi_slot,0
key_default,view_ramped_tiles,0
key_default,show_field_line,0
;key_default,psfile,''
key_default,ps,0
key_default,locate_pt,0
key_default,continuous_div,0
key_default,shadow1,180.
key_default,shadow2,180.
key_default,ps,0

xpix=480
ypix=640
if camera_view eq 'IR-A port' then begin
   xpix=320
   ypix=256
endif
if locate_pt then plot=1
port_names=['A','B','C','D','E','F','G','H','J','K']
loadct,45,/sil
;if strlen(psfile) gt 0 then ps=1 else ps=0
if ps then begin
;  ps,file=psfile,/psfont
  pso
endif
;print,'zrot,xshift,yshift=',zrot,xshift,yshift
;
; Decide if the common block is already loaded
;
if (not plot) and (not show_roi) and (not set_roi_slot) then begin
  if type_of(_roi) ne 0 then begin
    if _cmod eq cmod and _ramp_version eq ramp_version then begin
;      xbox=_xbox & ybox=_ybox & nseg=_nseg & lseg=_lseg & Xseg=_Xseg & Yseg=_Yseg & ColorSeg=_ColorSeg & LabelSeg=_LabelSeg & IndexSeg=_IndexSeg
;      Rpixel=_Rpixel & Zpixel=_Zpixel & Phipixel=_Phipixel & roi=_roi
;      return
    endif
  endif
endif

p_pos=!p.position
if plot then begin
  if not ps then x,title='C-Mod Vacuum Vessel',mag=1.5,window=1
  create_XYZ_Coord,[-1.5,1.5],[-1.5,1.5],[-.7,.7],xrot=60,zrot=170,ps=ps
endif
;
; Optionally read and plot C-Mod Tiles
;
if cmod and not continuous_div then Get_CMod_Tile_Surface,nseg_CMOD,lseg_CMOD,XYZseg_CMOD,ColorSeg_CMOD,LabelSeg_CMOD,IndexSeg_CMOD,plot=plot
if cmod and continuous_div then Get_CMod_continuous_div_Surface,nseg_CMOD,lseg_CMOD,XYZseg_CMOD,ColorSeg_CMOD,LabelSeg_CMOD,IndexSeg_CMOD,plot=plot
; XYZseg_CMOD is a [3,120,714] array of XYZ coordinates of tile
; corners such that XYZseg_CMOD(3,i,j) are the XYZ coordinate vectors,
; nseg is the number of coordinate vectors,
; lseg(n) is the number of elements describing the path of the nth segment
; colorseg (intarr(nseg)) is the plot color of that segment
; labelseg=strarr(nseg)  &; text label that segment
; indexseg=intarr(nseg)  &; index that segment 

;
; Read and optionally plot J-divertor tile gaps and hardware
;
if view_ramped_tiles then begin
;   Get_J_Divertor_Surface,ramp_version=ramp_version,nseg,lseg,XYZseg,ColorSeg,LabelSeg,IndexSeg,RZPhiseg=RZPhiseg,plot=plot
   Get_J_Divertor_Surface,ramp_version=ramp_version,nseg_CMOD,lseg_CMOD,XYZseg_CMOD,ColorSeg_CMOD,LabelSeg_CMOD,IndexSeg_CMOD,RZPhiseg=RZPhiseg,plot=plot
   cmod=0
;
; A-port periscope is roughly at R=.76, Z=+.44, phi=0 degrees, i.e., on x=0 plane
;
;   rzphi2xyz,0.62913,-0.41641,270,xyz_spot
endif else begin
   cmod=1
endelse
;   print, 'enter coordinates [R(m),Z(m)] of the A port periscope pinhole '
;   read,Rpin,Zpin
rzphi2xyz,Rpin,Zpin,phipin,xyz_eye
print, 'enter the viewing angles (degs) at the pinhole (alpha (rel to the R_maj of the pinhole), beta (rel to horizontal), and the full angle of the f-o-v'
read, alpha, beta,f_o_v_angle
;
; search the XYZseg_CMOD array for a segment point that makes the line
; between that segment point and the viewing point (xyz_eye) have
; angles close to the input viewing angles (alpha & beta) 
for i=0,nseg_CMOD-1 do begin
   for j=0,lseg_cmod(i)-1 do begin
       v=(XYZseg_CMOD(*,j,i)-xyz_eye)
;       if sqrt(total(v*v)) gt 0.1 and sqrt(total(v(0:1)*v(0:1))) gt 0.1 then begin
       if sqrt(total(v*v)) gt 0.1 and sqrt(total(v(0:1)*v(0:1))) gt 0.1 and sqrt(total(v*v)) lt 1.5*xyz_eye(1) then begin 
; note that the last condition is to keep the xyz_spot on the same
; side of the tokamak as the view
;       if sqrt(total(v*v)) gt 0.1 and sqrt(total(v(0:1)*v(0:1))) gt 0.1 and sqrt(total(v*v)) lt 1.4 then begin 
;       if sqrt(total(v*v)) gt 0.1 and sqrt(total(v(0:1)*v(0:1))) gt 0.1 then begin
           alpha_test=sign(XYZseg_CMOD(0,j,i))*acos(total(-xyz_eye(0:1)*v(0:1))/sqrt(total(v(0:1)*v(0:1)))/sqrt(total(xyz_eye(0:1)*xyz_eye(0:1))) < 1.)/!dtor
           beta_test=asin(v(2)/sqrt(total(v*v)) < 1.)/!dtor
;           if (abs(alpha-alpha_test) lt 5. and abs(beta-beta_test) lt 5.) then begin
;               print,XYZseg_CMOD(*,j,i),i,j,alpha_test,beta_test
               if (abs(alpha-alpha_test) lt 1. and abs(beta-beta_test) lt .3) then begin
                   xyz_spot=XYZseg_CMOD(*,j,i)
; Center of view is to be the point that approximately aligns with the desired alpha and beta
                   goto,get_out
               endif
;           endif
       endif
   endfor
endfor
; if we are here, then we were not able to find the vector from
; xyz_eye to a segment point that aligned with the view as defined by
; xyz_eye, alpha, and beta.
print, ' unable to find a line from xyz_eye to a segment point that aligns with alpha,beta - aborting'
goto, finish
get_out:
index1=j
index2=i
v_sight=xyz_spot-xyz_eye
v_sight_length=sqrt(total(v_sight*v_sight))
; see how close the xyz_spot is to the desired alpha and beta
alpha_act=asin(v_sight(0)/(sqrt(total(v_sight(0:1)*v_sight(0:1)))))/!dtor
beta_act=asin(v_sight(2)/v_sight_length)/!dtor
; calculate the angel between the central ray (v_sight) and the axis
; of the endoscope ([0,Rpin,Zpin]-[0,0,Zpin]
;v_endo=-[0,Rpin,Zpin]+[0,0,Zpin]
v_endo=-xyz_eye+[0,0,xyz_eye(2)]
gamma=acos(total(v_endo*v_sight)/v_sight_length/sqrt(total(v_endo*v_endo)) < 1.)
print,'angle between central ray and the endoscope axis is ',gamma/!dtor
;stop
;
; Draw line from center of ramped tiles to periscope aperture
;
xyz_sightline=[[xyz_spot],[xyz_eye]]
if plot then plots,xyz_sightline(0,*),xyz_sightline(1,*),xyz_sightline(2,*) ,/t3d,/data,color=2,thi=3
;
; If CMod and view_ramped_tiles, then concatenate arrays
;
if cmod and view_ramped_tiles then begin
  _XYZseg=fltarr(3,max([lseg_CMOD,lseg]),nseg_CMOD+nseg)
  _XYZseg(*,0:max(lseg_cmod)-1,0:nseg_CMOD-1)=XYZseg_CMOD
  _XYZseg(*,0:max(lseg)-1,nseg_CMOD:nseg_CMod+nseg-1)=XYZseg
  XYZSeg=_XYZseg
  _RZPhiseg=fltarr(3,max([lseg_CMOD,lseg]),nseg_CMOD+nseg)
  _RZPhiseg(*,0:max(lseg)-1,nseg_CMOD:nseg_CMod+nseg-1)=RZPhiseg
  RZPhiSeg=_RZPhiseg
  nseg=nseg_CMOD+nseg
  lseg=[lseg_CMOD,lseg]
  ColorSeg=[ColorSeg_CMOD,ColorSeg]
  LabelSeg=[LabelSeg_CMOD,LabelSeg]
  IndexSeg=[IndexSeg_CMOD,IndexSeg]
endif else begin
   XYZSeg=XYZseg_CMOD
   nseg=nseg_CMOD
   lseg=lseg_CMOD
   ColorSeg=ColorSeg_CMOD
   LabelSeg=LabelSeg_CMOD
   IndexSeg=IndexSeg_CMOD
endelse
if camera_view eq 'X-pt' then begin
   phi_cal=316.5
   cal_r=[0.51,0.583,0.5825,0.5095]
   cal_z=[-0.355,-0.357,-0.431,-0.429]
   for i=0,n_elements(cal_r)-1 do begin
      SOL_LC,1150717011,1.0,reform(cal_R(i)),reform(cal_Z(i)),0.5*!pi,0.01,LC,error,phipos=phi_cal,plot=0,/limiters,cr=cr,cz=cz,cphi=cphi
; now make segments of this field line in the same format as the
; machine feature segments. Remember:
; XYZseg_CMOD is a [3,120,714] array of XYZ coordinates of tile
; corners such that XYZseg_CMOD(3,i,j) are the XYZ coordinate vectors,
; nseg is the number of coordinate vectors,
; lseg(n) is the number of elements describing the path of the nth segment
; colorseg (intarr(nseg)) is the plot color of that segment
; labelseg=strarr(nseg)  &; text label that segment
; indexseg=intarr(nseg)  &; index that segment
      _CPhi=CPhi(0)+(CPhi(n_elements(CPhi)-1)-CPhi(0))*findgen(1000)/999.0
      _CR=(interpol(CR,CPhi,_CPhi))[0:1] ; plot only the first two points on the fieldline
      _CZ=(interpol(CZ,CPhi,_CPhi))[0:1] ; plot only the first two points on the fieldline
      _CPhi=_Cphi(0:1)                   ; plot only the first two points on the fieldline
      rzphi2xyz,_CR,_CZ,_Cphi/!dtor,xyzseg_cal
      nseg_cal=1
      lseg_cal=n_elements(_CPhi)
      colorseg_cal=intarr(nseg_cal)+2 ; assign them a red color
      indexseg_cal=nseg+indgen(nseg_cal)
      labelseg_cal=strarr(nseg_cal)+'X-pt view corner'
      _XYZseg=fltarr(3,max([lseg,lseg_cal]),nseg+nseg_cal)
      _XYZseg(*,0:max(lseg)-1,0:nseg-1)=XYZseg
      _XYZseg(*,0:max(lseg_cal)-1,nseg)=XYZseg_cal
      XYZSeg=_XYZseg
      nseg=nseg+nseg_cal
      lseg=[lseg,lseg_cal]
      ColorSeg=[ColorSeg,ColorSeg_cal]
      LabelSeg=[LabelSeg,LabelSeg_cal]
      IndexSeg=[IndexSeg,IndexSeg_cal]
      if plot then plots,xyzseg_cal(0,*),xyzseg_cal(1,*),xyzseg_cal(2,*),/t3d,/data,color=3,thi=3
   endfor
endif
if show_field_line then begin
   fl_shot=1L
;;   print,' enter shot #, time, R [m], and Z [m] of the field-line launch point'
;;   read, fl_shot,fl_time,fl_R,fl_Z
   fl_shot=1150717011
   fl_time=1.0
;   fl_r=[.53, .53, .60, .60] ;[0.52,0.57,0.57,0.565]
   n_fl_R=20
   n_fl_Z=20
   fl_R=0.53+0.07*(indgen(n_fl_R)/19.)
   fl_Z=-0.479+0.133*(indgen(n_fl_Z)/19.)
;   fl_z=[-.479, -.346, -.479, -.346] ;[-0.365,-0.367,-0.421,-0.47]
   for i=0,n_elements(fl_z)-1 do begin
      for ii=0,n_elements(fl_r)-1 do begin
         SOL_LC,fl_shot,fl_time,reform(fl_R(ii)),reform(fl_Z(i)),0.5*!pi,0.01,LC,error,phipos=phipin,plot=0,/limiters,cr=cr,cz=cz,cphi=cphi
; now make segments of this field line in the same format as the
; machine feature segments. Remember:
; XYZseg_CMOD is a [3,120,714] array of XYZ coordinates of tile
; corners such that XYZseg_CMOD(3,i,j) are the XYZ coordinate vectors,
; nseg is the number of coordinate vectors,
; lseg(n) is the number of elements describing the path of the nth segment
; colorseg (intarr(nseg)) is the plot color of that segment
; labelseg=strarr(nseg)  &; text label that segment
; indexseg=intarr(nseg)  &; index that segment
         if cphi(1) gt cphi(0) then begin
            SOL_LC,fl_shot,fl_time,reform(fl_R(ii)),$
                   reform(fl_Z(i)),-0.5*!pi,0.01,LC,error,phipos=phipin,$
                   plot=0,/limiters,cr=cr,cz=cz,cphi=cphi
;            fl_sign=-1.
            if cphi(1) gt cphi(0) then begin
;               print,'neither sol_lc direction gives a fieldline
;               going in the correct direction, so skipping
;               ',fl_R(ii),fl_Z(i)
               print, 'at R,Z='+sval(fl_R(ii),l=5)+','+sval(fl_z(i),l=6)+$
                      ' neither sol_lc direction gives a fieldline going in the correct direction, assuming fieldline is toroidal'
               _CPhi=CPhi(0)+(((alpha/abs(alpha))*!pi/1.5)-CPhi(0))*findgen(100)/99.0
               _CR=fltarr(100)+fl_R(ii)
               _CZ=fltarr(100)+fl_z(i)
               goto,jump_fl    
            endif
         endif
;         print,fl_R(ii),fl_Z(i)
         _CPhi=CPhi(0)+(CPhi(n_elements(CPhi)-1)-CPhi(0))*findgen(1000)/999.0
         _CPhi=_Cphi(0:locate(_Cphi,_Cphi(0)+(alpha/abs(alpha))*!pi/1.5))
         _CR=interpol(CR,CPhi,_CPhi)
         _CZ=interpol(CZ,CPhi,_CPhi)
         jump_fl:
         rzphi2xyz,_CR,_CZ,_Cphi/!dtor,xyzseg_fl
         nseg_fl=1
         lseg_fl=n_elements(_CPhi)
         colorseg_fl=intarr(nseg_fl)+4
         indexseg_fl=nseg+indgen(nseg_fl)
         labelseg_fl=strarr(nseg_fl)+'fieldline'+strtrim(i*n_elements(fl_r)+ii+1,2)
         _XYZseg=fltarr(3,max([lseg,lseg_fl]),nseg+nseg_fl)
         _XYZseg(*,0:max(lseg)-1,0:nseg-1)=XYZseg
         _XYZseg(*,0:max(lseg_fl)-1,nseg)=XYZseg_fl
         XYZSeg=_XYZseg
         nseg=nseg+nseg_fl
         lseg=[lseg,lseg_fl]
         ColorSeg=[ColorSeg,ColorSeg_fl]
         LabelSeg=[LabelSeg,LabelSeg_fl]
         IndexSeg=[IndexSeg,IndexSeg_fl]
         if plot then plots,xyzseg_fl(0,*),xyzseg_fl(1,*),xyzseg_fl(2,*),/t3d,/data,color=4,thi=3
;         jump_fl:
      endfor
   endfor
endif
;stop 
;
; find those segments in the XYZseg array that are within the desired cone of the view, ie
; those vectors [XYZseg-xyz_eye] whose dot product with
; [xyz_spot-xyz_eye] yields an angle le 0.5 of f_o_v_angle
;

nseg_sub=0
XYZseg_sub=fltarr(3,n_elements(XYZseg(0,*,0)),nseg) & ;XYZ coordinate vectors
lseg_sub=intarr(nseg)
colorseg_sub=intarr(nseg)
labelseg_sub=strarr(nseg)
indexseg_sub=intarr(nseg)
;stop
for i=0,nseg-1 do begin
   lcount=0
   for j=0,lseg(i)-1 do begin
      v=(XYZseg(*,j,i)-xyz_eye)
      ang=acos((total(v*v_sight)/v_sight_length/sqrt(total(v*v))) <1.)/!dtor
      if ang le f_o_v_angle/2. then begin
; If you are at this point, then this segment is to be included with
; those in the prescibed f-o-v
; now try to eliminate those segments that are shadowed by other
; segments with this view
          xyz2rzphi,XYZseg(*,j,i),a,b,c
          xyz2rzphi,xyz_spot,aa,bb,c_spot
; c is in degrees
          jt1=min([abs(c_spot-c),abs(c_spot+360.-c)],njt)
          if ((jt1 le shadow1 and njt eq 0) or (jt1 le shadow2 and njt eq 1)) then begin 
              XYZseg_sub(*,lcount,nseg_sub)=XYZseg(*,j,i)
              lcount=lcount+1
;              if j eq 0 then print,i
          endif
      endif
   endfor
;   if lcount gt 1 then begin
   if lcount ge 1 then begin
      lseg_sub(nseg_sub)=lcount
      colorseg_sub(nseg_sub)=colorseg(i)
      labelseg_sub(nseg_sub)=labelseg(i)
      indexseg_sub(nseg_sub)=indexseg(i)
;      print,nseg_sub, '  ',labelseg(i),'  ', colorseg(i)
      nseg_sub=nseg_sub+1
   endif
endfor
;
; trim the sub arrays
;
;nseg_sub=nseg_sub-1
XYZseg_sub=XYZseg_sub(*,0:max(lseg_sub)-1,0:nseg_sub-1)
; make the RZPhiseg_sub array and
; find the indices in the sub array that coorespond to the center of
; the view using the indices of the center from the large array 
RZPhiseg_sub=XYZseg_sub
for i=0,nseg_sub-1 do begin
   for j=0,lseg_sub(i)-1 do begin
      xyz2rzphi,XYZseg_sub(*,j,i),a,b,c
      RZPhiseg_sub(0,j,i)=a
      RZPhiseg_sub(1,j,i)=b
      RZPhiseg_sub(2,j,i)=c
      if XYZseg_sub(0,j,i) eq XYZseg(0,index1,index2) and XYZseg_sub(1,j,i) eq XYZseg(1,index1,index2) and $
XYZseg_sub(2,j,i) eq XYZseg(2,index1,index2) then begin
         spot_index1=j
         spot_index2=i
;         print,' the spot indices for the sub array are ',spot_index1,spot_index2
      endif
   endfor
endfor
lseg_sub=lseg_sub(0:nseg_sub-1)
ColorSeg_sub=ColorSeg_sub(0:nseg_sub-1)
;stop
;jtt=where(RZPhiseg_sub(0,0,*) gt 0.48,njtt)
;if njtt gt 0 then ColorSeg_sub(jtt)=2
jt2=where((RZPhiseg_sub(2,0,*) MOD 36) eq 0 and (RZPhiseg_sub(2,1,*) MOD 36) eq 0,njt)
if njt gt 0 then ColorSeg_sub(jt2)=1+nint((RZPhiseg_sub(2,0,jt2) MOD 360)/36)
LabelSeg_sub=LabelSeg_sub(0:nseg_sub-1)
IndexSeg_sub=IndexSeg_sub(0:nseg_sub-1)

;stop
;
; The image plane should be located a focal distance away from xyz_spot, along the sight line
;
fd=0.3 &;focal distance
xyz_slope=xyz_sightline(*,1)-xyz_sightline(*,0)
distance=sqrt(total((xyz_sightline(*,1)-xyz_sightline(*,0))^2))
image_plane_origin=xyz_eye-fd*xyz_slope/distance
; the minus sign moves the image (along the xyz_sightline) to a point IN FRONT
; OF xyz_eye, since this will keep the pinhole image of the view
; upright and not inverted/mirrored
; if you want to move the image plane behind the pinhole then comment
; out the next line
;image_plane_origin=xyz_eye+fd*xyz_slope/distance

;
; Build an image plane in a local x-y space.
;
scale=0.05/mag
dx=1.0 & dy=ypix*dx/xpix
xbox=[-dx,dx,dx,-dx]-xshift & ybox=[-dy,-dy,dy,dy]-yshift & z_image_plane=[0,0,0,0]
;
; Define the box in a local 3D space
;
_xyz_image_plane=fltarr(3,4)  
_xyz_image_plane(0,*)=xbox*scale
_xyz_image_plane(1,*)=ybox*scale
_xyz_image_plane(2,*)=z_image_plane
;
; Define a set of local coordinate axes to align the image plane with the view.
;  Select the local Z axis to be along the sightline
Zaxis=-xyz_spot+xyz_eye
;
;  Select the local X axis to be aligned with major radius coordinate to the image plane origin
Xaxis=image_plane_origin
Xaxis(*)=Xaxis(*)-[0,0,image_plane_origin(2)]
;
; Plot the desired local Z and X axes in C-Mod space
;
;if plot then begin
;  plots,[image_plane_origin(0),image_plane_origin(0)+Zaxis(0)],[image_plane_origin(1),image_plane_origin(1)+Zaxis(1)],$
;        [image_plane_origin(2),image_plane_origin(2)+Zaxis(2)],/t3d,/data,color=1
;  plots,[image_plane_origin(0),image_plane_origin(0)+Xaxis(0)],[image_plane_origin(1),image_plane_origin(1)+Xaxis(1)],$
;        [image_plane_origin(2),image_plane_origin(2)+Xaxis(2)],/t3d,/data,color=3
;
;  press_return
;endif
;
; Build a coordinate transformation centered on the image plane origin and aligned with new Z and X axes. Make Z-axis alignment exact.
;  Use it to transform the local x-y image plane to C-Mod coordinates.

xyz_image_plane=rotate_XYZ_axes(_xyz_image_plane,-1,XYZorigin=image_plane_origin,Xaxis=Xaxis,Zaxis=Zaxis,fixed='Z',$
                 xprime=xprime,yprime=yprime,zprime=zprime) 
;
; Plot xprime,yprime,zprime coordinate axes
;
if plot then begin
  plots,[image_plane_origin(0),image_plane_origin(0)+xprime(0)],[image_plane_origin(1),image_plane_origin(1)+xprime(1)],$
        [image_plane_origin(2),image_plane_origin(2)+xprime(2)],/t3d,/data,color=3
  plots,[image_plane_origin(0),image_plane_origin(0)+yprime(0)],[image_plane_origin(1),image_plane_origin(1)+yprime(1)],$
        [image_plane_origin(2),image_plane_origin(2)+yprime(2)],/t3d,/data,color=4
  plots,[image_plane_origin(0),image_plane_origin(0)+zprime(0)],[image_plane_origin(1),image_plane_origin(1)+zprime(1)],$
        [image_plane_origin(2),image_plane_origin(2)+zprime(2)],/t3d,/data,color=1
;
; Plot Image plane in C-Mod space
;
  plots,[reform(xyz_image_plane(0,*)),xyz_image_plane(0,0)],[reform(xyz_image_plane(1,*)),xyz_image_plane(1,0)],$
    [reform(xyz_image_plane(2,*)),xyz_image_plane(2,0)],/t3d,/data,color=4
;  stop
  press_return
endif
;
; Trace rays from XYZseg_sub (the gaps in the tiles) to XYZ_eye, intersecting with the image plane
;
;xyz_image=XYZseg
xyz_image=XYZseg_sub
for n=0,nseg_sub-1 do begin
  xyzRay1=reform(XYZseg_sub(*,0:lseg_sub(n)-1,n))
  xyzRay2=xyzRay1
  for m=0,n_elements(xyzRay2(0,*))-1 do xyzRay2(*,m)=xyz_eye
   ray_intersect,xyzRay1,xyzRay2,xyz_image_plane,xyzPoint,/allow_negative_rays
   xyz_image(*,0:lseg_sub(n)-1,n)=xyzPoint
   if showrays then begin
;
; Optional: plot rays
;
    for p=0,n_elements(xyzRay2(0,*))-1 do begin
      plots,[xyzRay1(0,p),xyzRay2(0,p)],[xyzRay1(1,p),xyzRay2(1,p)],[xyzRay1(2,p),xyzRay2(2,p)],/t3d,/data,color=4
    endfor
    press_return
  endif
endfor
if showrays then press_return
;
; Plot image on 3D image plane
;   
if plot then begin
  for n=0,nseg_sub-1 do begin
    color=colorseg_sub(n)
    plots,xyz_image(0,0:lseg_sub(n)-1,n),xyz_image(1,0:lseg_sub(n)-1,n),xyz_image(2,0:lseg_sub(n)-1,n),/t3d,/data,color=color
  endfor
  press_return
endif
;
; Now apply the forward transformation to map the polygons drawn on the 3D image plane to local x-y space
;
image_2D=XYZ_image
for n=0,nseg_sub-1 do begin
   XYZ_2D=rotate_XYZ_axes(XYZ_image(*,0:lseg_sub(n)-1,n),1)  
; Optionally rotate view about Z-axis
   image_2D(*,0:lseg_sub(n)-1,n)=rotate_XYZ(XYZ_2D,[0.0,0.0,0.0],zrot=zrot*!dtor)  
endfor
if plot then begin
   window,0,xsize=700,ysize=700,title='field-of-view in image plane'
   pos=[0.1,0.1,0.95,0.95]
;   plot,[min(image_2D(0,*,*)),max(image_2D(0,*,*))],[min(image_2D(1,*,*)),max(image_2D(1,*,*))],col=1,/nodata,pos=pos
   plot,[min(image_2D(1,*,*)),max(image_2D(1,*,*))],[min(image_2D(1,*,*)),max(image_2D(1,*,*))],col=1,/nodata,pos=pos
   for n=0,nseg_sub-1 do plots,image_2D(0,0:lseg_sub(n)-1,n),image_2D(1,0:lseg_sub(n)-1,n),col=colorseg_sub(n)
;   for nn=0,njt-1 do plots,image_2D(0,0:lseg_sub(jt2(nn))-1,jt2(nn)),image_2D(1,0:lseg_sub(jt2(nn))-1,jt2(nn)),col=colorseg_sub(jt2(nn)),thi=2
   xyouts,.15,.9,/norm,'(R,Z) of aperture=('+sval(rpin,l=5)+','+sval(zpin,l=5)+') m',col=1,charsiz=1.
   xyouts,.15,.85,/norm,'alpha, beta angles for view center='+sval(alpha_act,l=5)+','+sval(beta_act,l=6)+' deg',col=1,charsiz=1.
   xyouts,.15,.8,/norm,'full angle of f-o-v='+sval(f_o_v_angle,l=3)+' deg',col=1,charsiz=1.
;   for nn=0,njt-1 do xyouts,.85,.95-colorseg_sub(jt2(nn))*.05,/norm,port_names(colorseg_sub(jt2(nn))-1)+'-port',charsiz=1.,col=colorseg_sub(jt2(nn))
      plots,[image_2D(0,spot_index1,spot_index2)-0.002,image_2D(0,spot_index1,spot_index2)+0.002],[image_2D(1,spot_index1,spot_index2),image_2D(1,spot_index1,spot_index2)],col=1,thi=5
      plots,[image_2D(0,spot_index1,spot_index2),image_2D(0,spot_index1,spot_index2)],[image_2D(1,spot_index1,spot_index2)-0.002,image_2D(1,spot_index1,spot_index2)+0.002],col=1,thi=5
   if locate_pt then begin
      pt_alpha=fltarr(10)
      pt_beta=fltarr(10)
      dist_to_pin=fltarr(10)
      npts=1
      next:
      print,'click on image point whose distances/coordinates you want [out of plot box to move on]'
      cursor,a,b,/data
      wait,0.2
;      print,a,b
      if a lt !x.crange(0) or a gt !x.crange(1) or b lt !y.crange(0) or b gt !y.crange(1) then goto,done_with_pts
; now find image_2d point closest to the click image point
      dum=min(abs(image_2D(0,*,*)-a)+abs(image_2D(1,*,*)-b),ndum)
      ind = ARRAY_INDICES(reform(image_2D(0,*,*)),ndum)
      print,image_2D(0,ind(0),ind(1)),image_2D(1,ind(0),ind(1))
      plots,[image_2D(0,ind(0),ind(1))-0.005,image_2D(0,ind(0),ind(1))+0.005],[image_2D(1,ind(0),ind(1)),image_2D(1,ind(0),ind(1))],col=npts+1,thi=3
      plots,[image_2D(0,ind(0),ind(1)),image_2D(0,ind(0),ind(1))],[image_2D(1,ind(0),ind(1))-0.005,image_2D(1,ind(0),ind(1))+0.005],col=npts+1,thi=3
      v_to_pt=xyzseg_sub(*,ind(0),ind(1))-xyz_eye
      dist_to_pin(npts-1)=sqrt(total(v_to_pt*v_to_pt))
      pt_alpha(npts-1)=asin(v_to_pt(0)/(sqrt(total(v_to_pt(0:1)*v_to_pt(0:1)))))/!dtor
; this is the angle that the ray from center of pinhole to the test
; point makes with the vertical plane that contains the central rayis
; defined by phi=0
      pt_beta(npts-1)=asin(v_to_pt(2)/dist_to_pin(npts-1))/!dtor
; this is the angle that the ray from center of pinhole to the test
; point makes with the horizontal plane that contains the pinhole
      xyouts,.3,0.1+npts*0.03,/norm,'(R,Z,phi)!d'+sval(npts,l=1)+'!N=('+sval(rzphiseg_sub(0,ind(0),ind(1)),l=5)+','+sval(rzphiseg_sub(1,ind(0),ind(1)),l=5)+','+sval(rzphiseg_sub(2,ind(0),ind(1)),l=4)+')',col=npts+1,align=0.5
      xyouts,.7,0.1+npts*0.03,/norm,'(alpha,beta,length)!d'+sval(npts,l=1)+'!N=('+sval(pt_alpha(npts-1),l=5)+','+sval(pt_beta(npts-1),l=5)+','+sval(dist_to_pin(npts-1),l=4)+')',col=npts+1,align=0.5
      npts=npts+1
      goto,next
   endif
   done_with_pts:
   press_return
endif

;
; Simplify to x and y vector segments
;
xseg=reform(image_2d(0,*,*))/scale
yseg=reform(image_2d(1,*,*))/scale
zseg=reform(image_2d(2,*,*))/scale
if not ps then window,0,xsize=600,ysize=600,title='specified field-of-view'
pos=[0.1,0.1,0.95,0.1+(!d.x_vsize/!d.y_vsize)*0.85]
plot,[min(xseg),max(xseg)],[min(yseg),max(yseg)],col=1,/nodata,pos=pos
for n=0,nseg_sub-1 do plots,xseg(0:lseg_sub(n)-1,n),yseg(0:lseg_sub(n)-1,n),noclip=0,col=colorseg_sub(n)
;for nn=0,njt-1 do plots,xseg(0:lseg_sub(jt2(nn))-1,jt2(nn)),yseg(0:lseg_sub(jt2(nn))-1,jt2(nn)),col=colorseg_sub(jt2(nn)),thi=2
; above plots the segments in planes at phis at the centers of ports
xyouts,.15,.9,/norm,'R,Z of aperture='+sval(rpin,l=5)+','+sval(zpin,l=5)+' m',col=1,charsiz=1.
xyouts,.15,.85,/norm,'alpha, beta angles for view center='+sval(alpha_act,l=6)+','+sval(beta_act,l=6)+' deg',col=1,charsiz=1.
xyouts,.15,.8,/norm,'full angle of f-o-v='+sval(f_o_v_angle,l=3)+' deg',col=1,charsiz=1.
;for nn=0,njt-1 do xyouts,.15,.05+colorseg_sub(jt2(nn))*.05,/norm,port_names(colorseg_sub(jt2(nn))-1)+'-port',charsiz=1.,col=colorseg_sub(jt2(nn))
press_return
if proj_onto_pixels and show_field_line then begin
; define the seg point at the corners of the X-pt camera image
   corn_ind=where(labelseg_sub eq 'X-pt view corner',mm)
   ll=0
   tl_corn=[xseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll)),yseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll))]
   ll=1
   tr_corn=[xseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll)),yseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll))]
   ll=2
   br_corn=[xseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll)),yseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll))]
   ll=3
   bl_corn=[xseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll)),yseg(0:lseg_sub(corn_ind(ll))-2,corn_ind(ll))]
; now find fieldline segments thst project into the camera image
   fl_seg_ind=where(strmid(labelseg_sub,0,9) eq 'fieldline',mm)
   if mm gt 0 then begin
      xfl=fltarr(max(lseg_sub(fl_seg_ind)),mm)
      yfl=xfl
      xpixfl=xfl
      ypixfl=xfl
      fl_label=strarr(mm)
      n_fl_seg=intarr(mm)
      fl_count=0
      for mmm=0,mm-1 do begin
; apply transformation from xseg,yseg] space into pixel space - make
; the bl corner the (0,0) pixel and the tr corner the (63,63) pixel
;   xpix=nint((xseg_test-bl_corn(0))/(tr_corn(0)-bl_corn(0))*63.)
;   ypix=nint((yseg_test-bl_corn(1))/(tr_corn(1)-bl_corn(1))*63.)
         xtest=(xseg(0:lseg_sub(fl_seg_ind(mmm))-1,fl_seg_ind(mmm))-bl_corn(0))/(tr_corn(0)-bl_corn(0))*63.
         ytest=(yseg(0:lseg_sub(fl_seg_ind(mmm))-1,fl_seg_ind(mmm))-bl_corn(0))/(tr_corn(0)-bl_corn(0))*63.
         lll=where((abs(xtest-32.) lt 36.) and (abs(ytest-32.) lt 36.),n_lll)
         if n_lll gt 5 then begin
            
            n_fl_seg(fl_count)=n_lll
            xpixfl(0:n_fl_seg(fl_count)-1,fl_count)=xtest(lll)
            ypixfl(0:n_fl_seg(fl_count)-1,fl_count)=ytest(lll)
            fl_label(fl_count)=labelseg_sub(fl_seg_ind(mmm))
            xfl(0:n_fl_seg(fl_count)-1,fl_count)=xseg(lll,fl_seg_ind(mmm))
            yfl(0:n_fl_seg(fl_count)-1,fl_count)=yseg(lll,fl_seg_ind(mmm))
            fl_count=fl_count+1
         endif
      endfor
   endif 
;   stop
; now create structure and write structure to saveset
   if fl_count gt 0 then begin
;      rz_ind=intarr(fl_count)
      r_ind=intarr(fl_count)
      z_ind=intarr(fl_count)
      for mmmm=0,fl_count-1 do begin
         z_ind(mmmm)=(int(strtrim(strmid(fl_label(mmmm),9,5),2))-1)/n_fl_R
         r_ind(mmmm)=(int(strtrim(strmid(fl_label(mmmm),9,5),2))-1) MOD n_fl_Z
;         rz_ind(mmmm)=int(strtrim(strmid(fl_label(mmmm),9,5),2))-1
         if mmmm eq 0 then begin
            _fl_map={fl_label:fl_label(0), fl_shot:fl_shot, fl_time:fl_time,$
               fl_R:fl_R(r_ind(0)), fl_Z:fl_Z(z_ind(0)), phipin:phipin,$
               n_fl_seg:n_fl_seg(0),xpixfl:xpixfl(*,0),ypixfl:ypixfl(*,0)}
            fl_map=replicate(_fl_map,fl_count)            
         endif else begin
            fl_map(mmmm).fl_label=fl_label(mmmm)
            fl_map(mmmm).fl_R=fl_R(r_ind(mmmm))
            fl_map(mmmm).fl_Z=fl_Z(z_ind(mmmm))
            fl_map(mmmm).n_fl_seg=n_fl_seg(mmmm)
            fl_map(mmmm).xpixfl=xpixfl(*,mmmm)
            fl_map(mmmm).ypixfl=ypixfl(*,mmmm)
         endelse
      endfor
      save,file='Xpt_view_fieldline_map_'+strtrim(fl_shot,2)+'.sav',fl_map,/verb
      print,'wrote the save set for Xpt_view_fieldline_map_'+strtrim(fl_shot,2)
   endif
stop
endif
;
; Repack Xseg, Yseg, RZPhiseg arrays
;
Xseg_data=[0.0]
Yseg_data=[0.0]
Rseg_data=[0.0]
Zseg_data=[0.0]
Phiseg_data=[0.0]
RZPhiseg_sub=XYZseg_sub*0.
for n=0,nseg_sub-1 do begin
;
; If looking at the J instrumented divertor, then the Segments with negative indexes are used to interpolate R,Z,Phi coordinates
;
   if view_ramped_tiles and indexseg(n) lt 0 then begin
      Xseg_data=[Xseg_data,xseg(0:lseg_sub(n)-1,n)]
      Yseg_data=[Yseg_data,yseg(0:lseg_sub(n)-1,n)]
      Rseg_data=[Rseg_data,reform(RZPhiseg_sub(0,0:lseg_sub(n)-1,n))]
      Zseg_data=[Zseg_data,reform(RZPhiseg_sub(1,0:lseg_sub(n)-1,n))]
      Phiseg_data=[Phiseg_data,reform(RZPhiseg_sub(2,0:lseg_sub(n)-1,n))]
   endif else begin
      Xseg_data=[Xseg_data,xseg(0:lseg_sub(n)-1,n)]
      Yseg_data=[Yseg_data,yseg(0:lseg_sub(n)-1,n)]
      Rseg_data=[Rseg_data,reform(RZPhiseg_sub(0,0:lseg_sub(n)-1,n))]
      Zseg_data=[Zseg_data,reform(RZPhiseg_sub(1,0:lseg_sub(n)-1,n))]
      Phiseg_data=[Phiseg_data,reform(RZPhiseg_sub(2,0:lseg_sub(n)-1,n))]
   endelse
endfor
Xseg_data=Xseg_data(1:*)
Yseg_data=Yseg_data(1:*)
Rseg_data=Rseg_data(1:*)
Zseg_data=Zseg_data(1:*)
Phiseg_data=Phiseg_data(1:*)
;
  xsize=xpix*Windowmag
  ysize=dy*xsize/dx
;
; Construct Delaunay triangles for intperpolation
;
TRIANGULATE, Xseg_data, Yseg_data, Triangles 
;
; Interpolate for R, Z and Phi at pixel locations
;
xpixelwidth=(max(xbox)-min(xbox))/xsize
xpixelcenter=min(xbox)+0.5*xpixelwidth+findgen(xsize)*xpixelwidth

ypixelwidth=(max(ybox)-min(ybox))/ysize
ypixelcenter=min(ybox)+0.5*ypixelwidth+findgen(ysize)*ypixelwidth
Rpixel = TRIGRID( Xseg_data, Yseg_data, Rseg_data, Triangles, xout=xpixelcenter, yout=ypixelcenter)
Zpixel = TRIGRID( Xseg_data, Yseg_data, Zseg_data, Triangles, xout=xpixelcenter, yout=ypixelcenter)
Phipixel = TRIGRID( Xseg_data, Yseg_data, Phiseg_data, Triangles, xout=xpixelcenter, yout=ypixelcenter)

if rzphiplot then begin
    if not ps then x,title='Endoscope View - R,Z Phi',window=3,xsize=xsize,ysize=ysize
    p_pos=!p.position
    !p.position=[0,0,1,1]
    plot,xbox,ybox,xstyle=5,ystyle=5,/nodata
    for n=0,nseg_sub-1 do begin
      if indexseg(n) lt 0 then oplot,xseg(0:lseg_sub(n)-1,n),yseg(0:lseg_sub(n)-1,n),color=colorseg_sub(n)
    endfor
    ii=where(Rpixel gt 0)
    levels=min(rpixel(ii))+(max(rpixel)-min(rpixel(ii)))*findgen(100)/100
    contour, Rpixel, xpixelcenter, ypixelcenter,title='Major Radius',xtitle='X (m)',ytitle='Y (m)' ,levels=levels,/overplot
    press_return

    plot,xbox,ybox,xstyle=5,ystyle=5,/nodata
    for n=0,nseg_sub-1 do begin
      if indexseg(n) lt 0 then oplot,xseg(0:lseg_sub(n)-1,n),yseg(0:lseg_sub(n)-1,n),color=colorseg_sub(n)
    endfor
    ii=where(Zpixel lt 0)
    levels=min(Zpixel(ii))+(max(Zpixel)-min(Zpixel(ii)))*findgen(50)/50
    contour, Zpixel, xpixelcenter, ypixelcenter,title='Z location',xtitle='X (m)',ytitle='Y (m)' ,levels=levels,/overplot
    press_return

    plot,xbox,ybox,xstyle=5,ystyle=5,/nodata
    for n=0,nseg_sub-1 do begin
      if indexseg(n) lt 0 then oplot,xseg(0:lseg_sub(n)-1,n),yseg(0:lseg_sub(n)-1,n),color=colorseg_sub(n)
    endfor
    ii=where(Phipixel gt 0)
    levels=min(Phipixel(ii))+(max(Phipixel)-min(Phipixel(ii)))*findgen(50)/50
    contour, Phipixel, xpixelcenter, ypixelcenter,title='Toroidal Angle',xtitle='X (m)',ytitle='Y (m)' ,levels=levels,/overplot
    press_return
 endif ; RZphiplot
!p.position=p_pos


;
; Setup ROI structure array
;
;   Tags are: {name:name, active:active, npoly:npoly, xrange:xrange, yrange:yrange, xpoly:xpoly, ypoly:ypoly, $
;              nindices:nindices, indices:indices, minimize:minimize}
;
;   Each region of interest is defined two different ways - by a closed polygon (pixel indices are contained in xpoly, ypoly),
;   and by an array of pixel indices.
;
;                name - string, name of ROI ('ramped tiles', 'divertor slot', 'tile42',...)
;                active - 0 = ignore this ROI
;                         1 = use this ROI to lock image
;                         2 = use this ROI to determine color scale of image displayed when using 'view_image' keyword option
;                npoly - number of active elements stored in xpoly(*), ypoly(*)
;                xrange - left, right x coordinates of image boundary in xpoly units (example: xrange=[-1.0,1.0])
;                yrange - bottom, top y coordinates of image boundary in ypoly units (example: yrange=[0,255])
;                xpoly(*) - x indices of polygon
;                ypoly(*) - y indices of polygon
;                nroi - number of active elements in indices(*)
;                nindices - numnber of elements in indices corresponding to region of interest
;                indices(*) - array of pixel indices in region of interest
;                minimize - 0 = maximize cross-correlation between image and this ROI
;                           1 = minimize cross-correlation between image and this ROI
  maxpoly=1000 & maxindices=long(xpix)*long(ypix)
;  maxpoly=1000 & maxindices=windowmag*320L*256L
  name='' & active=0 & npoly=0L & xrange=[0.0,1.0] & yrange=[0.0,1.0] & xpoly=fltarr(maxpoly) & ypoly=fltarr(maxpoly) 
  nindices= 0L & indices=lonarr(maxindices) & minimize=0 & nRZPhi=0L & _R=fltarr(maxpoly) & _Z=_R & _Phi=_R
  _ROI={name:name, active:active, npoly:npoly, xrange:xrange, yrange:yrange, xpoly:xpoly, ypoly:ypoly,$
        nindices:nindices, indices:indices, minimize:minimize, nRZPhi: nRZPhi, R:_R, Z:_Z, Phi:_Phi}
;  ROI=replicate(_ROI,200)
  ROI=replicate(_ROI,2)
  goto, skip_ROIs

  nROI=0
;
;
;==================================================
; ROI - Ramped Tiles
;==================================================;
; Compute xpoly and ypoly
;
  ROI(nROI).name='Ramped Tiles'
  ROI(nROI).xrange=[min(xbox),max(xbox)]
  ROI(nROI).yrange=[min(ybox),max(ybox)]
  xpoly=[0.0]
  ypoly=[0.0]
  _R=[0.0]
  _Z=[0.0]
  _Phi=[0.0]
  for n=0,nseg-1 do begin
    if (indexseg(n) eq 1020) then begin
      xpoly=[xpoly,xseg(0:lseg(n)-1,n)] 
      ypoly=[ypoly,yseg(0:lseg(n)-1,n)] 
      _R=[_R,reform(RZPhiseg(0,0:lseg(n)-1,n))] 
      _Z=[_Z,reform(RZPhiseg(1,0:lseg(n)-1,n))] 
      _Phi=[_Phi,reform(RZPhiseg(2,0:lseg(n)-1,n))] 
    endif
    if (indexseg(n) eq 1040) then begin
      xpoly=[xpoly,reverse(xseg(0:lseg(n)-1,n))] 
      ypoly=[ypoly,reverse(yseg(0:lseg(n)-1,n))] 
      _R=[_R,reverse(reform(RZPhiseg(0,0:lseg(n)-1,n)))] 
      _Z=[_Z,reverse(reform(RZPhiseg(1,0:lseg(n)-1,n)))] 
      _Phi=[_Phi,reverse(reform(RZPhiseg(2,0:lseg(n)-1,n)))] 
    endif
  endfor
  xpoly=xpoly(1:*)
  ypoly=ypoly(1:*)
  _R=_R(1:*)
  _Z=_Z(1:*)
  _Phi=_Phi(1:*)
  n=n_elements(xpoly)
  ROI(nROI).npoly=n
  ROI(nROI).xpoly(0:n-1)=xpoly
  ROI(nROI).ypoly(0:n-1)=ypoly
  ROI(nROI).nRZPhi=n
  ROI(nROI).R(0:n-1)=_R
  ROI(nROI).Z(0:n-1)=_Z
  ROI(nROI).Phi(0:n-1)=_Phi
;
; Determine pixel indices (Note: it is important to map data to the center of the pixels)
;
  xpixelwidth=(max(xbox)-min(xbox))/xpix
  xroi=xpix*(xpoly-min(xbox))/(max(xbox)-min(xbox))-0.5*xpixelwidth & xroi=xroi < (xpix-1) & xroi=xroi > 0
  ypixelwidth=(max(ybox)-min(ybox))/ypix
  yroi=ypix*(ypoly-min(ybox))/(max(ybox)-min(ybox))-0.5*ypixelwidth & yroi=yroi < (ypix-1) & yroi=yroi > 0

  xroi=nint(xroi) & yroi=nint(yroi)
  indices=POLYFILLV( xroi, yroi, xpix, ypix)
  n=n_elements(indices)
  ROI(nROI).nindices=n
  ROI(nROI).indices(0:n-1)=indices
  nROI=nROI+1

  Phi0_ramp=min(_Phi) & Z0_ramp=min(_Z)
  Phi1_ramp=max(_Phi) & Z1_ramp=max(_Z)
  Phi0_ramp=min(_Phi) & Z0_ramp=min(_Z)
  Phi1_ramp=max(_Phi) & Z1_ramp=max(_Z)
;
if ramp_version lt 2 then goto,skip_full_ramp_tiles
;==================================================
; ROI - Full Ramped Tiles - including ramp-down and plateau tiles
;==================================================;
; Compute xpoly and ypoly
;
  ROI(nROI).name='Full Ramped Tiles'
  ROI(nROI).xrange=[min(xbox),max(xbox)]
  ROI(nROI).yrange=[min(ybox),max(ybox)]
  xpoly=[0.0]
  ypoly=[0.0]
  _R=[0.0]
  _Z=[0.0]
  _Phi=[0.0]
  for n=0,nseg-1 do begin
    if (indexseg(n) eq 1015) then begin
      xpoly=[xpoly,xseg(0:lseg(n)-1,n)] 
      ypoly=[ypoly,yseg(0:lseg(n)-1,n)] 
      _R=[_R,reform(RZPhiseg(0,0:lseg(n)-1,n))] 
      _Z=[_Z,reform(RZPhiseg(1,0:lseg(n)-1,n))] 
      _Phi=[_Phi,reform(RZPhiseg(2,0:lseg(n)-1,n))] 
    endif
    if (indexseg(n) eq 1045) then begin
      xpoly=[xpoly,reverse(xseg(0:lseg(n)-1,n))] 
      ypoly=[ypoly,reverse(yseg(0:lseg(n)-1,n))] 
      _R=[_R,reverse(reform(RZPhiseg(0,0:lseg(n)-1,n)))] 
      _Z=[_Z,reverse(reform(RZPhiseg(1,0:lseg(n)-1,n)))] 
      _Phi=[_Phi,reverse(reform(RZPhiseg(2,0:lseg(n)-1,n)))] 
    endif
  endfor
  xpoly=xpoly(1:*)
  ypoly=ypoly(1:*)
  _R=_R(1:*)
  _Z=_Z(1:*)
  _Phi=_Phi(1:*)
  n=n_elements(xpoly)
  ROI(nROI).npoly=n
  ROI(nROI).xpoly(0:n-1)=xpoly
  ROI(nROI).ypoly(0:n-1)=ypoly
  ROI(nROI).nRZPhi=n
  ROI(nROI).R(0:n-1)=_R
  ROI(nROI).Z(0:n-1)=_Z
  ROI(nROI).Phi(0:n-1)=_Phi
;
; Determine pixel indices (Note: it is important to map data to the center of the pixels)
;
  xpixelwidth=(max(xbox)-min(xbox))/xpix
  xroi=xpix*(xpoly-min(xbox))/(max(xbox)-min(xbox))-0.5*xpixelwidth & xroi=xroi < (xpix-1) & xroi=xroi > 0
  ypixelwidth=(max(ybox)-min(ybox))/ypix
  yroi=ypix*(ypoly-min(ybox))/(max(ybox)-min(ybox))-0.5*ypixelwidth & yroi=yroi < (ypix-1) & yroi=yroi > 0

  xroi=nint(xroi) & yroi=nint(yroi)
  indices=POLYFILLV( xroi, yroi, xpix, ypix)
  n=n_elements(indices)
  ROI(nROI).nindices=n
  ROI(nROI).indices(0:n-1)=indices
  nROI=nROI+1

  Phi0_ramp=min(_Phi) & Z0_ramp=min(_Z)
  Phi1_ramp=max(_Phi) & Z1_ramp=max(_Z)
  Phi0_ramp=min(_Phi) & Z0_ramp=min(_Z)
  Phi1_ramp=max(_Phi) & Z1_ramp=max(_Z)
skip_full_ramp_tiles:
;
;==================================================
; ROI - Ramped tile region for surface normal computation
;         Make this a little bigger than the ramped tiles
;==================================================
;
  Phi0=Phi0_ramp-0.5 & Phi1=Phi1_ramp+0.5 & Z0=z0_ramp-0.01 & Z1=z1_ramp+0.01
  set=phipixel ge Phi0 and phipixel le Phi1 and zpixel ge Z0 and zpixel le Z1
  indices=where(set)

  roi(nROI).name='Ramped Tiles Surface Normals'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;==================================================
; ROI - Individual Ramped Tiles
;==================================================;
;
  Tile_list=['Ramped Tile:73','Ramped Tile:74','Ramped Tile:75','Ramped Tile:76','Ramped Tile:77/78A','Ramped Tile:77/78B',$
             'Ramped Tile:79A','Ramped Tile:79B','Ramped Tile:80A','Ramped Tile:80B']
  for ii=0,n_elements(Tile_list)-1 do begin  
    ROI(nROI).name=Tile_List(ii)
    ROI(nROI).xrange=[min(xbox),max(xbox)]
    ROI(nROI).yrange=[min(ybox),max(ybox)]
    found=0
    for n=0,nseg-1 do begin
      if (labelseg(n) eq Tile_List(ii)) then begin
        found=1
        xpoly=xseg(0:lseg(n)-1,n)
        ypoly=yseg(0:lseg(n)-1,n)
        _R=reform(RZPhiseg(0,0:lseg(n)-1,n))
        _Z=reform(RZPhiseg(1,0:lseg(n)-1,n)) 
        _Phi=reform(RZPhiseg(2,0:lseg(n)-1,n))
      endif
    endfor
    if not found then message,Tile_List(ii)+' not found!'
    n=n_elements(xpoly)
    ROI(nROI).npoly=n
    ROI(nROI).xpoly(0:n-1)=xpoly
    ROI(nROI).ypoly(0:n-1)=ypoly
    ROI(nROI).nRZPhi=n
    ROI(nROI).R(0:n-1)=_R
    ROI(nROI).Z(0:n-1)=_Z
    ROI(nROI).Phi(0:n-1)=_Phi
  ;
  ; Determine pixel indices (Note: it is important to map data to the center of the pixels)
  ;
    xpixelwidth=(max(xbox)-min(xbox))/xpix
    xroi=xpix*(xpoly-min(xbox))/(max(xbox)-min(xbox))-0.5*xpixelwidth & xroi=xroi < (xpix-1) & xroi=xroi > 0
    ypixelwidth=(max(ybox)-min(ybox))/ypix
    yroi=ypix*(ypoly-min(ybox))/(max(ybox)-min(ybox))-0.5*ypixelwidth & yroi=yroi < (ypix-1) & yroi=yroi > 0

    xroi=nint(xroi) & yroi=nint(yroi)
    indices=POLYFILLV( xroi, yroi, xpix, ypix)
    n=n_elements(indices)
    ROI(nROI).nindices=n
    ROI(nROI).indices(0:n-1)=indices
    nROI=nROI+1
endfor
;
;==================================================
; ROI - ODIV sensors
;==================================================;
;
iList=where(indexseg ge 3000 and indexseg lt 4000,count)
;if (count ne 34) and (count ne 41) then message,'Only '+sval(count)+' ODIV sensors found'
for jj=0,count-1 do begin
  ii=iList(jj)
  ROI(nROI).name=labelseg(ii)
  ROI(nROI).xrange=[min(xbox),max(xbox)]
  ROI(nROI).yrange=[min(ybox),max(ybox)]
  n=lseg(ii)
  ROI(nROI).npoly=n
  xpoly=xseg(0:lseg(ii)-1,ii)
  ypoly=yseg(0:lseg(ii)-1,ii)
  ROI(nROI).xpoly(0:n-1)=xpoly
  ROI(nROI).ypoly(0:n-1)=ypoly
  ROI(nROI).nRZPhi=n
  _R=reform(RZPhiseg(0,0:lseg(ii)-1,ii))
  _Z=reform(RZPhiseg(1,0:lseg(ii)-1,ii)) 
  _Phi=reform(RZPhiseg(2,0:lseg(ii)-1,ii))
  ROI(nROI).R(0:n-1)=_R
  ROI(nROI).Z(0:n-1)=_Z
  ROI(nROI).Phi(0:n-1)=_Phi
;
; Determine pixel indices (Note: it is important to map data to the center of the pixels)
;
  xpixelwidth=(max(xbox)-min(xbox))/xpix
  xroi=xpix*(xpoly-min(xbox))/(max(xbox)-min(xbox))-0.5*xpixelwidth & xroi=xroi < (xpix-1) & xroi=xroi > 0
  ypixelwidth=(max(ybox)-min(ybox))/ypix
  yroi=ypix*(ypoly-min(ybox))/(max(ybox)-min(ybox))-0.5*ypixelwidth & yroi=yroi < (ypix-1) & yroi=yroi > 0

  xroi=nint(xroi) & yroi=nint(yroi)
  indices=POLYFILLV( xroi, yroi, xpix, ypix)
  n=n_elements(indices)
  ROI(nROI).nindices=n
  ROI(nROI).indices(0:n-1)=indices
  nROI=nROI+1
endfor
;
;
;==================================================
; ROI - full array of individual Tiles
;==================================================
;
; Set up coordinates of tile boundaries
;
  dphi_tile=2.98
  phi_tile=findgen(13)*dphi_tile+270.0-5*dphi_tile
  phi_tile(12)=phi_tile(12)-.1
  z_tile=[-0.5739,-0.54613,-0.51764,-0.48907,-0.4613,-0.43711,-0.41728,-0.39647,-0.37536,-0.36215]
;
; Load active pixel indices and determine enclosing polygon line segments from active pixels
;
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  for jrow=1,9 do begin
    for icolumn=1,12 do begin
      phi0=phi_tile(icolumn-1) & phi1=phi_tile(icolumn)
      z0=z_tile(jrow-1) & z1=z_tile(jrow)
      active=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
      indices=where(active,nindices)
      if nindices gt 0 then begin
         roi(nROI).nindices=nindices & roi(nROI).indices(0:nindices-1)=indices
         roi(nROI).name='TILE:'+sval(jrow)+':'+sval(icolumn)
;
         image=bytarr(xpix,ypix) & image(indices)=1  
         poly_from_pixels,image,xrange,yrange,xpoly,ypoly 
         roi(nROI).xrange=xrange
         roi(nROI).yrange=yrange
         n=n_elements(xpoly)
         roi(nROI).npoly=n
         roi(nROI).xpoly(0:n-1)=xpoly
         roi(nROI).ypoly(0:n-1)=ypoly
         roi(nROI).nRZPhi=4
         roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
         roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
         nROI=nROI+1
      endif
    endfor
  endfor
;
;==================================================
; ROI - Instrumented Ramped Tiles
;==================================================;
;
; Load active pixel indices and determine enclosing polygon line segments from active pixels
;
goto,skip1
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
;
; Double Ramped Tiles 73, 74, 75, 76
;
  for jrow=1,4 do begin
    phi0=phi_tile(4) & phi1=phi_tile(6)
    z0=z_tile(jrow-1) & z1=z_tile(jrow)
    active=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
    indices=where(active,nindices)
    if nindices gt 0 then begin
       roi(nROI).nindices=nindices & roi(nROI).indices(0:nindices-1)=indices
       roi(nROI).name='Ramped Tile:'+sval(jrow+72)
;
       image=bytarr(xpix,ypix) & image(indices)=1  
       poly_from_pixels,image,xrange,yrange,xpoly,ypoly 
       roi(nROI).xrange=xrange
       roi(nROI).yrange=yrange
       n=n_elements(xpoly)
       roi(nROI).npoly=n
       roi(nROI).xpoly(0:n-1)=xpoly
       roi(nROI).ypoly(0:n-1)=ypoly
       roi(nROI).nRZPhi=4
       roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
       roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
       nROI=nROI+1
    endif
  endfor
;
; Ramped Tiles 77/78, 79, 80 columns A and B
;
  for icolumn=5,6 do begin
    AB='B' & if icolumn eq 6 then AB='A'
;
; B tiles are thicker than A tiles
;
; 77/78 - nose tile
;
    phi0=phi_tile(icolumn-1) & phi1=phi_tile(icolumn)
    z0=z_tile(4) & z1=z_tile(6)
    active=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
    indices=where(active,nindices)
    if nindices gt 0 then begin
       roi(nROI).nindices=nindices & roi(nROI).indices(0:nindices-1)=indices
       roi(nROI).name='Ramped Tile:77/78'+AB
;
       image=bytarr(xpix,ypix) & image(indices)=1  
       poly_from_pixels,image,xrange,yrange,xpoly,ypoly 
       roi(nROI).xrange=xrange
       roi(nROI).yrange=yrange
       n=n_elements(xpoly)
       roi(nROI).npoly=n
       roi(nROI).xpoly(0:n-1)=xpoly
       roi(nROI).ypoly(0:n-1)=ypoly
       roi(nROI).nRZPhi=4
       roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
       roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
       nROI=nROI+1
    endif
;
; Tiles 79 and 80
;
    for jrow=7,8 do begin
      phi0=phi_tile(icolumn-1) & phi1=phi_tile(icolumn)
      z0=z_tile(jrow-1) & z1=z_tile(jrow)
      active=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
      indices=where(active,nindices)
      if nindices gt 0 then begin
         roi(nROI).nindices=nindices & roi(nROI).indices(0:nindices-1)=indices
         roi(nROI).name='Ramped Tile:'+sval(jrow+72)+AB
  ;
         image=bytarr(xpix,ypix) & image(indices)=1  
         poly_from_pixels,image,xrange,yrange,xpoly,ypoly 
         roi(nROI).xrange=xrange
         roi(nROI).yrange=yrange
         n=n_elements(xpoly)
         roi(nROI).npoly=n
         roi(nROI).xpoly(0:n-1)=xpoly
         roi(nROI).ypoly(0:n-1)=ypoly
         roi(nROI).nRZPhi=4
         roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
         roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
         nROI=nROI+1
      endif
    endfor
  endfor
skip1:
;
;==================================================
; ROI - Lower Ramped Tiles (nose and below)
;==================================================;
;
; Cut off ramped tiles ROI to nose tiles and below
;
  phi0=phi_tile(4) & phi1=phi_tile(6)
  z0=z_tile(0) & z1=z_tile(6)
  set=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)

  roi(nROI).name='Lower Ramped Tiles'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Active Image - subregion of full IR image. Sets SUBFRAME region in Spectroscopy tree, displayed in IR_SUBVIDEO
;====================================================================================================
;
  phi0=phi_tile(4)-1.0 & phi1=phi_tile(6)
  z0=z_tile(0) & z1=z_tile(7)
  set=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)

  roi(nROI).name='Active Image'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Active Image 2 - subregion of full IR image (version 2). Sets SUBFRAME region in Spectroscopy tree, displayed in IR_SUBVIDEO
;       This is a larger region than Active Image, which accommodates background ROIs on either side of the ramped tiles.
;====================================================================================================
;
  phi0=phi_tile(3) & phi1=phi_tile(7)
  z0=z_tile(0) & z1=z_tile(7)
  set=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)

  roi(nROI).name='Active Image 2'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Active Image 3 - subregion of full IR image (version 3). Sets SUBFRAME region in Spectroscopy tree, displayed in IR_SUBVIDEO
;       This is a larger region than Active Image, which accommodates background ROIs on either side of the ramped tiles.
;====================================================================================================
;
  phi0=phi_tile(3) & phi1=phi_tile(7)
  z0=z_tile(0) & z1=0.5*(z_tile(7)+z_tile(8))+0.001
  set=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)
  roi(nROI).name='Active Image 3'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Heat Flux Region - region of IR image used for temperature and heat flux analysis
;====================================================================================================
;
  Phi0=Phi0_ramp & Phi1=Phi1_ramp & Z0=z_tile(0)+.01 & Z1=-0.402
  set=phipixel ge Phi0 and phipixel le Phi1 and zpixel ge Z0 and zpixel le Z1
  indices=where(set)

  roi(nROI).name='Heat Flux Region'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Heat Flux Region 2 - region of IR image used for temperature and heat flux analysis
;====================================================================================================
;
  Phi0=Phi0_ramp & Phi1=Phi1_ramp & Z0=z_tile(0)+.002 & Z1=-0.397
  set=phipixel ge Phi0+.3 and phipixel le Phi1-.3 and zpixel ge Z0 and zpixel le Z1
  crack=phipixel ge 270-.4 and phipixel le 270+.2 and zpixel gt z_tile(4)-.002
  set=set and not crack
  HFR2_set=set
  Z0_HFR2=Z0
  Z1_HFR2=Z1
  indices=where(set)

  roi(nROI).name='Heat Flux Region 2'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Heat Flux Region Slices - 6 regions of different toroidal anagles, used for temperature and heat flux analysis
;       These are based on Heat Flux Region 2 ROI
;====================================================================================================
;
  bins=indgen(7)*(272.6-267.4)/6.+267.4
  phi_bin=bins(0:5)+(272.6-267.4)/12.
  for jj=0,5 do begin
    Phi0=bins(jj) & Phi1=bins(jj+1)
    Z0=Z0_HFR2
    Z1=Z1_HFR2
    set=HFR2_set and phipixel ge Phi0 and phipixel le Phi1
    indices=where(set)
    roi(nROI).name='Heat Flux Region 2 Slice:'+sval(jj)
    n=n_elements(indices)
    roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
    xrange=[min(xbox),max(xbox)]
    yrange=[min(ybox),max(ybox)]
    poly_from_pixels,set,xrange,yrange,xpoly,ypoly
    roi(nROI).xrange=xrange
    roi(nROI).yrange=yrange
    n=n_elements(xpoly)
    roi(nROI).npoly=n
    roi(nROI).xpoly(0:n-1)=xpoly
    roi(nROI).ypoly(0:n-1)=ypoly
    roi(nROI).nRZPhi=4
    roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
    roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
    nROI=nROI+1
  endfor
;
;
;====================================================================================================
; ROI - Heat Flux Region 3 - region of IR image used for temperature and heat flux analysis
;====================================================================================================
;
  Phi0=Phi0_ramp & Phi1=Phi1_ramp & Z0=z_tile(0)+.002 & Z1=0.5*(z_tile(7)+z_tile(8))
  set=phipixel ge Phi0+.3 and phipixel le Phi1-.3 and zpixel ge Z0 and zpixel le Z1
  crack=phipixel ge 270-.4 and phipixel le 270+.2 and zpixel gt z_tile(4)-.002
  set=set and not crack
  HFR3_set=set
  Z0_HFR3=Z0
  Z1_HFR3=Z1
  indices=where(set)
  roi(nROI).name='Heat Flux Region 3'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Heat Flux Region Slices - 6 regions of different toroidal anagles, used for temperature and heat flux analysis
;       These are based on Heat Flux Region 3 ROI
;====================================================================================================
;
  bins=indgen(7)*(272.6-267.4)/6.+267.4
  phi_bin=bins(0:5)+(272.6-267.4)/12.
  for jj=0,5 do begin
    Phi0=bins(jj) & Phi1=bins(jj+1)
    Z0=Z0_HFR3
    Z1=Z1_HFR3
    set=HFR3_set and phipixel ge Phi0 and phipixel le Phi1
    indices=where(set)
    roi(nROI).name='Heat Flux Region 3 Slice:'+sval(jj)
    n=n_elements(indices)
    roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
    xrange=[min(xbox),max(xbox)]
    yrange=[min(ybox),max(ybox)]
    poly_from_pixels,set,xrange,yrange,xpoly,ypoly
    roi(nROI).xrange=xrange
    roi(nROI).yrange=yrange
    n=n_elements(xpoly)
    roi(nROI).npoly=n
    roi(nROI).xpoly(0:n-1)=xpoly
    roi(nROI).ypoly(0:n-1)=ypoly
    roi(nROI).nRZPhi=4
    roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
    roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
    nROI=nROI+1
  endfor
;
;====================================================================================================
; ROI - TSurface - used to define pixels that are included in the TSurface spline fit
;====================================================================================================
;
; Use slices 0 through 2 as defined above - the left half of the ramped tiles
;
  Phi0=bins(0) & Phi1=bins(2)
  Z0=Z0_HFR3
  Z1=Z1_HFR3

  set=HFR3_set and phipixel ge Phi0 and phipixel le Phi1
  indices=where(set)
  roi(nROI).name='TSurface'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;====================================================================================================
; ROI - Heat Flux Region Backgrounds, used for temperature and heat flux analysis
;====================================================================================================
;
;  Heat Flux Region 2 Background:0
;
  Phi0=Phi_Tile(3) & Phi1=Phi_Tile(4)-1.0
  Z0=z_tile(0) & Z1=z_tile(7)
  set=(phipixel ge Phi0) and (phipixel le Phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)
  roi(nROI).name='Heat Flux Region 2 Background:0'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;  Heat Flux Region 2 Background:1
;
  Phi0=Phi_Tile(6)+1.0 & Phi1=Phi_Tile(7)
  Z0=z_tile(0) & Z1=z_tile(7)
  set=(phipixel ge Phi0) and (phipixel le Phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)
  roi(nROI).name='Heat Flux Region 2 Background:1'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;  Heat Flux Region 3 Background:0
;
  Phi0=Phi_Tile(3) & Phi1=Phi_Tile(4)-1.0
  Z0=z_tile(0) & Z1=0.5*(z_tile(7)+z_tile(8))+.001
  set=(phipixel ge Phi0) and (phipixel le Phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)
  roi(nROI).name='Heat Flux Region 3 Background:0'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;  Heat Flux Region 3 Background:1
;
  Phi0=Phi_Tile(6)+1.0 & Phi1=Phi_Tile(7)
  Z0=z_tile(0) & Z1=0.5*(z_tile(7)+z_tile(8))+.001
  set=(phipixel ge Phi0) and (phipixel le Phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)
  roi(nROI).name='Heat Flux Region 3 Background:1'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;  Heat Flux Region 4 Background
;
  Phi0=Phi_Tile(6)+0.8 & Phi1=Phi_Tile(7)-1.2
  Z0=z_tile(0) & Z1=0.5*(z_tile(7)+z_tile(8))+.001
  set=(phipixel ge Phi0) and (phipixel le Phi1) and (zpixel ge z0) and (zpixel le z1)
  indices=where(set)
  roi(nROI).name='Heat Flux Region 4 Background'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  poly_from_pixels,set,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
; Plot Reference Tile image and optionally set Slot ROI
;
if set_roi_slot then begin
  if not ps then x,title='A-port IR Periscope View',window=2,xsize=xsize*2,ysize=ysize*2,colortable=3
  ref_frame=get_ref_frame(shot)
;  restore,file='/usr/local/cmod/codes/spectroscopy/ir/FLIR/ref_frame_1090730027_last.sav'
  tvscl,congrid(ref_frame(30:349,30:285),xsize*2,ysize*2,/interp,/center)
  p_pos=!p.position
  !p.position=[0,0,1,1]
  plot,xbox,ybox,xstyle=5,ystyle=5,/nodata,/noerase
  xyouts,/normal,.1,.9,'zrot = '+sval(zrot),charsize=1.2
  oplot_roi,roi,'Ramped Tiles',color=200
  oplot_roi,roi,'Lower Ramped Tiles',color=180
  for i=0,nROI-1 do begin
     name=roi(i).name
     if strpos(name,'TILE') gt -1 then oplot_roi,roi,name,color=160
  endfor
     print,'Click on divertor slot ROI (click to far left to exit)'
     yy=[0] & xx=[0]
again:
     cursor,x,y,/device,/down
     print,x,y
     if x lt 20 then goto,done
     xx=[xx,x] & yy=[yy,y]
     plots,xx(1:*),yy(1:*),/device,color=255
     goto,again
done:
     xx=xx(1:*)/2 & yy=yy(1:*)/2
     indices_roi_slot= POLYFILLV( xx, yy, xpix, ypix)
     save,file='/home/labombard/edge/modelling/geometry/indices_roi_slot.dat',indices_roi_slot
endif
;
;==================================================
; ROI - Divertor Slot
;==================================================;
;
  restore,'/home/labombard/edge/modelling/geometry/indices_roi_slot.dat'
  indices=indices_roi_slot
  roi(nROI).name='Divertor Slot'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
;==================================================
; ROI - Background (to scale image color table)
;==================================================;
;
;  Do this on the basis of the R,Z Phi coordinates
;
; Set1: tile columns 7 through 10, between phi=273 and phi=285 and tile types 73 through 79
;
  set1=phipixel gt phi_tile(6)+1.0 and phipixel lt phi_tile(10) and zpixel gt z_tile(0) and zpixel lt z_tile(7)
;
; Set2: tile columns 11 through 12 above tile 77
  set2=phipixel gt phi_tile(10) and phipixel lt phi_tile(12) and zpixel gt z_tile(5) and zpixel lt z_tile(9)

; Set3: tile columns 3 through 8 above tile 77
  set3=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(8) and zpixel gt z_tile(6) and zpixel lt z_tile(9)

; Set4: tile columns 3 through 4 above tile 73
  set4=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(4)-0.5 and zpixel gt z_tile(0) and zpixel lt z_tile(9)

  indices=where(set1 or set2 or set3 or set4)
  roi(nROI).name='Background'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
;====================================================================================================
; ROI - Camera Alignment Tiles - ramp_verion = 1 - for shots 1090701001 through 1100601000  
;       This ROI is used to stabilize the IR image for times during a plasma shot
;====================================================================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
; Set1: tile columns 2.5 through 12 above upper half of upper nose tiles
;
  set1=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(12) and zpixel gt 0.5*(z_tile(5)+z_tile(6)) and zpixel lt z_tile(9)
;
; Set2: tile columns 2 through 12 above tile 79
  set2=phipixel gt phi_tile(1)+2.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(7)-0.005 and zpixel lt z_tile(9)

; Set3: tile columns 1 through 12 above tile 80
  set3=phipixel gt phi_tile(0)+2.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(8)-0.005 and zpixel lt z_tile(9)

  indices=where(set1 or set2 or set3)
  roi(nROI).name='Camera Alignment Tiles'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
;
;====================================================================================================
; ROI - Camera Alignment Tiles2 - for all shots (ramp_verion=1 or ramp_verion=2)
;       This ROI is used to stabilize the IR image for times after a plasma shot
;====================================================================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
; Set1: all tile columns 2.5 and greater and tile rows 1 through 9
  set1=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(0) and zpixel lt z_tile(9)

; Set2: tile columns 2 through 12 above tile 79
  set2=phipixel gt phi_tile(1)+2.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(7)-0.005 and zpixel lt z_tile(9)

; Set3: tile columns 3 through 12 above nose tiles
  set3=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(6)-0.005 and zpixel lt z_tile(9)

; Set4: tile columns 1 through 12 above tile 80
  set4=phipixel gt phi_tile(0)+2.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(8)-0.005 and zpixel lt z_tile(9)

  indices=where(set1 or set2 or set3 or set4)
  roi(nROI).name='Camera Alignment Tiles2'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
;====================================================================================================
; ROI - Camera Alignment TilesB - for shots 1100601001 to present (ramp_verson=2)
;       This ROI is used to stabilize the IR image for times during a plasma shot
;====================================================================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
; Set1: all tile columns 2.5 through 4 and tile rows 2 through 9
  set1=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(4) and zpixel gt z_tile(1) and zpixel lt z_tile(9)

; Set2: tile columns 2 through 12 above tile 79
  set2=phipixel gt phi_tile(1)+2.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(7)-0.005 and zpixel lt z_tile(9)

; Set3: tile columns 3 through 12 above nose tiles
  set3=phipixel gt phi_tile(2)+1.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(6)-0.005 and zpixel lt z_tile(9)

; Set4: tile columns 1 through 12 above tile 80
  set4=phipixel gt phi_tile(0)+2.0 and phipixel lt phi_tile(12) and zpixel gt z_tile(8)-0.005 and zpixel lt z_tile(9)

  indices=where(set1 or set2 or set3 or set4)
  roi(nROI).name='Camera Alignment TilesB'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
;==================================================
; ROI - Divertor Tiles (all visible tiles on J outer divertor)
;==================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
;
; Set1: tile columns 4 through 11, all z
;
  set1=phipixel gt phi_tile(3) and phipixel lt phi_tile(12) and zpixel gt z_tile(0) and zpixel lt z_tile(9)

  indices=where(set1)
  roi(nROI).name='Divertor Tiles'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
;==================================================
; ROI - Divertor Tiles Above Nose (all visible tiles on J outer divertor above the nose)
;==================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
;
; Set1: tile columns 3 through 11
;
  phi0=phi_tile(2)+1.0 & phi1=phi_tile(12)
  z0=z_tile(6)-.002 & z1=z_tile(9)
  set1=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)

  indices=where(set1)
  roi(nROI).name='Divertor Tiles Above Nose'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;
;==================================================
; ROI - Divertor Tiles, Nose and Above 
;==================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
;
; Set1: tile columns 3 through 11
;
  phi0=phi_tile(2)+1.0 & phi1=phi_tile(12)
  z0=z_tile(4) & z1=z_tile(9)
  set1=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)

  indices=where(set1)
  roi(nROI).name='Divertor Tiles, Nose and Above'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;==================================================
; ROI - Divertor Tiles Vertical Face
;==================================================
;
;  Do this on the basis of the R,Z Phi coordinates
;
;
; Set1: tile columns 3 through 11
;
  phi0=phi_tile(2)+1.5 & phi1=phi_tile(12)+3.0
  z0=z_tile(0) & z1=z_tile(4)
  set1=(phipixel ge phi0) and (phipixel le phi1) and (zpixel ge z0) and (zpixel le z1)

  indices=where(set1)
  roi(nROI).name='Divertor Tiles Vertical Face'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  roi(nROI).nRZPhi=4
  roi(nROI).Phi(0:3)=[Phi0,Phi1,Phi1,Phi0]
  roi(nROI).z(0:3)=[Z0,Z0,Z1,Z1]
  nROI=nROI+1
;
;==================================================
; ROI - All but vertical face and slot
;==================================================
;
  indices=where(set1)
  data=bytarr(xpix,ypix)+1 & data(indices)=0
  data(indices_roi_slot)=0
  indices=where(data)
  roi(nROI).name='All but vertical face and slot'
  n=n_elements(indices)
  roi(nROI).nindices=n & roi(nROI).indices(0:n-1)=indices
  xrange=[min(xbox),max(xbox)]
  yrange=[min(ybox),max(ybox)]
  data=bytarr(xpix,ypix) & data(indices)=1
  poly_from_pixels,data,xrange,yrange,xpoly,ypoly
  roi(nROI).xrange=xrange
  roi(nROI).yrange=yrange
  n=n_elements(xpoly)
  roi(nROI).npoly=n
  roi(nROI).xpoly(0:n-1)=xpoly
  roi(nROI).ypoly(0:n-1)=ypoly
  nROI=nROI+1
;
skip_ROIs:
;
_xbox=xbox & _ybox=ybox & _nseg=nseg_sub & _lseg=lseg_sub & _Xseg=Xseg & _Yseg=Yseg & _ColorSeg=ColorSeg_sub & _LabelSeg=LabelSeg_sub & _IndexSeg=IndexSeg_sub
_Rpixel=Rpixel & _Zpixel=Zpixel & _Phipixel=Phipixel & _roi=roi & _cmod=cmod & _ramp_version=ramp_version

gar=' '
if strmid(camera_view,0,1) eq ' ' or camera_view eq 'X-pt' then goto,finish 
if not ps then x,title='Endoscope View - Stick Figure Alignment',window=2,xsize=xsize*2,ysize=ysize*2,colortable=3
if camera_view eq 'DIV' then begin
   restore,file='/usr/local/cmod/idl/div2_registration.sav'
   ref_frame=div2(0:479,0:639)
endif
if camera_view eq 'IR-A port' then begin
   ref_frame=get_ref_frame(shot)
   ref_frame=ref_frame(30:349,30:285)
endif
tvscl,congrid(ref_frame,xsize*2,ysize*2,/interp,/center)
p_pos=!p.position
!p.position=[0,0,1,1]
plot,xbox,ybox,xstyle=5,ystyle=5,/nodata,/noerase
;  xyouts,/normal,.1,.9,'zrot = '+sval(zrot),charsize=1.2,color=255
;  oplot_periscope_view,cmod=cmod,shot=shot,col=0
oplot_periscope_view,cmod=cmod,shot=shot,col=255
!p.position=p_pos
read,gar


;if show_roi then begin
;  if not ps then x,title='A-port IR Periscope View - Regions of Interest',window=3,xsize=xsize*4,ysize=ysize*4,colortable=43,xpos=700,ypos=100
;  ref_frame=get_ref_frame(shot)
;  restore,file='/usr/local/cmod/codes/spectroscopy/ir/FLIR/ref_frame_1090730027_last.sav'
;  f=bytscl(congrid(ref_frame(30:349,30:285),xsize*4,ysize*4,/interp,/center),top=255-16)+16
;  tv,f


;  oplot_roi,roi,'Lower Ramped Tiles',color=180,thick=2
;  for i=0,nROI-1 do begin
;     name=roi(i).name
;     if strpos(name,'TILE:') gt -1 then oplot_roi,roi,name,color=3,thick=2
;  endfor
;  oplot_roi,roi,'Divertor Slot',color=50,thick=2
;  oplot_roi,roi,'Ramped Tiles',color=0,thick=2
;  oplot_roi,roi,'Full Ramped Tiles',color=0,thick=2
;  oplot_roi,roi,'Ramped Tiles Test',color=14,thick=2
;  oplot_roi,roi,'Background',color=2,thick=2
;  oplot_roi,roi,'Camera Alignment Tiles',color=4,thick=2
;  oplot_roi,roi,'Camera Alignment Tiles2',color=5,thick=2
;  for ii=73,76 do begin
;    oplot_roi,roi,'Ramped Tile:'+sval(ii),color=6,thick=2
;  endfor
;  oplot_roi,roi,'Ramped Tile:77/78A',color=7,thick=2
;  oplot_roi,roi,'Ramped Tile:77/78B',color=8,thick=2
;  for ii=79,80 do begin
;    oplot_roi,roi,'Ramped Tile:'+sval(ii)+'B',color=12,thick=2
;    oplot_roi,roi,'Ramped Tile:'+sval(ii)+'A',color=15,thick=2
;  endfor
;  for ii=73,76 do begin
;    oplot_roi,roi,'Ramped Tile:'+sval(ii),color=6,thick=2
;  endfor
;  oplot_roi,roi,'Ramped Tile:77/78A',color=7,thick=2
;  oplot_roi,roi,'Ramped Tile:77/78B',color=8,thick=2

;  oplot_roi,roi,'Heat Flux Region 3',color=2,thick=2

; ilist=[0,3,4]
;  for ii=0,n_elements(ilist)-1 do begin
;    oplot_roi,roi,'Heat Flux Region 3 Slice:'+sval(ilist(ii)),color=9+ii,thick=2
;  endfor
;;  oplot_roi,roi,'Active Image 3',color=4,thick=2
;  if ramp_version eq 1 then begin
;    oplot_roi,roi,'Heat Flux Region 3 Background:0',color=1,thick=2
;    oplot_roi,roi,'Heat Flux Region 3 Background:1',color=1,thick=2
;  endif
;  if ramp_version eq 2 then begin
;    oplot_roi,roi,'Heat Flux Region 4 Background',color=1,thick=2
;  endif
;;  oplot_roi,roi,'Ramped Tiles',color=1,thick=2
;endif

finish:
if ps then begin
;   ps,/close
    psc
endif
stop
return
end

