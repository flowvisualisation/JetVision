;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function field_interp, U, V, nx, ny, x, y, dx, dy, xp, yp

;
;
;
;
; PURPOSE:  Linearly interpolate the vector
;           fields  U and V at point xp and yp
;
;
;
; ARGUMENTS:
;
;   U,V:    2-D vectors to be interpolated
;   nx,ny:  number of points in the two directions
;   x,y:    1-D abscissas and ordinates
;   dx,dy:  mesh sizes
;   xp,yp:  location where U and V should be interpolated
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; -----------------------------------------
; q is the two element vector with the
; interpolated values of U and V
; at xp and yp
; -----------------------------------------

  q = fltarr(2)

; ----------------------------------------
;  Locate i0 and j0, indexes of
;  cell to which xp,yp belong to
; ----------------------------------------

  i0 = (WHERE(ABS(xp - x) LT 0.51*dx))[0]
  j0 = (WHERE(ABS(yp - y) LT 0.51*dy))[0]

  IF i0 EQ -1 THEN BEGIN
    i0 = 0 + (xp GT x(nx-1))*(nx-1)
  ENDIF

  IF j0 EQ -1 THEN BEGIN
    j0 = 0 + (yp GT y(ny-1))*(ny-1)
  ENDIF

;
; interpolate U along x
;

  scrh1 = interpol(U(*,j0),x,xp)
  scrh2 = interpol(U(i0,*),y,yp)
  q(0)  = scrh1 + scrh2 - U(i0,j0)

;print,"q0 :",q(0)
;
; interpolate V along x
;

  scrh1 = interpol(V(*,j0),x,xp)
  scrh2 = interpol(V(i0,*),y,yp)
  q(1)  = scrh1 + scrh2 - V(i0,j0)
;print,"q1 :",q(1)

  return, q
end


;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; NAME:    field_line
;
; AUTHOR:  A. Mignone
;          mignone@to.astro.it
;
; PURPOSE: Given a 2-D vector field (U,V), it computes
;          the field line passing through x0, y0
;
;
; SYNTAX: field_line, U, V, x, y, x0, y0, xfield, yfield
;
; REQUIRES: interpol
;
; Arguments (INPUT):
;
;   U: 2-D vector giving the 1st component of the field
;   V: 2-D vector giving the 2nd component of the field
;
;   x: 1-D vector with abscissas
;   y: 1-D vector with ordinates
;
;   x0,y0: point coordinates through which the
;          field line goes.
;
; Arguments (OUTPUT):
;
;   xfield: 1-D vector containing the abscissa of the
;           field line
;
;   yfield: 1-D vector containing the ordinates of the
;           field line
;
;
; Last Modified: Nov 5, 2005
;
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pro field_line, U,V, x, y, x0, y0, xfield, yfield

  max_tol = 1.e-5

  sa = size(U)
  nx = sa[1]
  ny = sa[2]

; -------------------
;  Define grid sizes
; -------------------

  dx = fltarr(nx)
  dy = fltarr(ny)

;
; Notice that this definition for the mesh
; spacing does not necessarily coincide with
; the real one, when stretched grid are used,
; i.e.  x(i+1/2) != x(i) + 0.5*dx(i).
;

  dx = x - shift(x,1)
  dy = y - shift(y,1)
  ;print,'dx = ',dx
  ;print,'dy = ',dy

  dx(0) = x(1) - x(0)
  dy(0) = y(1) - y(0)


;  for i = 1, nx-1 do dx(i) = x(i) - x(i - 1)
;  for j = 1, ny-1 do dy(j) = y(j) - y(j - 1)

; grid ranges:

  xbeg = x(0)    - 0.5*dx(0)
  xend = x(nx-1) + 0.5*dx(nx-1)

  ybeg = y(0)    - 0.5*dy(0)
  yend = y(ny-1) + 0.5*dy(ny-1)

; -------------------------------
;  check whether the point x0,y0
;  falls  outside grid ranges
; -------------------------------

  inside_domain = (x0 LT xend) AND $
                  (x0 GT xbeg) AND $
                  (y0 LT yend) AND $
                  (y0 GT ybeg)

  IF (inside_domain EQ 0) THEN BEGIN
    print," Point ",x0,y0," outside grid range"
    print,xbeg, xend
    print,ybeg, yend
    print,min(x),max(x)
    print,min(y),max(y)
    return
  ENDIF

; ----------------------------------------------------------
;  xln_fwd, yln_fwd  = coordinates for forward  integration
;  xln_bck, yln_bck  = coordinates for backward integration
; ----------------------------------------------------------

  MAX_STEPS = 20000
  xln_fwd = fltarr(MAX_STEPS)
  yln_fwd = fltarr(MAX_STEPS)
  xln_bck = fltarr(MAX_STEPS)
  yln_bck = fltarr(MAX_STEPS)

; ----------------------------------
;       Initial conditions
; ----------------------------------

  xln_fwd(0) = x0
  yln_fwd(0) = y0
  xln_bck(0) = x0
  yln_bck(0) = y0

  rhs = fltarr(2)

; ---------------------------------------------
;            Integrate Forward
; ---------------------------------------------

  inside_domain = 1
  k             = 0
  WHILE inside_domain DO BEGIN

    ; -------------------------------
    ;   Take one step, using  RK II
    ; -------------------------------

    R1 = field_interp(U, V, nx, ny, x, y, dx, dy, $
                      xln_fwd(k), yln_fwd(k))

   ; --  Adjust step size so that dl \sim dx --

    dl = min([dx,dy])/sqrt(R1(0)^2 + R1(1)^2 + 1.d-14)

    xscrh = xln_fwd(k) + 0.5*dl*R1(0)
    yscrh = yln_fwd(k) + 0.5*dl*R1(1)

    R2 = field_interp(U, V, nx, ny, x, y, dx, dy, $
                      xscrh, yscrh)

    xln_one = xln_fwd(k) + dl*R2(0)
    yln_one = yln_fwd(k) + dl*R2(1)

    xln_fwd(k + 1) = xln_one
    yln_fwd(k + 1) = yln_one

   ; ----------------------------------------------------
   ;  Check to see whether we're still inside the domain
   ; ----------------------------------------------------

    inside_domain = (xln_fwd(k+1) GT xbeg) AND $
                    (xln_fwd(k+1) LT xend) AND $
                    (yln_fwd(k+1) GT ybeg) AND $
                    (yln_fwd(k+1) LT yend)

   ; ----------------------------------------------------
   ;     or exit anyway if k exceed maximum iteration #
   ; ----------------------------------------------------

    inside_domain = inside_domain AND (k LE MAX_STEPS-3)

    k = k+1
    IF (k GE MAX_STEPS - 1) THEN break
  ENDWHILE
  k_fwd = k

; ---------------------------------------------
;            Integrate Backward
; ---------------------------------------------

  inside_domain = 1
  k             = 0
  WHILE inside_domain DO BEGIN

    ; -------------------------------
    ;   Take one step, using  RK IV
    ; -------------------------------

    R1 = field_interp(U, V, nx, ny, x, y, dx, dy, $
                      xln_bck(k), yln_bck(k))

   ; --  Adjust step size so that dl \sim dx --

    dl = -min([dx,dy])/sqrt(R1(0)^2 + R1(1)^2 + 1.d-14)

    xscrh = xln_bck(k) + 0.5*dl*R1(0)
    yscrh = yln_bck(k) + 0.5*dl*R1(1)

    R2 = field_interp(U, V, nx, ny, x, y, dx, dy, $
                      xscrh, yscrh)

    xln_one = xln_bck(k) + dl*R2(0)
    yln_one = yln_bck(k) + dl*R2(1)

    xln_bck(k + 1) = xln_one
    yln_bck(k + 1) = yln_one

   ; ----------------------------------------------------
   ;  Check to see whether we're still inside the domain
   ; ----------------------------------------------------

    inside_domain = (xln_bck(k+1) GT xbeg) AND $
                    (xln_bck(k+1) LT xend) AND $
                    (yln_bck(k+1) GT ybeg) AND $
                    (yln_bck(k+1) LT yend)

   ; ----------------------------------------------------
   ;     or exit anyway if k exceed maximum iteration #
   ; ----------------------------------------------------

    inside_domain = inside_domain AND (k LE MAX_STEPS-3)

    k = k+1
    IF (k GE MAX_STEPS - 1) THEN break
  ENDWHILE

  k_bck = k
  print,"Fwd integration steps: ",strcompress(string(k_fwd,format='(i5)')),$
        "; bck integration steps: ",strcompress(string(k_bck,format='(i5)'))

  print, "Field line range  (", $
         string(xln_bck(k_bck),format=('(f6.2)')), ",",$
         string(yln_bck(k_bck),format=('(f6.2)')), ") - (",$
         string(xln_fwd(k_fwd),format=('(f6.2)')), ",",$
         string(yln_fwd(k_fwd),format=('(f6.2)')),")"

; --------------------------------------------
;         return arrays
; --------------------------------------------


  xfield = [reverse(xln_bck(0:k_bck)), xln_fwd(0:k_fwd)]
  yfield = [reverse(yln_bck(0:k_bck)), yln_fwd(0:k_fwd)]

;print,k_fwd, k_bck
end

