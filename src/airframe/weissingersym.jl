
weissingerSymmetric(args...) = weissingerSymmetric(atleast1d(args[1:8])...,(args[9:end])...)

function weissingerSymmetric(
    xrle::SVector{Nelem, Float64},
    xrte::SVector{Nelem, Float64},
    xtle::SVector{Nelem, Float64},
    xtte::SVector{Nelem, Float64},
    yrle::SVector{Nelem, Float64},
    ytle::SVector{Nelem, Float64},
    rooti::SVector{Nelem, Float64},
    tipi::SVector{Nelem, Float64},
    alpha::Float64,
    xref::Float64,
    Sref::Float64,
    cref::Float64,
    bref::Float64) where Nelem
    """ Simple VLM for planar wings with single chordwise panels (Weissinger).
    The total lift, drag, and ptching moment coefficients are computed, 
    along with the local lift (Cl c / cref) and lift coefficient (Cl) between the root and tip.
    In this version of the code, the geometry is planar, composed of multiple coplanar trapezoidal 
    elements. Also, the geometry is assumed symmetric and only the right side (y>0) needs to be input.
        
    Inputs:
    xrle, xrte, xtle, xtte = x (~streamwise) position of leading and trailing edges at root and tip
    yrle, yrte = y location of element root and tip (we assume that yrte = yrle and ytte = ytle)
    rooti, tipi = element incidence at root and tip (in degrees)
    alpha = system angle of attack (angle of freestream from x axis in degrees)
    xref, Sref, bref, cref = configuration reference values for CL, CD, Cm, and e (moment is about xref) 

    Outputs:
    xctl = x coordinates of control points
    yctl = y coordinates of control points
    chord = local wing chord at panel lateral centers
    Clcbar = Cl*chord/cref  (proportional to local lift, l(y))
    Cl = section lift coefficient, l(y)/(q*chord)
    CL, Cm, CDi, e = lift, pitching moment, induced drag coefficients, and span efficiency

    Examples (all with wing root leading edge at x=0, but this is not necessary):
    Sample call: AR 10, rectangular wing (1 element)
    xctl,yctl,gamma,chord,Clcbar,Cl,CL,Cm,CDi,e = weissingerSymmetric(0.,1.,0.,1.,0.,5.,0.,0.,10.,0.,10.,1.,10.)

    Sample call: swept tapered wing (1 element)
    xctl,yctl,gamma,chord,Clcbar,Cl,CL,Cm,CDi,e = weissingerSymmetric(0.,1.,1.,1.5,0.,5.,0.,-1.,10.,0.,7.5,0.75,10.)

    Sample call: rectangular wing + tail
    xctl,yctl,gamma,chord,Clcbar,Cl,CL,Cm,CDi,e = weissingerSymmetric([0.,5.],[1.,5.5],[0.,5.],[1.,5.5],[0.,0],[5.,2.5],[0.,-8.],[0.,-8.],10.,0.25,10.,1.,10.)

    Sample call: swept wing (2 elements)
    xctl,yctl,gamma,chord,Clcbar,Cl,CL,Cm,CDi,e = weissingerSymmetric([0.,0.4],[1.,1.4],[0.4,1.],[1.4,2.],[0.,2.],[2.,5.],[0.,0.],[0.,0.],10.,0.,10.,1.0,10.)

    Sample call: canard (2 elements)
    xctl,yctl,gamma,chord,Clcbar,Cl,CL,Cm,CDi,e = weissingerSymmetric([-5.,0.],[-4.5,1.],[-5.,0.],[-4.5,1.],[0.,0.],[2.5,5.],[5.,0.],[5.,0.],10.,-1.2,10.,1.,10.)
    """

    # Element geometry:
    # Auto-paneling: Creates 100 panels over the semi-span and aligns all panels
    ymax = maximum(ytle) 
    Nelempan=round.(Int, 100*(ytle-yrle)./ymax)

    croot = (xrte-xrle)
    ctip  = (xtte-xtle) 
    sspan = ytle-yrle 
    area  = (croot+ctip).*sspan/2 
    dyelem = sspan ./ Nelempan
    xrqc = xrle+croot/4 
    xtqc = xtle+ctip/4 
    tans = (xtqc-xrqc)./sspan

    # Configuration totals (note that Sref, bref, cref are used for CL, Cm, CD, e, Clcbar)
    Stotal = 2*sum(area)
    btotal = 2*ymax
    xwake = 100*btotal/2
    Npantotal = sum(Nelempan)

    # Preallocate
    dy, yvc, ybvr, ybvt, xvc, xbvr, xbvt, xctl, yctl, eta, chord, paneli = (zeros(Npantotal) for _=1:12)
    AIC = zeros(Npantotal,Npantotal)
    rhs = zeros(Npantotal)
    gamma = zeros(Npantotal)
    Cl= zeros(Npantotal)
    Clcbar = zeros(Npantotal)

    # Panel geometry (all panels in a 1-D array)
    ipanel = 0
    for ielem = 1:Nelem
        for iy = 1:Nelempan[ielem]
            ipanel = ipanel+1
            dy[ipanel]    = dyelem[ielem]
            yvc[ipanel]   = yrle[ielem] + iy*dyelem[ielem] - dyelem[ielem]/2 
            ybvr[ipanel]  = yvc[ipanel] - dyelem[ielem]/2 
            ybvt[ipanel]  = yvc[ipanel] + dyelem[ielem]/2
            xvc[ipanel]   = xrqc[ielem] + (yvc[ipanel] - yrle[ielem])*tans[ielem] 
            xbvr[ipanel]  = xrqc[ielem] + (ybvr[ipanel] - yrle[ielem])*tans[ielem]
            xbvt[ipanel]  = xrqc[ielem] + (ybvt[ipanel] - yrle[ielem])*tans[ielem]
            eta[ipanel]   = (yvc[ipanel]- yrle[ielem])/sspan[ielem] 
            chord[ipanel] = croot[ielem]*(1-eta[ipanel]) + ctip[ielem]*eta[ipanel] 
            xctl[ipanel]  = xvc[ipanel] + 0.5*chord[ipanel] 
            yctl[ipanel]  = yvc[ipanel]
            paneli[ipanel] = rooti[ielem]*(1-eta[ipanel])+tipi[ielem]*eta[ipanel]
        end
    end

    # Compute AIC matrix: downwash at the control point of panel i
    # due to unit vortex strength of horseshoe vortex j and its image
    for i=1:Npantotal
        for j=1:Npantotal
        # Bound vortex    
        Gx = xbvt[j]-xbvr[j]
        Gy = ybvt[j]-ybvr[j]
        Rx = xbvr[j]-xctl[i]
        Ry = ybvr[j]-yctl[i]
        Ryi = ybvr[j]+yctl[i]
        AIC[i,j] = VzVortex(Gx,Gy,Rx,Ry) + VzVortex(Gx,Gy,Rx,Ryi) 
            
        # Trailing vortex root
        Gx = xbvr[j]-xwake 
        Gy = 0
        Rx = xwake-xctl[i]
        Ry = ybvr[j]-yctl[i]
        Ryi = ybvr[j] + yctl[i]
        AIC[i,j] = AIC[i,j] + VzVortex(Gx,Gy,Rx,Ry) + VzVortex(Gx,Gy,Rx,Ryi)
        
        # Trailing vortex tip
        Gx = xwake - xbvt[j]
        Gy = 0
        Rx = xbvt[j]-xctl[i]
        Ry = ybvt[j]-yctl[i]
        Ryi = ybvt[j]+yctl[i]
        AIC[i,j] = AIC[i,j] + VzVortex(Gx,Gy,Rx,Ry) + VzVortex(Gx,Gy,Rx,Ryi)        
        end
    end

    # Compute freestream flow through control points:  U_\infty sin alpha_i
    rhs .= -sind.(alpha .+ paneli)

    # Solve the system
    gamma .= (AIC \ rhs)

    # Lift and drag based on Trefftz plane downwash
    CL = 0.
    Cm = 0.
    CDi = 0.
    for i=1:Npantotal
        Cl[i] = 2 * gamma[i] / chord[i] # = rho U gamma / qc =   2 gamma / (U c) = 2 gamma / c when U=1
        Clcbar[i] = Cl[i] * chord[i] / cref # = 2 * gamma / cavg
        CL = CL + 4 * gamma[i] * dy[i] / Sref # = 4 * gamma * sspan / n / (span * cavg) = 2 * gamma / cavg / n 
        Cm = Cm + 4 * gamma[i] * dy[i] / Sref * (xref-xvc[i]) / cref
        wi = 0.
        for j=1:Npantotal
            if (abs(ybvt[j]-yctl[i]) > 1e-12) && (abs(ybvr[j]-yctl[i]) > 1e-12)
                wi = wi + gamma[j] * (1/(ybvt[j]-yctl[i]) - 1/(ybvr[j]-yctl[i]) )/(4*pi)
                wi = wi + gamma[j] * (1/(ybvt[j]+yctl[i]) - 1/(ybvr[j]+yctl[i]) )/(4*pi)
            end
        end
        CDi = CDi + 4 * gamma[i] * dy[i] / Sref * wi
    end
    e = CL^2/(pi*bref^2/Sref) / CDi

    return xctl,yctl,gamma,chord,Clcbar,Cl,CL,Cm,CDi,e
end