##### Velocity due to horseshoe vortices ####
function VzVortex(Gx, Gy, Rx, Ry)
    # x is positive downstream, y is positive from centerline to right tip
    # Gx, Gy are the lengths of a vortex filament in the x and y directions
    # Rx, Ry is a vector from a vortex root to a control point. Rx = xroot - xctl
    R   = hypot(Rx,Ry)
    G2   = Gx*Gx + Gy*Gy
    GR   = Gx*Rx + Gy*Ry;
    GXRz = Gx*Ry - Gy*Rx
    E1   = GXRz*GXRz
    E    = hypot(Gx + Rx, Gy + Ry);
    V1   = (GR*E-(G2 + GR)*R)
    V2   = 4.0*pi*E1*E * R
    if (abs(GXRz) < 1e-12) && (R^2 < G2) && (GR<0) # Avoid infinite if query point is on vortex
        Vz = 0.
    else
        Vz = V1*GXRz / (V2+eps())
    end
end
VzVortexSemiInf(Rx, Ry) = (abs(Ry) < 1e-12) ? 0. : 1/(4π*Ry) * (1 - Rx / hypot(Rx,Ry))
symmetric_index(i, Npanhalfwing) = (Npanhalfwing-i+1>0) ? (Npanhalfwing-i+1) : (i-Npanhalfwing) # full span to semi-span index

function influence!(target::Union{Matrix{T},Tuple{Matrix{T},Matrix{T}}} where T, xtarget, ytarget, xr, xt, yr, yt)
    # If target is a Matrix, it's the influence coefficient matrix AIC
    #   AIC is a 4-blocs matrix if separated between left and right side
    #   AIC = [left->left (=right->right), right->left; left->right, right->right].
    #   left->left  = right->right
    #   left->right = (right->left)'
    # If target is a tuple, it is the tuple (Γ,w), where Γ is the vortex strength and w the induced velocity
    
    Npanhalfwing = length(xtarget)
    
    if target isa Tuple 
        @. target[2] = 0. # set w to 0
    end
    
    for i=1:Npanhalfwing
        iright = Npanhalfwing + i
        ileft  = Npanhalfwing - i + 1
        for j=1:Npanhalfwing
            jright = Npanhalfwing + j
            jleft  = Npanhalfwing - j + 1

            # Bound vortex size G and distance R to control point
            Gx = xt[j]-xr[j]
            Gy = yt[j]-yr[j]
            Rx = xr[j]-xtarget[i]
            Ry = yr[j]-ytarget[i]
            Rysym = yr[j]+ytarget[i]

            # Induced normal velocity
            Vzrr = VzVortex(Gx,Gy,Rx,Ry) + VzVortexSemiInf(Rx, Ry) - VzVortexSemiInf(Rx+Gx, Ry+Gy) # right -> right
            Vzlr = VzVortex(Gx,Gy,Rx,Rysym) + VzVortexSemiInf(Rx, Rysym) - VzVortexSemiInf(Rx+Gx, Rysym+Gy) # left -> right
            
            if target isa Matrix
                target[iright, jright] = -Vzrr # right -> right
                target[ileft , jleft ] = -Vzrr # left -> left
                target[iright, jleft ] = -Vzlr # left -> right
                target[ileft , jright] = -Vzlr # right -> left
            elseif target isa Tuple
                for i=1:size(target[2],2)
                    target[2][iright, i] += -Vzrr*target[1][jright,i]-Vzlr*target[1][jleft, i] # (right -> right) + (left -> right)
                    target[2][ileft,  i] += -Vzrr*target[1][jleft, i]-Vzlr*target[1][jright,i] # (left -> left) + (right -> left)
                end
            end
        end
    end
end



#### Velocity due to environment ####
# ωb = [ca*p-sa*r,q,sa*p+ca*r] (omega in body frame)
# V = V0 + [Lx,Ly,0.] x ωb
Vbody_x(ca,sa,Lx,Ly,p,q,r) = SVector(
     ca + Ly * ( ca*r + p*sa), # _
    -sa + Ly * (-sa*r + p*ca), # α
     Ly * sa, # p
     0., # q
     Ly*ca # r
)
Vbody_y(ca,sa,Lx,Ly,p,q,r) = SVector(
    - Lx*( ca*r + p*sa), # _
    - Lx*(-sa*r + p*ca), # α
    - Lx*sa, # p
    0., # q
    - Lx*ca, # r
)
Vbody_z(ca,sa,Lx,Ly,p,q,r) = SVector(
     sa + Lx*q - Ly*(ca*p - r*sa), # _
     ca - Ly*(-sa*p - r*ca), # α
    -Ly * ca, # p
     Lx, # q
     Ly * sa # r
)

#### Main function ####
const tpl{Nelem,T} = Union{NTuple{Nelem,T},StaticVector{Nelem,T}}
function weissinger!(
    results::AbstractWeissingerResults,
    xrle ::tpl{Nelem,<:Real},
    xrte ::tpl{Nelem,<:Real},
    xtle ::tpl{Nelem,<:Real},
    xtte ::tpl{Nelem,<:Real},
    yrle ::tpl{Nelem,<:Real},
    ytle ::tpl{Nelem,<:Real},
    rooti::tpl{Nelem,<:Real},
    tipi ::tpl{Nelem,<:Real},
    controls::tpl{Nelem,<:Real},
    alpha::Real,
    omega::tpl{3,<:Real},
    airfoil::Tuple, # Tuple should contain subelements of Airfoil class
    signsymcontrol::tpl{Nelem,Int64},
    ref  ::tpl{5, <:Real};
    symmetric::Union{Val{true},Val{false}}=Val{false}(),
    Δy::Union{NTuple{Nelem},Nothing}=nothing,
    VX::Union{Function,Nothing}=nothing,VY::Union{Function,Nothing}=nothing,VZ::Union{Function,Nothing}=nothing,
    reuseGeom::Union{Val{true},Val{false}}=Val{false}()) where {Nelem}
    """ Simple VLM for planar wings with single chordwise panels (Weissinger).
        The total lift, drag, and ptching moment coefficients are computed, 
        along with the local lift (Cl c / cref) and lift coefficient (Cl) between the root and tip.
        In this version of the code, the geometry is planar, composed of multiple coplanar trapezoidal 
        elements. Also, the geometry is assumed symmetric and only the right side (y>0) needs to be input.
            
        Inputs:
        xrle, xrte, xtle, xtte = x (~streamwise) position of leading and trailing edges at root and tip
        yrle , yrte            = y location of element root and tip (we assume that yrte = yrle and ytte = ytle)
        rooti, tipi            = element incidence at root and tip (in degrees)
        controls               = control surface deflection (in degrees). Array of length Nelem (on control per wing element).
        alpha                  = system angle of attack (angle of freestream from x axis in degrees)
        omega                  = [p,q,r] *non-dimensionalized* rotation rates in *stability* frame.
        ref                    = (xref , Sref, bref, cref, Vref) configuration reference values for CL, CD, Cm, and e (moment is about xref)
        airfoil                = Airfoil instance for each wing element (array of size nelem).
        signsymcontrol         = array of size Nelen, 0 if control is symmetric, 1 if a sign flip in delta is necessary
        symmetric              = Boolean type allowing only solve for positive y. Note that assymetric coefficients and derivatives (q,r,...) will all be 0.
                                 Currently this relies on views unfortunately does not provide performance gains.
        Δy                     = Optional vector of length Nelem containing vectors of panel size
        VX,VY,VZ               = Optional function to add flow velociy at a given location. VX/Y/Z(x,y,ca,sa) returns an SVector similar to body_x/y/z
        reuseGeom              = Boolean type indicating that the geometry has not changed: the mesh and AIC factorization can be reused
        
        Outputs:
        xctl = x coordinates of control points
        yctl = y coordinates of control points
        chord = local wing chord at panel lateral centers
        Clcbar = Cl*chord/cref  (proportional to local lift, l(y))
        Cl = section lift coefficient, l(y)/(q*chord)
        CL, Cm, CDff, e = lift, pitching moment, induced drag coefficients, and span efficiency

        Examples (all with wing root leading edge at x=0, but this is not necessary):
        Sample call: AR 10, rectangular wing (1 element)
        results = weissinger(0.,1.,0.,1.,0.,5.,0.,0.,10.,0.,0.,0.,0.,10.,1.,10.)

        Sample call: swept tapered wing (1 element)
        results = weissinger(0.,1.,1.,1.5,0.,5.,0.,-1.,10.,0.,0.,0.,0.,7.5,0.75,
        10.)

        Sample call: rectangular wing + tail
        results = weissinger([0.,5.],[1.,5.5],[0.,5.],[1.,5.5],[0.,0],[5.,2.5],[0.,-8.],[0.,-8.],10.,0.,0.,0.,0.25,10.,1.,10.)

        Sample call: swept wing (2 elements)
        results = weissinger([0.,0.4],[1.,1.4],[0.4,1.],[1.4,2.],[0.,2.],[2.,5.],[0.,0.],[0.,0.],10.,0.,0.,0.,0.,10.,1.0,10.)

        Sample call: canard (2 elements)
        results = weissinger([-5.,0.],[-4.5,1.],[-5.,0.],[-4.5,1.],[0.,0.],[2.5,5.],[5.,0.],[5.,0.],10.,0.,0.,0.,-1.2,10.,1.,10.)
    """
    # Unpack input
    xref,Sref,cref,bref,Vref = ref
    p,q,r = omega[1]*2/bref, omega[2]*2/cref, omega[3]*2/bref # make dimensional
    T = promote_type(
        eltype(xrle), eltype(xrte), eltype(xtle), eltype(xtte), eltype(yrle),
        eltype(ytle), eltype(rooti), eltype(tipi), typeof(alpha), eltype(omega), eltype(ref)
        )

    # Assert basics
    Nelemcheck = reduce(&, (Nelem == length(input)) for input=(xrle, xtle, xrte, xtte, yrle, ytle, rooti, tipi, airfoil, controls,))
    @assert Nelemcheck """Inconsistent number of wing elements throughout input:
    (xrle:$(length(xrle)), xtle:$(length(xtle)), xrte:$(length(xrte)), xtte:$(length(xtte)), yrle:$(length(yrle)), ytle:$(length(ytle)), rooti:$(length(rooti)), tipi:$(length(tipi)), airfoil:$(length(airfoil)), controls:$(length(controls)),)"""

    # Extract arrays from preallocated structure
    results_arrays = get_results(results,T)
    (; CL, Cm, Cl, CD, CDp, CDi, CDff, gamma, w, cl, clcbar, AIC) = results_arrays
    @assert length(CL) == (length(deriv) + Nelem) "WeissingerResults structure is not initialized with enough control entries."

    # Extract grid arrays 
    (; xc, yc, xr, yr, xt, yt, xctl, xrc, xtc, chord, paneli) = results_arrays
    # Zero out results
    for A=(CL,Cm,Cl,CD,CDff, gamma,w,cl,clcbar)
        A .= 0.
    end
    for A = CDp; A .= 0.; end
    for A = CDi; A .= 0.; end
    
    # Element geometry:
    Nelempan = results.Nelempan
    Npanhalfwing = sum(Nelempan) 
    Npantotal = 2 * Npanhalfwing # include mirrored side

    # Panel geometry (all panels in a 1-D array)
    # for y:
    #   r: root
    #   t: tip
    #   c: center
    # for x:
    #   c: bound vortex center
    #   r: bound vortex root
    #   t: bound vortex tip
    #   rc: root trailing vortex center
    #   tc: tip trailing vortex center
    #   ctl: control point
    if reuseGeom isa Val{false}
        for ielem = 1:Nelem
            croot = (xrte[ielem]-xrle[ielem])
            ctip  = (xtte[ielem]-xtle[ielem])
            xrqc  = xrle[ielem] + croot/4 
            xtqc  = xtle[ielem] + ctip/4 
            sspan = ytle[ielem]-yrle[ielem] 
            tans  = (xtqc-xrqc)./sspan
            
            y = yrle[ielem]
            for iy = 1:Nelempan[ielem]
                if Δy isa Nothing
                    dyelem = sspan / Nelempan[ielem]
                else
                    dyelem = Δy[ielem][iy]
                end
                ipanel = results_arrays.ielem[ielem]+iy
                
                # Root of bound vortex
                yr[ipanel] = y
                xr[ipanel] = xrqc + (yr[ipanel] - yrle[ielem]) * tans

                # Center of bound vortex
                yc[ipanel] = yr[ipanel] + dyelem/2 
                xc[ipanel] = xrqc + (yc[ipanel] - yrle[ielem]) * tans
                
                # Tip of bound vortex
                yt[ipanel] = yr[ipanel] + dyelem
                xt[ipanel] = xrqc + (yt[ipanel] - yrle[ielem]) * tans
                
                # Center of trailing root leg
                eta_r   = (yr[ipanel]- yrle[ielem])/sspan
                chord_r = croot * (1-eta_r) + ctip * eta_r 
                xrc[ipanel] = xr[ipanel] + chord_r * 3/8
                
                # Center of trailing tip leg
                eta_t   = (yt[ipanel]- yrle[ielem])/sspan
                chord_t = croot * (1-eta_t) + ctip * eta_t 
                xtc[ipanel] = xt[ipanel] + chord_t * 3/8
                
                # Control point on panel
                eta = (yc[ipanel]- yrle[ielem])/sspan
                chord[ipanel] = croot * (1-eta) + ctip * eta
                Re = Vref * chord[ipanel] / envt.nu
                cla = dcl_α(airfoil[ielem], 0., 0., Re) # 0-lift lift curve slope
                ctl_chord = cla / 4π # location of control point to obtain provided lift curve slope
                xctl[ipanel] = xc[ipanel] + ctl_chord*chord[ipanel]

                # Local effective twist: geometric + camber + control
                tw_geom = (rooti[ielem]*(1 - eta) + tipi[ielem] * eta)
                α0 = liftcoefficient(airfoil[ielem], 0., 0., Re)/ cla
                tw_camber = α0 * 180/pi
                paneli[ipanel] = tw_geom + tw_camber

                y+= dyelem
            end
        end
        
        # Compute AIC matrix: downwash at the control point of panel i
        # due to unit vortex strength of horseshoe vortex j
        influence!(AIC, xctl, yc, xr, xt, yr, yt)
    end

    # Compute flow through control points
    isym(i) = symmetric_index(i, Npanhalfwing)
    i = 0
    sa, ca = sincosd(alpha)
    for e=-Nelem:Nelem
        (e==0) && continue
        ielem = abs(e)
        snc = (e<0) ? signsymcontrol[ielem] : 1
        δ = controls[ielem]

        for _=1:Nelempan[ielem]
            i+= 1
            i_ = isym(i)
            Lx = xctl[i_]-xref
            Ly = yc[i_] * ((i<Npanhalfwing+1) ? -1 : 1)
            Vz = Vbody_z(ca,sa,Lx,Ly,p,q,r)
            Vx = Vbody_x(ca,sa,Lx,Ly,p,q,r)
            
            if !isa(VZ,Nothing)
                Vz += VZ(xc[i_],yc[i_],sa,ca)
            end
            if !isa(VX,Nothing)
                Vx += VX(xc[i_],yc[i_],sa,ca)
            end

            Re =  Vref * chord[i_] / envt.nu
            clα = dcl_α(airfoil[ielem], 0., 0., Re) # 0-lift lift curve slope
            clδ = dcl_δ(airfoil[ielem], 0., 0., Re) * snc # derivative at 0 lift
            tw_controls = clδ * δ / clα
            stw, ctw = sincosd(paneli[i_] + tw_controls)

            # Normal velocity at control point & derivatives wrt/ (α,p,q,r)
            for j = eachindex(Vz) 
                w[i,j] = (stw * Vx[j] + ctw * Vz[j])
            end

            # Derivative with respect to controls through "effective twist"
            w[i,ielem+5] = clδ / clα * (ctw * Vx[1] - stw * Vz[1])
        end
    end

    # Solve the system
    if symmetric isa Val{false}
        linear_solve!(AIC, gamma, w, results, reuseGeom)
    else
        for i=1:Npanhalfwing
            @. AIC[1:Npanhalfwing,i] += AIC[1:Npanhalfwing,2*Npanhalfwing-i+1]
        end
        A = view(AIC, 1:Npanhalfwing, 1:Npanhalfwing) # views will create a small number of allocations but it won't be as bad as LU of a full system
        g = view(gamma, 1:Npanhalfwing, :)
        w_= view(w, 1:Npanhalfwing, :)
        linear_solve!(A, g, w_, results)
        for i=1:Npanhalfwing
            @. gamma[2*Npanhalfwing-i+1, :] = g[i,:]
        end
    end

    ## Aerodynamic forces with Kutta-Joukovski theorem
    # Forces on central, left and right vortex with Kutta Joukovski theorem
    vtx_center_root = (
        (xc , yc, xr, yr,  1), # central
        (xrc, yr, xr, yr, -1), # root
        (xtc, yt, xt, yt,  1), # tip
        )
        
        for kvtx=eachindex(vtx_center_root)
            xvc, yvc, xvr, yvr, flipsign = vtx_center_root[kvtx]
            
            # Downwash at central vortex
            if kvtx == 1
                influence!((gamma,w),xvc,yvc, xr, xt, yr, yt)
                
                # influence!(AIC,xvc,yvc, xr, xt, yr, yt)
                # dualmul!(w, AIC, gamma, results)
                
                # mul!(w, AIC, gamma)
            else
                w .= 0.
            end
            
        # Forces and moments
        i = 0
        sa, ca = sincosd(alpha)
        for e=-Nelem:Nelem
            (e==0) && continue
            ielem = abs(e)
            δ = controls[ielem]

            for iy=1:Nelempan[ielem]

                i += 1
                leftwing = ((i<Npanhalfwing+1) ? -1 : 1)
                i_ = isym(i)
                dyelem = (yt[i_] - yr[i_])
                Lx = (xvc[i_]- xref)
                Ly = yvc[i_] * leftwing
                Gx = 2 * (xvc[i_] - xvr[i_]) * leftwing * flipsign
                Gy = 2 * (yvc[i_] - yvr[i_]) * flipsign
                Rei = Vref * chord[i_] / envt.nu

                # Local velocity at central vortex (in body frame)
                # ωb = [ca*p-sa*r,q,sa*p+ca*r]
                # V = V0 + [Lx,Ly,0.] x ω + [0,0,w]
                Vx = Vbody_x(ca,sa,Lx,Ly,p,q,r)
                Vy = Vbody_y(ca,sa,Lx,Ly,p,q,r)
                Vz = Vbody_z(ca,sa,Lx,Ly,p,q,r)

                if !isa(VZ,Nothing)
                    Vz += VZ(xvc[i_],yvc[i_],sa,ca)
                end
                if !isa(VX,Nothing)
                    Vx += VX(xvc[i_],yvc[i_],sa,ca)
                end
                
                # Local force on vortex (body frame)
                fx =-Gy * (Vz[1] - w[i,1]) * gamma[i,1]
                fy = Gx * (Vz[1] - w[i,1]) * gamma[i,1] 
                fz = (Gy * Vx[1] - Gx * Vy[1]) * gamma[i,1]
                
                # Local lift
                L = (fz * ca - fx * sa)
                CL[1] += L
                clcbar[i] += L / cref / dyelem
                cl[i,1] += L / (1/2 * chord[i_] * dyelem)

                # Induced drag
                CDi[ielem][1] += (fz * sa + fx * ca)

                # Local parasitic drag
                if kvtx == 1 # only do once
                    cd = 1/2 * dragcoefficient(airfoil[ielem], cl[i,1], δ, Rei) * chord[i_] * dyelem
                    Vn = hypot(Vx[1], Vy[1], Vz[1]-w[i,1])
                    dfx = cd * Vx[1] / Vn
                    dfy = cd * Vy[1] / Vn
                    dfz = cd * (Vz[1]-w[i,1]) / Vn
                    fx += dfx; fy += dfy; fz += dfz
                    CDp[ielem][1] += (dfz * sa + dfx * ca)
                end 

                # Drag
                CD[1] += (fz * sa + fx * ca)

                # Pitching moment due to aerodynamic force (z component in body frame)
                Cm[1] -= fz * Lx
                
                # Pitching moment due to local camber and control surface deflection
                # Here we assume that 2d sections only create a pitching moment, 
                # eventhough they may locally be in a sideslip, which would rotate the axis of rotation.
                # Neglected compared to actual roll due to control surfaces.
                if kvtx == 1
                    Cm[1] += 1/2 * momentcoefficient(airfoil[ielem], cl[i,1], δ, Rei) * chord[i_]^2 * dyelem
                end

                # Rolling moment due to lift force (z component in wind frame)
                Cl[1] += fz * Ly * ca - (fx * Ly - fy * Lx) * sa

                # Compute derivatives with respect to state variables
                for k=2:5
                    # Change in local force on vortex (body frame)
                    fx_= -Gy * ((Vz[1] - w[i,1]) * gamma[i,k] + (Vz[k] - w[i,k]) * gamma[i,1]) 
                    fy_=  Gx * ((Vz[1] - w[i,1]) * gamma[i,k] + (Vz[k] - w[i,k]) * gamma[i,1])
                    fz_=  ((Gy * Vx[k] - Gx * Vy[k]) * gamma[i,1] +
                            (Gy * Vx[1]  - Gx * Vy[1]) * gamma[i,k])
                    
                    # Change in induced drag
                    CDi[ielem][k] += (fz_ * sa + fx_ * ca)
                    
                    # Add change in local parasitic drag direction
                    if kvtx == 1
                        d1_Vn = - (Vx[k]*Vx[1] + Vy[k]*Vy[1] + Vz[k]*Vz[1]) / Vn^3
                        dfx_ = cd * (Vx[1] * d1_Vn + Vx[k] / Vn)
                        dfy_ = cd * (Vy[1] * d1_Vn + Vy[k] / Vn)
                        dfz_ = cd * (Vz[1] * d1_Vn + Vz[k] / Vn)
                        fx_ += dfx_; fy_ += dfy_; fz_ += dfz_
                        CDp[ielem][k] += (dfz_ * sa + dfx_ * ca)
                    end

                    # Change in local lift
                    L = (fz_* ca - fx_* sa) - ((k==2) ? (fz * sa + fx * ca) : 0.)
                    CL[k]   += L
                    cl[i,k] += L / (1/2 * chord[i_] * hypot(Gx,Gy))
                    
                    # Change in local parasitic drag magnitude
                    if kvtx == 1
                        dcd = 1/2 * cl[i,k] * dcd_cl(airfoil[ielem], cl[i,1], δ, Rei) * chord[i_] * dyelem
                        dfx_= dcd * Vx[1] / Vn
                        dfy_= dcd * Vy[1] / Vn
                        dfz_= dcd * Vz[1] / Vn
                        fx_ += dfx_; fy_ += dfy_; fz_ += dfz_
                        CDp[ielem][k] += (dfz_ * sa + dfx_ * ca) + ((k==2) ? (fz * ca - fx * sa) : 0.)
                    end
                    
                    # Change in drag
                    CD[k] += (fz_ * sa + fx_ * ca) + ((k==2) ? (fz * ca - fx * sa) : 0.)
                    
                    # Pitching moment due to aerodynamic force (z component in body frame)
                    Cm[k] -= fz_* Lx

                    # Change in Pitching moment due to local camber and control surface deflection
                    if kvtx == 1
                        Cm[k] += 1/2 * cl[i,k] * dcm_cl(airfoil[ielem], cl[i,1], δ, Rei) * chord[i_]^2 * dyelem
                    end
                    # Change in pitching moment direction
                    # Cm[1] += (dot(controls, clδ)*180/pi) * chord[i_] * dyelem
                    # Neglected
                    
                    # Rolling moment due to lift force (z component in wind frame)
                    Cl[k] += fz_* Ly * ca - (fx_* Ly - fy_* Lx) * sa + 
                            ((k==2) ? - fz * Ly * sa - (fx * Ly - fy * Lx) * ca : 0.)
                end
        

                # Compute derivatives with respect to controls
                for kc=1:Nelem
                    k = 5 + kc
                    # Change in local force on vortex (body frame), driven by changes in gamma and w
                    fx_= -Gy * ((Vz[1] - w[i,1]) * gamma[i,k] - w[i,k] * gamma[i,1])
                    fy_=  Gx * ((Vz[1] - w[i,1]) * gamma[i,k] - w[i,k] * gamma[i,1])
                    fz_=  (Gy * Vx[1]  - Gx * Vy[1]) * gamma[i,k]
                
                    # Local lift
                    L = (fz_* ca - fx_* sa)
                    CL[k] += L
                    cl[i,k] += L / (1/2 * chord[i_] * hypot(Gx,Gy))
                
                    # Induced drag
                    CDi[ielem][k] += (fz_ * sa + fx_ * ca)

                    # Local parasitic drag change with controls
                    if kvtx == 1
                        dcd  = 1/2 * (dcd_cl(airfoil[ielem], cl[i,1], δ, Rei) * cl[i,k]
                                + dcd_δ(airfoil[ielem], cl[i,1], δ, Rei)) * chord[i_] * dyelem
                        dfx_ = dcd * Vx[1] / Vn
                        dfy_ = dcd * Vy[1] / Vn
                        dfz_ = dcd * (Vz[1] - w[i,1]) / Vn
                        fx_ += dfx_; fy_ += dfy_; fz_ += dfz_
                        CDp[ielem][k] += (dfz_ * sa + dfx_ * ca)
                    end

                    # Change in drag
                    CD[k] += (fz_ * sa + fx_ * ca)
                
                    # Pitching moment due to aerodynamic force (z component in body frame)
                    Cm[k] -= fz_* Lx

                    # Pitching moment due to added camber
                    if kvtx == 1
                        Cm[k] += 1/2 * (dcm_cl(airfoil[ielem], cl[i,1], δ, Rei) * cl[i,k]
                                      + dcm_δ(airfoil[ielem], cl[i,1], δ, Rei)) * chord[i_]^2 * dyelem
                    end
                    
                    # Rolling moment due to lift force (z component in wind frame)
                    Cl[k] += fz_* Ly * ca - (fx_* Ly - fy_* Lx) * sa
                end
            end
        end
    end

    CL ./= Sref / 2
    CD ./= Sref / 2
    for (cdp,cdi)=zip(CDp,CDi)
        cdp ./= Sref / 2
        cdi ./= Sref / 2
    end
    Cm ./= Sref * cref / 2
    Cl ./= Sref * bref / 2
    
    for val=(CL,Cm,Cl,CD)
        val[3] *= 2 / bref # p
        val[4] *= 2 / cref # q
        val[5] *= 2 / bref # r
    end
    for val=CDp
        val[3] *= 2 / bref # p
        val[4] *= 2 / cref # q
        val[5] *= 2 / bref # r
    end
    for val=CDi
        val[3] *= 2 / bref # p
        val[4] *= 2 / cref # q
        val[5] *= 2 / bref # r
    end

    # Section drag in Trefftz plane
    CDff[1] = 0.
    for i=1:Npanhalfwing 
            iright = Npanhalfwing + i
            ileft  = Npanhalfwing - i + 1
            wright = 0.
            wleft  = 0.
        for j=1:Npanhalfwing
            jright = Npanhalfwing + j
            jleft  = Npanhalfwing - j + 1

            di    = (1/(yt[j]-yc[i]) - 1/(yr[j]-yc[i]) )/(4*pi)
            dsym  = (1/(yt[j]+yc[i]) - 1/(yr[j]+yc[i]) )/(4*pi)
            
            wright += gamma[jright, 1] * di + gamma[jleft,  1] * dsym
            wleft  += gamma[jleft,  1] * di + gamma[jright, 1] * dsym
        end
        CDff[1] += 2 * (gamma[iright, 1] * wright + gamma[ileft, 1] * wleft) * (yt[i]-yr[i])
    end
    CDff[1] = CDff[1] / Sref

    return (CL, Cm, Cl, CD, CDff, cl)

end


## Custom version of mul! for dual numbers (mul! is not pure Julia, it wrapps BLAS)
dualmul!(w::AbstractArray{T},AIC::AbstractArray{T},
gamma::AbstractArray{T}, wres::AbstractWeissingerResults) where T<:Number = mul!(w,AIC,gamma)

function dualmul!(w::AbstractArray{ForwardDiff.Dual{T,V,N}},
            AIC::AbstractArray{ForwardDiff.Dual{T,V,N}},
            gamma::AbstractArray{ForwardDiff.Dual{T,V,N}},
            wres::AbstractWeissingerResults) where {T,V,N} 
    # w = AIC * gamma
    # dw/dx = AIC * dgamma/dx + dAIC/dx * gamma

    # Use primal matrices as buffers
    gammap = wres.primal.gamma
    wp = wres.primal.w
    AICp = wres.primal.AIC

    # 1/ Value
    map!(a->a.value, AICp, AIC) # copy dual AIC into primal matrix
    map!(a->a.value, gammap, gamma) # copy dual gamma into primal vector
    map!(a->a.value, wp, w) # copy dual w into primal w
    mul!(wp, AICp, gammap) # perform multiplication with BLAS
    map!((w, wdual)->FD.Dual{T}(w,wdual.partials), w, wp, w) # set w dual value

    # 2/ Set dw/dx to AIC * dgamma/dx for partials
    for ipl=1:N
        map!(a->a.partials[ipl], gammap, gamma) # copy dual dgamma/dx into primal buffer
        mul!(wp, AICp, gammap) # perform multiplication with BLAS, wp contains dw/dx
        map!((wdual,dwdx)-> FD.Dual{T}(wdual.value, 
                FD.Partials{N,V}(setindex(wdual.partials.values,dwdx,ipl))), 
                # FD.Partials{N,V}(ntuple(i -> (i==ipl) ? w : wdual.partials[i], N))), 
                w, w, wp) # set dw/dx
    end

    # 3/ Add dAIC/dx * gamma
    map!(a->a.value, gammap, gamma) # copy dual gamma into primal vector
    for ipl=1:N
        map!(a->a.partials[ipl], AICp, AIC) # copy dual dAIC/dx into primal buffer
        mul!(wp, AICp, gammap) # perform multiplication with BLAS, wp contains dAIC/dx * gamma
        map!((wdual,dwdx)-> FD.Dual{T}(wdual.value, 
                # FD.Partials{N,V}(ntuple(i -> wdual.partials[i] + ((i==ipl) ? w : 0.), N))), 
                FD.Partials{N,V}(setindex(wdual.partials.values,wdual.partials.values[ipl]+dwdx,ipl))), 
                w, w, wp) # add to dw/dx
    end
end

## Custom version of linear solver for dual numbers (getrs! & getrf! are from LAPACK)
function linear_solve!(AIC::AbstractArray{T}, gamma::AbstractArray{T}, 
                        w::AbstractArray{T}, wres::AbstractWeissingerResults,
                        reuseGeom::Union{Val{false},Val{true}}) where T<:Number
    if reuseGeom isa Val{false}
        (_,ipiv,_) = LAPACK.getrf!(AIC);
        @. wres.ipiv = ipiv
    end
    
    gamma .= w
    LAPACK.getrs!('N', AIC, wres.ipiv, gamma)
end

function linear_solve!(AIC::AbstractArray{ForwardDiff.Dual{T,V,N}}, 
                        gamma::AbstractArray{ForwardDiff.Dual{T,V,N}},
                        w::AbstractArray{ForwardDiff.Dual{T,V,N}}, 
                        wres::AbstractWeissingerResults,
                        reuseGeom::Union{Val{false},Val{true}}) where {T,V,N}
    # (RHS means right-hand-side.)
    
    # Use primal matrices as buffers
    gammap = wres.primal.gamma
    wp     = wres.primal.w
    AICp   = wres.primal.AIC
    
    # If symmetric cases, create views
    if AIC isa SubArray
        gammap = view(gammap, gamma.indices...)
        wp     = view(wp,     w.indices)
        AICp   = view(AICp,   AIC.indices)
    end
    
    # 1/ LU decomposition of AIC
    if reuseGeom isa Val{false}
        map!(a->a.value, AICp, AIC); # copy dual value into primal matrix
        (_,ipiv,_) = LAPACK.getrf!(AICp); # lu decomposition, main source of runtime & allocations
        @. wres.ipiv = ipiv
    end

    # 2/ Solve for value of gamma
    map!(a->a.value, gammap, w); # copy RHS into primal gamma    
    LAPACK.getrs!('N', AICp, wres.ipiv, gammap) # solve w/ LU
    
    # 3/ Store LU(AIC) and gamma in dual.value
    map!((g, gdual)->FD.Dual{T}(g,gdual.partials), gamma, gammap, gamma) # set gamma dual value
    map!((g, gdual)->FD.Dual{T}(g,gdual.partials), AIC, AICp, AIC) # set AIC dual value to LU

    # 4/ Solve linear system for each partial
    for ipl=1:N
        # 4.1/ Create right hand side, store in wp : RHS = - dAIC/dx * Gamma + dw/dx
        map!(a->a.partials[ipl], AICp, AIC) # dAICdx stored in AICp
        map!(a->a.partials[ipl], wp, w) # dwdx stored in wp
        mul!(wp, AICp, gammap, -1., 1.) # wp = - AICp * gammap + wp

        # 4.2/ Solve linear system
        map!(a->a.value, AICp, AIC); # bring LU(AIC) back to AICp
        LAPACK.getrs!('N', AICp, wres.ipiv, wp) # solve w/ LU (dgamma/dx stored in wp)
        
        # 4.3/ Update dGamma/dx in partial
        map!((gdual,g)-> FD.Dual{T}(gdual.value, 
                FD.Partials{N,V}(setindex(gdual.partials.values,g,ipl))), 
                gamma, gamma, wp) # set partial value
    end
end

function weissinger(
    xrle ::SVector{Nelem, T},
    xrte ::SVector{Nelem, T},
    xtle ::SVector{Nelem, T},
    xtte ::SVector{Nelem, T},
    yrle ::SVector{Nelem, T},
    ytle ::SVector{Nelem, T},
    rooti::SVector{Nelem, T},
    tipi ::SVector{Nelem, T},
    alpha::T,
    omega::SVector{3,T},
    ref  ::SVector{5,T},) where {Nelem,T}

    # Element geometry: Creates 100 panels over the semi-span and aligns all panels
    ymax = maximum(ytle) 
    Nelempan = round.(Int, 100*(ytle-yrle)./ymax)
    
    wres=WeissingerResults(Float64, Nelempan);
    airfoil = fill(ThinAirfoil(), Nelem)
    signsymcontrol = fill(0, Nelem)
    controls = fill(0., Nelem)
    (CL, Cm, Cl, CD, CDff, cl) = weissinger!(wres, xrle, xrte, xtle, xtte, yrle, ytle, rooti, tipi, controls, alpha, omega, airfoil, signsymcontrol, ref);

    # controls::AbstractVector{<:Real},
    # alpha::Real,
    # omega::Union{AbstractVector{<:Real},NTuple{3, <:Real},StaticVector{3,<:Real}},
    # airfoil::AbstractVector,
    # signsymcontrol::AbstractVector{<:Real},
    # ref  ::Union{AbstractVector{<:Real},NTuple{4, <:Real},StaticVector{4,<:Real}}, 
    return  (wres, CL, Cm, Cl, CD, CDff, cl)
end

function quad_wing_tail_yspacing!(y, n1, n2, span_1, span_2)
    N = n1+n2
    k = ((span_2+span_1)-span_1*N/n1)*(N)/n2
    map!(i-> i/N*((span_1+span_2)+k*(i/N-1)), y,0:n1+n2)
end

function quad_wing_tail_yspacing!(Δy, span_1, span_2)
    n1 = length(Δy[1]) 
    n2 = length(Δy[2])
    y_ = 0.
    N = n1 + n2
    k = ((span_2+span_1)-span_1*N/n1)*(N)/n2
    for i=1:n1
        y = i/N*((span_1+span_2)+k*(i/N-1))
        Δy[1][i] = y-y_
        Δy[3][i] = Δy[1][i]
        y_ = y
    end
    for i=1:n2
        y = (n1+i)/N*((span_1+span_2)+k*((n1+i)/N-1))
        Δy[2][i] = y-y_
        y_ = y
    end
end

function cos_wing_tail_yspacing!(y, n1, n2, n3, span_1, span_2, span_3)
    c1=c2=c3=0.
    for i=0:n1-1
        y[i+1] = span_1*(sin(π/2*i/n1).+c1)/(1+c1)
    end
    for i=0:n2-1
        y[n1+i+1] = (span_1+span_2/2) + span_2*(sin.(π/2*(-1+2*i/n2)).+c2)/(2+c2)
    end
    for i=0:n3
        y[n1+n2+i+1] = (span_1+span_2+span_3/2) + span_3*(sin.(π/2*(-1+2*i/n3)).+c3)/(2+c3)
    end
end

function wing_tail_yspacing!(Δy, span_1, span_2)
    """
    Computes Δy spacing for 2 lifting surfaces behind each other, where the vortices must be aligned.
    The surfaces are in the following order [inboard wing section, outboard wing section]. 
    The tail is behind the inboard wingsection.
    Δy is the input to weissinger! and has therefore 3 components: [inboard wing, outboardwing, tail]
    """
    n1 = length(Δy[1])
    n2 = length(Δy[2])

    slope_1 = (span_1/n1-span_2/n2) * 6 / (-(1-2*n2)*(n2-1)*n1/n2 + (n1+1)*(1+2*n1))
    slope_2 = -n1*slope_1/n2
    ff(n,b,l,N) = b/N + l/6*(N-1)*(N+1) - l * n * (n-1)/2

    map!(n->ff(n, span_1, slope_1, n1), Δy[1], 1:n1)
    @. Δy[3] = Δy[1]
    map!(n->ff(n2+n1-n+1,span_2, slope_2, n2), Δy[2], n1+1:n2+n1)
end

# function cos_wing_tail_yspacing!(y, n1, n2, span_1, span_2)
#     y[1:n1+1] = LinRange(0., span_1, n1+1) # equal panel spacing on inboard section
#     θ₁ = acos(1/(n1+1))
#     # θ₁ = acos(span_1 / (span_1+span_2))
#     θ  = LinRange(θ₁, 0.0, n2+1)
#     @. y[n1+2:n1+n2+1] = y[n1+1] + (cos(θ[2:end])-cos(θ[1])) * (span_1+span_2)
# end

function wing_tail_yspacing!(y, n1, n2, span_1, span_2)
    """
    Computes Δy spacing for 2 lifting surfaces behind each other, where the vortices must be aligned.
    The surfaces are in the following order [inboard wing section, outboard wing section]. 
    The tail is behind the inboard wingsection.
    Δy is the input to weissinger! and has therefore 3 components: [inboard wing, outboardwing, tail]
    """
    slope_1 = (span_1/n1-span_2/n2) * 6 / (-(1-2*n2)*(n2-1)*n1/n2 + (n1+1)*(1+2*n1))
    slope_2 = -n1*slope_1/n2
    ff(n,b,l,N) = b/N + l/6*(N-1)*(N+1) - l * n * (n-1)/2

    y[1] = 0.
    for n=1:n1
        y[n+1] = y[n] + ff(n, span_1, slope_1, n1)
    end
    for n=n1+1:n1+n2
        y[n+1] = y[n] + ff(n2+n1-n+1,span_2, slope_2, n2)
    end
    return y
end
wing_tail_yspacing(n1, n2, span_1, span_2) = wing_tail_yspacing!(zeros(1+n1+n2), n1, n2, span_1, span_2)

function striptheory(f2d, coef::Vector{T}, V::Real, airfoil, controls, wres; 
                        elems=eachindex(wres.Nelempan),
                        rightwing::Union{Val{true},Val{false}}=Val{true}()) where T
    Npanhalfwing = sum(wres.Nelempan)
    cl = get_array(wres, :cl, T)
    chord = get_array(wres, :chord, T)

    # Right wing only
    if rightwing isa Val{true}
        for elt=elems
            afl = airfoil[elt]
            ctl = controls[elt]
            ini = (elt==1) ? 0 : sum(wres.Nelempan[i] for i=1:elt-1)
            for i=ini+1:ini+wres.Nelempan[elt]
                Rei = V*chord[i]/envt.nu
                coef[i] = f2d(afl,cl[Npanhalfwing+i],ctl,Rei)
            end
        end
    
    # Both wings
    else
        Ne = length(wres.Nelempan)
        (Ne != length(elems)) && raise("Both wings only implemented for full vortex system")
        ini = 0
        for e=-Ne:Ne
            if e==0
                continue
            end
            elt = abs(e)
            afl = airfoil[elt]
            ctl = controls[elt]
            for i=ini+1:ini+wres.Nelempan[elt]
                Rei = V*chord[symmetric_index(i, Npanhalfwing)]/envt.nu
                coef[i] = f2d(afl,cl[i],ctl,Rei)
            end
            ini += wres.Nelempan[elt]
        end
    end
end

function elementwise_cl(wres::AbstractWeissingerResults, elt::Int64, Δα::Real=0.)
    """
    Computes lift coefficient of a lifting element from the lifting distribution.
    It is possible to get the lift coefficient at a different angle of attack.
    Δα in radians!
    """
    (; Nelempan, yt, yr, chord, ielem, cl) = wres
    lift(j) = (cl[j] + Δα*cl[j,2])
    Npanhalfwing = sum(Nelempan)
    S = 0
    CLright = 0
    CLleft = 0
    for i=ielem[elt]+1:ielem[elt]+Nelempan[elt]
        S  += 2*(yt[i]-yr[i])*chord[i] # assuming symmetric lifting surface
        CLright += lift(Npanhalfwing+i)   * (yt[i]-yr[i]) * chord[i]
        CLleft  += lift(Npanhalfwing-i+1) * (yt[i]-yr[i]) * chord[i]
    end
    return SA[(CLright+CLleft)/S, S]
end

function elementwise_cl(wres::AbstractWeissingerResults, elt::Tuple, Δα::Real=0.)
    CL = 0.
    S  = 0.
    for e=elt
        CL_, S_ = elementwise_cl(wres, e, Δα)
        CL += CL_*S_
        S  += S_
    end
    return CL/S, S
end