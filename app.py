## WEBPAGE
from flask import Flask, render_template, request,jsonify
## MATH LIBRARIES
from scipy.special import erfc
from numpy import exp,sqrt,linspace,abs,minimum,expm1

app=Flask(__name__)

## constants
kappa=10**(-7)      ## thermal diffusivity                                  [s/m^2]
rho_MD=1000         ## maximum density of water                             [kg/m^3]
c_P=4210            ## heat capacitance of water                            [J-*C/kg]
T_MD=4              ## temperature of max density                           [*C]
I_S=500             ## midday spring solar irradiance at high lattitude     [W-m^-2] -- MAKE PARAMETER

'''def graph(z0,H,T_H,t,kw,kb,p,npts=400):
    ## making sure values are treated right
    kw=float(kw)    ## k_w:     attenuation coefficient for white ice           slider to change based on range (in paper) 
    kb=float(kb)    ## k_b:     attenuation coefficient for black ice ~[0.1]
    p=float(p)      ## p:       proportionality                                 on [0,1] w/dP~0.05
    H=float(H)
    T_H=float(T_H)
    z0=float(z0)
    npts=int(npts)
    
    clarity=0.015
    L=0.5  ## L:       ice thickness
    alpha=0.3 if p>0 else 0.1   ## alpha:   surface albedo
    
    I_0=(1-alpha)*I_S*exp(-1*L*(kw*p+kb*(1-p)))     ## light intensity
    
    Z = abs(z)     # fine since z is negative
    H = -z0        # depth to bottom (your choice)
    B = T_H / H    # slope so T0(H)=TH
    T_0 = B * Z    # linear baseline (A=0)

    
    z=linspace(z0,0.0,npts)
    
    H=-z0 #fixes the T-> 300 issue by making H the depth to the bottom of the lake
    
    #T_0=-T_H*z/H
    
    t=t*3600
    t = max(t, 1e-9)
    
    ##lambda=e2,A=0,B=-T_H/H
    
    #coeff=I_0/(clarity*rho_MD*c_P)
    #T=T_0+coeff*np.exp(z/clarity)*t ##nondiffusive solution
    #clarity=0.015
    Z=abs(z)
    
    B = T_H / H        # slope
    T_0 = B * Z        # linear baseline (A = 0)

    
    L_diff=sqrt(kappa*t)
    e2=clarity
    C=I_0*e2/(rho_MD*c_P)
    term1=(C/(kappa*e2**2))*erfc(Z/(2*L_diff))
    term2a=-1*(C/(2*kappa*e2**2))*exp(kappa*t*e2**2)
    term2bi=exp(-1*e2*Z)*erfc(Z/(2*L_diff)-e2*L_diff)
    term2bii=exp(e2*Z)*erfc(Z/(2*L_diff)+e2*L_diff)
    term2=term2a*(term2bi+term2bii)
    term3 = (C/(kappa*e2**2))*exp(-e2*Z)*expm1((L_diff*e2)**2)
    T=term1+term2+term3+T_0 if t !=0 else T_0
    
    dT = T - T_0
    print("dT at surface:", dT[-1])      # if z array ends at 0
    print("dT at bottom:", dT[0])
    print("max|dT|:", float(abs(dT).max()))
    
    return z.tolist(), T.tolist()'''
    
def graph(z0,H,T_H,t,kw,kb,p,npts=400):
    kw=float(kw); kb=float(kb); p=float(p)
    H=float(H); T_H=float(T_H); z0=float(z0)
    npts=int(npts)

    # parameters
    clarity = 0.5
    L = 0.5
    alpha = 0.3 if p>0 else 0.1

    I_0 = (1-alpha)*I_S*exp(-L*(kw*p + kb*(1-p)))

    # enforce sign conventions
    T_H = abs(T_H)
    z0  = -abs(z0)

    # depth grid (z <= 0)
    z = linspace(z0, 0.0, npts)
    Z = abs(z)              # depth below ice (>=0)

    # baseline: A=0, linear BZ (this is the one that "rejoins")
    Hdepth = -z0            # lake depth
    B = T_H / Hdepth
    T_0 = B * Z

    # time
    t = max(float(t)*3600.0, 1e-9)

    # heating solution (Eq 12a with A=0)
    L_diff = sqrt(kappa*t)
    e2 = clarity
    C  = I_0*e2/(rho_MD*c_P)

    term1 = (C/(kappa*e2**2))*erfc(Z/(2*L_diff))

    term2a  = -(C/(2*kappa*e2**2))*exp(kappa*t*e2**2)
    term2bi = exp(-e2*Z)*erfc(Z/(2*L_diff) - e2*L_diff)
    term2bii= exp(+e2*Z)*erfc(Z/(2*L_diff) + e2*L_diff)
    term2   = term2a*(term2bi + term2bii)

    term3 = (C/(kappa*e2**2))*exp(-e2*Z)*expm1((L_diff*e2)**2)

    T = term1 + term2 + term3 + T_0

    # debug (safe)
    dT = T - T_0
    print("dT at surface:", float(dT[-1]))
    print("dT at bottom:", float(dT[0]))
    print("max|dT|:", float(abs(dT).max()))

    return z.tolist(), T.tolist()

    

@app.route("/", methods=["GET"])
def index():
    ##  inputs variables and gets ready to assign default values
    defaults={
        "z0":10.0,
        "H":1.0,
        "T_H":10.0,
        "t":2.0,
        "kw":6.9,
        "kb":2.2,
        "p":0.0
    }
    return render_template("index.html",**defaults)

def safe_float(val, default):
    try:
        return float(val)
    except:
        return default


@app.route("/api/profile",methods=["GET"])
def api_profile():
    try:
        z0 = safe_float(request.args.get("z0", ""), 10)
        H  = safe_float(request.args.get("H", ""), 1)
        T_H= safe_float(request.args.get("T_H",""), 10)
        t  = safe_float(request.args.get("t", ""), 2)
        kw = safe_float (request.args.get("kw", ""), 6.9)
        kb = safe_float(request.args.get("kb", ""), 2.2)
        p  = safe_float(request.args.get("p", ""), 0)
    except ValueError:
        z0,H,T_H,t,kw,kb,p=10.0,1.0,10.0,2.0,6.9,2.2,0.0
    
    z,T=graph(z0,H,T_H,t,kw,kb,p,npts=400)
    return jsonify({"z":z,"T":T})

if __name__=="__main__":
    app.run(debug=True)
