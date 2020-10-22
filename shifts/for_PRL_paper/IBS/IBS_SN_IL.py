import numpy as np
from numpy import sqrt, sin, cos, fabs, cbrt, pi,log

#Symmetric elliptic integral of the second kind:
def R_D(x, y, z) :
      
    C1 = 3./14.
    C2 = 1./6.
    C3 = 9./22.
    C4 = 3./26. 
    C5 = 9./88.
    C6 = 9./52.
    ERRTOL = 0.05
    
    xt = x
    yt = y
    zt = z
    delx = 1.
    dely = 1.
    delz = 1.
    asum = 0.0
    fac = 1.0
    
    while max(fabs(delx), fabs(dely), fabs(delz)) > ERRTOL :
        sqrtx = sqrt(xt)
        sqrty = sqrt(yt)
        sqrtz = sqrt(zt)
        alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        asum += fac/(sqrtz*(zt+alamb))
        fac = 0.25*fac
        xt = 0.25*(xt+alamb)
        yt = 0.25*(yt+alamb)
        zt = 0.25*(zt+alamb)
        ave = 0.2*(xt+yt+3.0*zt)
        delx = (ave-xt)/ave
        dely = (ave-yt)/ave
        delz = (ave-zt)/ave

    
    ea = delx*dely
    eb = delz*delz
    ec = ea - eb
    ed = ea - 6.0*eb
    ee = ed + ec + ec
    
    return 3.0*asum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))

    #calculation of local IBS rates, associated with particular ring azimuth point 
def IBS_lrates(gam,emx,sgmx,emy,sgmy,sgms,sgmp, betx, Dx, Phix, bety) :
    ax = betx/emx  
    ay = bety/emy  
    a_s = ax*(Dx*Dx/(betx*betx)+Phix*Phix)+ 1./(sgmp*sgmp) 
    a1 = (ax + gam*gam*a_s)/2.   
    a2 = (ax - gam*gam*a_s)/2.
    l1=ay 
    l2=(a1 + sqrt(a2*a2+gam*gam*ax*ax*Phix*Phix))  
    l3=(a1 - sqrt(a2*a2+gam*gam*ax*ax*Phix*Phix)) 
    x,y,z=1/l1,1/l2,1/l3
    T_ave = x+y+z
    h = 1./cbrt(x*y*z)
    X = h*x
    Y = h*y
    Z = h*z
    h12 = sqrt(h)
    r1 = R_D(Y, Z, X)
    r2 = R_D(Z, X, Y)
    r3 = 3.-r1-r2
    R1 = h12*X*r1      
    R2 = h12*Y*r2    
    R3 = h12*Z*r3      
    a2xp = sqrt(a2*a2+gam*gam*ax*ax*Phix*Phix)      
    a32 = 3.*a2/a2xp
    Sp = gam*gam*(2.*R1-R2*(1.-a32)-R3*(1.+a32))/2.    
    Sx = (2.*R1-R2*(1.+a32)-R3*(1.-a32))/2.           
    Sxp = (R3-R2)*3.*gam*gam*Phix*Phix*ax/a2xp        
    Psiy = -2.*R1 + R2 + R3                           
    rate_sgmp2 = Sp/sgmx/sgmy                       
    rate_emy = bety*Psiy/sgmx/sgmy 
    rate_emx = betx*(Sx+(Dx*Dx/(betx*betx)+Phix*Phix)*Sp+ Sxp)/sgmx/sgmy       
    return np.array([rate_emx, rate_emy, rate_sgmp2]),T_ave

def CalcRates(Aibs,re,gam,emx,emy,sgms,sgmp,dz_arr, bet_x, D_x, Phi_x, bet_y,showLcMinMaxMean = False):
    beta = sqrt(1-1/gam**2)
    coef = re/gam**2/beta**2
    allrates = np.array([0,0,0])
    Lcs = []
    for dz,betx,Dx,Phix,bety in zip(dz_arr,bet_x,D_x,Phi_x,bet_y):
        sgmx = sqrt(emx*betx+Dx**2*sgmp**2) 
        sgmy = sqrt(emy*bety)
        dallrates,T_ave = IBS_lrates(gam,emx,sgmx,emy,sgmy,sgms,sgmp, betx, Dx, Phix, bety)
        rmax = min(sgmx,sgmy,gam*sgms)
        rmin = coef/T_ave
        Lc = log(rmax/rmin)
        if showLcMinMaxMean:
            Lcs.append(Lc)
        allrates = allrates+dz*Lc*dallrates 
    if showLcMinMaxMean:
        print('Lcmin,Lcmax,Lcmean = {:.2e}, {:.2e}, {:.2e}'.format(min(Lcs),max(Lcs),np.mean(Lcs)))
    return Aibs/sgms*allrates

def CalcRelRates(Aibs,re,gam,emx,emy,sgms,sgmp,dz_arr, bet_x, D_x, Phi_x, bet_y,showLcMinMaxMean = False):
    rates = CalcRates(Aibs,re,gam,emx,emy,sgms,sgmp,dz_arr, bet_x, D_x, Phi_x, bet_y,showLcMinMaxMean)
    ems = [emx,emy,sgmp**2]
    return np.array([i/j for i,j in zip(rates,ems)])