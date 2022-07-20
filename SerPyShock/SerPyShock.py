import sys, os
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

# To make it innto pip etc
# https://gist.github.com/jgieseler/1840eab2a3345041aebb08d5da0c959b

# For the loaders is going to be somethign liike
# pip install git+https://github.com/jgieseler/solo-mag-loader


# Licenses
# https://choosealicense.com/

#realpython.com


class ShockParameters:
    pass


def get_time_indices(tstart, tend, tm):
    """
    Function that computes indices corresponding to start and end time in a given stream of datetimes

    Parameters
    ----------
    tstart : 'datetime'
        Starting time.
    tend : 'datetime'
        End time.
    tm : 'DatetimeIndex'
        Stream of time in which the search of indices is desired

    Returns
    -------
    its : 'int'
        index corresponding to tstart.
    ite : 'int'
        Index corresponding to tend.
    """
    its = min(min(np.where(tm > tstart)))
    ite = max(max(np.where(tm < tend)))
    return its, ite
    

def select_subS(F, tF, its, ite, dims):
    """
    Function that selects a time interval for a measured quantity

    Parameters
    ----------
    F : 'numpy array'
        numpy nd array with n_timestamps x n_components dimensions.
    tF : 'DatetimeIndex'
        Stream of times over which the measurment is obtained.
    its : 'int'
        Index corresponding to starting time
    ite : 'int'
        Index corresponding to end time
    dims : 'int'
        Number of components of the measured quantity.

    Returns
    -------
    t_sS : 'DatetimeIndex'
        Interval time stream
    ssF : 'numpy array'
        Array containing measurement for give interval

    """
    t_sS = tF[its:ite]
    ssF = np.squeeze(np.zeros([t_sS.shape[0],dims]))
    if dims > 1:
        for id in range(dims):
            ssF[:,id] = F[its:ite,id]
    else:
        ssF[:] = np.squeeze(F[its:ite])
    return t_sS, ssF


def calc_MC(Bu,Bd,frame):
    """
    Routine that computes shock normal vector and theta_Bn angle in a given frame
    using the Magnetic Coplanarity method

    Parameters
    ----------
    Bu : 'numpy array'
        Upstream mean magnetic field (3 components).
    Bd : 'numpy array'
        Downstream mean magnetic field (3 components).
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'

    Returns
    -------
    n : 'numpy array'
        Shock normal vector.
    tbn : 'float64'
        Shock normal angle.
    """
    DB = Bd - Bu
    BdXBu = np.cross(Bd,Bu)
    BdXBuXDB = np.cross(BdXBu,DB)
    n = BdXBuXDB/np.linalg.norm(BdXBuXDB)
    if (frame == 'RTN'):
        if(n[0] < 0):
            n = -n
    elif(frame == 'GSE'):
        if(n[0] > 0):
            n = -n
    else:
        print('frame not valid. Implemented: RTN,GSE')
    
    tbn = np.arccos( np.dot(Bu,n)/np.linalg.norm(Bu))*180/np.pi
    if tbn > 90:
        tbn = 180 -tbn
        
    return n,tbn


def calc_MX1(Bu,Bd,Vu, Vd,frame):
    """
    Routine that computes shock normal vector and theta_Bn angle in a given frame
    using the Mixed Mode 1 method

    Parameters
    ----------
    Bu : 'numpy array'
        Upstream mean magnetic field (3 components).
    Bd : 'numpy array'
        Downstream mean magnetic field (3 components).
    Vu : 'numpy array'
        Upstream mean bulk flow speed (3 components).
    Vd : 'numpy array'
        Downstream mean bulk flow speed (3 components).
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'

    Returns
    -------
    n : 'numpy array'
        Shock normal vector.
    tbn : 'float64'
        Shock normal angle.
    """
    
    DB = Bd - Bu 
    DV = Vd - Vu 
    BuXDV    = np.cross(Bu, DV) 
    BuXDVXDB = np.cross(BuXDV,DB) 
    n = BuXDVXDB/(np.linalg.norm(BuXDVXDB)) 
    if (frame == 'RTN'):
        if(n[0] < 0):
            n = -n
    elif(frame == 'GSE'):
        if(n[0] > 0):
            n = -n
    else:
        print('frame not valid. Implemented: RTN,GSE')
    
    tbn = np.arccos( np.dot(Bu,n)/np.linalg.norm(Bu))*180/np.pi
    if tbn > 90:
        tbn = 180 -tbn
        
    return n,tbn


def calc_MX2(Bu,Bd,Vu, Vd,frame):
    """
    Routine that computes shock normal vector and theta_Bn angle in a given frame
    using the Mixed Mode 2 method

    Parameters
    ----------
    Bu : 'numpy array'
        Upstream mean magnetic field (3 components).
    Bd : 'numpy array'
        Downstream mean magnetic field (3 components).
    Vu : 'numpy array'
        Upstream mean bulk flow speed (3 components).
    Vd : 'numpy array'
        Downstream mean bulk flow speed (3 components).
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'

    Returns
    -------
    n : 'numpy array'
        Shock normal vector.
    tbn : 'float64'
        Shock normal angle.
    """
    
    DB = Bd - Bu 
    DV = Vd - Vu 
    BdXDV    = np.cross(Bd, DV) 
    BdXDVXDB = np.cross(BdXDV,DB) 
    n = BdXDVXDB/(np.linalg.norm(BdXDVXDB)) 
    if (frame == 'RTN'):
        if(n[0] < 0):
            n = -n
    elif(frame == 'GSE'):
        if(n[0] > 0):
            n = -n
    else:
        print('frame not valid. Implemented: RTN,GSE')
    
    tbn = np.arccos( np.dot(Bu,n)/np.linalg.norm(Bu))*180/np.pi
    if tbn > 90:
        tbn = 180 -tbn
        
    return n,tbn




def calc_MX3(Bu,Bd,Vu, Vd,frame):
    
    """
    Routine that computes shock normal vector and theta_Bn angle in a given frame
    using the Mixed Mode 3 method

    Parameters
    ----------
    Bu : 'numpy array'
        Upstream mean magnetic field (3 components).
    Bd : 'numpy array'
        Downstream mean magnetic field (3 components).
    Vu : 'numpy array'
        Upstream mean bulk flow speed (3 components).
    Vd : 'numpy array'
        Downstream mean bulk flow speed (3 components).
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'

    Returns
    -------
    n : 'numpy array'
        Shock normal vector.
    tbn : 'float64'
        Shock normal angle.
    """
        
    DB = Bd - Bu 
    DV = Vd - Vu 
    DBXDV    = np.cross(DB, DV) 
    DBXDVXDB = np.cross(DBXDV,DB) 
    n = DBXDVXDB/(np.linalg.norm(DBXDVXDB)) 
    if (frame == 'RTN'):
        if(n[0] < 0):
            n = -n
    elif(frame == 'GSE'):
        if(n[0] > 0):
            n = -n
    else:
        print('frame not valid. Implemented: RTN,GSE')
    
    tbn = np.arccos( np.dot(Bu,n)/np.linalg.norm(Bu))*180/np.pi
    if tbn > 90:
        tbn = 180 -tbn
        
    return n,tbn




def MX_stats(tB, B, tV, V, shock_time, up_shk, dw_shk, min_up_dur, max_up_dur, min_dw_dur, max_dw_dur, tcad, frame, method='cross'):

    """
    Routine computing shock normal vector, theta_bn and magnetic compression ratio 
    for an ensemble of upstream/downstream averaging windows, that are systematically changed.
    
    Parameters
    ----------
    tB : 'DatetimeIndex'
        Stream of times for magnetic field measurements.
    B : 'numpy array'
        Magnetic field measurements with n_timestamps x 3 dimensions.
    tV : 'DatetimeIndex'
        Stream of times for bulk flow speed measurements.
    V : 'numpy array'
        Bulk flow speed measurements with n_timestamps x 3 dimensions.
    shock_time : 'datetime'
        Time of the shock crossing
    up_shk : 'datetime'
        Beginning of the shock upstream.
    dw_shk : 'datetime'
        Beginning of the shock downstream.
    min_up_dur : 'timedelta'
        Duration of the smallest upstream averaging window
    max_up_dur : 'timedelta'
        Duration of the largest upstream averaging window
    min_dw_dur : 'timedelta'
        Duration of the smallest downstream averaging window.
    max_dw_dur : 'timedelta'
        Duration of the largest downstream averaging window.
    tcad : 'timedelta'
        Time cadence at which windows are enlarged.
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'
    method : 'str', optional
        Method by which different windows are considered. The default is 'cross'.

    Returns
    -------
    n : 'ShockParameters'
        Object with one attribute per technique (MC, Mx1, MX2, MX3).
        Each attribute is an array of shock normal vectors with dimensions
        n_windows_combination x 3
    tbn : 'ShockParameters'
        Object with one attribute per technique (MC, Mx1, MX2, MX3).
        Each attribute is an array of shock normal vectors with dimensions
        n_windows_combination x 3
    rB : 'numpy array'
        Array containing magnetic compression ratio computed per each window choice
    ex : 'ShockParameters'
        Values of thetabn and magnetic compression ratio obtained witht the 
        smallest and largest possible choice of upstream/downstream windows.
    """
    n   = ShockParameters()
    tbn = ShockParameters()
    rB  = ShockParameters()
    ex  = ShockParameters()
    
    tuf = shock_time - up_shk # Duration. Start of upstream
    tu1 = tuf + min_up_dur    # Duration of first window
    tu2 = tuf + max_up_dur    # Duration of last window 
    sldu = int(np.floor((tu2-tu1)/tcad)) 
    
    tdi = dw_shk - shock_time  # Duration. Start of downstream 
    td1 = tdi + min_dw_dur     # Duration of first window
    td2 = tdi + max_dw_dur     # Duration of last window 
    sldd = int(np.floor((td2-td1)/tcad))
    #disp([sldu,sldd]) 
    iw = 0
    
    n.MC   = np.zeros([sldu*sldd,3])
    n.MX1  = np.zeros([sldu*sldd,3])
    n.MX2  = np.zeros([sldu*sldd,3])
    n.MX3  = np.zeros([sldu*sldd,3])
    tbn.MC = np.zeros([sldu*sldd,1])
    tbn.MX1= np.zeros([sldu*sldd,1])    
    tbn.MX2= np.zeros([sldu*sldd,1])
    tbn.MX3= np.zeros([sldu*sldd,1])
    rB     = np.zeros([sldu*sldd,1]) 
    ex.tbn_susd = []
    ex.rB_susd = []
    ex.tbn_luld = []
    ex.rB_luld = []
    for i in range(sldu):
        stut = shock_time- (tu1 + (i)*tcad)
        enut = shock_time - tuf
        Bitsu, Biteu = get_time_indices(stut, enut, tB)
        tsb, sBu     = select_subS(B, tB, Bitsu, Biteu, 3)

        Vitsu, Viteu = get_time_indices(stut, enut, tV)
        tsv, sVu     = select_subS(V, tV, Vitsu, Viteu, 3)


        Bu = np.nanmean(sBu,axis=0)
        Vu = np.nanmean(sVu,axis=0)   
        
        
        for j in range(sldd):
            stdt = shock_time + tdi 
            endt = shock_time + (td1 + (j)*tcad)
            
            Bitsd, Bited = get_time_indices(stdt, endt, tB)
            tsb, sBd     = select_subS(B, tB, Bitsd, Bited, 3)
            Vitsd, Vited = get_time_indices(stdt, endt, tV)
            tsv, sVd     = select_subS(V, tV, Vitsd, Vited, 3)

            Bd = np.nanmean(sBd,axis=0)
            Vd = np.nanmean(sVd,axis=0)   
            
            rB[iw] = np.linalg.norm(Bd)/np.linalg.norm(Bu)
                           
            n.MC[iw,:], tbn.MC[iw] = calc_MC(Bu,Bd,frame)
            n.MX1[iw,:], tbn.MX1[iw] = calc_MX1(Bu,Bd,Vu,Vd,frame)
            n.MX2[iw,:], tbn.MX2[iw] = calc_MX2(Bu,Bd,Vu,Vd,frame)
            n.MX3[iw,:], tbn.MX3[iw] = calc_MX3(Bu,Bd,Vu,Vd,frame)            
            
            
            #Special cases
            if (i == 0 and j == 0):
                ex.tbn_susd = tbn.MX3[iw]  # Smallest Upstream Smallest Downstream
                ex.rB_susd  = rB[iw]   # Smallest Upstream Smallest Downstream

            iw = iw +1

        print('Upstream windows = ' + str(i) + ' / ' + str(sldu))
        
    ex.tbn_luld = tbn.MX3[iw-1]  # Smallest Upstream Smallest Downstream
    ex.rB_luld  = rB[iw-1]   # Smallest Upstream Smallest Downstream
    
    return n, tbn, rB, ex


def Vsh_stats(n, tP, V, Rho, shock_time, up_shk, dw_shk, min_up_dur, max_up_dur, min_dw_dur, max_dw_dur, tcad, method='cross'):
    """
    Routine that computes shock speed along the shock normal direction for an ensemble
    of upstream/downstream averaging windows, that are systematically changed.

    Parameters
    ----------
    n : 'numpy array'
        Shock normal vector.
    tP : 'DatetimeIndex'
        Stream of times for plasma measurements.
    V : 'numpy array'
        Bulk flow speed measurements with dimensions n_timestamps x 3.
    Rho : 'numpy array'
        Plasma density measurements with n_timestamps x 3 dimensions.
    shock_time : 'datetime'
        Time of the shock crossing
    up_shk : 'datetime'
        Beginning of the shock upstream.
    dw_shk : 'datetime'
        Beginning of the shock downstream.
    min_up_dur : 'timedelta'
        Duration of the smallest upstream averaging window
    max_up_dur : 'timedelta'
        Duration of the largest upstream averaging window
    min_dw_dur : 'timedelta'
        Duration of the smallest downstream averaging window.
    max_dw_dur : 'timedelta'
        Duration of the largest downstream averaging window.
    tcad : 'timedelta'
        Time cadence at which windows are enlarged.
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'
    method : 'str', optional
        Method by which different windows are considered. The default is 'cross'.
        
    Returns
    -------
    vsh : 'numpy array'
        Array containing shock speed computed per each window choice.
    ex : 'ShockParameters'
        Object containing values of shock speed obtained with the 
        smallest and largest possible choice of upstream/downstream windows.
    """
    vsh  = ShockParameters()
    ex  = ShockParameters()
    
    tuf = shock_time - up_shk # Duration. Start of upstream
    tu1 = tuf + min_up_dur    # Duration of first window
    tu2 = tuf + max_up_dur    # Duration of last window 
    sldu = int(np.floor((tu2-tu1)/tcad)) 
    
    tdi = dw_shk - shock_time  # Duration. Start of downstream 
    td1 = tdi + min_dw_dur     # Duration of first window
    td2 = tdi + max_dw_dur     # Duration of last window 
    sldd = int(np.floor((td2-td1)/tcad))
    #disp([sldu,sldd]) 
    iw = 0
    
    vsh   = np.zeros([sldu*sldd,1])
    ex.vsh_min = []
    ex.vsh_max = []
    
    for i in range(sldu):
        stut = shock_time- (tu1 + (i)*tcad)
        enut = shock_time - tuf
        #Bitsu, Biteu = get_time_indices(stut, enut, tB)
        #tsb, sBu     = select_subS(B, tB, Bitsu, Biteu, 3)

        Vitsu, Viteu   = get_time_indices(stut, enut, tP)
        tsv, sVu       = select_subS(V, tP, Vitsu, Viteu, 3)
        tsr, sRhou     = select_subS(Rho, tP, Vitsu, Viteu, 1)

        Rhou = np.nanmean(sRhou,axis=0)
        Vu = np.nanmean(sVu,axis=0)   
        
        
        for j in range(sldd):
            stdt = shock_time + tdi 
            endt = shock_time + (td1 + (j)*tcad)
            
            Vitsd, Vited   = get_time_indices(stdt, endt, tP)
            tsv, sVd       = select_subS(V, tP, Vitsd, Vited, 3)
            tsr, sRhod     = select_subS(Rho, tP, Vitsd, Vited, 1)
            
            Rhod = np.nanmean(sRhod,axis=0)
            Vd = np.nanmean(sVd,axis=0)   
            
            DRHO = Rhod - Rhou;
            rdVdn  = Rhod*np.dot(Vd,n);
            ruVun  = Rhou*np.dot(Vu,n);                
            vsh[iw]  = (rdVdn - ruVun)/DRHO;
            
            #Special cases
            if (i == 0 and j == 0):
                ex.vshmin = vsh[iw]  # Smallest Upstream Smallest Downstream
                
            iw = iw +1

        print('Upstream windows = ' + str(i) + ' / ' + str(sldu))
        
    ex.vshmax = vsh[iw-1]  # Smallest Upstream Smallest Downstream
    
    return vsh, ex

def rgas_stats(trho, Rho, shock_time, up_shk, dw_shk, min_up_dur, max_up_dur, min_dw_dur, max_dw_dur, tcad, method='cross'):
    """
    Routine computing the shock gas compression ratio for an ensemble
    of upstream/downstream averaging windows, that are systematically changed.

    Parameters
    ----------
    trho : 'DatetimeIndex'
        Stream of times for plasma measurements.
    Rho : 'numpy array'
        Plasma density measurements with n_timestamps x 3 dimensions.
    shock_time : 'datetime'
        Time of the shock crossing
    up_shk : 'datetime'
        Beginning of the shock upstream.
    dw_shk : 'datetime'
        Beginning of the shock downstream.
    min_up_dur : 'timedelta'
        Duration of the smallest upstream averaging window
    max_up_dur : 'timedelta'
        Duration of the largest upstream averaging window
    min_dw_dur : 'timedelta'
        Duration of the smallest downstream averaging window.
    max_dw_dur : 'timedelta'
        Duration of the largest downstream averaging window.
    tcad : 'timedelta'
        Time cadence at which windows are enlarged.
    frame : 'str'
        Frame in which the calculation is done. implemented: 'RTN', 'GSE'
    method : 'str', optional
        Method by which different windows are considered. The default is 'cross'.

    Returns
    -------
    r : 'numpy array'
        Array containing shock gas compression ratio computed per each window choice.
    ex : 'ShockParameters'
        Object containing values of gas compression ratio obtained with the 
        smallest and largest possible choice of upstream/downstream windows.

    """
    
    r  = ShockParameters()
    ex  = ShockParameters()
    
    tuf = shock_time - up_shk # Duration. Start of upstream
    tu1 = tuf + min_up_dur    # Duration of first window
    tu2 = tuf + max_up_dur    # Duration of last window 
    sldu = int(np.floor((tu2-tu1)/tcad)) 
    
    tdi = dw_shk - shock_time  # Duration. Start of downstream 
    td1 = tdi + min_dw_dur     # Duration of first window
    td2 = tdi + max_dw_dur     # Duration of last window 
    sldd = int(np.floor((td2-td1)/tcad))
    #disp([sldu,sldd]) 
    iw = 0
    
    r   = np.zeros([sldu*sldd,1])
    ex.r_min = []
    ex.r_max = []
    
    for i in range(sldu):
        stut = shock_time- (tu1 + (i)*tcad)
        enut = shock_time - tuf
        #Bitsu, Biteu = get_time_indices(stut, enut, tB)
        #tsb, sBu     = select_subS(B, tB, Bitsu, Biteu, 3)

        Vitsu, Viteu   = get_time_indices(stut, enut, trho)
        tsr, sRhou     = select_subS(Rho, trho, Vitsu, Viteu, 1)

        Rhou = np.nanmean(sRhou,axis=0)
        
        
        for j in range(sldd):
            stdt = shock_time + tdi 
            endt = shock_time + (td1 + (j)*tcad)
            
            Vitsd, Vited   = get_time_indices(stdt, endt,  trho)
            tsr, sRhod     = select_subS(Rho,  trho, Vitsd, Vited, 1)
            
            Rhod = np.nanmean(sRhod,axis=0)
            
                         
            r[iw]  = Rhod/Rhou;
            
            #Special cases
            if (i == 0 and j == 0):
                ex.r_min = r[iw]  # Smallest Upstream Smallest Downstream
                
            iw = iw +1

        print('Upstream windows = ' + str(i) + ' / ' + str(sldu))
        
    ex.r_max = r[iw-1]  # Smallest Upstream Smallest Downstream
    
    return r, ex

#% Functions to rankinise
def compute_Rankine_condition(Bx1, By1):
    """
    Routine that computes Rankine-Hugoniot relations

    Parameters
    ----------
    Bx1 : 'float64'
        Bx vlue upstream.
    By1 : 'float64'
        By value downstream.

    Returns
    -------
    Bx1,Bx2,By1,By2,Vx1,Vx2,Vy1,Vy2,rho1,rho2: 'float64'
        Magnetic field, bulk flow speed and density values upstream and dowsntream

    """
    
    gam = 5/3
    r   = 3.8
    theta1 = 0 # Angle between bulk flow speed and shock normal (0 is NIF)
    rho1 = 1
    Vx1  = 200
    Vy1  = 0
    P1   = 10 
    #rm   = (gam+1)/(gam-1)
    Bm1  = np.sqrt(Bx1**2+By1**2)
    va1 = np.sqrt((Bm1**2)/rho1)
    vs1 = np.sqrt(gam*P1/rho1)
    v1sq  = (2*r*vs1**2)/((gam +1) - (gam-1)*r) 
    rho2 = r*rho1
    Bx2  = Bx1
    nom = v1sq - np.cos(theta1)**2*va1**2
    den = v1sq - r*np.cos(theta1)**2*va1**2
    By2 = By1*(r*nom/den)
    Vx2 = Vx1/r
    Vy2 = Vy1*nom/den
    
    return Bx1,Bx2,By1,By2,Vx1,Vx2,Vy1,Vy2,rho1,rho2

def create_ranki_stream(Bx1,Bx2,By1,By2,Vx1,Vx2,Vy1,Vy2,rho1,rho2,nts, mode, noise_level = 0.):
    """
    Create Stream of data compliant with Rankine-Hugoniot jump conditions

    Parameters
    ----------
    Bx1 : 'float64'
        Value of Bx upstream.
    Bx2 : 'float64'
        Value of Bx downstream.
    By1 : 'float64'
        Value of By upstream.
    By2 : 'float64'
        Value of By downstream.
    Vx1 : 'float64'
        Value of Vx upstream.
    Vx2 : 'float64'
        Value of Vx downstream.
    Vy1 : 'float64'
        Value of Vy upstream.
    Vy2 : 'float64'
        Value of Vy downstream.
    rho1 : 'float64'
        Value of plasma density upstream.
    rho2 : 'float64'
        Value of plasma density downstream.
    nts : 'int'
        Number of timestamps needed.
    mode : 'str'
        Choose if other signals need to be superimposed to the RH compliant one.
        Implemented: 'clean': no superimpositions, 'white_noise': include white noise
    noise_level : 'float64', optional
        Level of white noise desired. The default is 0..

    Returns
    -------
    B : 'numpy array'
        Array containing magnetic field synthetic measurements, dimensions nts x 4 (components+magnitude)
    V : 'numpy array'
        Array containing bulk flow speed synthetic measurements, dimensions nts x 4 (components+magnitude).
    rho : 'numpy array'
        Array containing plasma density synthetic measurements, dimensions nts x 1.

    """
    B = np.zeros([nts,4])
    V = np.zeros([nts,4])
    rho = np.zeros([nts,1])
    nts2 = int(nts/2)
    if(mode == 'clean'):
        for i in range(nts2):
            B[i,0] = Bx1
            B[i,1] = By1
            B[i,2] = 0.0
            B[i,3] = np.sqrt(Bx1**2+By1**2)
            
            B[i+nts2,0] = Bx2
            B[i+nts2,1] = By2
            B[i+nts2,2] = 0.0
            B[i+nts2,3] = np.sqrt(Bx2**2+By2**2)
            
            V[i,0] = Vx1
            V[i,1] = Vy1
            V[i,2] = 0.0
            V[i,3] = np.sqrt(Vx1**2+Vy1**2)
            
            V[i+nts2,0] = Vx2
            V[i+nts2,1] = Vy2
            V[i+nts2,2] = 0.0
            V[i+nts2,3] = np.sqrt(Vx2**2+Vy2**2)
            
            rho[i] = rho1
            rho[i+nts2] = rho2
            
    if(mode == 'white_noise'):
        for i in range(nts2):
            B[i,0] = Bx1 + np.random.normal(0. , Bx1*noise_level)
            B[i,1] = By1 + np.random.normal(0. , By1*noise_level)
            B[i,2] = 0.0 + np.random.normal(0. , noise_level)
            B[i,3] = np.sqrt(B[i,0]**2+B[i,1]**2+B[i,2]**2)
            
            B[i+nts2,0] = Bx2 + np.random.normal(0. , Bx2*noise_level)
            B[i+nts2,1] = By2 + np.random.normal(0. , By2*noise_level)
            B[i+nts2,2] = 0.0 + np.random.normal(0. , noise_level)
            B[i+nts2,3] = np.sqrt(B[i+nts2,0]**2+B[i+nts2,1]**2+B[i+nts2,2]**2)
            
            V[i,0] = Vx1 + np.random.normal(0. , Vx1*noise_level)
            V[i,1] = Vy1 + np.random.normal(0. , Vy1*noise_level)
            V[i,2] = 0.0 + np.random.normal(0. , noise_level)
            V[i,3] = np.sqrt(V[i,0]**2+V[i,1]**2+V[i,2]**2)
            
            V[i+nts2,0] = Vx2 + np.random.normal(0. , Vx2*noise_level)
            V[i+nts2,1] = Vy2 + np.random.normal(0. , Vy2*noise_level)
            V[i+nts2,2] = 0.0 + np.random.normal(0. , noise_level)
            V[i+nts2,3] = np.sqrt(V[i+nts2,0]**2+V[i+nts2,1]**2+V[i+nts2,2]**2)
            
            rho[i] = rho1 + np.random.normal(0. , rho1*noise_level)
            rho[i+nts2] = rho2 + np.random.normal(0. , rho2*noise_level)


    return B,V,rho


