# -*- coding: utf-8 -*-
"""
The Rankinisator - shock overview

Show shock diagnostics for three clean streams and let user define theta bns. 

@author: trotta
"""
#%% Init and import
import sys

sys.path.append(r'..\\SerPyShock')

                
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import math
import SerPyShock 



#%% Generate synthetic stream
tot_dur = dt.timedelta( hours   = 24)
res_dur = dt.timedelta( seconds = 1 )
nts = int(tot_dur/res_dur)
shock_time = dt.datetime(1991,6,1,17,15,00)
start_time = shock_time - tot_dur/2

t = []

for i in range(nts):
    t.append(start_time + i*res_dur)


teta = 25
Bx1,Bx2,By1,By2,Vx1,Vx2,Vy1,Vy2,rho1,rho2 = SerPyShock.compute_Rankine_condition(np.cos(math.radians(teta)),np.sin(math.radians(teta)))
B_s2,V_s2,rho_s2 = SerPyShock.create_ranki_stream(Bx1,Bx2,By1,By2,Vx1,Vx2,Vy1,Vy2,rho1,rho2,nts, 'clean')


#%% Prototype selector (This code could be much better)

#st_ups  = shock_time - dt.timedelta(hours = 6)
#end_ups = shock_time - dt.timedelta(minutes = 1)
tm = np.array(t)

# Select overview
st_ovi_t  = shock_time  - dt.timedelta(minutes = 10)
end_ovi_t = shock_time  + dt.timedelta(minutes = 10)
Bbitso, Bbiteo = SerPyShock.get_time_indices(st_ovi_t, end_ovi_t, tm)
toBb, sBbo     = SerPyShock.select_subS(B_s2, tm, Bbitso, Bbiteo, 4)
toBb, sVbo     = SerPyShock.select_subS(V_s2, tm, Bbitso, Bbiteo, 4)
toBb, rhoo     = SerPyShock.select_subS(rho_s2, tm, Bbitso, Bbiteo, 1)
print('Overview selected');


# Select Upstream
st_ups_t  = shock_time - (dt.timedelta(hours = 0) + dt.timedelta(minutes = 1))
end_ups_t = shock_time - (dt.timedelta(hours = 0) + dt.timedelta(minutes = 1))
Bitsu, Biteu = SerPyShock.get_time_indices(st_ups_t, end_ups_t, tm)
tuB, sBu     = SerPyShock.select_subS(B_s2, tm, Bitsu, Biteu, 4)
tuV, sVu     = SerPyShock.select_subS(V_s2, tm, Bitsu, Biteu, 4)
tuR, srhou     = SerPyShock.select_subS(rho_s2, tm, Bitsu, Biteu, 1)
print('Upstream selected');

# Select Upstream
st_dw_t  = shock_time + (dt.timedelta(hours = 0) + dt.timedelta(minutes = 1))
end_dw_t = shock_time + (dt.timedelta(hours = 0) + dt.timedelta(minutes = 1))
Bitsd, Bited = SerPyShock.get_time_indices(st_dw_t, end_dw_t, tm)
tdB, sBd     = SerPyShock.select_subS(B_s2, tm, Bitsd, Bited, 4)
tdV, sVd     = SerPyShock.select_subS(V_s2, tm, Bitsd, Bited, 4)
tdR, srhod     = SerPyShock.select_subS(rho_s2, tm, Bitsd, Bited, 1)
print('Downstream selected');
#%% Statistics parameters

up_shk = end_ups_t;
dw_shk = st_dw_t;
min_up_dur = dt.timedelta(minutes=2)
max_up_dur = dt.timedelta(minutes=6)

min_dw_dur =  dt.timedelta(minutes=1)
max_dw_dur =  dt.timedelta(minutes=8)

tcad =  dt.timedelta(seconds=8)


tuf = shock_time - up_shk # Duration. Start of upstream
tu1 = tuf + min_up_dur    # Duration of first window
tu2 = tuf + max_up_dur    # Duration of last window 
sldu = int(np.floor((tu2-tu1)/tcad)) 

tdi = dw_shk - shock_time  # Duration. Start of downstream 
td1 = tdi + min_dw_dur     # Duration of first window
td2 = tdi + max_dw_dur     # Duration of last window 
sldd = int(np.floor((td2-td1)/tcad))


stut_min = shock_time- (tu1 + (0)*tcad)
stut_max = shock_time- (tu1 + (sldu)*tcad)

enut = shock_time - tuf

stdt = shock_time + tdi 
endt_min = shock_time + (td1 + (0)*tcad)
endt_max = shock_time + (td1 + (sldd)*tcad)


#% Plot Overview of Ranki Stream
plt.clf()
plt.close('all')
plt.ion()
fig = plt.figure(figsize = (8,12))
fsax = 12 #fontsize ax 


gs = fig.add_gridspec(3,1)
#------------B MAGNITUDE---------------
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot( tm[Bbitso:Bbiteo], B_s2[Bbitso:Bbiteo,0], color = 'b', label = '$B_R$' )
ax1.plot( tm[Bbitso:Bbiteo], B_s2[Bbitso:Bbiteo,1], color = 'r', label = '$B_T$' )
ax1.plot( tm[Bbitso:Bbiteo], B_s2[Bbitso:Bbiteo,2], color = 'g', label = '$B_N$' )
ax1.plot( tm[Bbitso:Bbiteo], B_s2[Bbitso:Bbiteo,3], color = [0,0,0], label = '$B$' )
ax1.set_xlabel(r'Time',fontsize=fsax)
ax1.set_ylabel(r'$B$',fontsize=fsax)
ax1.axvspan(stut_min,enut, facecolor = 'b', alpha = 0.2)
ax1.axvspan(stut_max,enut, facecolor = 'b', alpha = 0.2)
#ax1.text(stut_min + dt.timedelta(seconds=10),.2,'Smallest Upstream Window', rotation='vertical')
#ax1.text(stut_max + dt.timedelta(seconds=10),.2,'Largest Upstream Window', rotation='vertical')
ax1.set_xlim(st_ovi_t,end_ovi_t)
ax1.axvspan(stdt,endt_min, facecolor   = 'r', alpha = 0.2)
ax1.axvspan(stdt,endt_max, facecolor   = 'r', alpha = 0.2)
#ax1.text(endt_min - dt.timedelta(seconds=30),.2,'Smallest Downstream Window', rotation='vertical')
#ax1.text(endt_max - dt.timedelta(seconds=30),.2,'Largest Downstream Window', rotation='vertical')
ax1.legend()



ax2 = fig.add_subplot(gs[1, 0])
ax2.plot( tm[Bbitso:Bbiteo], V_s2[Bbitso:Bbiteo,0], color = 'b', label = '$U_R$' )
ax2.plot( tm[Bbitso:Bbiteo], V_s2[Bbitso:Bbiteo,1], color = 'r', label = '$U_T$' )
ax2.plot( tm[Bbitso:Bbiteo], V_s2[Bbitso:Bbiteo,2], color = 'g', label = '$U_N$' )
ax2.plot( tm[Bbitso:Bbiteo], V_s2[Bbitso:Bbiteo,3], color = [0,0,0], label = '$U$' )
ax2.set_xlabel(r'Time',fontsize=fsax)
ax2.set_ylabel(r'$U$',fontsize=fsax)
ax2.axvspan(stut_min,enut, facecolor = 'b', alpha = 0.2)
ax2.axvspan(stut_max,enut, facecolor = 'b', alpha = 0.2)
#ax1.text(stut_min + dt.timedelta(seconds=10),.2,'Smallest Upstream Window', rotation='vertical')
#ax1.text(stut_max + dt.timedelta(seconds=10),.2,'Largest Upstream Window', rotation='vertical')
ax2.set_xlim(st_ovi_t,end_ovi_t)
ax2.axvspan(stdt,endt_min, facecolor   = 'r', alpha = 0.2)
ax2.axvspan(stdt,endt_max, facecolor   = 'r', alpha = 0.2)
#ax1.text(endt_min - dt.timedelta(seconds=30),.2,'Smallest Downstream Window', rotation='vertical')
#ax1.text(endt_max - dt.timedelta(seconds=30),.2,'Largest Downstream Window', rotation='vertical')
ax2.legend()



ax3 = fig.add_subplot(gs[2, 0])
ax3.plot( tm[Bbitso:Bbiteo], rho_s2[Bbitso:Bbiteo], color = 'b' )
ax3.set_xlabel(r'Time',fontsize=fsax)
ax3.set_ylabel(r'$n$',fontsize=fsax)
ax3.axvspan(st_ups_t,end_ups_t, facecolor = 'b', alpha = 0.2)
ax3.axvspan(st_dw_t,end_dw_t, facecolor   = 'r', alpha = 0.2)
ax3.axvspan(stut_min,enut, facecolor = 'b', alpha = 0.2)
ax3.axvspan(stut_max,enut, facecolor = 'b', alpha = 0.2)
ax3.text(stut_min + dt.timedelta(seconds=10),1.2,'Smallest Upstream Window', rotation='vertical')
ax3.text(stut_max + dt.timedelta(seconds=10),1.2,'Largest Upstream Window', rotation='vertical')
ax3.set_xlim(st_ovi_t,end_ovi_t)
ax3.axvspan(stdt,endt_min, facecolor   = 'r', alpha = 0.2)
ax3.axvspan(stdt,endt_max, facecolor   = 'r', alpha = 0.2)
ax3.text(endt_min - dt.timedelta(seconds=30),1.2,'Smallest Downstream Window', rotation='vertical')
ax3.text(endt_max - dt.timedelta(seconds=30),1.2,'Largest Downstream Window', rotation='vertical')

#plt.plot(t,B[:,1])
#plt.plot(t[its:ite],B[its:ite,1])
plt.tight_layout()
#plt.savefig('Fig1_oview_Ranki.jpg')
plt.show()
#%% Diagnostics

up_shk = end_ups_t;
dw_shk = st_dw_t;
min_up_dur = dt.timedelta(minutes=1)
max_up_dur = dt.timedelta(minutes=6)

min_dw_dur =  dt.timedelta(minutes=1)
max_dw_dur =  dt.timedelta(minutes=8)

tcad =  dt.timedelta(seconds=20)

n2, tbn2, rB2, ex2 = SerPyShock.MX_stats( tm, B_s2[:,0:3], tm, V_s2[:,0:3],  shock_time, up_shk, dw_shk, min_up_dur, max_up_dur, min_dw_dur, max_dw_dur, tcad, 'RTN');

#%% Plot Diagnostics

nbns = 5000
hist2,bin_edges = np.histogram(tbn2.MX3, bins = nbns, range=(0,90))

dtbn = bin_edges[1] - bin_edges[0]
PDF2 = hist2/(dtbn*tbn2.MC.shape[0])

#%% Plot
plt.clf()
plt.close('all')
plt.ion()
fig = plt.figure(figsize = (6,6))
fsax = 12 #fontsize ax 
CB = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
plt.plot(bin_edges[:-1],PDF2, color = CB[1], label = '$25^\circ$')
#plt.title('a')
plt.xlabel(r'$\theta_{Bn} [^\circ]$',fontsize=fsax)
plt.ylabel(r'PDF',fontsize=fsax)
plt.legend(fontsize=9,title = r'Nominal $ \theta_{Bn}$')

plt.xlim(0,90)
#plt.ylim(0,1.02)
#plt.savefig('Fig2_6shocks.jpg')
plt.show()



