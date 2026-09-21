"""
Reproducible computational implementation for the Diabetic Kidney Latent Reserve Model.

This script reproduces the theoretical model analyses reported in the accompanying
manuscript: reference and intervention simulations, equilibrium/Jacobian analysis,
continuation, a two-parameter risk surface, Latin hypercube uncertainty analysis,
PRCC, local elasticity, solver robustness, basin analysis, figures, and CSV outputs.

The model is a literature-parameterized theoretical proof-of-concept and is not a
clinically calibrated diagnostic, prognostic, or treatment tool.
"""

import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import root, brentq
from scipy.stats import qmc, rankdata, pearsonr
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

# Repository-relative paths. This script can be run from any working directory.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = REPO_ROOT / 'data' / 'outputs'
FIG_DIR = REPO_ROOT / 'figures'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

MASTER_CSV = OUTPUT_DIR / 'renal_risk_surface_results.csv'

# ---------------- MODEL ----------------
# Dimensionless states, all constrained to [0,1]
# G glycaemic burden; A albuminuric injury; C cystatin-C discordance burden;
# F frailty burden; M metabolic-inflammatory stress; R renal reserve potential.
p0=dict(
    ug=0.23, kg=0.34,
    ag=0.22, am=0.08, af=0.05, ka=0.28, era=0.25,
    ca=0.18, cm=0.10, kc=0.32, erc=0.20,
    fg=0.08, fc=0.12, fa=0.09, kf=0.24, erf=0.16,
    mg=0.10, mf=0.11, km=0.30, erm=0.18,
    rho=0.20, dg=0.10, da=0.22, dc=0.16, df=0.18, dm=0.12,
    ur=0.00
)
y0=np.array([0.28,0.16,0.10,0.14,0.12,0.74],float)
state_names=['Glycaemic burden G','Albuminuric injury A','Cystatin-C discordance C','Frailty burden F','Metabolic stress M','Renal reserve R']


def rhs(t,y,p):
    G,A,C,F,M,R=y
    dG=p['ug']*(1-G)-p['kg']*G
    dA=(p['ag']*G+p['am']*M+p['af']*F)*(1-A)-p['ka']*A-p['era']*R*A
    dC=(p['ca']*A+p['cm']*M)*(1-C)-p['kc']*C-p['erc']*R*C
    dF=(p['fg']*G+p['fc']*C+p['fa']*A)*(1-F)-p['kf']*F-p['erf']*R*F
    dM=(p['mg']*G+p['mf']*F)*(1-M)-p['km']*M-p['erm']*R*M
    damage=p['dg']*G+p['da']*A+p['dc']*C+p['df']*F+p['dm']*M
    dR=(p['rho']+p['ur'])*(1-R)-damage*R
    return np.array([dG,dA,dC,dF,dM,dR])


def simulate(p=None,years=12,yinit=None,rtol=1e-9,atol=1e-11,npts=601):
    if p is None: p=p0
    if yinit is None: yinit=y0
    t=np.linspace(0,years,npts)
    sol=solve_ivp(lambda tt,yy: rhs(tt,yy,p),(0,years),yinit,t_eval=t,method='RK45',rtol=rtol,atol=atol)
    if not sol.success: raise RuntimeError(sol.message)
    return sol

weights=np.array([0.15,0.25,0.20,0.18,0.12])
def derived(sol):
    G,A,C,F,M,R=sol.y
    reserve=R-(weights[0]*G+weights[1]*A+weights[2]*C+weights[3]*F+weights[4]*M)
    surface=0.14*G+0.24*A+0.18*C+0.18*F+0.12*M+0.10*A*F+0.08*C*F-0.50*R
    hazard=1/(1+np.exp(-8*(surface+0.10)))
    return reserve,surface,hazard

def regime(q):
    if q>=0.45: return 'Preserved'
    if q>=0.25: return 'Constrained'
    return 'Depleted'

def metric_row(sol):
    q,s,h=derived(sol)
    vals=sol.y[:,-1]
    return dict(final_G=vals[0],final_A=vals[1],final_C=vals[2],final_F=vals[3],final_M=vals[4],final_R=vals[5],
                reserve_potential=q[-1],risk_surface=s[-1],hazard=h[-1],regime=regime(q[-1]),minimum_state=float(sol.y.min()),maximum_state=float(sol.y.max()))

# Equilibrium
long=simulate(years=100,npts=2001)
xguess=long.y[:,-1]
eq=root(lambda x: rhs(0,x,p0),xguess)
xeq=eq.x

def num_jac(x,p,eps=1e-6):
    n=len(x); J=np.zeros((n,n))
    for j in range(n):
        d=np.zeros(n); d[j]=eps
        J[:,j]=(rhs(0,x+d,p)-rhs(0,x-d,p))/(2*eps)
    return J
J=num_jac(xeq,p0); eig=np.linalg.eigvals(J)

# Main scenarios
scenarios={
    'Reference':{},
    'Glycaemic input reduced 30%':{'ug':p0['ug']*0.70},
    'Albuminuric propagation reduced 35%':{'ag':p0['ag']*0.65},
    'Frailty coupling reduced 30%':{'fg':p0['fg']*0.70,'fc':p0['fc']*0.70,'fa':p0['fa']*0.70},
    'Reserve restoration added':{'ur':0.08},
    'Combined multi-domain strategy':{'ug':p0['ug']*0.70,'ag':p0['ag']*0.65,'fg':p0['fg']*0.70,'fc':p0['fc']*0.70,'fa':p0['fa']*0.70,'ur':0.08}
}
scenario_rows=[]; sols={}
for label,mods in scenarios.items():
    p=p0.copy();p.update(mods); sol=simulate(p); sols[label]=sol
    scenario_rows.append({'section':'scenario','item':label,**metric_row(sol)})
scenario_df=pd.DataFrame(scenario_rows)
ref=scenario_df.iloc[0]
scenario_df['hazard_reduction_pct']=100*(ref['hazard']-scenario_df['hazard'])/ref['hazard']
scenario_df['reserve_gain_pct']=100*(scenario_df['reserve_potential']-ref['reserve_potential'])/abs(ref['reserve_potential'])

# Continuation in glycaemic input and albuminuric propagation
ugs=np.linspace(0.08,0.90,165); cont=[]
for ug in ugs:
    sol=simulate({**p0,'ug':float(ug)},npts=401)
    m=metric_row(sol); cont.append({'ug':ug,**m})
cont_df=pd.DataFrame(cont)
crit_q25=brentq(lambda u: metric_row(simulate({**p0,'ug':u},npts=401))['reserve_potential']-0.25,0.08,0.90)
crit_h60=brentq(lambda u: metric_row(simulate({**p0,'ug':u},npts=401))['hazard']-0.60,0.08,0.90)

# 2D sensitivity surface in ug and ag
UG=np.linspace(0.10,0.55,46); AG=np.linspace(0.10,0.34,46); Z=np.zeros((len(AG),len(UG))); H=np.zeros_like(Z)
for i,ag in enumerate(AG):
    for j,ug in enumerate(UG):
        m=metric_row(simulate({**p0,'ug':float(ug),'ag':float(ag)},years=10,npts=301,rtol=3e-8,atol=1e-10))
        Z[i,j]=m['reserve_potential']; H[i,j]=m['hazard']

# Global uncertainty and PRCC
ranges={
'ug':(0.14,0.36),'kg':(0.26,0.45),'ag':(0.14,0.32),'ka':(0.20,0.38),'ca':(0.10,0.27),'kc':(0.24,0.42),
'fg':(0.04,0.14),'kf':(0.16,0.34),'rho':(0.12,0.30),'da':(0.14,0.30),'dc':(0.10,0.24),'df':(0.10,0.26)
}
N=2500; names=list(ranges)
sampler=qmc.LatinHypercube(d=len(names),seed=20260727); u=sampler.random(N); X=np.zeros_like(u)
for j,n in enumerate(names):
    lo,hi=ranges[n];X[:,j]=lo+u[:,j]*(hi-lo)
Y=[]; failures=0
for row in X:
    pp=p0.copy();
    for j,n in enumerate(names): pp[n]=float(row[j])
    try:
        m=metric_row(simulate(pp,years=10,npts=251,rtol=5e-7,atol=1e-9));Y.append([m['reserve_potential'],m['hazard'],m['final_R']])
    except Exception:
        failures+=1;Y.append([np.nan,np.nan,np.nan])
Y=np.array(Y); valid=np.isfinite(Y).all(axis=1); Xv=X[valid];Yv=Y[valid]

def prcc(X,y):
    Xr=np.apply_along_axis(rankdata,0,X); yr=rankdata(y); vals=[]
    for j in range(X.shape[1]):
        idx=[k for k in range(X.shape[1]) if k!=j]
        rx=Xr[:,j]-LinearRegression().fit(Xr[:,idx],Xr[:,j]).predict(Xr[:,idx])
        ry=yr-LinearRegression().fit(Xr[:,idx],yr).predict(Xr[:,idx])
        vals.append(pearsonr(rx,ry)[0])
    return np.array(vals)
pr_q=prcc(Xv,Yv[:,0]);pr_h=prcc(Xv,Yv[:,1])
prcc_df=pd.DataFrame({'parameter':names,'PRCC_reserve':pr_q,'PRCC_hazard':pr_h})
prcc_df=prcc_df.reindex(prcc_df.PRCC_reserve.abs().sort_values(ascending=False).index)

# local elasticity
base=metric_row(simulate())
elastic=[]
for n in names:
    p1=p0.copy();p2=p0.copy();p1[n]*=1.01;p2[n]*=.99
    m1=metric_row(simulate(p1));m2=metric_row(simulate(p2))
    e_q=(math.log(abs(m1['reserve_potential']))-math.log(abs(m2['reserve_potential'])))/(math.log(p1[n])-math.log(p2[n]))
    e_h=(math.log(m1['hazard'])-math.log(m2['hazard']))/(math.log(p1[n])-math.log(p2[n]))
    elastic.append({'parameter':n,'elasticity_reserve':e_q,'elasticity_hazard':e_h})
elastic_df=pd.DataFrame(elastic).sort_values('elasticity_reserve',key=np.abs,ascending=False)

# solver robustness
rob=[]
for rtol,atol in [(1e-6,1e-8),(1e-8,1e-10),(1e-10,1e-12)]:
    m=metric_row(simulate(rtol=rtol,atol=atol));rob.append({'rtol':rtol,'atol':atol,**m})
rob_df=pd.DataFrame(rob)

# Initial-state basin grid A0,F0
A0s=np.linspace(.05,.55,31);F0s=np.linspace(.05,.55,31); basin=[]
for a0 in A0s:
    for f0 in F0s:
        yi=y0.copy();yi[1]=a0;yi[3]=f0
        m=metric_row(simulate(yinit=yi,years=12,npts=301,rtol=3e-8,atol=1e-10))
        basin.append({'A0':a0,'F0':f0,'reserve_potential':m['reserve_potential'],'hazard':m['hazard'],'regime':m['regime']})
basin_df=pd.DataFrame(basin)

# Figures
plt.figure(figsize=(7.4,4.8))
for label in ['Reference','Glycaemic input reduced 30%','Albuminuric propagation reduced 35%','Combined multi-domain strategy']:
    sol=sols[label]; q,_,_=derived(sol);plt.plot(sol.t,q,label=label)
plt.axhline(.45,ls='--',lw=.8);plt.axhline(.25,ls='--',lw=.8)
plt.xlabel('Time (model years)');plt.ylabel('Dimensionless renal reserve potential');plt.legend(fontsize=8);plt.tight_layout()
fig1=FIG_DIR / 'figure1_reserve_trajectories.png';plt.savefig(fig1,dpi=220);plt.close()

plt.figure(figsize=(7.2,4.7));plt.plot(cont_df.ug,cont_df.reserve_potential)
plt.axhline(.25,ls='--',lw=.8);plt.axvline(crit_q25,ls='--',lw=.8)
plt.xlabel('Dimensionless glycaemic input u_g');plt.ylabel('Reserve potential at year 12');plt.tight_layout()
fig2=FIG_DIR / 'figure2_continuation.png';plt.savefig(fig2,dpi=220);plt.close()

plt.figure(figsize=(7.2,4.8));cs=plt.contourf(UG,AG,Z,levels=18);plt.colorbar(cs,label='Reserve potential at year 10')
plt.contour(UG,AG,Z,levels=[.25,.45],linewidths=1)
plt.xlabel('Glycaemic input u_g');plt.ylabel('Albuminuric propagation a_g');plt.tight_layout()
fig3=FIG_DIR / 'figure3_risk_surface.png';plt.savefig(fig3,dpi=220);plt.close()

plt.figure(figsize=(7.2,5));plotdf=prcc_df.sort_values('PRCC_reserve');plt.barh(plotdf.parameter,plotdf.PRCC_reserve)
plt.axvline(0,lw=.8);plt.xlabel('PRCC with renal reserve potential');plt.tight_layout()
fig4=FIG_DIR / 'figure4_prcc.png';plt.savefig(fig4,dpi=220);plt.close()

pivot=basin_df.pivot(index='F0',columns='A0',values='reserve_potential')
plt.figure(figsize=(7.2,5));cs=plt.contourf(pivot.columns,pivot.index,pivot.values,levels=18);plt.colorbar(cs,label='Reserve potential at year 12')
plt.contour(pivot.columns,pivot.index,pivot.values,levels=[.25,.45],linewidths=1)
plt.xlabel('Initial albuminuric injury A(0)');plt.ylabel('Initial frailty burden F(0)');plt.tight_layout()
fig5=FIG_DIR / 'figure5_basin.png';plt.savefig(fig5,dpi=220);plt.close()

# combined results CSV long format
frames=[]
for name,df in [('scenario',scenario_df),('continuation',cont_df),('prcc',prcc_df),('elasticity',elastic_df),('robustness',rob_df),('basin',basin_df)]:
    d=df.copy();d.insert(0,'table',name);frames.append(d.astype(object))
master=pd.concat(frames,ignore_index=True,sort=False)
master.to_csv(MASTER_CSV,index=False)

# Save standalone tables
scenario_df.to_csv(OUTPUT_DIR / 'scenario_results.csv',index=False)
prcc_df.to_csv(OUTPUT_DIR / 'global_sensitivity_prcc.csv',index=False)
elastic_df.to_csv(OUTPUT_DIR / 'local_elasticity.csv',index=False)
cont_df.to_csv(OUTPUT_DIR / 'continuation_results.csv',index=False)
basin_df.to_csv(OUTPUT_DIR / 'basin_results.csv',index=False)


# Console summary
print('Analysis complete.')
print(f'Repository root: {REPO_ROOT}')
print(f'Outputs: {OUTPUT_DIR}')
print(f'Figures: {FIG_DIR}')
print(f'Q=0.25 continuation boundary: u_g={crit_q25:.6f}')
print(f'H=0.60 continuation boundary: u_g={crit_h60:.6f}')
print(f'Dominant Jacobian real part: {float(np.max(eig.real)):.6f}')
print(f'Valid Latin hypercube runs: {len(Yv)} / {N}')
