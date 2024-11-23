#%%
import numpy as np
import sympy as sp
import matplotlib.pyplot as mpl

sp.init_printing( use_latex='mathjax' )

t = sp.var('t')
f = sp.Function('f')
diffeq = sp.Eq(sp.harmonic(f(t)), sp.diff(f(t), t))

def euler(q, dt, ti, tf, qi, *argv):
    dfdt = sp.solve(q, sp.diff(f(t), t))
    dfdt_eval = lambda x: sp.N(dfdt[0].subs(f(t), x))

    steps = int((tf-ti)/dt)
    f_approx = np.empty(steps+1)
    f_approx[0] = qi

    for i in range(steps):
        f_approx[i+1] = f_approx[i] + (dfdt_eval(f_approx[i]) * dt)
    
    return(f_approx)

#%%
eul = euler(diffeq, 0.00001, 0, 10, 1)
mpl.plot(np.arange(0, 10.00001, 0.00001), eul)

print(eul)