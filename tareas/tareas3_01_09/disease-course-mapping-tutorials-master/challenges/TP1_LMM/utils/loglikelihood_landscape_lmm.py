import copy
from itertools import product
import inspect

from joblib import Parallel, delayed
import numpy as np
from statsmodels.regression.mixed_linear_model import MixedLMParams #MixedLM

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm #DivergingNorm
#from matplotlib.ticker import AutoMinorLocator

def cov_re_to_re_params(cov_re, k_fe=2, k_re=2):
    #assert len(cov_re) == k_re*(k_re+1) // 2

    # <!> we need to pass use_sqrt=False manually
    return MixedLMParams.from_packed(
        np.array(cov_re), k_fe=k_fe, k_re=k_re,
        use_sqrt=False, has_fe=False
    )

def ll_params(mlm_params, ll, profile_fe=True):
    return ll(mlm_params, profile_fe=profile_fe)

def ll_point(sx, sy, corr, model):

    # not sqrt
    assert model.k_re == 2
    cov_re = np.array([sx**2, corr*sx*sy, sy**2])

    # <!> we need to pass use_sqrt=False manually
    re_params = cov_re_to_re_params(cov_re, k_fe=model.k_fe, k_re=model.k_re)

    return ll_params(re_params, ll=model.loglike, profile_fe=True)

def ll_landscape(model, N=51, n_jobs=1,
                 #sx_min=.01, sx_max=50,
                 #sy_min=1e-2, sy_max=1.,
                 sx_min=.1, sx_max=20,
                 sy_min=.1, sy_max=20,
                 corr_min=-1+1e-3, corr_max=1-1e-3,
                 **fix):

    # 1 fixed params (2D)
    assert len(fix) == 1

    # X, Y
    leg = []
    M = []

    if 'sx' in fix:
        vx = fix['sx']**2
    else:
        leg.append('sx')
        M.append( np.linspace(sx_min, sx_max, N) )

    if 'sy' in fix:
        vy = fix['sy']**2
    else:
        leg.append('sy')
        M.append( np.linspace(sy_min, sy_max, N) )

    if 'corr' not in fix:
        leg.append('corr')
        M.append( np.linspace(corr_min, corr_max, N) )

    X,Y = np.meshgrid(M[0], M[1])
    #Z = np.zeros_like(X)

    COV = X*Y*next(iter(fix.values())) # element-wise # <!> *std
    X2 = X*X
    Y2 = Y*Y

    if 'corr' in fix:
        def xyp_to_re_params(x2,y2,cov):
            return cov_re_to_re_params([x2, cov, y2])
    elif 'sx' in fix:
        def xyp_to_re_params(x2,y2,cov):
            return cov_re_to_re_params([vx, cov, x2])
    else: # 'sy' in fix:
        def xyp_to_re_params(x2,y2,cov):
            return cov_re_to_re_params([x2, cov, vy])

    ll = model.loglike

    Z = Parallel(n_jobs=n_jobs)(
        delayed(ll_params)(xyp_to_re_params(X2[i,j], Y2[i,j], COV[i,j]), ll=ll)
        for i,j in product(range(X.shape[0]),range(X.shape[1]))
    )

    # reshape
    Z = np.array(Z).reshape(N,N)

    return X,Y,Z,{'x':leg[0],'y':leg[1],**fix}

def plt_ll_landscape_(R, vmid=None, vmin=1600, vmax=2600, levels=None):

    lbls = {
        'sx': 'std of random intercepts (unscaled)',
        'sy': 'std of random slopes (unscaled)',
        'corr': 'correlation btw random effects'
    }

    # extract
    X,Y,Z,leg = R

    leg_p = copy.copy(leg)

    # vmid and levels
    if vmid is None:
        vmid = .5*(vmin+vmax)
    if levels is None:
        levels = [vmin, vmid, vmax]

    plt.subplots(figsize=(14,6))
    c = plt.pcolormesh(X, Y, Z,
                   cmap='RdBu_r',
                   shading='gouraud', #snap=True,
                   antialiased=True,
                   zorder=-1., # grid over
                   norm=TwoSlopeNorm(vmid, vmin, vmax),
                  )
    plt.xlabel(lbls[leg_p.pop('x')])
    plt.ylabel(lbls[leg_p.pop('y')])

    contours = plt.contour(X, Y, Z, levels,
                           colors='black', alpha=.8, linewidths=1)
    cs = plt.clabel(contours, inline=True,
               #inline_spacing=300,
               rightside_up=True, fmt='%.1f',
               use_clabeltext=True, fontsize=10)

    ## move level labels up a bit
    #for csi in cs:
    #    csi.set_x(csi._x*1.1)
    #    csi.set_y(csi._y*1.1)

    plt.grid(which='major')
    plt.gca().tick_params(which='minor', top=True, right=True, bottom=True, left=True)

    if leg['y'] == 'corr':
        plt.ylim(-1,1)

    plt.colorbar(c, extend='both')

    fix_lbl = ", ".join([f"{k} = {v:.2f}" for k,v in leg_p.items()])
    plt.title('Log-likelihood of model wrt random effects var-covar matrix\n'
              f'Fix: {fix_lbl}')
    plt.show()


# wrapper
def plt_ll_landscape(model, **kws):

    plotting_kws = {
        kw: kws.pop(kw)
        for kw in inspect.signature(plt_ll_landscape_).parameters.keys()
        if kw != 'R' and kw in kws # result
    }

    R = ll_landscape(model, **kws)

    # plot
    plt_ll_landscape_(R, **plotting_kws)

    # return res
    return R
