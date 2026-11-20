"""Smoke tests using synthetic signals. No real subjects included."""
import importlib.util
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('nf_300_recompute', R / 'scripts' / 'recompute_full_300hz.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def test_welch_high_beta_concentrates_power():
    fs=300
    t=np.arange(fs*20)/fs
    x=np.sin(2*np.pi*23*t) + 0.1*np.sin(2*np.pi*10*t)
    p=m.psd_features(x,fs)
    assert p['highbeta_power']>p['smr_power']
    assert 0.0 < p['highbeta_power'] / p['total_power'] < 1.0

def test_exact_assignment_counts():
    a=np.array([0.,1.,2.,3.,4.,5.,6.,7.])
    b=np.array([-1.,-2.,-3.,-4.,-5.,-6.,-7.,-8.])
    res=m.exact(a,b)
    assert res['allocations']==12870
    assert np.isfinite(res['contrast']) and np.isfinite(res['p']) and 0<=res['p']<=1
    assert 0<=res['extreme']<=res['allocations']

def test_pac_on_synthetic_signal():
    fs=300
    t=np.arange(fs*22)/fs
    rng=np.random.default_rng(2026)
    x=np.sin(2*np.pi*2*t)+0.3*np.sin(2*np.pi*13*t)+0.01*rng.normal(size=len(t))
    z=m.pac_features(x,fs)
    assert np.isfinite(z)
