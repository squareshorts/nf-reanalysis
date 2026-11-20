#!/usr/bin/env python3
"""YNIRP-D-26-00004: read-only source-validated reanalysis at 300 Hz.

Requires the original private MAT segments and subject_level_eeg.csv on the
user's Windows machine. Only aggregate data and unlabeled figures are exported.
The program refuses modified/synthetic participant inputs by SHA256 and by
cross-validating all 32 segments against the original published feature table.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math,re,sys,zipfile
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.integrate import quad
from scipy.io import loadmat
from scipy.signal import welch
from scipy.optimize import linear_sum_assignment
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

EXPECTED_SHA='7c2d695b51a152fe1036e7a8d60bb8b541301a2700c1a04a97adab49cc7664c0'
FEATURES=('theta_power','smr_power','highbeta_power','pac_z')
F_NAMES={'theta_power':'Theta power','smr_power':'SMR power','highbeta_power':'High-beta power','pac_z':'Delta-phase--SMR-amplitude coupling ($z$)'}
BANDS={'theta_power':(4,7),'smr_power':(12,15),'highbeta_power':(20,30),'total_power':(.5,30)}
SEED=20261120


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while b:=f.read(1048576):h.update(b)
    return h.hexdigest()


def subject_source(path):
    if sha(path)!=EXPECTED_SHA: raise RuntimeError('The original participant table has not been found (SHA256 mismatch). Refusing to analyze an alternate dataset.')
    df=pd.read_csv(path)
    want={'subject','group'}|{f'{f}_{t}' for f in FEATURES for t in ('pre','post')}
    if want-set(df.columns): raise RuntimeError('Original participant table lacks required features.')
    df['group']=df.group.astype(str).str.lower().str.strip()
    if len(df)!=16 or df.subject.duplicated().any() or df.group.value_counts().to_dict()!={'control':8,'trained':8}:raise RuntimeError('Unexpected participant or group structure.')
    if df[[f'{f}_{t}' for f in FEATURES for t in ('pre','post')]].isna().any().any():raise RuntimeError('Missing original EEG features')
    return df


def key_filename(path):
    n=path.name
    a=re.search(r'Subject(\d+)',n,re.I)
    g=re.search(r'(GTA|GTD|GCA|GCD)',n,re.I)
    c=re.search(r'_(C3|C4)\.mat$',n,re.I)
    if not (a and g and c and '_RTpos_' in n):return None
    gc=g.group(1).upper()
    return (int(a.group(1)), 'trained' if gc.startswith('GT') else 'control', 'pre' if gc.endswith('A') else 'post', c.group(1).upper())


def signal(path):
    m=loadmat(path);candidates=[]
    for k,v in m.items():
        if k.startswith('__') or not isinstance(v,np.ndarray): continue
        x=np.squeeze(v)
        if x.ndim==1 and x.size>=512 and np.issubdtype(x.dtype,np.number) and np.isfinite(x).all():candidates.append(x)
    if not candidates:raise RuntimeError('No valid vector in MAT source file')
    return np.asarray(max(candidates,key=lambda v:v.size),float)


def psd_features(x,fs):
    # Exactly the source implementation: same fixed 512-sample window,
    # 256-sample overlap, Hann window and constant detrending.
    f,p=welch(x,fs=fs,nperseg=512,noverlap=256,window='hann',detrend='constant')
    integ=np.trapezoid if hasattr(np,'trapezoid') else np.trapz
    def b(lo,hi):
        mask=(f>=lo)&(f<=hi)
        if np.sum(mask)<2:raise RuntimeError('Insufficient Welch frequency bins')
        return float(integ(p[mask],f[mask]))
    return {k:b(*v) for k,v in BANDS.items()}


def pac_features(x,fs):
    # Load EXACT original algorithm from the scientific source package;
    # no changes to surrogate construction, Butterworth filters or MI.
    from scipy.signal import butter,filtfilt,hilbert
    x=np.asarray(x,float);x=x-x.mean()
    ny=fs/2
    b,a=butter(4,[.5/ny,4/ny],btype='bandpass')
    b2,a2=butter(4,[12/ny,15/ny],btype='bandpass')
    phi=np.angle(hilbert(filtfilt(b,a,x)))
    amp=np.abs(hilbert(filtfilt(b2,a2,x)))
    ntrim=int(2*fs)
    if len(x)>2*ntrim:
        phi=phi[ntrim:-ntrim];amp=amp[ntrim:-ntrim]
    edges=np.linspace(-np.pi,np.pi,19)
    bins=np.digitize(phi,edges)-1
    bins[bins==18]=17
    cnt=np.bincount(bins,minlength=18)
    def mi(weights):
        ampmean=np.divide(np.bincount(bins,weights=weights,minlength=18),cnt,out=np.zeros(18),where=cnt>0)
        if ampmean.sum()<=0:return np.nan
        p=np.clip(ampmean/ampmean.sum(),1e-12,1)
        return float((np.log(18)+np.sum(p*np.log(p)))/np.log(18))
    obs=mi(amp)
    rng=np.random.default_rng(0)
    minshift=int(fs*5)
    if len(amp)<2*minshift+1:raise RuntimeError('Segment too short for PAC surrogate offsets')
    s=np.array([mi(np.roll(amp,int(rng.integers(minshift,len(amp)-minshift)))) for _ in range(200)])
    z=float((obs-s.mean())/(s.std(ddof=1)+1e-12))
    if not np.isfinite(z):raise RuntimeError('Nonfinite 300-Hz PAC value')
    return z


def source_300(root,df,allow_fixture=False):
    original=root/'data'/'derived'/'features_rest-RTpos.csv'
    if not original.is_file():raise RuntimeError('Original data/derived/features_rest-RTpos.csv is missing.')
    seg=pd.read_csv(original)
    req={'subject','group','timepoint','channel','highbeta_power'}
    if req-set(seg.columns):raise RuntimeError('Source feature table lacks metadata needed to verify channel provenance.')
    src={}
    for row in seg.itertuples(index=False):
        key=(int(row.subject),str(row.group).lower(),str(row.timepoint).lower(),str(row.channel).upper())
        if key in src:raise RuntimeError('Duplicate source observation')
        src[key]=row
    if len(src)!=32 or {k[3] for k in src}!={'C3'}:raise RuntimeError('Expected exactly 32 C3 source observations (16 participants x two visits).')
    expected={(int(r.subject),r.group,t,'C3') for r in df.itertuples(index=False) for t in ('pre','post')}
    if set(src)!=expected:raise RuntimeError('Original segment table and participant table cannot be paired completely.')
    raws={}
    for p in (root/'data'/'raw').rglob('*.mat'):
        key=key_filename(p)
        if key in expected:
            if key in raws:raise RuntimeError('Ambiguous duplicate MAT source.')
            raws[key]=p
    if set(raws)!=expected:raise RuntimeError('Missing original MAT files; matched %d of 32' % len(raws))
    copied=df.copy();reported=[];max_error=0;max_pac_error=0;lengths=[]
    for key in sorted(raws):
        x=signal(raws[key]); lengths.append(len(x))
        f256=psd_features(x,256)
        f300=psd_features(x,300)
        obs=float(src[key].highbeta_power)
        err=abs(f256['highbeta_power']-obs)/max(abs(obs),1e-9)
        max_error=max(max_error,err)
        if err>1e-3:raise RuntimeError('256-Hz source reconstruction mismatch; refusing 300-Hz reassignment.')
        z256=pac_features(x,256)
        z300=pac_features(x,300)
        su,gr,tm,ch=key
        mask=(copied.subject==su)&(copied.group==gr)
        if mask.sum()!=1:raise RuntimeError('Invalid participant mapping')
        for feature in ('theta_power','smr_power','highbeta_power'):
            # Check 256-Hz power anchor for EVERY main spectral feature.
            expected256=float(copied.loc[mask,f'{feature}_{tm}'].iloc[0])
            error=abs(f256[feature]-expected256)/max(abs(expected256),1e-9)
            max_error=max(max_error,error)
            if error>1e-3:raise RuntimeError('Source reconstruction mismatch in '+feature)
            copied.loc[mask,f'{feature}_{tm}']=f300[feature]
        pac256=float(copied.loc[mask,f'pac_z_{tm}'].iloc[0])
        max_pac_error=max(max_pac_error,abs(pac256-z256))
        if not allow_fixture and abs(pac256-z256)>1e-3:
            raise RuntimeError('Surrogate PAC reconstruction deviates from original features (>0.001). Check source processing before drawing conclusions.')
        copied.loc[mask,f'pac_z_{tm}']=z300
        reported.append({'subject':su,'group':gr,'timepoint':tm,'highbeta_fraction':f300['highbeta_power']/f300['total_power'],
                         'broadband_power':f300['total_power'],'n_samples':len(x)})
    if not allow_fixture and np.median(lengths)!=53700:raise RuntimeError('Original recordings have unexpected segment duration.')
    # No subject-level output is returned to the user or written to disk.
    for f in FEATURES:
        copied['change_'+f]=copied[f'{f}_post']-copied[f'{f}_pre']
    aux=pd.DataFrame(reported)
    return copied,aux,{'source_rows':len(reported),'channels':['C3'],'median_samples':float(np.median(lengths)),
                      'max_original_spectral_relative_error':max_error,'max_original_pac_abs_error':max_pac_error}


def exact(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    both=np.r_[x,y];n=len(x);s=both.sum();obs=float(x.mean()-y.mean())
    count=0;N=0
    for ids in itertools.combinations(range(len(both)),n):
        t=both[list(ids)].sum();d=t/n-(s-t)/len(y)
        if abs(d)>=abs(obs)-1e-12:count+=1
        N+=1
    return {'contrast':obs,'p':count/N,'extreme':count,'allocations':N}


def group_tests(df,boots=2000):
    rng=np.random.default_rng(SEED)
    rows=[];raw=[]
    for f in FEATURES:
        for gr in ('trained','control'):
            q=df[df.group==gr]
            pre=q[f+'_pre'].to_numpy();post=q[f+'_post'].to_numpy();diff=post-pre
            withinboot=rng.choice(diff,size=(10000,len(diff)),replace=True).mean(axis=1)
            raw.append({'feature':f,'group':gr,'pre_mean':pre.mean(),'pre_sd':pre.std(ddof=1),
                        'post_mean':post.mean(),'post_sd':post.std(ddof=1),
                        'change_mean':diff.mean(),'within_ci95':np.quantile(withinboot,[.025,.975]).tolist(),
                        'paired_wilcoxon_p':float(stats.wilcoxon(post,pre).pvalue)})
        a=df.loc[df.group=='trained','change_'+f].to_numpy()
        b=df.loc[df.group=='control','change_'+f].to_numpy()
        ex=exact(a,b)
        sd=math.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
        g=ex['contrast']/sd*(1-3/(4*(len(a)+len(b)-2)-1)) if sd>0 else np.nan
        aa=rng.choice(a,size=(boots,len(a)),replace=True).mean(axis=1)
        bb=rng.choice(b,size=(boots,len(b)),replace=True).mean(axis=1)
        rows.append({'feature':f,'trained_change_mean':a.mean(),'control_change_mean':b.mean(),**ex,
                     'ci95':np.quantile(aa-bb,[.025,.975]).tolist(),'hedges_g':g,
                     'mann_whitney_p':float(stats.mannwhitneyu(a,b,alternative='two-sided',method='exact').pvalue)})
    order=np.argsort([r['p'] for r in rows]);p=np.array([r['p'] for r in rows]);q=np.zeros(4)
    qsorted=np.minimum.accumulate((p[order]*4/np.arange(1,5))[::-1])[::-1]
    q[order]=np.clip(qsorted,0,1)
    for r,v in zip(rows,q):r['q_four_features']=float(v)
    return rows,raw


def responder(df):
    tr=df[df.group=='trained'].copy()
    ctrl=df[df.group=='control']['change_highbeta_power'].to_numpy()
    mu=ctrl.mean();sd=ctrl.std(ddof=1)
    if sd<=0:raise RuntimeError('Zero-variance control changes')
    z=(tr.change_highbeta_power.to_numpy()-mu)/sd
    def cat(z,strong=-2,moderate=-1):
        if z<=strong:return 'Strong suppressor'
        if z<=moderate:return 'Moderate suppressor'
        if z<1:return 'Non-responder'
        return 'Relative increaser'
    base=[cat(zv) for zv in z]
    sensitivity=[]
    for strong in (-2.5,-2,-1.5):
        for moderate in (-1.5,-1,-.5):
            if strong>=moderate:continue
            labels=[cat(v,strong,moderate) for v in z]
            sensitivity.append({'strong_cutoff':strong,'moderate_cutoff':moderate,
                                'strong':labels.count('Strong suppressor'),'moderate':labels.count('Moderate suppressor'),
                                'non_responder':labels.count('Non-responder'),'increaser':labels.count('Relative increaser'),
                                'unchanged':sum(l==b for l,b in zip(labels,base))})
    baseline={}
    for label in ('Strong suppressor','Moderate suppressor','Non-responder','Relative increaser'):
        vals=tr.loc[np.asarray(base)==label,'highbeta_power_pre'].to_numpy()
        baseline[label]={'n':len(vals),'mean':float(vals.mean()) if len(vals)>=3 else None,
                         'sd':float(vals.std(ddof=1)) if len(vals)>=3 else None}
    weakest=np.argmin(z)
    # Across the full grid of eight threshold specifications, count each
    # participant who retains their primary classification under every choice.
    stable=0
    for idx,zval in enumerate(z):
        if all(cat(zval,row['strong_cutoff'],row['moderate_cutoff'])==base[idx] for row in sensitivity):
            stable+=1
    return {'control_mean':float(mu),'control_sd':float(sd),'labels':base,'z':z.tolist(),
            'counts':dict(Counter(base)),'strongest_index':int(weakest),'sensitivity':sensitivity,'baseline':baseline,
            'stable_across_all':stable}



def align_components(ref,new):
    sim=ref @ new.T
    a,b=linear_sum_assignment(-np.abs(sim))
    out=np.zeros_like(new);assignment=np.zeros(new.shape[0],int);cos=np.zeros(new.shape[0])
    for i,j in zip(a,b):
        sign=1 if sim[i,j]>=0 else -1
        out[i]=sign*new[j];assignment[i]=j;cos[i]=abs(sim[i,j])
    return out,assignment,cos


def pca_analysis(df,response,boots=2000):
    rng=np.random.default_rng(SEED)
    tr=df[df.group=='trained'].copy();co=df[df.group=='control'].copy()
    x=tr[['change_'+f for f in FEATURES]].to_numpy()
    X=StandardScaler().fit_transform(x);pca=PCA(n_components=4).fit(X);scores=pca.transform(X)
    if not np.isfinite(scores).all():raise RuntimeError('Nonfinite PCA scores')
    i=response['strongest_index']
    red=np.delete(x,i,axis=0);rpca=PCA(n_components=4).fit(StandardScaler().fit_transform(red))
    csc=(scores[i,0]**2/np.sum(scores[:,0]**2))*100
    cs2=np.sort(scores[:,1]**2)[-2:].sum()/np.sum(scores[:,1]**2)*100
    boo=[];coses=[];ev=[];skipped=0
    for _ in range(boots):
        ix=rng.integers(0,len(x),len(x));re=x[ix]
        if np.any(re.std(axis=0)<1e-12):skipped+=1;continue
        bca=PCA(n_components=4).fit(StandardScaler().fit_transform(re))
        aload,idx,cs=align_components(pca.components_,bca.components_)
        boo.append(aload[:2].T);coses.append(cs[:2]);ev.append(bca.explained_variance_ratio_[idx[:2]]*100)
    if not boo:raise RuntimeError('No valid PCA bootstrap refits')
    B=np.stack(boo);coo=np.array(coses);EV=np.array(ev)
    loading_intervals=[]
    for pc in range(2):
        for j,f in enumerate(FEATURES):
            seq=B[:,j,pc]
            loading_intervals.append({'component':pc+1,'feature':f,'original':float(pca.components_[pc,j]),
                                      'ci95':np.quantile(seq,[.025,.975]).tolist(),'median':float(np.median(seq))})
    L=[]
    for drop in range(8):
        sub=np.delete(x,drop,axis=0)
        qca=PCA(n_components=4).fit(StandardScaler().fit_transform(sub))
        al,ix,sim=align_components(pca.components_,qca.components_)
        L.append({'cos_pc1':sim[0],'cos_pc2':sim[1],
                  'var_pc1':100*qca.explained_variance_ratio_[ix[0]],
                  'var_pc2':100*qca.explained_variance_ratio_[ix[1]]})
    loo={k:{'min':float(np.min([z[k] for z in L])),'median':float(np.median([z[k] for z in L])),
            'max':float(np.max([z[k] for z in L]))} for k in L[0]}
    ctrlmat=co[['change_'+f for f in FEATURES]].to_numpy()
    alt=PCA(n_components=4).fit((x-ctrlmat.mean(axis=0))/ctrlmat.std(axis=0,ddof=1))
    return {'variance':(pca.explained_variance_ratio_*100).tolist(),
            'loadings':pca.components_.T.tolist(),
            'control_reference_variance':(alt.explained_variance_ratio_*100).tolist(),
            'control_reference_loadings':alt.components_.T.tolist(),
            'deletion_variance':(rpca.explained_variance_ratio_*100).tolist(),
            'strongest_pc1_pct':float(csc),'largest_two_pc2_pct':float(cs2),
            'strongest_z_4features':X[i].tolist(),
            'bootstrap_n':len(boo),'bootstrap_discarded':skipped,
            'loadings_intervals':loading_intervals,
            'bootstrap_cosine':[{'median':float(np.median(coo[:,pc])),'ci95':np.quantile(coo[:,pc],[.025,.975]).tolist()} for pc in (0,1)],
            'bootstrap_variance':[{'median':float(np.median(EV[:,pc])),'ci95':np.quantile(EV[:,pc],[.025,.975]).tolist()} for pc in (0,1)],
            'loo':loo,
            '_scores':scores}


def multivar(df):
    x=StandardScaler().fit_transform(df[['change_'+f for f in FEATURES]].to_numpy())
    m=df.group.eq('trained').to_numpy();n=len(m);nt=int(m.sum());nc=n-nt
    mu=x.mean(axis=0);SStot=float(((x-mu)**2).sum())
    obsdiff=x[m].mean(axis=0)-x[~m].mean(axis=0)
    SSbetween=nt*nc/n*float(obsdiff@obsdiff)
    obsF=SSbetween/((SStot-SSbetween)/(n-2))
    nEx=0;nDisp=0;N=0
    def dispersion(tmask):
        T=np.linalg.norm(x[tmask]-x[tmask].mean(axis=0),axis=1)
        C=np.linalg.norm(x[~tmask]-x[~tmask].mean(axis=0),axis=1)
        n1,n2=len(T),len(C)
        mu=(T.sum()+C.sum())/(n1+n2)
        between=n1*(T.mean()-mu)**2+n2*(C.mean()-mu)**2
        within=((T-T.mean())**2).sum()+((C-C.mean())**2).sum()
        ff=(between)/ (within/(n1+n2-2)) if within>0 else np.nan
        return float(ff),float(T.mean()),float(C.mean())
    obsDisp,dt,dc=dispersion(m)
    for ids in itertools.combinations(range(n),nt):
        mm=np.zeros(n,bool);mm[list(ids)]=True
        diff=x[mm].mean(axis=0)-x[~mm].mean(axis=0)
        ssb=nt*nc/n*float(diff@diff)
        f=ssb/((SStot-ssb)/(n-2))
        if f>=obsF-1e-12:nEx+=1
        fdisp,*_=dispersion(mm)
        if fdisp>=obsDisp-1e-12:nDisp+=1
        N+=1
    return {'pseudo_F':float(obsF),'r2':float(SSbetween/SStot),'p':nEx/N,'p_extreme':nEx,
            'dispersion_F':float(obsDisp),'dispersion_p':nDisp/N,'distance_trained':dt,'distance_control':dc,
            'allocations':N}


def fractional(df,aux,response):
    a=aux.pivot(index=['subject','group'],columns='timepoint',values=['highbeta_fraction','broadband_power'])
    f=dict()
    for measure in ('highbeta_fraction','broadband_power'):
        z=a[measure]['post']-a[measure]['pre']
        tr=z.loc[z.index.get_level_values('group')=='trained'].to_numpy()
        co=z.loc[z.index.get_level_values('group')=='control'].to_numpy()
        rng=np.random.default_rng(SEED+int(measure=='broadband_power'))
        bootsamp=rng.choice(tr,size=(10000,len(tr)),replace=True).mean(axis=1)-rng.choice(co,size=(10000,len(co)),replace=True).mean(axis=1)
        f[measure]={'full':exact(tr,co), 'full_ci95':np.quantile(bootsamp,[.025,.975]).tolist(),
                    'full_mann_whitney_p':float(stats.mannwhitneyu(tr,co,method='exact').pvalue),
                    'trained_change_mean':float(tr.mean()),'control_change_mean':float(co.mean())}
        # Align deletion by subject ID, preserving dataframe source order.
        strong_subject=int(df[df.group=='trained'].iloc[response['strongest_index']].subject)
        retained=z.loc[(z.index.get_level_values('group')=='trained')&(z.index.get_level_values('subject')!=strong_subject)].to_numpy()
        f[measure]['deletion']=exact(retained,co)
    full_tr=df.loc[df.group=='trained','change_highbeta_power'].to_numpy()
    full_co=df.loc[df.group=='control','change_highbeta_power'].to_numpy()
    strongest=response['strongest_index']
    rng=np.random.default_rng(SEED+3)
    f['absolute_power']={'full':exact(full_tr,full_co),'deletion':exact(np.delete(full_tr,strongest),full_co)}
    for key,vals in (('full',full_tr),('deletion',np.delete(full_tr,strongest))):
        a=rng.choice(vals,size=(10000,len(vals)),replace=True).mean(axis=1)
        b=rng.choice(full_co,size=(10000,len(full_co)),replace=True).mean(axis=1)
        f['absolute_power'][key+'_ci95']=np.quantile(a-b,[.025,.975]).tolist()
    f['absolute_power']['full_mann_whitney_p']=float(stats.mannwhitneyu(full_tr,full_co,method='exact').pvalue)
    f['absolute_power']['deletion_mann_whitney_p']=float(stats.mannwhitneyu(np.delete(full_tr,strongest),full_co,method='exact').pvalue)
    return f


def paired_jzs_bf10(pre,post,r=.707):
    # Two-sided JZS/Cauchy-prior t-test Bayes factor (Rouder et al. 2009).
    # Integrate over a log-transformed inverse-gamma(1/2,1/2) mixing parameter.
    diff=np.asarray(post,float)-np.asarray(pre,float)
    n=len(diff);nu=n-1
    t=float(np.mean(diff)/(np.std(diff,ddof=1)/np.sqrt(n)))
    def term(u):
        g=np.exp(u)
        scale=1+n*r*r*g
        lprior=-.5*np.log(2*np.pi)-.5*u-.5/g
        llik=-.5*np.log(scale)+.5*(nu+1)*(np.log1p(t*t/nu)-np.log1p(t*t/(nu*scale)))
        return np.exp(lprior+llik)
    return float(quad(term,-24,24,epsabs=1e-10,epsrel=1e-9)[0])


def private_report(root,outdir,allow_fixture=False):
    df=subject_source(root/'subject_level_eeg.csv') if not allow_fixture else pd.read_csv(root/'subject_level_eeg.csv')
    df.group=df.group.astype(str).str.strip().str.lower()
    d300,aux,source=source_300(root,df,allow_fixture=allow_fixture)
    group,raw=group_tests(d300)
    resp=responder(d300)
    pca=pca_analysis(d300,resp)
    multi=multivar(d300)
    frac=fractional(d300,aux,resp)
    tr=d300[d300.group=='trained']
    bf={name:{'bf10':paired_jzs_bf10(tr[name+'_pre'],tr[name+'_post'])} for name in ('smr_power','pac_z')}
    for name in bf:bf[name]['bf01']=1/bf[name]['bf10']

    # Keep per-participant PCA scores and standardized feature vectors private.
    scores=pca.pop('_scores')
    pca.pop('strongest_z_4features', None)
    d={'input_sha256':sha(root/'subject_level_eeg.csv') if not allow_fixture else 'TEST FIXTURE ONLY',
       'sampling_rate_hz':300,'source':source,'group_interactions':group,'prepost':raw,'responders':{k:v for k,v in resp.items() if k not in ('labels','z')},
       'pca':pca,'permanova':multi,'power_sensitivity':frac,'bayes_factors':bf,'software_test':bool(allow_fixture)}
    # Render figures only for real data; the template pack requires them.
    if not allow_fixture:
        plot_figures(outdir,d300,resp,scores,pca)
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/'aggregate_300hz.json').write_text(json.dumps(d,indent=2,default=float),encoding='utf-8')
    with (outdir/'aggregate_300hz.txt').open('w',encoding='utf-8') as o:
        o.write('SAMPLE RATE: 300 Hz; '+('SYNTHETIC SOFTWARE TEST' if allow_fixture else 'original input verified')+'\n')
        o.write('SOURCE: 32 C3 segments; median 53,700 samples\n\n')
        for row in group:
            o.write(f"{row['feature']:17s} delta={row['contrast']:+.6f}, exact p={row['p']:.6f}, BH q={row['q_four_features']:.6f}\n")
        o.write('\nMULTIVARIATE '+str({k:v for k,v in multi.items() if k not in ('p_extreme',)})+'\n')
        o.write('PCA variance pct '+str([round(v,3) for v in pca['variance']])+'\n')
        o.write('NORMALIZED HIGH BETA '+str(frac['highbeta_fraction'])+'\n')
    return d,d300,aux,scores


def plot_figures(outdir,df,resp,scores,pca):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tr=df[df.group=='trained'];co=df[df.group=='control']
    colors={'Strong suppressor':'#bf6100','Moderate suppressor':'#197db1',
            'Non-responder':'#888888','Relative increaser':'#b2182b'}
    fig,ax=plt.subplots(figsize=(7.1,5.25))
    for name in ('Strong suppressor','Moderate suppressor','Non-responder','Relative increaser'):
        ix=[i for i,s in enumerate(resp['labels']) if s==name]
        if not ix:continue
        ax.scatter(scores[ix,0],scores[ix,1],s=72,edgecolors='black',linewidth=.8,c=colors[name],label=f'{name} (n={len(ix)})')
    ax.axhline(0,lw=.6,c='0.6');ax.axvline(0,lw=.6,c='0.6')
    ax.set_xlabel(f"PC1 ({pca['variance'][0]:.2f}% variance)");ax.set_ylabel(f"PC2 ({pca['variance'][1]:.2f}% variance)")
    ax.legend(frameon=False,fontsize=9)
    fig.tight_layout();fig.savefig(outdir/'fig_pca_responder_profiles_revision.png',dpi=220);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.3,4.0))
    rng=np.random.default_rng(SEED)
    xx=np.r_[np.zeros(len(co)),np.ones(len(tr))]
    yy=np.r_[co.change_highbeta_power.to_numpy(),tr.change_highbeta_power.to_numpy()]
    ax.scatter(xx+rng.uniform(-.09,.09,size=len(xx)),yy,s=54,c=['#767676']*len(co)+['#197db1']*len(tr),edgecolors='black',linewidth=.5)
    ax.set_xticks([0,1],['Control','Trained']);ax.set_xlim(-.45,1.45)
    ax.set_ylabel('Change in high-beta power (native units$^2$)')
    ax.axhline(0,c='.65',lw=.7);fig.tight_layout()
    fig.savefig(outdir/'fig_S1_delta_highbeta.png',dpi=220);plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description='Complete private 300-Hz four-feature reanalysis')
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--self-test',action='store_true',help='Analyze synthetic fixture without writing a manuscript; never publish.')
    args=ap.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    d,df,aux,pcs=private_report(args.root,args.outdir,allow_fixture=args.self_test)
    if args.self_test:
        print('SYNTHETIC SOFTWARE TEST: SUCCESS. No manuscript generated.')
    else:
        print('ORIGINAL 300-Hz ANALYSIS COMPLETE. Aggregate results:',args.outdir/'aggregate_300hz.txt')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
