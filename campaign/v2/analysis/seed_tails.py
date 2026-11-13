#!/usr/bin/env python3
"""Analyze independent seed-survival batches without editing raw outputs.

The predeclared primary analysis requires all sixteen confirmatory 500k shards.
The original eight 2M shards remain a separate supplementary fixed ensemble.
Runtime-truncated shards and partial collections of completed shards are excluded: elapsed-time stopping can select against expensive long survivors.
Legacy shards are reported separately with nominal intervals.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
OLD_GLOB="data/runs/E8a_long/out/p0.99_a0.5_j2_q0.0012954_L8192_T1e6_k*.txt"
NEW_GLOB="data/v2/cpu/tails/tail_p0.99_q0.0012954_L8192_shard*.txt"
CONFIRM_GLOB="data/v2/cpu/confirm8m/tails/confirm8m_p0.99_q0.0012954_L8192_shard*.txt"
FIXED_POPULATIONS={"original_fixed_16M", "confirmatory_8M"}
Z95=1.959963984540054


def binomial_cdf(k,n,p):
    if k<0:return 0.
    if k>=n or p<=0:return 1.
    if p>=1:return 0.
    if n*p<500:
        term=math.exp(n*math.log1p(-p));total=term
        ratio=p/(1-p)
        for j in range(k):
            term*=((n-j)/(j+1))*ratio;total+=term
    else:
        logpmf=(math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
                +k*math.log(p)+(n-k)*math.log1p(-p))
        term=math.exp(logpmf) if logpmf>-745 else 0.;total=term
        ratio=(1-p)/p
        for j in range(k,0,-1):
            term*=j/(n-j+1)*ratio;total+=term
    return min(1.,max(0.,total))


def solve_cdf(k,n,target):
    lo,hi=0.,1.
    for _ in range(70):
        mid=(lo+hi)/2
        if binomial_cdf(k,n,mid)>target:lo=mid
        else:hi=mid
    return (lo+hi)/2


def binomial_interval(k,n,alpha=.05):
    """Exact equal-tailed intervals in sparse tails; Wilson for dense counts."""
    if not 0<=k<=n or n<=0:raise ValueError("Invalid binomial counts")
    if k>n/2:
        lo,hi,method=binomial_interval(n-k,n,alpha)
        return 1-hi,1-lo,method
    if k==0:
        return 0.,-math.expm1(math.log(alpha/2)/n),"Clopper-Pearson"
    if k<=1000:
        return (solve_cdf(k-1,n,1-alpha/2),solve_cdf(k,n,alpha/2),
                "Clopper-Pearson")
    p=k/n;z=Z95;denom=1+z*z/n
    centre=(p+z*z/(2*n))/denom
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/denom
    return max(0.,centre-half),min(1.,centre+half),"Wilson"


def read_batch(path):
    text=path.read_text()
    match=lambda pattern:re.search(pattern,text,re.M)
    seed_line=match(r"^# Lx \d+ Ly \d+ T \d+ seed (\d+)")
    configuration=match(r"^# reps_max (\d+) target (\d+) groups (\d+)")
    if not seed_line or not configuration:raise ValueError(f"Missing replica header: {path}")
    seed=int(seed_line.group(1))
    footer=match(r"^# done replicas (\d+) survivors_T (\d+) escaped (\d+)")
    completed=bool(footer)
    rows=np.loadtxt(path,comments="#")
    if rows.ndim!=2 or rows.shape[1]<5:raise ValueError(f"Bad replica rows: {path}")
    groups=int(configuration.group(3))
    times=np.unique(rows[:,1]).astype(np.int64)
    values=np.zeros((len(times),3),dtype=np.int64)
    for i,t in enumerate(times):
        block=rows[rows[:,1]==t]
        if len(block)!=groups or set(block[:,0])!=set(range(groups)):
            raise ValueError(f"Missing or duplicate replica groups: {path}")
        selected=block[:,2:5]
        if not np.array_equal(selected,np.rint(selected)):
            raise ValueError(f"Noninteger replica counts: {path}")
        values[i]=np.rint(selected.sum(axis=0)).astype(np.int64)
    if len(set(values[:,0]))!=1:raise ValueError(f"Changing denominator: {path}")
    if np.any(values<0) or np.any(values[:,2]>values[:,1]) or np.any(values[:,1]>values[:,0]):
        raise ValueError(f"Impossible survival/escape counts: {path}")
    if np.any(np.diff(values[:,1])>0):
        raise ValueError(f"Nonmonotone survival: {path}")
    if np.any(np.diff(values[:,2])<0):
        raise ValueError(f"Nonmonotone cumulative escapes: {path}")
    if completed and tuple(map(int,footer.groups()))!=tuple(map(int,values[-1])):
        raise ValueError(f"Completion footer disagrees with integer counts: {path}")
    declared=int(configuration.group(1));target=int(configuration.group(2))
    return dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                seed=seed,low_key=seed&0xffffffff,declared_replicas=declared,
                adaptive_target=target,completed=completed,actual_replicas=int(values[0,0]),
                groups=groups,
                reached_fixed_max=int(values[0,0])==declared,times=times,values=values)


def pool_batches(batches):
    times=batches[0]["times"].copy()
    for batch in batches[1:]:times=np.intersect1d(times,batch["times"])
    if not len(times):raise ValueError("No shared output times")
    values=np.zeros((len(times),3),dtype=np.int64)
    for batch in batches:
        index=np.searchsorted(batch["times"],times)
        values+=batch["values"][index]
    return times,values


def survival_rows(name,batches):
    times,values=pool_batches(batches)
    rows=[]
    for t,(n,surv,esc) in zip(times,values):
        n,surv,esc=map(int,(n,surv,esc))
        lo,hi,method=binomial_interval(surv,n)
        if esc:
            lo=binomial_interval(surv-esc,n)[0]
            method+="; escape bounds"
        rows.append(dict(population=name,time=int(t),n=n,survivors=surv,escaped=esc,
                         alive_unescaped=surv-esc,survival=surv/n,
                         ci_low=lo,ci_high=hi,ci_method=(method if name in FIXED_POPULATIONS else "nominal "+method),
                         survival_lower_estimate=(surv-esc)/n,
                         survival_upper_estimate=surv/n))
    return rows


def exponent_rows(name,rows,min_survivors):
    result=[]
    times=np.array([r["time"] for r in rows])
    positive=np.where(times>0)[0]
    for j in positive:
        if times[j]<100:continue
        i=int(positive[np.argmin(abs(np.log(times[positive]/(times[j]/4))))])
        ratio=times[j]/times[i]
        if not 3.5<ratio<4.5:continue
        a,b=rows[i],rows[j]
        supported=(b["survivors"]>=min_survivors and not a["escaped"] and not b["escaped"])
        row=dict(population=name,t1=a["time"],t2=b["time"],time_ratio=ratio,
                 survivors_t1=a["survivors"],survivors_t2=b["survivors"],
                 escaped_t2=b["escaped"],supported=supported,
                 minimum_late_survivors=min_survivors,delta_eff=float("nan"),
                 ci_low=float("nan"),ci_high=float("nan"))
        if supported:
            r=b["survivors"]/a["survivors"]
            lo,hi,_=binomial_interval(b["survivors"],a["survivors"])
            row.update(delta_eff=-math.log(r)/math.log(ratio),
                       ci_low=-math.log(hi)/math.log(ratio),
                       ci_high=-math.log(lo)/math.log(ratio) if lo else float("inf"))
        result.append(row)
    return result


def batch_diagnostics(name,batches,rng):
    last=min(int(b["times"][-1]) for b in batches)
    ns=[];ss=[];ee=[];table=[]
    for b in batches:
        i=int(np.searchsorted(b["times"],last))
        n,s,e=map(int,b["values"][i]);ns.append(n);ss.append(s);ee.append(e)
        lo,hi,method=binomial_interval(s,n)
        table.append(dict(population=name,path=b["path"],seed=b["seed"],time=last,
                          n=n,survivors=s,escaped=e,survival=s/n,
                          ci_low=lo,ci_high=hi,ci_method=(method if name in FIXED_POPULATIONS else "nominal "+method),
                          adaptive_target=b["adaptive_target"]))
    ns=np.array(ns);ss=np.array(ss);total=int(ss.sum());weights=ns/ns.sum()
    diag=dict(time=last,batches=len(batches),replicas=int(ns.sum()),survivors=total,
              escaped=int(sum(ee)),conditional_test="Multivariate hypergeometric allocation of fixed total survivors across fixed-size independent batches",
              monte_carlo_repetitions=50000)
    if total and len(ns)>1:
        expected=total*weights
        stat=float(np.sum((ss-expected)**2/expected))
        simulated=rng.multivariate_hypergeometric(ns,total,size=50000)
        null=np.sum((simulated-expected)**2/expected,axis=1)
        diag.update(pearson_statistic=stat,
                    conditional_monte_carlo_p=(1+int(np.sum(null>=stat-1e-12)))/(1+len(null)),
                    minimum_expected_survivors=float(expected.min()))
    else:
        diag.update(pearson_statistic=None,conditional_monte_carlo_p=None,
                    minimum_expected_survivors=0)
    return table,diag


def fisher_two_sided(k1,n1,k2,n2):
    """Conditional comparison; stable for the sparse endpoints studied here."""
    K=k1+k2
    if K>10000:return None
    lower=max(0,K-n2);upper=min(K,n1)
    if lower!=0:return None
    prob0=math.prod((n2-i)/(n1+n2-i) for i in range(K))
    probs=[prob0]
    for x in range(upper):
        probs.append(probs[-1]*(n1-x)/(x+1)*(K-x)/(n2-K+x+1))
    probs=np.array(probs);probs/=probs.sum()
    return float(probs[probs<=probs[k1]+1e-15].sum())


def write_csv(path,rows):
    if not rows:return
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def selftest():
    assert abs(binomial_cdf(1,2,.5)-.75)<1e-14
    lo,hi,method=binomial_interval(0,100)
    assert lo==0 and abs(hi-(1-.025**.01))<1e-14
    lo,hi,_=binomial_interval(3,3000000)
    assert 2.0e-7<lo<2.1e-7 and 2.9e-6<hi<3.0e-6
    for k,n in [(1,10),(3,3000000),(1000,16000000)]:
        lo,hi,_=binomial_interval(k,n)
        assert abs(binomial_cdf(k-1,n,lo)-.975)<2e-7
        assert abs(binomial_cdf(k,n,hi)-.025)<2e-7
    assert abs(fisher_two_sided(1,10,1,10)-1)<1e-12
    import tempfile
    with tempfile.TemporaryDirectory(dir=ROOT/"campaign/v2") as directory:
        p=Path(directory)/"batch.txt"
        source=("# Lx 64 Ly 64 T 1 seed 1 threads 1\n# reps_max 4 target 0 groups 2\n"
                "0 0 2 2 0\n0 1 2 1 0\n1 0 2 2 0\n1 1 2 0 0\n"
                "# done replicas 4 survivors_T 1 escaped 0\n")
        p.write_text(source);assert read_batch(p)["completed"]
        p.write_text(source.replace("survivors_T 1","survivors_T 2"))
        try:read_batch(p)
        except ValueError:pass
        else:raise AssertionError("inconsistent completed footer accepted")
        p.write_text(source.replace('1 1 2 0 0','0 1 2 0 0'))
        try:read_batch(p)
        except ValueError:pass
        else:raise AssertionError("duplicate replica group accepted")
    print("seed_tails selftest PASS")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selftest",action="store_true")
    parser.add_argument("--min-survivors",type=int,default=20)
    parser.add_argument("--seed",type=int,default=20261103)
    parser.add_argument("--output",type=Path,default=ROOT/"campaign"/"v2"/"results"/"seed_tails")
    args=parser.parse_args()
    if args.selftest:selftest();return
    from .io import read_header
    def validate_physics(path,batch):
        h=read_header(path)
        expected=dict(Lx=8192,Ly=8192,T=1_000_000,p=.99,q=.0012954,alpha=.5,
                      init='disk:1:0.75',parameter_source='header')
        for key,value in expected.items():
            actual=h.get(key)
            valid=(isinstance(actual,(int,float)) and abs(actual-value)<1e-14
                   if isinstance(value,float) else actual==value)
            if not valid:raise ValueError(f'Tail physical header mismatch {key}: {path}')
        expected_times={0,1_000_000}|{t for d in range(11) for i in range(20)
                                     if (t:=math.floor(10.**(d+i/20)+.5))<=1_000_000}
        if set(batch['times'])!=expected_times:
            raise ValueError(f'Incomplete declared tail observation schedule: {path}')
    all_old=[read_batch(p) for p in sorted(ROOT.glob(OLD_GLOB))]
    for batch in all_old:validate_physics(ROOT/batch['path'],batch)
    def load_campaign(pattern, manifest_path, replicas, shards):
        manifest=json.loads(manifest_path.read_text())
        expected={Path(t["output"]).name:t for t in manifest["tasks"]["tails"]}
        assert len(expected)==shards
        batches=[];invalid=[]
        for path in sorted(ROOT.glob(pattern)):
            if path.name not in expected:raise ValueError(f"Unplanned tail file: {path}")
            try:
                batch=read_batch(path)
            except (ValueError,AttributeError,IndexError) as error:
                if "# done replicas" in path.read_text():raise
                invalid.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,reason=str(error)))
                continue
            if (batch["seed"]!=expected[path.name]["seed"] or batch["declared_replicas"]!=replicas
                    or batch["adaptive_target"]!=0 or batch["groups"]!=16):
                raise ValueError(f"Tail metadata differs from fixed manifest: {path}")
            validate_physics(path,batch)
            batches.append(batch)
        completed=[b for b in batches if b["completed"] and b["reached_fixed_max"]]
        return batches,completed,invalid,len(completed)==shards
    all_new,new,invalid_new,campaign_complete=load_campaign(NEW_GLOB,
        ROOT/"campaign/v2/hpc/cpu/manifest.json",2_000_000,8)
    all_confirm,confirm,invalid_confirm,confirm_complete=load_campaign(CONFIRM_GLOB,
        ROOT/"campaign/v2/hpc/cpu/confirm8m/manifest.json",500_000,16)
    old=[b for b in all_old if b["completed"] and b["reached_fixed_max"]]
    keys=[b["low_key"] for b in all_old+all_new+all_confirm]
    if len(keys)!=len(set(keys)):
        raise ValueError("RNG low-key overlap: inspect before inference")
    populations={}
    if old:populations["legacy_3M"]=old
    if campaign_complete:populations["original_fixed_16M"]=new
    if confirm_complete:populations["confirmatory_8M"]=confirm
    rng=np.random.default_rng(args.seed)
    all_rows=[];exponents=[];batch_rows=[];diagnostics={}
    for name,batches in populations.items():
        rows=survival_rows(name,batches);all_rows+=rows
        exponents+=exponent_rows(name,rows,args.min_survivors)
        table,diag=batch_diagnostics(name,batches,rng);batch_rows+=table;diagnostics[name]=diag
    comparison=None
    if campaign_complete and confirm_complete:
        a=diagnostics["original_fixed_16M"];b=diagnostics["confirmatory_8M"]
        comparison=dict(time=a["time"],original_survivors=a["survivors"],original_n=a["replicas"],
                        confirmatory_survivors=b["survivors"],confirmatory_n=b["replicas"],
                        fisher_two_sided_p=fisher_two_sided(a["survivors"],a["replicas"],b["survivors"],b["replicas"]),
                        caveat="Independent fixed ensembles, reported separately. Confirmatory8M remains primary regardless of outcomes or completion order.")
    serial=lambda b:{k:v for k,v in b.items() if k not in ["times","values"]}
    result=dict(created_sources=[serial(b) for b in all_old+all_new+all_confirm],
                excluded_incomplete_new=[b["path"] for b in all_new if b["path"] not in {x["path"] for x in new}],
                excluded_incomplete_confirmatory=[b["path"] for b in all_confirm if b["path"] not in {x["path"] for x in confirm}],
                primary="confirmatory_8M",primary_available=confirm_complete,
                interim_figure_population="confirmatory_8M" if confirm_complete else "legacy_3M",
                new_completed_shards=len(new),new_planned_shards=8,invalid_or_empty_new=invalid_new,
                new_campaign_complete=campaign_complete,
                confirmatory_completed_shards=len(confirm),confirmatory_planned_shards=16,
                invalid_or_empty_confirmatory=invalid_confirm,confirmatory_campaign_complete=confirm_complete,
                inference_status=("Completed predeclared confirmatory fixed-N campaign" if confirm_complete else
                                  "Confirmatory campaign incomplete. No inference from a wall-clock-selected subset."),
                all_observed_rng_low_words_disjoint=True,
                min_late_survivors_for_exponent=args.min_survivors,
                exponent_intervals="Conditional binomial retention from t1 to t2; overlapping windows are correlated, no global exponent fit.",
                survival_intervals="Clopper-Pearson for sparse counts<=1000 or complementary counts<=1000; Wilson otherwise.",
                escape_handling="Escapes reported separately; survival CI brackets escaped-failure and escaped-survival interpretations. No local exponent if escapes are present.",
                stopping_rule="No pooled inference from either new campaign until all its declared fixed-N shards complete. Confirmatory8M is always primary; original16M is supplementary. Legacy intervals are nominal because target stopping was configured.",
                legacy_caveat="Legacy target=5 was configured, but every included old shard reached375000 with fewer than5 endpoint survivors.",
                finite_time_only="No claim about P_infinity or generalized-voter exclusion is made.",
                batch_diagnostics=diagnostics,independent_fixed_ensemble_comparison=comparison)
    args.output.mkdir(parents=True,exist_ok=True)
    write_csv(args.output/"survival.csv",all_rows)
    write_csv(args.output/"local_exponents.csv",exponents)
    write_csv(args.output/"independent_batches.csv",batch_rows)
    (args.output/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    for name,diag in diagnostics.items():
        row=next(r for r in all_rows if r["population"]==name and r["time"]==diag["time"])
        print(f"{name}: N={row['n']} S(T)={row['survivors']} E(T)={row['escaped']} "
              f"P(T)={row['survival']:.8g}, 95% CI [{row['ci_low']:.8g},{row['ci_high']:.8g}]")
    print(f"Original fixed16M shards {len(new)}/8; confirmatory8M shards {len(confirm)}/16; incomplete ensembles excluded.")


if __name__=="__main__":main()

