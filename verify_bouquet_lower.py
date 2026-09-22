#!/usr/bin/env python3
import json
from collections import deque
from pathlib import Path


def construct(k,t,q):
    m=t-1; L=m*k; adj=[]
    def vertex():
        adj.append(set()); return len(adj)-1
    def edge(a,b):
        adj[a].add(b); adj[b].add(a)
    x=vertex(); hubs=[x]; paths=[]
    for a in range(q):
        y=vertex(); hubs.append(y)
        for b in range(m):
            p=[x]+[vertex() for _ in range(2*k)]+[y]; paths.append(p)
            for u,v in zip(p,p[1:]): edge(u,v)
            for _ in range(L): edge(p[1],vertex())
    return adj,hubs,paths


def distances(adj,root,cutoff=None):
    dist={root:0}; queue=deque([root])
    while queue:
        u=queue.popleft()
        if cutoff is not None and dist[u]==cutoff: continue
        for v in adj[u]:
            if v not in dist: dist[v]=dist[u]+1; queue.append(v)
    return dist


def verify(k,t,q):
    adj,hubs,paths=construct(k,t,q); n=len(adj); m=t-1
    balls=[distances(adj,v,k) for v in range(n)]
    sizes=list(map(len,balls))
    choices=[max(ball,key=lambda w:(sizes[w],w)) for ball in balls]
    for p in paths:
        for i in range(1,k+1):
            assert choices[p[k+i]]==p[i],(k,t,q,i)
            assert all(sizes[p[i]] > sizes[w]
                       for w in balls[p[k+i]] if w != p[i]), (k,t,q,i,'tie')
            expected=1+q*m*k-(q*m-2)*i+m*k*(1+(q*m-1)*(i<=k-2))
            assert sizes[p[i]]==expected
        for i in range(k+1,2*k+1):
            d=2*k+1-i
            assert sizes[p[i]]==1+m*k-(m-2)*d
    assert set().union(*(set(balls[h]) for h in hubs))==set(range(n))
    for i,h in enumerate(hubs):
        dist=distances(adj,h)
        for h2 in hubs[i+1:]: assert dist[h2]>2*k
    return dict(k=k,t=t,q=q,n=n,selected=len(set(choices)),optimum=q+1,
                certified_selected=q*k*m,ratio=len(set(choices))/(q+1))

if __name__=='__main__':
    cases=[(k,t,q) for k in range(2,9) for t in range(3,9) for q in (2,3,5)]
    cases += [(10,10,10),(20,3,10)]
    results=[verify(*case) for case in cases]
    Path(__file__).with_name('bouquet_lower_verification.json').write_text(json.dumps({'instances':len(results),'failures':0,'results':results},indent=2))
    print(json.dumps({'instances':len(results),'failures':0,'last':results[-2:]}))
