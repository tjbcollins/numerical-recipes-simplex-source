      subroutine amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      parameter (nmax=20,alpha=1.0,beta=0.5,gamma=2.0,itmax=500)
      dimension p(mp,np),y(mp),pr(nmax),prr(nmax),pbar(nmax)
      mpts=ndim+1
      iter=0
1     ilo=1
      if(y(1).gt.y(2))then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 11 i=1,mpts
        if(y(i).lt.y(ilo)) ilo=i
        if(y(i).gt.y(ihi))then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi))then
          if(i.ne.ihi) inhi=i
        endif
11    continue
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if(rtol.lt.ftol)return
      if(iter.eq.itmax) pause 'amoeba exceeding maximum iterations.'
      iter=iter+1
      do 12 j=1,ndim
        pbar(j)=0.
12    continue
      do 14 i=1,mpts
        if(i.ne.ihi)then
          do 13 j=1,ndim
            pbar(j)=pbar(j)+p(i,j)
13        continue
        endif
14    continue
      do 15 j=1,ndim
        pbar(j)=pbar(j)/ndim
        pr(j)=(1.+alpha)*pbar(j)-alpha*p(ihi,j)
15    continue
      ypr=funk(pr)
      if(ypr.le.y(ilo))then
        do 16 j=1,ndim
          prr(j)=gamma*pr(j)+(1.-gamma)*pbar(j)
16      continue
        yprr=funk(prr)
        if(yprr.lt.y(ilo))then
          do 17 j=1,ndim
            p(ihi,j)=prr(j)
17        continue
          y(ihi)=yprr
        else
          do 18 j=1,ndim
            p(ihi,j)=pr(j)
18        continue
          y(ihi)=ypr
        endif
      else if(ypr.ge.y(inhi))then
        if(ypr.lt.y(ihi))then
          do 19 j=1,ndim
            p(ihi,j)=pr(j)
19        continue
          y(ihi)=ypr
        endif
        do 21 j=1,ndim
          prr(j)=beta*p(ihi,j)+(1.-beta)*pbar(j)
21      continue
        yprr=funk(prr)
        if(yprr.lt.y(ihi))then
          do 22 j=1,ndim
            p(ihi,j)=prr(j)
22        continue
          y(ihi)=yprr
        else
          do 24 i=1,mpts
            if(i.ne.ilo)then
              do 23 j=1,ndim
                pr(j)=0.5*(p(i,j)+p(ilo,j))
                p(i,j)=pr(j)
23            continue
              y(i)=funk(pr)
            endif
24        continue
        endif
      else
        do 25 j=1,ndim
          p(ihi,j)=pr(j)
25      continue
        y(ihi)=ypr
      endif
      go to 1
      end
