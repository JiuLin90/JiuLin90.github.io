# A function to find out if partition kappa dominates or equal to partition lam
def dominate(kappa,lam):
     if len(kappa)>len(lam):
         return False
     else:
         s1=0;
         s2=0;
         for i in range(len(kappa)):
             s1=s1+kappa[i];
             s2=s2+lam[i];
             if s1<s2:
                 return False
         return True

# A function to find out if mu and lam has at most two distinct elements
def checkmu(lam,mu):
    m1=len(lam);
    m2=len(mu);
    if (m2<m1-1)|(m2>m1):
        return False
    i1=0;
    i2=0;
    c=0;
    while (i2<m2)&(i1<m1):
        if lam[i1]==mu[i2]:
            i1=i1+1;
            i2=i2+1;
        elif lam[i1]<mu[i2]:
            i2=i2+1;
        else:
            i1=i1+1;
            c=c+1;
    c=c+m1-i1;
    return c<=2

# A function to compute the transition matrix from Jack polynomials to 
# monomials of degree k.  It takes two optional arguments:
# alpha: default is 2 
# normalization: 'J','C','P', or 'Q', default is 'J'
def JacktoM(k,*arg):
    if len(arg)==0:
        alpha=2;
        norm='J';
    else:
       alpha=arg[0];
       if len(arg)==1:
           norm='J';
       else:
           norm=arg[1];
    ai = 2/alpha;
    A=Partitions(k).list();
    n=len(A);
    if (norm=='C')|(norm=='Q'):
        ck=[prod(flatten(A[i].upper_hook_lengths(alpha))) for i in range(n)];
# M is a transition matrix from J-normalization of Jack polynomials to monomials.
# It has the attractive property that all elements of this transition matrix are integers.
    M=identity_matrix(ZZ,n);
    for i in range(n):
        M[i,i]=A[i].hook_product(alpha);    # Use internal function to compute the diagonal elements of M
        A[i]=list(A[i]);
    kappa_l=[len(A[i]) for i in range(n)];
    rho=[sum(A[i][j]*(A[i][j]-1-ai*j) for j in range(kappa_l[i])) for i in range(n)];     
    M[:,-1]=factorial(k);               # last column of M is k!
# Create a matrix with coefficients a_lam(mu)
    n1= round(n*n^(2/5)+6);
    a=zero_vector(n1);
    alist=zero_vector(n1);
    count=zero_vector(n-1);
    for jj in range(1,n-1):
        count[jj]=count[jj-1];  
        lam=A[jj];
        m2=kappa_l[jj];
#       lam_c=Partition(lam).to_exp();
        for kk in range(jj):
            mu=A[kk];
            if checkmu(lam,mu):
                s=m2-1;
                if kappa_l[kk]==m2:
                    while lam[s]==mu[s]:
                        s=s-1;
                    t=mu[s];
                else:
                    t=0;
                while lam[s]==mu[s-1]:
                    s=s-1;
                t=lam[s]-t;
                r=s-1;
                while (mu[r]!=lam[r]+t)&((lam[r]==mu[r])|(lam[r-1]!=mu[r])):
                    r=r-1;
                w=lam.count(lam[r]);
                if lam[r]==lam[s]:
                    w=w*(w-1)/2;
                else:
                    w=w*lam.count(lam[s]);
# Store all the nonzero entries a_lam(mu)
                a[count[jj]]=(lam[r]-lam[s]+2*t)*w
                alist[count[jj]]=kk;
                count[jj]=count[jj]+1;  
# Compute the coefficients of M using recursion. 
    for ii in range(n-1):
        kappa=A[ii];
        for jj in range(ii+1,n-1):
            lam=A[jj];
            if dominate(kappa,lam):
                M[ii,jj]=sum(a[kk]*M[ii,alist[kk]] for kk in range(count[jj-1],count[jj]) if alist[kk]>=ii)*ai/(rho[ii]-rho[jj]); 
    if norm!='J':    
        M=M.base_extend(QQ);
        if norm=='P':
            for ii in range(n):
                M[ii,ii:]=M[ii,ii:]/M[ii,ii];
        elif norm=='Q':
            for ii in range(n):
                M[ii,ii:]=1/ck[ii]*M[ii,ii:];
        else:
            c0=alpha**k*factorial(k);
            for ii in range(n):
                M[ii,ii:]=c0/M[ii,ii]/ck[ii]*M[ii,ii:];
    return M

# Compute the transitition matrix from Jack polynomials to power-sum symmetric polynomials
def JacktoP(k,*arg):
    if len(arg)==0:
        alpha=2;
        norm='J';
    else:
       alpha=arg[0];
       if len(arg)==1:
           norm='J';
       else:
           norm=arg[1];
    m=SymmetricFunctions(QQ).m();
    p=SymmetricFunctions(QQ).p();
    A=JacktoM(k,alpha,norm)*m.transition_matrix(p,k);
    return A

# Compute the transitition matrix from power-sum symmetric polynomials to Jack polynomials.
# Note that we do not use the inverse of JacktoP(k,alpha,norm) to obtain this.
def PtoJack(k,*arg):
    if len(arg)==0:
        alpha=2;
        norm='J';
    else:
       alpha=arg[0];
       if len(arg)==1:
           norm='J';
       else:
           norm=arg[1];
    B=JacktoP(k,alpha,'J').transpose();
    n=B.nrows();
    for i in range(n):
        B[i,:]=B[i,:]/B[i,0];
    if norm!='C':
        B=alpha**k*factorial(k)*B;
        p=Partitions(k).list();
        if norm=='P':
            for ii in range(n):
                B[:,ii]=B[:,ii]/prod(flatten(p[ii].upper_hook_lengths(alpha)));
        elif norm=='J':
            for ii in range(n):
                B[:,ii]=B[:,ii]/prod(flatten(p[ii].upper_hook_lengths(alpha)))/p[ii].hook_product(alpha);
        else:
            for ii in range(n):
                B[:,ii]=B[:,ii]/p[ii].hook_product(alpha);
    return B

# Compute the transition matrix from monomial symmetric polynomials to Jack polynomials
def MtoJack(k,*arg):
    if len(arg)==0:
        alpha=2;
        norm='J';
    else:
       alpha=arg[0];
       if len(arg)==1:
           norm='J';
       else:
           norm=arg[1];
    m=SymmetricFunctions(QQ).m();
    p=SymmetricFunctions(QQ).p();
    C=m.transition_matrix(p,k);
    B=(JacktoM(k,alpha,'J')*C).transpose();
    n=B.nrows();
    for i in range(n):
        B[i,:]=B[i,:]/B[i,0];
    if norm!='C':
        B=alpha**k*factorial(k)*B;
        p=Partitions(k).list();
        if norm=='P':
            for ii in range(n):
                B[:,ii]=B[:,ii]/prod(flatten(p[ii].upper_hook_lengths(alpha)));
        elif norm=='J':
            for ii in range(n):
                B[:,ii]=B[:,ii]/prod(flatten(p[ii].upper_hook_lengths(alpha)))/p[ii].hook_product(alpha);
        else:
            for ii in range(n):
                B[:,ii]=B[:,ii]/p[ii].hook_product(alpha);
    return C*B
