def Calmi(p):                         #For a paritition with m_1 1's, m_2 2's and etc, return (m_1)!(m_2)!...(m_l)!. This is the reciprocal of the leading coefficient of M-polynomial.#
    t=Partition(p).to_exp()           #exponential form of the partition {m_1,m_2,...}#
    re=prod(factorial(i) for i in t)  #re=(m_1)!(m_2)!...(m_l)!#
    return re;                        #return the product.#

def listexp(list1,list2):                       #Given two lists of same length, say {a_1,a_2,...,a_n} and {b_1,b_2,...,b_n} return (a_1)^(b_1)(a_2)^(b_2)...(a_n)^(b_n)#
    n=len(list1);       
    re=[list1[i]^list2[i] for i in range(n)];
    return prod(i for i in re); 

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

def MZonal(partition,variables):              #Computing M-polynomial by (3.2)#
    re=0;                                     #return#
    m=len(partition);                         #m=length of partition#
    n=len(variables);                         #n=number of variables#
    if n>=m:                                   #If n<m, index is 1 and return 0#
        perm=Permutations(n,m).list();        #By (3.2), for the summation part, we choose all posible combination of the variables for the summand#
        for i in perm:                        #For each summmand,#
            temp=[variables[j-1] for j in i]; #pick the right combination of the variables#
            re=re+listexp(temp,partition);    #Get the term (y_1)^(\lambda_1)(y_2)^(\lambda_2)...(y_m)^(\lambda_m)#
        re=re/Calmi(partition);                   #Divide the sum by the leading coefficient.#
    return re;

#Given partitions lam<kappa, return all nonzero c_{mu,lam} with kappa>=mu>=lam#
def Lcoeffi(kappa,lam,*arg):                  
    if len(arg)==0:
        norm='J';
    else:
        norm=arg[0];     
    A=kappa.dominated_partitions();  #only consider partitions <=k#
    k=sum(kappa);
    if len(lam)<k:
        count=0;
        while A[count]>=lam:
            count=count+1;
        A=A[:count];    
    m=len(A);                #length of this partial list#
    M=list(zero_vector(m));
    M[0]=kappa.hook_product(2);
    if len(lam)==k:
        M[-1]=factorial(k);
    for i in range(m):
        A[i]=list(A[i]);
    kappa_l=[len(A[i]) for i in range(m)];
    rho=[sum(A[i][j]*(A[i][j]-j-1) for j in range(kappa_l[i])) for i in range(m)];     
# Create a matrix with coefficients a_lam(mu)
    m1=m-ZZ(len(lam)==k);
# Compute the coefficients of M using recursion. 
    for jj in range(1,m1):
        lam1=A[jj];
        if dominate(kappa,lam1):
            m2=kappa_l[jj];
            c=0;
            for kk in range(jj):
                mu=A[kk];
                if checkmu(lam1,mu):
                    s=m2-1;
                    if kappa_l[kk]==m2:
                        while lam1[s]==mu[s]:
                            s=s-1;
                        t=mu[s];
                    else:
                        t=0;
                    while lam1[s]==mu[s-1]:
                        s=s-1;
                    t=lam1[s]-t;
                    r=s-1;
                    while (mu[r]!=lam1[r]+t)&((lam1[r]==mu[r])|(lam1[r-1]!=mu[r])):
                        r=r-1;
                    w=lam1.count(lam1[r]);
                    if lam1[r]==lam1[s]:
                        w=w*(w-1)/2;
                    else:
                        w=w*lam1.count(lam1[s]);
                    c = c+(lam1[r]-lam1[s]+2*t)*w*M[kk];
            M[jj]=c/(rho[0]-rho[jj]); 
    if norm!='J':    
        if norm=='P':
            M=[M[i]/M[0] for i in range(m)];
        elif norm=='Q':
            c0=prod(flatten(kappa.upper_hook_lengths(2)));
            M=[M[i]/c0 for i in range(m)];
        else:
            c0=2**k*factorial(k)/prod(flatten(kappa.upper_hook_lengths(2)))/M[0];
            M=[M[i]*c0 for i in range(m)];
    return M

def Coeffi(k,l,*arg): 
    if l>k:
        return 0;
    else:
        if len(arg)==0:
            norm='J';
        else:
            norm=arg[0];     
        return Lcoeffi(k,l,norm)[-1];    

def FLcoeffi(kappa,*arg):                      #Full List of coefficients c_{kappa,lambda} for all lambda<=kappa#
    if len(arg)==0:
        norm='J';
    else:
        norm=arg[0];
    k = sum(kappa);
    return Lcoeffi(kappa,Partition([1]*k),norm)

def CZonal(k,v,*arg):                                                     #Given partition k and variables v, compute C-polynomial C_{k}(v)#
    if len(arg)==0:
        norm='J';
    else:
        norm=arg[0];     
    partiallist=k.dominated_partitions();                                 #only consider partitions <=k#
    coefftable=FLcoeffi(k,norm);                                          #list of all coefficients c_{k,l} for l<=k#
    Mtable=[MZonal(list(t),v) for t in partiallist];                      #list of all corresponding M_{l}(v)#
    re=sum(coefftable[t]*Mtable[t] for t in range(len(partiallist)));     #(3.3)#
    return re; 


# A function to compute the transition matrix from zonal polynomials to 
# monomials of degree k.  It takes an optional argument:
# normalization: 'J','C','P', or 'Q', default is 'J'
def ZonaltoM(k,*arg):
    if len(arg)==0:
        norm='J';
    else:
        norm=arg[0];     
    A=Partitions(k).list();
    n=len(A);
    if (norm=='C')|(norm=='Q'):
        ck=[prod(flatten(A[i].upper_hook_lengths(2))) for i in range(n)]; 
# M is a transition matrix from J-normalization of zonal polynomials to monomials.
# It has the attractive property that all elements of this transition matrix are integers.
    M=identity_matrix(ZZ,n);
    for i in range(n):
        M[i,i]=A[i].hook_product(2);    # Use internal function to compute the diagonal elements of M
        A[i]=list(A[i]);
    kappa_l=[len(A[i]) for i in range(n)];
    rho=[sum(A[i][j]*(A[i][j]-j-1) for j in range(kappa_l[i])) for i in range(n)];     
    M[:,-1]=factorial(k);              # last column of M is k!
# Create a matrix with coefficients a_lam(mu)
    n1= round(n*n^(2/5)+6); 
    a=zero_vector(n1);
    alist=zero_vector(n1);
    count=zero_vector(n-1);
    for jj in range(1,n-1):
        count[jj]=count[jj-1]; 
        lam=A[jj];
        m2=kappa_l[jj];
        lam_c=Partition(lam).to_exp();
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
                if lam[r]==lam[s]:
                    w=lam_c[lam[r]-1];
                    w=w*(w-1)/2;
                else:
                    w=lam_c[lam[r]-1]*lam_c[lam[s]-1];
# Store all the nonzero entries a_lam(mu)
                a[count[jj]]=(lam[r]-lam[s]+2*t)*w;
                alist[count[jj]]=kk;
                count[jj]=count[jj]+1; 
# Compute the coefficients of M using recursion. 
    for ii in range(n-1):
        kappa=A[ii];
        for jj in range(ii+1,n-1):
            lam=A[jj];
            if dominate(kappa,lam):
                M[ii,jj]=sum(a[kk]*M[ii,alist[kk]] for kk in range(count[jj-1],count[jj]) if alist[kk]>=ii)/(rho[ii]-rho[jj]);
    if norm!='J':    
        M=M.base_extend(QQ);
        if norm=='P':
            for ii in range(n):
                M[ii,ii:]=M[ii,ii:]/M[ii,ii];
        elif norm=='Q':
            for ii in range(n):
                M[ii,ii:]=M[ii,ii:]/ck[ii];
        else:
            c0=2**k*factorial(k);
            for ii in range(n):
                M[ii,ii:]=c0/M[ii,ii]/ck[ii]*M[ii,ii:]; 
    return M
