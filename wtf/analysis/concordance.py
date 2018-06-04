from numpy import *
import matplotlib.pyplot as plt
################################
## METHODE DE CONCORDANCE ######
################################
## cf Joel Attia

def mu(p,M):
    #p : vect proba d'apparition des elements
    #M : matrice de concordance
    sum1 = sum([i**2 for i in p])
    sum2 = 0
    nb_elts = len(p)
    for j in range(nb_elts) :
        for i in range(j):
            sum2+=(M[i,j]*p[i]*p[j])
    return sum1+2*sum2,sum2

def sigma(p,M):
    sum1 = sum([i**2 for i in p])
    res = mu(p,M)
    m = res[0]
    sum2 = res[1]
    return sum1+2*sum2-m**2

def val_obs_Cmsigma(seq,m,M):
    #seq : sequence
    #m : time step
    #M : matrice de concordance
    N=len(seq)
    Sm=seq[:N-m]
    Tm=seq[m:]
    # s = 0
    # for i in range(N-m):
    #     s+=M[Sm[i],Tm[i]]
    # print s
    return sum([M[Sm[i],Tm[i]] for i in range(N-m)])*1./(N-m)#sum([M[i,j] for i in Sm for j in Tm])

def Uobs(cmsigm,p,N,m,M):
    return (cmsigm-mu(p,M)[0])*1./(sigma(p,M)/sqrt(N-m))

def lim_curve(u,N,p,M,m):
    return mu(p,M)[0]+u*(sigma(p,M)*1./sqrt(N-m))

def plot_concord(seq,M,vect_m,u,fileout=-1):
    elts = list(set(seq))
    N = len(seq)
    p=[list(seq).count(e)*1./N for e in elts]
    
    valobs=[val_obs_Cmsigma(seq,m,M) for m in vect_m]
    limc=[lim_curve(u,N,p,M,m) for m in vect_m]
    plt.figure()
    plt.plot(vect_m,valobs)
    plt.plot(vect_m,limc)
    if fileout==-1: plt.show()
    else : plt.savefig(fileout)
    


#######################################################
## tests ##
#######################################################
# seq = array([0,0,0,1,0,2,0,1,0,0,0,1,0,2,0,1,0,0,0,1,0,1,0,2,0,1,0,0,0,0])
# M = array([[1,0,0],[0,1,0],[0,0,1]])
# seq = array([0,0,0,1,0,1,0,1,1,1,0,1,0,1,0,1,0,1,0,0,0,1,0,0,1,1,1,0,0,0])
# seq = random.randint(0,2,10000)
# M = array([[1,0],[0,1]])
# plot_concord(seq,M,vect_m=range(0,10,1),u=1.64)

import sys
import get_callseq as gcs
import get_projname as gpj
from get_seq_with_silences import get_seq_with_silences

proj_dir = sys.argv[1]
projname= gpj.get_projname(proj_dir)
plotsdir = "./plots_concordance/"
callseq,nb_birds = gcs.get_callseq(proj_dir)
bin_choice = float(sys.argv[2])
callseq = callseq[:10000]
seq = get_seq_with_silences(callseq,bin_choice)[1]
u05 = 1.64
f = plotsdir+"concord_"+projname+".pdf"
v = range(0,6000,1)
M = identity(nb_birds+1)
#print seq
print M
plot_concord(seq,M,vect_m=v,u=u05,fileout=f)
