
#######################################################
### extract bursts from a call sequence ###############
#######################################################
import sys
sys.append('../')
from utils.seq_analysis import get_nbcalls_bin_all
import time
from numpy import *
from matplotlib import pyplot as plt

def extract_bursts(callseq,nbsecsilence,nbcall_mini):
    # callseq : call sequence start end numbird
    # bin_choice : sampling freq. ex : bin_choice=0.01 : sequence cut into 10ms bins
    # nb_sec_silence : number of second of silence required to separate 2 bursts. (ex:nb_sec_silence=10sec : if 10 seconds without a call, create 2 bursts)

    ## return :
    # list of bursts with silences [[0,1,0,0,-1,0],[0,0,1,2,-1]]
    # bursts limits in '%H:%M:%S' format : [(00:25:01,00:28:06),(01:56:23,02:50:54)]

    d = diff(callseq[:,0])
    l = len(d)
    burst=[]
    prec=0
    for i,di in enumerate(d):
        if di>nbsecsilence:
            b = callseq[prec:i]
            if len(b)>nbcall_mini :
                burst.append(b)
                fb=open(proj_dir+"/bursts/burst_"+str(nbsec)+"_sec_"+str(len(burst)-1)+".txt","w")
                for c in b :fb.write('%f %f %d\n'%(c[0],c[1],c[2]))
                fb.close()
                print 'burst : ',callseq[prec,0],callseq[i-1,0],len(b)
            prec = i
    b = callseq[prec:]

    if len(b)>nbcall_mini : 
        burst.append(b)
        fb=open(proj_dir+"/bursts/burst_"+str(nbsec)+"_sec_"+str(len(burst)-1)+".txt","w")
        for c in b :
            fb.write('%f %f %d\n'%(c[0],c[1],c[2]))
        fb.close()

    #plt.hist(d[d>10],120)
    #plt.show()

def extract_bursts_withfreq(callseq,bin_choice,freq_min):
    ## extract burst using instantaneous call rate
    ## freq_min = threshold for call rate in call/min
    nbcalls = get_nbcalls_bin_all(callseq,bin_choice)
#    print nbcalls,freq_min,bin_choice
    freq=freq_min*bin_choice*1./60
    offset=callseq[0,0]
    start=0+offset
    burst=[]
    burst_opened=False
    nbburst=0
    for i,n in enumerate(nbcalls):
        if n>=freq : ## burst
            if burst_opened==False: ## if no current burst : create one
                nbburst+=1
                start = offset+i*bin_choice                            
                burst_opened=True
        else : ## no burst
            if burst_opened==True: ## if one burst open : close it
                c1=callseq[callseq[:,0]>=start]
                b=c1[c1[:,0]<=offset+i*bin_choice]
                burst.append(b)
                burst_opened=False
                #print 'burst=',start,offset+i*bin_choice
            #fb=open(proj_dir+"/bursts/burst_"+str(nbsec)+"_sec_"+str(len(burst)-1)+".txt","w")
            #for c in b :fb.write('%f %f %d\n'%(c[0],c[1],c[2]))
            #fb.close()
    if burst_opened==True :
        b = callseq[callseq[:,0]>=start]
        burst.append(b)
        #print 'burst : ',start,'end'
    return nbburst,burst
            


# ######################################################
# ### Tests ###
# ######################################################
#callseq = [0,0,0,1,0,1,0,1,1,1,1,1,0,0,0,-1,-1,-1,-1,0,0,1,0,1,-1,0,-1,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,1,1,0,1,0]

# callseq=array([[ 154.369887,  154.458322,    0.      ],
#                [ 156.090975,  156.174875,    3.      ],
#                [ 165.227029,  165.333605,    0.      ],
#                [ 170.088707,  170.202086,    2.      ],
#                [ 190.174875,  190.265578,    0.      ],
#                [ 198.349478,  198.428844,    3.      ],
#                [ 251.306395,  251.390295,    0.      ],
#                [ 252.054694,  252.122721,    0.      ]])
# callseq=array([[0,1,0],
#                [101.6,110.7,0.],
#                [110.75,110.76,0.],
#                [111.2,111.3,0.],
#                [121.1,121.2,0.],
#                [122.6,122.7,0.],
#                [303,304,0.],
#                [1201.,1202.,3.],
#                [1202.,1203.,4.],
#                [1204.,1205.,4.],
#                [1245,1246,5.],
#                [1246.8,1246.9,2],
#                [1500,1520,2]])

# extract_bursts_withfreq(callseq,30,5) #5 calls/min minimum
#print lim_b,bursts
# print callseq
#new_seq = extract_sequence_bursts(callseq,lim_b)
#print new_se
#bursts,lim_bursts = extract_bursts(callseq,bin_choice=0.1,nb_sec_silence=nbsec)
    # for l in lim_bursts:
    #     print str(l),len(bursts[l])
    # new_seq = extract_sequence_bursts(callseq,lim_bursts,proj_dir,nbsec)
    # f=open(proj_dir+"/callseq_bursts.txt","w")
    # for i in new_seq:f.write("%f %f %d\n"%(i[0],i[1],i[2]))
    # f.close()




