from numpy import *
def import_table(fname,delimit='\t',h=False):
    ## fname : file with columns
    ## return 
    from numpy import fromfile,genfromtxt,loadtxt
    f=open(fname,'r')
    #l=f.readline()
    #print l
    #nc=len(l.split())
    # l1=f.readlines()
    # for l2 in l1:
    #     print len(l2.split(',')),l2.count(',')
    #f.close()
    # print l
    
    
    g=loadtxt(fname,delimiter=delimit,dtype=str)
    #print g
    #print nc
    #g=reshape(g,(nc,-1))
    if h==True:g=g[1:]
    return g
