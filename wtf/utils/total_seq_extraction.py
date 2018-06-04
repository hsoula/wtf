__author__ = 'hsoula'

from numpy import *
import cPickle
import copy
import sys
try:
    import sndfile
except:
    print "pb with sndfile"
    pass
import matplotlib.pyplot as plt
from time import gmtime, strftime
import time
sys.path.append('../../')
from wtf.zf import wt_project,wt_decimate,wt_parse,wt_waves
from wtf import cl
from wtf.utils import extract_calls

def get_current_stage(wav_dir,projname):
    import os
    stage=-1
    for x in os.listdir(wav_dir):
        if x.endswith('.dump') :
            if int(x[:-5].split("-")[-1])>stage:
              stage=int(x[:-5].split("-")[-1])
    return stage

#print 'testing: wt_project'
print "Date:",strftime("%Y-%m-%d %H:%M:%S", gmtime())

print "######### PARAM ##############"
f_fftw=1024 #500# 128
f_div=100
Hz_LowPass=1000.0
dB = float(sys.argv[2]) # threshold for sound extraction from envelope
samplerate = float(sys.argv[3])
chunk_size=220500#int(10*samplerate) # 10*44100
f_pad=1012 #500+f_fftw/2
peak_pc = 0.05
thr_parse = 662 #8410#1765
div_parse = 45
thr_dec=1000 # decimate threshold
cut_noise_freq = 100 # not used
bool_keep_calls=False ## calls are not individually extracted and kept in the "calls" directory

## classif parameters :
nt=10
fmin=500
fmax=20000
nnt=2
nnf=2
f_fftw_classif = 1024 #f_fftw
file_train = "" 
classif_dump = "/Volumes/INRIA3/dippers/classifier_cincles.dump"#Nora_pairs_WR/classifier_small_silentbox/classif_keep_xms_100.dump"##"/Users/marie/Documents/recherche/wtf/cl/classifier/classif_ste_keep_xms_35pcnoncall.dump"
choice_classifier='rf'
ntr = -1
type_build_classif = 'keep_xms'#'reduce'

nb_calls_to_extract = -1
nb_calls_to_get = -1
wav_dir= sys.argv[1]
outdir= wav_dir+"/calls"

print "wav_dir=",wav_dir
print "out dir=",outdir
print "f_fftw=",f_fftw,"\nf_div=",f_div,"\nHz_LowPass=",Hz_LowPass,"\ndB=",dB,"\nsamplerate=",samplerate,"\nchunk_size=",chunk_size,"\nfpad=",f_pad,"\npeak_pc=",peak_pc
print "\ncl_parse parameters:"
print "thr_parse =",thr_parse,"\ndiv_parse=",div_parse,"\nthr_dec=",thr_dec
print "\nClassif parameters:"
print "nt=",nt,"fmin=",fmin,"fmax=",fmax,"nnt=",nnt,"nnf=",nnf,"f_fftw_classif=",f_fftw_classif,"file_classif=",classif_dump,"ntr=",ntr
print "nb_calls_to_extract=",nb_calls_to_extract,"nb_calls_to_get =",nb_calls_to_get 
print "classifier=",choice_classifier

t0= time.clock()
wt_pro=wt_project(wav_dir,f_fftw,f_div,Hz_LowPass,dB,chunk_size)
print "NAME:",wt_pro.name
print "\nChannels :"
for i in range(len(wt_pro.files)) : print "Channel",i,":",wt_pro.files[i] 
stage=-1

stage= get_current_stage(wav_dir,wt_pro.name)
print "stage",stage
wt_pro.stage = stage
print "\ncurrent stage=",stage
t0= time.time()
if stage<0:
    print "\n***** Step 1 : get sounds"
    wt_pro.get_sounds() # keep double calls
    print "nb_sounds:",len(wt_pro.s_sounds)
    wt_pro.stage+=1
    wt_pro.dump()

if stage<1:
    print "\n***** Step 2 : remove double calls (""parse"")"
    wt_pro.load()
    wt_pro.get_parse_sounds(-1,thr_parse,div_parse) # remove double calls
    print "nb_sounds:",len(wt_pro.s_sounds)
    wt_pro.stage+=1
    wt_pro.dump()



if stage<2:
    print "\n***** Step 3 : get calls"
    wt_pro.load()
    calls=wt_pro.get_calls(f_pad,peak_pc,nb_calls_to_get,0)
    print "nb_calls:",len(wt_pro.calls)
    
    wt_pro.stage+=1
    wt_pro.dump()



if stage<3:
    print "\n***** Step 4 : decimate : group calls"
    wt_pro.load()
    print "Calls before decimate"
    print "NUMBER OF CALLS:",len(wt_pro.calls)
    if len(wt_pro.calls)!=0:
        cl_dec=wt_decimate(wt_pro.calls)
        dec_calls=cl_dec.decimate(thr_dec)
        dec_calls=array(sorted(dec_calls,key=lambda x:x[1]))
    
        wt_pro.calls = dec_calls[:]
        print "nb_calls:",len(wt_pro.calls)

    wt_pro.stage+=1
    wt_pro.dump()





if stage<4:
    print "\n***** Step 5 : extract calls"
    wt_pro.load()
    cl_wav=wt_waves([wt_pro.wav_dir + '/'+f for f in wt_pro.files],wt_pro.calls)
    ## create call list.txt in calls directory
    wav_l=cl_wav.extract(cut_noise_freq,nb_calls_to_extract,outdir,keep_calls=bool_keep_calls) #f_fftw_classif,f_div,f_pad,dB,peak_pc,
    print "after extract calls"
    print "nb_calls:",len(wt_pro.calls)

    wt_pro.stage+=1
    wt_pro.dump()


if stage<5 :
    print "\n***** Step 6 : remove ""non-calls" ""
    v=0
    file_to_classify=outdir+'/call_list.txt'
    wt_pro.load()
    
    if len(wt_pro.calls)!=0:
        cl_classif = cl.cl_classifier(samplerate,fmin,fmax,f_fftw_classif,f_div,choice_cl=choice_classifier,type_build=type_build_classif,nt=-1,nnt=-1,nnf=-1,nbtree = 10,nb_ms=-1,thrpc_noncall=-1,verbose=0)
        cl_bin = cl_classif.classify(file_to_classify,classif_dump,verbose=v)
        wt_pro.res_classif=cl_bin ## list 0 or 1 if call or not

        calls = []
        nums=[]
        lc = len(wt_pro.calls)
        for i in range(lc):#nb_calls_to_extract):
            if int(cl_bin[i]) == 1 : #or int(cl_bin[i]) == 0:
                calls.append(wt_pro.calls[i])
                nums.append(i)
        calls=array(calls)
        wt_pro.calls = calls[:]
        wt_pro.numcalls = nums[:]
    f=open(wav_dir+'/'+'final_call_list.txt','w')
    for i,c in enumerate(wt_pro.calls):
        fnameinit='call-%d-%.3lf-%d'%(nums[i],c[1]*1./samplerate,c[0])+'.wav'
        f.write('%s %.3lf %.3lf %d\n'%(fnameinit,c[1]*1./samplerate,c[2]*1./samplerate,c[0]))
    f.close()
    print "nb calls:",len(wt_pro.calls)
    wt_pro.wav_dir=wav_dir

    wt_pro.stage+=1
    wt_pro.dump()

if stage<6 :

    print "\n***** Step 7 : save call sequence and print time course"
    wd = "_".join(wav_dir.split("/"))
    wt_pro.load()
    print wt_pro.name
    if wav_dir.endswith('/'):
        f_callseq = wav_dir+"/call_seq_"+wav_dir.split('/')[-2]+".txt"
    else:
        f_callseq = wav_dir+"/call_seq_"+wav_dir.split('/')[-1]+".txt"

    print "call sequence saved in",f_callseq
    div = 10# in frames
    
    #wt_pro.time_course(div)
    wt_pro.save_callseq(f_callseq,wt_pro.numcalls)
    print "coucou",f_callseq, wav_dir
    if len(wt_pro.calls)!=0:
        extract_calls(f_callseq,wav_dir,100,samplerate=samplerate)

    #time.sleep(1000)
    print round(time.time()-t0,2),"sec"
    wt_pro.stage+=1

    # for x in wt_pro.s_sounds:
    #     print x[0],x[1]/44100.0,x[2]/44100.0,(x[2]-x[1])/44100.*1000.
