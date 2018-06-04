import sys
sys.path.append("../")
import pysndfile as sndfile
import import_table
import numpy as np
import os


def extract_call_from_wavfile(start,end,sfile,samplerate,outfile, normalize=True,channels=1):
    ## start of call in sec
    ## end of call in sec
    ## sfile : sndfile already open
    
    start=int(start*samplerate)
    end=int(end*samplerate)
    
    sfile.seek(start,0) # go to the start of the call in the file
    wave=sfile.readf_float(end-start)
    #print end,start,sfile.info.frames
    outinfo = sndfile.SF_INFO(samplerate=samplerate, channels = channels,format = (sndfile.SF_FORMAT_WAV|sndfile.SF_FORMAT_PCM_16),sections = 1,seekable = 1)
    o_sfile = sndfile.open(outfile,mode="w",info=outinfo)
    if normalize==True:
        o_sfile.write_double(wave/max(abs(wave)))
    else:
        o_sfile.write_double(wave)
    o_sfile.close()



def extract_calls(calllist_file,workdir,nb_calls_per_bird=-1,shuffle=False,samplerate=44100, colstart=0, colend=1,colid=2,colname=3):
    
    calllist=np.array(import_table(calllist_file,delimit=' ',h=False))
    if shuffle==True:
        np.random.shuffle(calllist)
    
    try:
        nbbirds=len(set(calllist[:,colid]))
    except: ## call list with only one call
        nbbirds=int(calllist[2])
        pass
    sfiles=[]
    for wavf in sorted(os.listdir(workdir)):
        if (wavf.endswith(".WAV") or wavf.endswith(".wav")) and not wavf.startswith("."):
            sfiles.append(sndfile.open(workdir+'/'+wavf,'r'))
    print sfiles
    l=len(calllist)
    nbcalls=[0]*nbbirds
    i=0
    while nbcalls!=[nb_calls_per_bird]*nbbirds and i<l:
        c = calllist[i]
        print c
        bird=int(c[colid])
        
        if nb_calls_per_bird==-1 or nbcalls[bird]<nb_calls_per_bird:
            if colname==3:
                nn = c[colname].split('-')[-1]
            if colname==0:nn=c[colname].split('-')[1]
            fname=workdir+"/calls/"+'call-%d-%.3lf-%d'%(int(nn),float(c[colstart]),int(bird))+'.wav'
            sfile=sfiles[int(c[colid])] # get the file of the bird we want to extract the call
            start=float(c[colstart])
            end=float(c[colend])
            extract_call_from_wavfile(start,end,sfile,samplerate,fname)
            nbcalls[bird]+=1

        i+=1

    for s in sfiles:
        s.close()

def test_extract_calls():
    #d = "/Volumes/INRIA3/dippers/CINCLUS_2015_SM2_nest/20150325_SM2_cesolet/ExtractSeq_2_CESOLET_20150325_060000/"
    #d = "/Volumes/INRIA3/dippers/CINCLUS_2015_SM2_nest/20150320_SM2_chailles/ExtractSeq_2_chaille_20150320_060000/"
    d="/Volumes/VERBATIMHD/REC_Pair_nest_2016_Dippers/Nichoir_CB26_chemin_billard/Data_kit6_20160430_CB26_Chemin_billard/ExtractSeq_7_KIT6_20160501_093329/"
    callseqname = d+'calls/call_list.txt' #'calls/call_list.txt'
    extract_calls(callseqname,d,nb_calls_per_bird=-1,shuffle=True,samplerate=96000,colstart=1, colend=2,colid=3,colname=0)#, colstart=0, colend=1,colid=2,colname=3)#




if __name__=="__main__":
    test_extract_calls()