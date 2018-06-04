from numpy import *
import cPickle
try:
    import sndfile
except:
    print "pb with sndfile"
    pass
import copy
import sys
import sys
sys.path.append("../")
from zf import zf_chunk,zf_record,zf_call
from utils.seq_analysis import get_projname,get_callseq


class wt_chunk(zf_chunk):
    # class to handle energy of calls when extracted
    def __init__(self,sfile,f_size,overlap=1024):
        super(wt_chunk,self).__init__(sfile,f_size,overlap) #init with zf_chunk class __init__

    def get_env(self,f_fftw,f_div,Hz_LowPass,dB,method=2):
        ## first call after the base get_env function 
        zf_chunk.get_env(self,f_fftw,f_div,Hz_LowPass,dB,method)
        ## for these calls we need to compute its energy
        self.e_sounds=[]
        for call in self.f_sounds:
            start=call[0]/f_div
            end=call[1]/f_div
            if call[1]==-1: ## call interrupted we compute the energy till the end
                end=-1
                ev=sum(self.env[start-1:self.f_size+1]) # f_size+1 to reach the end of the file / start=1
            else :
                ev=sum(self.env[start-1:end+1])    
            
            self.e_sounds.append([call[0],call[1],ev])
        return self.env
    
class wt_record(zf_record):
    ## class of record to handle energy 
    def __init__(self,fname,f_size,timemax=-1):
        super(wt_record,self).__init__(fname,f_size,timemax)
        self.e_sounds=[]
        self.last_e=0.0

    def _get_next_chunk(self):
        self.c_curr=wt_chunk(self.sfile,self.f_size)
        if self.c_curr.frames==0:
            return 0
        self.nchunk+=1
        return 1

    def _get_calls(self,f_fftw,f_div,Hz_LowPass,dB,method=2):
        self.c_env=self.c_curr.get_env(f_fftw,f_div,Hz_LowPass,dB,method)
        if self.last_time!=-1: ## call started last time : we merge the two
            if len(self.c_curr.f_sounds)>0: # if more than one sound in the 
                self.c_curr.f_sounds[0][0]=self.last_time-self.offset
                # we copy the f_sounds in e_sounds 
                self.c_curr.e_sounds[0][0]=self.last_time-self.offset
                self.c_curr.e_sounds[0][2] += self.last_e
                self.last_time=-1
                self.last_e=0.0
            else: ## empty calls the call is longer than the chunk
                print 'warning: call longer than chunk: consider increasing chunk size'
                ## but nothing to do
        for sound in self.c_curr.e_sounds:
            f_start=sound[0]
            f_end=sound[1]
            e_all=sound[2]
            if f_end!=-1:
                self.c_sounds.append([f_start+self.offset+f_fftw/2,f_end+self.offset+f_fftw/2,e_all])
            else:
                self.last_time=f_start + self.offset
                self.last_e=e_all
        self.offset+= self.f_size
        
        return self.c_sounds

    def flush(self):
        # empty the c_sounds
        # get last interrupted call
        s=self.c_sounds[-1]
        if s[1]==-1: # if the end of the sound = -1 : not finished => keep it
            self.c_sounds=[s]
        else:
            self.c_sounds=[]

class wt_waves(object):
    ## class to handle wave extraction/selection/learning
    def __init__(self,files,calls):
        # calls : list of calls (1 call = (num bird, start, end, energy))
        self.wave_files=files
        self.calls=array(calls[:])
    
    def extract(self,low_lim,cmax,outdir,keep_calls=True): #fft_w,f_div,f_pad,dB,peak_pc,
        # create files containing one call.
        # filenames = "call-n-start-bird.wav", n=number of the call (after shuffle), start=start of the call in the initial wav file, in second, bird = number of the bird (micro)
        # cmax : max nb of calls
        # outdir : directory where the created files will be stored
        self.select=copy.copy(self.calls)
        #random.shuffle(self.select)
        if cmax!=-1:self.select=self.select[:cmax,:]
        fi=open(outdir+'/call_list.txt','w')
        self.sfiles=[sndfile.open(f,'r') for f in self.wave_files]
        self.files=[]
        for n,ca in enumerate(self.select):
            bird=int(ca[0])
            f_start=int(ca[1])
            f_end=int(ca[2])
            wavfilebird=self.wave_files[bird]
            samplerate=self.sfiles[bird].info.samplerate#float(wav_sfile.info.samplerate)
            fnameinit='call-%d-%.3lf-%d'%(n,f_start*1./samplerate,bird)+'.wav'
            fi.write('%s %.3lf %.3lf %s %s\n'%(fnameinit,f_start*1./samplerate,f_end*1./samplerate,bird,wavfilebird))
            fname=outdir+'/'+fnameinit
            self.files.append(fnameinit)
            if keep_calls==True:
                sfile=self.sfiles[bird] # get the file of the bird we want to extract the call
                sfile.seek(f_start,0) # go to the start of the call in the file
                wave_init=sfile.readf_double(f_end-f_start)
                wave=concatenate(([0.]*3000,wave_init,[0.]*3000)) # add 3000 frames of silence before and after the sound for audacity ###a virer + dans le classifier
                outinfo = sndfile.SF_INFO(samplerate=samplerate, channels = 1,format = (sndfile.SF_FORMAT_WAV|sndfile.SF_FORMAT_PCM_16),sections = 1,seekable = 1)
                o_sfile = sndfile.open(fname,mode="w",info=outinfo)
                o_sfile.write_double(wave/max(abs(wave)))
                o_sfile.close()
        [f.close() for f in self.sfiles]
        fi.close()
        return self.files




class wt_parse(object):
    # class to parse multi-channels call list and decide which was the one used
    def __init__(self,call_list,waves):
        self.waves=waves[:]
        self.nchannels=len(self.waves)
        self.c_list=[list(e) for e in call_list]
        
        self.ncalls=len(self.c_list)
        self.resort()
    
    def resort(self):
        ## we need to assure its properly sorted
        ## so we resort
        self.c_list.sort(key=lambda x:x[1]) ## sort using the second column : start of calls
    
    def get_close(self,f_thr,div,chunk_size):
        # thr : threshold to determine if 2 calls are actually one call recorded by 2 microphones
        # if false, remove every overlapping call
        # self.clist[:,3]= call energy
        i=0
        
        while i<len(self.c_list) :
            sound1 = self.c_list[i]
            start1 = sound1[1]
            len1 = sound1[2]-sound1[1]
            channel1 = sound1[0]
            j=i+1
            while j<len(self.c_list) :
                sound2 = self.c_list[j]
                start2 = sound2[1]
                channel2 = sound2[0]
                
                if start2 < start1 + len1 and channel1 != channel2 :
                    len2=sound2[2]-sound2[1]
                    e1 = sound1[3]*1./len1
                    e2 = sound2[3]*1./len2
                    if e2 < e1:
                        self.c_list.remove(sound2)
                    
                    elif start2==start1 and e1<=e2 :
                        self.c_list.remove(sound1)
                        break
                    else : j+=1
                else : j+=1
    
            i+=1



class wt_decimate(object):
    ## class object to decimate overlapping calls on the same channel
    def __init__(self,calls):
        self.calls=array(calls[:])
        self.nchannels=unique(self.calls[:,0]) # list of bird numbers
    
    def decimate(self,THR):
        # THR
        self.decimated=[]
        
        for nch in self.nchannels:
            va=self.calls[self.calls[:,0]==nch,:] # calls from channel nch
            va=self.dec_end(va)
            self.decimated=append(self.decimated,va)
        self.decimated=reshape(self.decimated,(-1,4)) # don't keep the energy of the call
        return self.decimated
    
    
    def dec_end(self,va):
        # va : list of calls for one channel
        ## cut depending on last
        v0=[]
        i=0
        cpt_dec = 0
        while i<len(va):
            #print i
            #print va
            xc=va[i,:] #ith call
            j=i+1
            done=0
            while done==0 and j<len(va):
                xn=va[j,:] # jth call
                if xn[1]<xc[2]: ## end is after the preceding start
                    cpt_dec += 1
                    #print xn[1]/44100.,"<",xc[2]/44100.
                    xc[2]=max(xn[2],xc[2]) ## merge them
                    j+=1
                    i+=1
                else:
                    done=1
            v0.append(xc)
            i+=1
        print "nb decimated = ",cpt_dec
        return array(v0)


class wt_multi(object):
    # wrapper class for handling multiple waves file simultaneously
    # create a list of records
    def __init__(self,nb_birds,filenames,f_size):
        self.nb_birds=nb_birds
        self.fnames=filenames[:]
        if len(self.fnames)!=self.nb_birds:
            raise Exception("sizes should be equal")
        self.records=[]
        for f in self.fnames: # add every wt_record object to the list self.records
            new_wtr=wt_record(f,f_size)
            self.records.append(new_wtr)
        self.n=[0]*self.nb_birds
        self.s_sounds = []

    def get_next_calls(self,f_fftw,f_div,Hz_LowPass,dB):
        # fill self.s_sounds list with current chunks, sorted by start time
        self.current_chunk=[]
        ndone=0
        for i,zr in enumerate(self.records): # read each record (wt_record object)
            self.n[i]=zr.next_chunk(f_fftw,f_div,Hz_LowPass,dB) # function next_chunk inherited from zf_record class

            if self.n[i]==-1: # there are other chunks
                ndone +=1

                for c in zr.c_sounds: # c_sounds contains sounds of all previous chunks
                    self.current_chunk.append([i,c[0],c[1],c[2]]) # i = number of the record (or bird)
        if ndone==0:
            ## no more chunks to parse
            return [x.nchunk for x in self.records] # number of chunks for each record
        # we get all the sounds for this chunk
        # sorts them by start time
        self.s_sounds=sorted(self.current_chunk,key=lambda x:x[1]) #x[1] = start time

        #        print self.s_sounds
        del self.current_chunk
        return -1 # done

    def get_waves(self): # return list of sounds (all records)
        return [f.c_curr.r for f in self.records]

    def save_wave(self,fname): # save wave in nb_bird different files
        for i,f in enumerate(self.records):
            f.save_wave(fname+'-%d'%i)


class wt_project(object):
    def __init__(self,wav_dir,f_fftw,f_div,Hz_LowPass,dB,chunk_size):
        # names
        self.wav_dir=wav_dir
        self.name = get_projname(wav_dir)

        import os
        self.files=[f for f in sorted(os.listdir(wav_dir)) if ((f.endswith('WAV') or f.endswith('wav')) and f.startswith('.')==False)]
        self.nbirds=len(self.files)
        
        ## get samplerate
        sfile0 = sndfile.open(self.wav_dir + '/'+self.files[0],'r')
        self.samplerate = float(sfile0.info.samplerate)
        sfile0.close()

        # spectrum parmaters
        self.f_fftw=f_fftw
        self.f_div=f_div
        self.Hz_LowPass=Hz_LowPass
        self.dB=dB
        self.chunk_size=chunk_size
        self.calls=[]
        self.numcalls=[]## call number after removing noise
        self.res_classif=[]## list of 0 or 1 if call or not
        # create multi objects

        # starting stage
        self.stage=0

    def get_sounds(self,verbose=0,nmax=-1,corr=True):
        # thr_parse : threshold for cl_parse (remove double calls). If thr=-1 : keep double calls
        wt_mul=wt_multi(self.nbirds,[self.wav_dir+ '/'+f for f in self.files],self.chunk_size)
        ttimes=[rec.frames for rec in wt_mul.records] # ? number of frame for each wav file 
        if nmax==-1:
            nmax=max(ttimes)/self.chunk_size # max number of frame for all records / chunk : init max size for chunks
        n=-1
        ntime=0 # counter
        self.chunk_sounds=[]
        
        while n==-1 and ntime<nmax: # n==-1 if there is a next chunk
            nb_sounds_init = len(wt_mul.s_sounds)
            n=wt_mul.get_next_calls(self.f_fftw,self.f_div,self.Hz_LowPass,self.dB) # fill wt_mult.sounds with all sounds from all records, sorted by start time
            nb_sounds_end = len(wt_mul.s_sounds)
            
            self.chunk_sounds.append((nb_sounds_init,nb_sounds_end))
            
            ntime +=1
            if verbose>0:
                print [len(f.c_sounds) for f in wt_mul.records]
                print round(float(ntime)/nmax*100.,2),"%"
            try:
                print 'position in the wav file:',wt_mul.s_sounds[-1][1]/self.samplerate/60 ,'min'
            except:pass
        self.s_sounds=wt_mul.s_sounds[:] 

    def get_parse_sounds(self,nmax=-1,thr_parse=-1,div_parse=0,corr=True):
        wt_mul=wt_multi(self.nbirds,[self.wav_dir + '/'+f for f in self.files],self.chunk_size)    
        ttimes=[rec.frames for rec in wt_mul.records] # ? number of frame for each wav file ?
        filename_maxcor = "./max_corr.txt"
        if nmax==-1:
            nmax=max(ttimes)/self.chunk_size # max number of frame for all records / chunk : init max size for chunks
        n=-1
        ntime=0 # counter
        sounds_parse = []
        max_cor = []
        #print len(self.chunk_sounds)
        #print nmax
        while n==-1 and ntime<nmax: # n==-1 if there is a next chunk
            #print ntime
            lim_chunksounds = self.chunk_sounds[ntime] #init call, end call for the current chunk
            #print lim_chunksounds
            sounds_chunk = copy.copy(self.s_sounds[lim_chunksounds[0]:lim_chunksounds[1]]) # get sounds of current chunk
            [f._get_next_chunk() for f in wt_mul.records]
            
            cl_s=wt_parse(sounds_chunk,wt_mul.get_waves())
            cl_s.get_close(thr_parse,div_parse,self.chunk_size)
            [sounds_parse.append(s) for s in cl_s.c_list]
            ntime +=1

        f = open(filename_maxcor,"w")
        tau = range(len(arange(0,thr_parse,div_parse)))
        max_cor.sort(key=lambda x:x[0])
        max_cor = array(max_cor)
        max_cor.flatten()
        list_tau_cor = array([x[0] for x in max_cor])
        for x in tau :
            ind = where(list_tau_cor==x)[0]
            list_corri = [max_cor[i][1] for i in ind]
            f.write(str(x)+"\t"+str(mean(list_corri))+"\n")
        f.close()
        self.s_sounds=sounds_parse[:]


    def dump(self): # save the current wt_project object
        fname=self.wav_dir+"/dump_stage-%d.dump"%self.stage
        f = open(fname,'wb')
        cPickle.dump(self.__dict__,f)
        f.close()

    def load(self): #load saved object
        fname=self.wav_dir+"/dump_stage-%d.dump"%self.stage
        f = open(fname,'rb')
        tmp_dict = cPickle.load(f)
        f.close()
        self.__dict__.update(tmp_dict) # update new attributes

    def get_calls(self,f_pad,thr,smax,verbose=0):
        # thr : threshold for extract_calls
        sfiles=[sndfile.open(self.wav_dir + '/'+f,'r') for f in self.files]
        if smax==-1 : smax=len(self.s_sounds)
        #print "start=",self.s_sounds[0][1],'end=',self.s_sounds[0][2]
        for i,s in enumerate(self.s_sounds[:smax]): ##  read every sound from s_sounds until smax
            sfile=sfiles[int(s[0])] # get file number (or bird number)
            start=s[1]
            end=s[2]
            
            call=zf_call(sfile,start,end,f_pad) # zf_call object
            call.get_env(self.f_fftw,self.f_div,self.Hz_LowPass) # compute envelope
            call.extract_call(thr)

            self.calls.append([int(s[0]),call.f_start,call.f_end,call.f_end-call.f_start])
            if verbose>0:
                print s[0],'%.4lf'%(call.f_start/float(sfile.info.samplerate)),'%.2lf'%(call.f_end/float(sfile.info.samplerate)),
                print '%.4lf'%((call.f_end-call.f_start)/float(sfile.info.samplerate)*1000)
        [f.close() for f in sfiles]
        #print "start=",self.calls[0][1],'end=',self.calls[0][2]

        return self.calls


    def time_course(self,div):
        timecourse(array(self.calls[:]),div)
        # calls=array(self.calls[:])
        # channels = set(calls[:,0])
        # nb_bird = len(channels)

        # start = calls[0][1]-44100
        # end = calls[len(calls)-1][2]+44100
        # x = arange(start,end,div)
        # y = [zeros((len(x)))+i*0.1-0.5 for i in range(nb_bird)]

        # for c in calls[:1000] :
        #     ch = int(c[0])
        #     st = int((c[1]-start)/div)
        #     dur = st+len(arange(c[1],c[2],div))
        #     for i in range(st,dur): y[ch][i]+=0.08

        # # x in seconds
        # range_x = 10
        # xsec = [xi*1./44100 for xi in x]
        # lines_plot = [zeros((len(x)))+i+1 for i in range(nb_bird)]
        # plt.plot(xsec,y[0],"r",xsec,y[1],"g",xsec,y[2],"b",xsec,y[3],"y")
        # plt.axis([start*1./44100., end*1./44100., -1, 1])
        # plt.xticks(range(0,int(xsec[len(xsec)-1]),range_x))
        # plt.title("Time course")
        # plt.xlabel("time (in sec)")
        # plt.ylabel("Birds")
        # plt.show()

#calls_test = array([[0,10.2,15.3],[1,15.8,16.2],[1,17.8,18.5],[2,25.2,26.3]])
        

    def save_callseq(self,fname,call_num):
        f=open(fname,'w')
        for i,c in enumerate(self.calls):
            f.write('%f %f %d %s\n' %((c[1]/self.samplerate,c[2]/self.samplerate,c[0],'call-'+str(call_num[i]))))
        f.close()
        
