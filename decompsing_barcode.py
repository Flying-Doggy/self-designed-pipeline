# decompsing_barcode

import collections
import multiprocessing
import argparse
import sys
import itertools
import gzip
import time
import os

class barcode_pools:
    def __init__( self, barcodes:list = [] , bar_len:list = [] ,  linkers:list = [] ):
        self.bar = barcodes
        self.links = linkers
        self.bar_l = bar_len
        pass

    def get_offset( self ):
        self.off_idx = [ 0 ]        # bar1_start , bar1_end/linker1_start , linker1_end/bar2_start ...
        for i in range( len(self.bar_l) ):
            self.off_idx.append( self.off_idx[-1] + self.bar_l[i] )
            if i < len(self.links):
                self.off_idx.append( self.off_idx[-1] + len(self.links[i])  )
        pass

    def read_match( self , query:str ):
        bar_res = []
        linker_res = []
        for i in range( len(self.off_idx) - 1 ):
            tmp_idx = i // 2
            tmp_part = query[ self.off_idx[i]:self.off_idx[i+1] ]
            if i%2 == 0:   # map to barcode
                for tmp_id,tmp_seq in self.bar[tmp_idx].items():
                    if mismatch( tmp_seq , tmp_part ):
                        bar_res.append( tmp_id )
                        break
                else:
                    bar_res.append( 'NO_BAR_MATCH' )
            else:
                linker_res.append( mismatch( self.links[tmp_idx] , tmp_part ))
        return bar_res, linker_res

    def __repr__( self ):
        return ", ".join( ['barcode%d-%dreads'%( i ,len(self.bar[i])) for i in range(len(self.bar))]) + '\n' + \
                ', '.join( [ 'linker%d:%s'%(i,self.links[i]) for i in range(len(self.links))])

class fastq:
    def __init__( self , ID:str , r1:str = '' , r2:str = '' , a1:str= '' , a2:str='', q1:str='' , q2:str='' ):
        # r- read , a- annotation , q- quality
        self.id = ID
        self.r1 , self.r2 = r1, r2
        self.a1 , self.a2 = a1, a2
        self.q1 , self.q2 = q1, q2
        pass

    def add_read( self, r:int , infos:list ):
        if r == 1:
            self.r1 , self.a1, self.q1 = infos
        else:
            self.r2 , self.a2, self.q2 = infos
        pass

    def cut_fastq( self , r:int , cut_point:int ):
        if r == 1:
            self.r1 = self.r1[cut_point:]
            self.q1 = self.q1[cut_point:]
        else:
            self.r2 = self.r2[cut_point:]
            self.q2 = self.q2[cut_point:]
        pass

    def print_seq( self , out = None):
        if out:
            print( '@%s/1\n%s\n%s\n%s'%(self.id, self.r1 , self.a1 , self.q1) ,file= out[0] )
            print( '@%s/2\n%s\n%s\n%s'%(self.id, self.r2 , self.a2 , self.q2) ,file= out[1] )
        else:
            print( '@%s/1\n%s\n%s\n%s'%(self.id, self.r1 , self.a1 , self.q1) )
            print( '@%s/2\n%s\n%s\n%s'%(self.id, self.r2 , self.a2 , self.q2) )
        pass

    def __repr__( self ):
        return "ID:%s \n read1:%s \n read2:%s"%( self.id, self.r1 , self.r2)

def is_fq_structured( barcode_pool , fq_object , r:int = 1 ):
    if r==1:
        bar_flag, link_flag = barcode_pool.read_match( fq_object.r1  )
    else:
        bar_flag, link_flag = barcode_pool.read_match( fq_object.r2 )
    return ( all(link_flag) and not any([i.startswith('N') for i in bar_flag]) ), bar_flag , link_flag

def multi_identify_barcodes( fq_list:list , bar_pool ,
                            out_seq_prefix:str='out_seq' , out_barcode:str='out_barcode',
                           threads:int = 1 , rm_flag:bool=True,
                           r:int = 1):
    
    chunk_size = (len(fq_list) // threads )+ 1

    time_p_start = time.time()
    process_list = []
    for i in range( threads ):
        tmp_p = multiprocessing.Process( target=process_chunk , args=( fq_list[i*chunk_size:(i+1)*chunk_size],
                                                                    bar_pool ,
                                                                    out_seq_prefix+'_'+str(i), out_barcode+'_'+str(i),
                                                                    r) )
        tmp_p.start()
        process_list.append( tmp_p )

    time_p_end = time.time()
    print(  'spend %f seconds to finish assigning processes '%( time_p_end - time_p_start) )


    for tmp_p in process_list:
        tmp_p.join()
    time_run_end = time.time()
    print(  'spend %f seconds to finish running processes '%( time_run_end - time_p_end ) )

    os.system('cat %s_*_R1.fq | pigz > %s_R1.fq.gz'%( out_seq_prefix , out_seq_prefix) )
    os.system('cat %s_*_R2.fq | pigz > %s_R2.fq.gz'%( out_seq_prefix , out_seq_prefix) )
    os.system('cat %s_*  > %s '%( out_barcode , out_barcode) )
    time_merge_end = time.time()
    print(  'spend %f seconds to finish merging files '%( time_merge_end - time_run_end ) )
    
    if rm_flag:
        os.system( 'rm -rf %s_*_R*.fq   %s_*'%( out_seq_prefix, out_barcode ))

def process_chunk( fq_list:list , barcode_pool , out_seq_prefix:str , out_barcode_name:str, r:int=1):
    out1 , out2 = open( out_seq_prefix+'_R1.fq' , 'w' ), open( out_seq_prefix+'_R2.fq' , 'w' )
    bar_out = open( out_barcode_name , 'w' )
    for tmp_fq in fq_list:
        tmp_flag, tmp_bar , tmp_link = is_fq_structured( barcode_pool= barcode_pool , fq_object= tmp_fq , r=r )
        if tmp_flag:
            tmp_fq.cut_fastq( r=r , cut_point= barcode_pool.off_idx[-1] )
            tmp_fq.print_seq( out = [out1 , out2 ] )
        print( tmp_fq.id + '\t' + ','.join(tmp_bar) + '\t' + ','.join( [ str(i) for i in tmp_link]) , file= bar_out )

    out1.close()
    out2.close()
    bar_out.close()
    return None

def read_fq_files( fq1:str , fq2:str ) -> dict:
    fq_dic = {}
    for tmp_file in [fq1 , fq2]:
        if tmp_file.endswith( '.gz'):
            with gzip.open( tmp_file , 'rt' ) as f:
                datas = f.readlines()
        else:
            with open( tmp_file ) as f:
                datas = f.readlines()

        for i in range( 0,len(datas),4):
            tmp_id = datas[i]
            tmp_id , tmp_r = tmp_id[1:-1].split('/')
            if fq_dic.get( tmp_id ) == None:
                fq_dic[tmp_id] = fastq( ID= tmp_id )
            fq_dic[tmp_id].add_read( int(tmp_r), [ j.strip() for j in datas[i+1:i+4] ] )
    
    return fq_dic

def mismatch( s1:str , s2:str , mis_lim:int=1 ) -> bool:
    if len(s1) != len(s2):
        return False
    mis_cnt = 0
    for i in range(len(s1)):
        mis_cnt += (s1[i]!=s2[i])
    return mis_cnt <= mis_lim

def read_barcode_file( file_name:str , linkers:list = []):
    with open(file=file_name) as f:
        datas = f.readlines()
        barcodes = [ ]
        barcode_len = [ ]
        for i in datas:
            if i.startswith('#'):
                barcodes.append( {} )
                barcode_len.append(0)
            else:
                tmp_id,tmp_seq = i.strip().split()
                barcodes[-1][tmp_id] = tmp_seq
                barcode_len[-1] = len(tmp_seq)
    tmp_barcode_pool = barcode_pools( barcodes=barcodes , bar_len=barcode_len, linkers=linkers )
    return tmp_barcode_pool


def get_args( ):
    parser = argparse.ArgumentParser( description="get the arguments of analyzing fq datas to barcode information")
    parser.add_argument( '-r1', '--read1' , type=str , help='path of read1 file' ,  dest='r1' , required=True )
    parser.add_argument( '-r2', '--read2' , type=str , help='path of read2 file' ,  dest='r2' , required=True )
    parser.add_argument( '-b', '--barcode' , type=str , help='path of barcode file' ,  dest='bc' , required=True )
    parser.add_argument( '-l', '--linker' , type=str , help='provide the sequence of linkers, seperated by comma' , dest='l' , required=True )
    parser.add_argument( '-s', '--selected' , type=int , help="provide the structured read number 1/2" ,dest='s' , default=1 )
    parser.add_argument( '-t' , '--threads' , type=int , help='the thread number of running decompsing function', dest='t' , default=1)
    parser.add_argument( '-o' , '--output' , type=str , help='prefix of output' , dest='o' , default='out')
    return parser.parse_args()


def main():
    args = get_args()

    read_bar_s = time.time()
    tmp_bar_pool = read_barcode_file( args.bc  , linkers= args.l.split(',') )
    tmp_bar_pool.get_offset()
    read_bar_e = time.time()
    print( 'spend %f seconds to finish reading barcode_pools '%( read_bar_e- read_bar_s))

    read_fq_s = time.time()
    tmp_fq = read_fq_files( args.r1, args.r2)
    read_fq_e = time.time()
    print( 'spend %f seconds to finish reading fastq files'%( read_fq_e- read_fq_s))


    multi_identify_barcodes( list(tmp_fq.values()) , tmp_bar_pool , args.o+'_trimmed'  , args.o+'_bar' ,  args.t , r=args.s )

    return None

main()