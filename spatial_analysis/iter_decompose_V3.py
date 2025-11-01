# I want to write a script similar to BarcodeIdentification.jar to analyze the barcode structure in samll memory
# The input and output files is similar to BarcodeIdentification.jar

# compared to the V2 script, the difference is the output of barcode format for each read
# {barcode}_{id} → {id}:CB:{barcode}
# barcode: L1 + B1 + L2 + B2 + .... -> B1 + B2 , skip the linkage of linker.

from functools import cache
import collections
import numpy as np
import time
import gzip
import argparse

class POOL:
    def __init__( self , candidates:dict = None , size:int = 0 ):
        self.barcodes = candidates
        self.l = size
        pass

    @cache
    def is_match( self, query:str = '' ):
        if len(query) != self.l:
            return 'NOT_MATCH'
        for tmp_id, tmp_seq_tup in self.barcodes.items():
            tmp_seq, tmp_mis = tmp_seq_tup
            if mismatch( query ,tmp_seq , tmp_mis):
                return tmp_id
        return 'NOT_FOUND'

class Barcode_Structure:
    def __init__( self , barcode_dic:dict = None , spacer:dict = None , pattern:dict = None ):
        self.bar = barcode_dic
        self.spa = spacer
        self.pat = pattern
        # merge the barcode and linkers information
        self.all_pat = dict( collections.ChainMap( self.spa , self.bar))    
        pass

    def get_pat_off_idx( self ):
        def default_pat_off_idx( ):
            return [ (0,0,0,0) ]
        self.pat_offset = collections.defaultdict( default_pat_off_idx )
        for tmp_pat_ID, tmp_order in self.pat.items():
            tmp_list = [] 
            prev_idx = 0                      
            for tmp_pat in  tmp_order:
                if isinstance(self.all_pat[ tmp_pat ] , int ):    # # judging the int spacer
                    tmp_list.append( (prev_idx , prev_idx + self.all_pat[ tmp_pat ] , self.all_pat[ tmp_pat ] , tmp_pat ) )
                elif isinstance(self.all_pat[ tmp_pat ] , str ):
                    tmp_list.append( (prev_idx , prev_idx + len( self.all_pat[ tmp_pat ] )  , self.all_pat[ tmp_pat ] , tmp_pat ) )
                else:
                    tmp_list.append( (prev_idx , prev_idx +self.all_pat[tmp_pat].l   , self.all_pat[ tmp_pat ] , tmp_pat) ) 
                prev_idx = tmp_list[-1][1]
            self.pat_offset[tmp_pat_ID] = tmp_list
        return None

    def is_match( self, query:str = '' , pat_mode:str = 'READ1'):
        match_res = []
        match_flag = 1
        for off_s,off_e,tmp_pat,tmp_pat_ID in self.pat_offset[pat_mode]:
            if isinstance( tmp_pat  , int ):    # # judging the int spacer
                continue
            elif isinstance( tmp_pat  , str ):    # judging the string linker
                # if mismatch(tmp_pat ,query[ off_s: off_e ]) <= 1:    
                #     match_res.append( tmp_pat_ID  )
                # else:
                #     match_flag = 0
                #     match_res.append( f'NOT_{tmp_pat_ID}'  )
                if mismatch(tmp_pat ,query[ off_s: off_e ]) > 1:
                    match_flag = 0
            else:                                           # # judging the barcode pool
                match_res.append( tmp_pat.is_match( query[off_s:off_e] ) )
                if match_res[-1] == 'NOT_FOUND':        # the barcode structure is not complete
                    match_flag = 0
        return match_res, match_flag


def parse_config_file(file_path: str) -> Barcode_Structure:
    """
    This function reads a configuration file and constructs a Barcode_Structure object.
    The configuration file contains information about barcode pools, spacers, and read patterns.

    Parameters:
    file_path (str): The path to the configuration file.

    Returns:
    Barcode_Structure: An instance of the Barcode_Structure class, containing the parsed information.
    """
    read_pattern = {}
    barcode_pool = collections.defaultdict(list)
    spacer = {}

    with open(file_path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('\n'):
                continue
            else:
                line = line.strip()
                if line.startswith('READ'):  # Read Pattern
                    tmp_r, tmp_pat = line.split('=')
                    read_pattern[tmp_r.strip()] = tmp_pat.strip().split('|')
                elif '=' in line:  # Other Structure
                    tmp_id, tmp_pat = line.split('=')
                    try:
                        spacer[tmp_id.strip()] = int(tmp_pat)
                    except:
                        spacer[tmp_id.strip()] = tmp_pat.strip()
                else:  # barcode pool
                    tmp_cat, tmp_name, tmp_seq, tmp_mis = line.strip().split()
                    barcode_pool[tmp_cat].append((tmp_name, tmp_seq, int(tmp_mis)))

    for tmp_cat in list(barcode_pool.keys()):
        barcode_pool[tmp_cat] = POOL(candidates={i[0]: i[1:] for i in barcode_pool[tmp_cat]},
                                     size=len(barcode_pool[tmp_cat][0][1]))
    res_bs = Barcode_Structure(barcode_dic=barcode_pool, pattern=read_pattern, spacer=spacer)
    res_bs.get_pat_off_idx()
    return res_bs

def mismatch(s1: str, s2: str, mis_lim: int = 1) -> bool:
    # 使用 numpy 数组比较字符串加速分析
    arr1 = np.array(list(s1))
    arr2 = np.array(list(s2))
    return len(arr1) == len(arr2) and np.sum(arr1!= arr2) <= mis_lim

def get_identifier( match_res:list = [] ):
    if len(match_res) == 0:
        return ''
    else:
        return '-'.join(match_res) 
    

def process_fastq(in_fq1:str , in_fq2:str , out_fq1:str , out_fq2:str , test_bar:Barcode_Structure , bar_out:str = 'bar_out.txt' ,
                  trim_flag:bool = False , only_structured:bool = False , bar_out_flag:bool = False , identifier_flag:bool = True):  
    in1 = None  
    in2 = None  
    try:  
        if in_fq1.endswith('.gz'):  
            in1 = gzip.open(in_fq1, 'rt')  
        else:  
            in1 = open(in_fq1, 'r')  
  
        if in_fq2.endswith('.gz'):  
            in2 = gzip.open(in_fq2, 'rt')  
        else:  
            in2 = open(in_fq2, 'r')  

        if bar_out_flag:
            bar_out_file = open( bar_out  , 'w' )
  
        read_cnt = 0  
        with gzip.open(out_fq1, 'wb') as out1, gzip.open(out_fq2, 'wb') as out2:  
            start_time = time.time()
            while True:  
                ID_1 = next(in1, None)
                if ID_1 == None:  # End of the fastq file
                    break
                else:
                    ID_1 = ID_1.strip()
                read_cnt += 1  
                # 从第一个文件读取fastq单元（三行） 
                seq_1, _1, qual_1 = (  
                    next(in1).strip() for _ in range(3)  
                )  
                # 从第二个文件读取fastq单元（四行）    
                ID_2, seq_2, _2, qual_2 = (  
                    next(in2).strip() for _ in range(4)  
                )  
                    
  
                match_1, flag_1 = test_bar.is_match(seq_1, pat_mode="READ1")  
                match_2, flag_2 = test_bar.is_match(seq_2, pat_mode="READ2")  

                tmp_identifier_1 = get_identifier(match_1) 
                tmp_identifier_2 = get_identifier(match_2)

                if trim_flag:
                    if flag_1:
                        seq_1 = seq_1[ test_bar.pat_offset["READ1"][-1][1]: ]
                        qual_1 = qual_1[ test_bar.pat_offset["READ1"][-1][1]: ]
                    if flag_2:
                        seq_2 = seq_2[ test_bar.pat_offset["READ2"][-1][1]: ]
                        qual_2 = qual_2[ test_bar.pat_offset["READ2"][-1][1]: ]

                if only_structured==False or (flag_1 and flag_2):
                    if identifier_flag:
                        ID_1 = ID_1[1:]
                        ID_2 = ID_2[1:]
                        outstring1 = f'@{ID_1} CB:{tmp_identifier_1+tmp_identifier_2}\n{seq_1}\n{_1}\n{qual_1}\n'
                        outstring2 = f'@{ID_2} CB:{tmp_identifier_1+tmp_identifier_2}\n{seq_2}\n{_2}\n{qual_2}\n'
                    else:
                        outstring1 = f'{ID_1}\n{seq_1}\n{_1}\n{qual_1}\n'
                        outstring2 = f'{ID_2}\n{seq_2}\n{_2}\n{qual_2}\n'
                    # 只有在 没有structured的约束 下  或者 两条reads都符合structured结构  才输出序列
                    out1.write( outstring1.encode('utf - 8') )  
                    out2.write( outstring2.encode('utf - 8') )   
                

                if bar_out_flag:
                    bar_out_file.write( '%s\t%s\n'%( ID_1.split('/')[0][1:] , '\t'.join(match_1 + match_2) ) )                    
  
                if read_cnt % 50_000_000 == 0:  
                    print('%f seconds to identify %d reads'%( time.time()-start_time, read_cnt ))  
            
  
    finally:    # 关闭已打开的文件
        in1.close()     
        in2.close() 
        if bar_out_flag:
            bar_out_file.close()

def get_args():  
    parser = argparse.ArgumentParser(description="Get the arguments for analyzing FQ data to barcode information")  
      
    # 输入参数  
    parser.add_argument('-i1', '--input1', type=str, help='Path of read1 file', dest='i1',  required=True)  
    parser.add_argument('-i2', '--input2', type=str, help='Path of read2 file', dest='i2',  required=True)  
    parser.add_argument('-c', '--config', type=str, help='Path of configuration file', dest='c' , required=True)  
      
    # 输出参数  
    parser.add_argument('-o1', '--output1', type=str, help='Name of output read1 file', dest='o1' , required=True)  
    parser.add_argument('-o2', '--output2', type=str, help='Name of output read2 file', dest='o2' ,  required=True)  # 提供默认值或保持 required=True  
    parser.add_argument('-ob', '--out_bar', type=str, help='Name of output barcode stats file', dest='ob' ,default='tmp_bar_output.txt')  
      
    # 布尔标志  
    parser.add_argument('--enable-trim', action='store_true', dest='trim', help='Enable trimming of output')  
    parser.add_argument('--disable-trim', action='store_false', dest='trim', help='Disable trimming of output')  
    parser.set_defaults(trim=False)  
      
    parser.add_argument('--enable-filter', action='store_true',  dest='filter' , help='Enable filtering of reads')  
    parser.add_argument('--disable-filter', action='store_false', dest='filter', help='Disable filtering of reads')  
    parser.set_defaults(filter=False)  
      
    parser.add_argument('--enable-report', action='store_true', dest='report', help='Output the report of barcode identification')  
    parser.add_argument('--disable-report', action='store_false', dest='report', help='Disable the report of barcode identification')  
    parser.set_defaults(report=False)  

    parser.add_argument('--enable-fq-ID', action='store_true', dest='fq_ID', help='Add barcode ID to the fastq ID field')   
    parser.add_argument('--disable-fq-ID', action='store_false', dest='fq_ID', help='Disable barcode ID to the fastq ID field')  
    parser.set_defaults(fq_ID=True)  
      
    return parser.parse_args()


def main():
    args = get_args()
    config_file = args.c
    barcode_structure = parse_config_file(file_path=config_file)
    
    process_fastq(in_fq1=args.i1, in_fq2=args.i2, out_fq1=args.o1, out_fq2=args.o2, test_bar=barcode_structure,
                 bar_out=args.ob, trim_flag=args.trim, only_structured=args.filter, bar_out_flag=args.report , identifier_flag= args.fq_ID )
    return None

main()
