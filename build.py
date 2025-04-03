import argparse
import pandas as pd
import os
import zipfile
import tarfile
import gzip
import shutil
import sys


headers=['SNP', 'A1', 'A2', 'Z']
parameters=['out','N','snp','a1','a2','p','frq','signedSumstats']

class InvalidSumstats(Exception):
    def __init__(self,columns,message="sumstats file doesn't have all the required headers"):
        missing=[]
        columns=columns[0].split()
        for header in headers:
            if header not in columns:
                missing.append(header)
        super().__init__(f"{message}, missing: {','.join(missing)}")

class MultipleFileAfterUncompression(Exception):
    def __init__(self,com_lis,message="Multiple files in inflated after uncompression please run the script again with specified file"):
        super().__init__(f"{message}, files: {com_lis}")


class MissingColumn(Exception):
    def __init__(self,columns,message="File doesn't have all the required columns"):
        super().__init__(f"{message}, missing: {','.join(columns)}")


class IncompleteParameterPass(Exception):
    def __init__(self,args,message="Not all required parameters have been passed"):
        missing=[]
        for parameter in parameters:
            if getattr(args,parameter)==None:
                missing.append(parameter)
        super().__init__(f"{message}, missing: {','.join(missing)}")

class MultipleStudiesDetected(Exception):
    def __init__(self,para,message="Multiple studies detected"):
        super().__init__(f"{message}, parameter {para[0]} detected multiple times at row {para[1]}, \
 set SupMultiStudy to 1 to suppress this error catch,SupMultiStudy=1")


def count_word(string, word):
    return string.lower().count(word.lower())

def chechMultiStudy(file_path):  
    names=[getattr(args,para) for para in parameters[2:]]
    with open(file_path,'r') as f:
        for line,row in enumerate(f.readlines(),start=1):
            for name in names:
                if(count_word(row, name)>1):
                    raise MultipleStudiesDetected((name,line))
                    
                    


def uncompress_file(file_path, output_dir='.'):
    os.makedirs(output_dir, exist_ok=True)
    uncompressed_files = []

    if file_path.endswith('.zip'):
        print(f'Extracting ZIP file: {file_path}')
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(output_dir)
            uncompressed_files = [os.path.join(output_dir, name) for name in zip_ref.namelist()]

    elif file_path.endswith(('.tar.gz', '.tgz', '.tar')):
        print(f'Extracting TAR file: {file_path}')
        mode = 'r:gz' if file_path.endswith(('.tar.gz', '.tgz')) else 'r:'
        with tarfile.open(file_path, mode) as tar_ref:
            tar_ref.extractall(output_dir)
            uncompressed_files = [os.path.join(output_dir, member.name) for member in tar_ref.getmembers() if member.isfile()]

    elif file_path.endswith('.gz'):
        print(f'Extracting GZ file: {file_path}')
        output_file = os.path.join(output_dir, os.path.basename(file_path)[:-3])
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        uncompressed_files = [output_file]

    else:
        raise ValueError(f'Unsupported file format: {file_path}')

    uncompressed_files = [f for f in uncompressed_files if not os.path.basename(f).startswith('._')]
    if len(uncompressed_files)>1:
        raise MultipleFileAfterUncompression(uncompressed_files)
    print(f'Extraction complete: {output_dir}')
    return uncompressed_files

def drop_irr(path,extension,header,output='temp'):
    os.makedirs(output, exist_ok=True)

    if extension=='csv':
        df=pd.read_csv(path)
    elif extension=='tsv':
        df=pd.read_csv(path,sep='\t')
    drop_list=list(set(df.columns)-set(header))
    df.drop(columns=drop_list,inplace=True)
    change_dic={}
    for para in parameters[2:]:
        change_dic.update({getattr(args,para):para})
    missing_cols = [col for col in change_dic.keys() if col not in df.columns]
    if missing_cols:
        raise MissingColumn(missing_cols)
    df.rename(columns=change_dic,inplace=True)
    path_output=os.path.join(output,f"{os.path.basename(path).split('.')[0]}.txt")
    df.to_csv(path_output,sep=' ',index=False)
    return path_output
    
    
def check_parameter():
    header_name=[]
    for parameter in parameters:
        if getattr(args,parameter)==None:
            raise IncompleteParameterPass(args)
        else:
            header_name.append(getattr(args,parameter))
    return header_name


parser = argparse.ArgumentParser(description="A script to preprocess file for ldsc file formatting \
                                support file formats are csv,tsv and there compressed versions")
parser.add_argument('--file', type=str, help='input file',required=True)
parser.add_argument('--out', type=str, help='output file')
parser.add_argument('--snp', type=str, help='file path')
parser.add_argument('--a1', type=str, help='Name of A1 column')
parser.add_argument('--a2', type=str, help='Name of A2 column')
parser.add_argument('--p', type=str, help='Name of p-value column')
parser.add_argument('--frq', type=str, help='Name of FRQ or MAF column')
parser.add_argument('--signedSumstats', type=str, help='Name of signed sumstat column, comma null value (e.g.,\
                        Z,0 or OR,1)')
parser.add_argument('--N', type=str, help='number of samples')
parser.add_argument('--nullV', type=str, help='null value of signed sumstat selected')
parser.add_argument('--Rscript', type=str, help='Rscript for twas analysis')
parser.add_argument('--sumstats', type=str, help='sumstats file')
parser.add_argument('--weights', type=str, help='weight pos file')
parser.add_argument('--weights_dir', type=str, help='weights_dir directory')
parser.add_argument('--ref_ld_chr', type=str, help='ref_ld_chr')
parser.add_argument('--chr', type=str, help='chromosome number if not passed, will run for all autosomal chromosomes')
parser.add_argument('--outF', type=str, help='output file')
parser.add_argument('--SupMultiStudy', type=int, help='suppress multi study check,pass 1 to suppress',default=0)
parser.add_argument('--datC', type=int,default=None ,help="internal parameter to convert dat to csv don't pass \
                    anything")
parser.add_argument('--outFolder', type=str, help='output folder')


args = parser.parse_args()
path=args.file

if(args.datC):
    df=pd.read_csv(path,sep='\t')
    print(f"output file:{os.path.basename(path).split('.')[0]}.csv")
    basePath=os.path.dirname(os.path.dirname(path))
    output_dir = os.path.join(basePath if basePath!='' else '..', "chromosomes")
    output_file = f"{os.path.basename(path).split('.')[0]}.csv"
    df.to_csv(os.path.join(output_dir, output_file))
    sys.exit(0)
if(not(os.path.exists(path))):
    raise FileNotFoundError(f"file {path} doesn't exist please check and re-attempt")

file=os.path.basename(path)
extension=file.split('.')[-1]

if extension=='sumstats':
    df=pd.read_csv(path,sep='\t')
    count=0
    for header in df.columns:
        if header in  headers:
            count+=1
    if count!=4:
        raise InvalidSumstats(df.columns)
    os.makedirs('temp', exist_ok=True)  
    shutil.copy(path,os.path.join('temp',file))
    makefile_content = f"""SHELL := /bin/zsh\nCONDA_PREFIX = $(shell conda info --base)\nbase ?= /Users/shehzailabbas/Desktop/fusion\n.SILENT:\n.ONESHELL:\n.EXPORT_ALL_VARIABLES:\nrun:\n\tsource $(CONDA_PREFIX)/etc/profile.d/conda.sh && conda activate ldsc && Rscript {args.Rscript} --sumstats {file} --weights {args.weights} --weights_dir {args.weights_dir} --ref_ld_chr {args.ref_ld_chr} --chr {args.chr} --out ../{args.outF} \n\tpython $(base)/build.py --file ../{args.outF} --datC 1"""

    with open('temp/Makefile', 'w') as f:
        f.write(makefile_content)
    sys.exit(0)

header_name=check_parameter()

if not(extension in ["csv",'tsv']):
    path=uncompress_file(path,output_dir="temp")[0]
    extension=os.path.basename(path).split('.')[-1]

if not(args.SupMultiStudy):
    chechMultiStudy(path)
out_path=os.path.basename(drop_irr(path,extension,header_name))


makefile_content = f"""
SHELL := /bin/zsh
CONDA_PREFIX = $(shell conda info --base)
base ?= /Users/shehzailabbas/Desktop/fusion
.SILENT:
.ONESHELL:
.EXPORT_ALL_VARIABLES:
ldsc:
	source $(CONDA_PREFIX)/etc/profile.d/conda.sh && conda activate ldsc && $(base)/ldsc/munge_sumstats.py \
--sumstats {out_path} \
--out output \
--snp snp \
--a1 a1 \
--a2 a2 \
--p p \
--frq frq \
--signed-sumstats signedSumstats,{args.nullV} \
--N {args.N} && gunzip output.sumstats.gz && Rscript {args.Rscript} \
--sumstats {args.sumstats} \
--weights {args.weights} \
--weights_dir {args.weights_dir} \
--ref_ld_chr {args.ref_ld_chr} \
--chr {args.chr} \
--out ../{args.outF} \n\tpython $(base)/build.py --file ../{args.outF} --datC 1
"""

with open('temp/Makefile', 'w') as f:
    f.write(makefile_content)

    