from __future__ import print_function
import os
import sys
import boto3
import argparse
import uuid
import datetime
import subprocess
from rdkit import Chem
import multiprocessing as mp
from multiprocessing import Pool,Process
from rdkit.Chem import Descriptors as desc

# establish an Amazon S3 connection; an IAM role is needed to access Amazon S3
s3=boto3.client('s3')
# change to execution dir
os.chdir("/data")

# execute the main RDKit descriptor worker
def smiles_desc(smiles_str):
    output_desc=[]
    m = Chem.MolFromSmiles(smiles_str)
    output_desc=smiles_str, \
           desc.MolWt(m), \
           desc.Ipc(m), \
           desc.TPSA(m), \
           desc.LabuteASA(m), \
           desc.NumHDonors(m), \
           desc.NumHAcceptors(m), \
           desc.MolLogP(m), \
           desc.HeavyAtomCount(m), \
           desc.NumRotatableBonds(m), \
           desc.RingCount(m), \
           desc.NumValenceElectrons(m)

    desc_shard=str(output_desc).strip('()')
    smiles_out(desc_shard,csv_header)

    return
# write out the results to the local filesystem
def smiles_out(single_des_out,csv_header):
    fd=open("%s_smiles_result.csv" %csv_header,"a")
    fd.write(single_des_out)
    fd.write('\n')
    fd.close
    
# upload to an s3 bucket
def s3_upload(csv_header):
    env_S3_OUT=os.environ['OUTPUT_SMILES_S3']
    s3.upload_file("%s_smiles_result.csv" %csv_header,"%s" %env_S3_OUT,"%s_smiles_result.csv" %csv_header)

# calculate stat function
def calc_perf(input_file,timedelta):
    calc_perf_stat=int(len(input_file)/timedelta.total_seconds())
    return calc_perf_stat

if __name__ == "__main__":
# accepted args for command line execution
    parser = argparse.ArgumentParser(description='Process SMILES')
    parser.add_argument("-i", dest="smiles", required=False, help='SMILES string')
    parser.add_argument("-f", dest="smiles_file", required=False, help="File for Batch SMILES processing")
    args = parser.parse_args()

    smiles_str=args.smiles
    smiles_files=args.smiles_file
    csv_header=str(uuid.uuid4())

# if defined on the command line
    if smiles_str:
        print (smiles_desc(smiles_str))
        smiles_out(str(smiles_desc(smiles_str)),csv_header)
        s3_upload(csv_header)

# if importing a file
# by default, parallelism expands to available cores exposed to the container
    elif smiles_files:
        infile=open(smiles_files,"r")
        pool = mp.Pool()
        start_calc=datetime.datetime.now()
        smiles_list=list(infile)
        pool.map(smiles_desc,smiles_list)
        pool.close()
        end_calc=datetime.datetime.now()
        delta_calc=end_calc-start_calc
        stat_calc=calc_perf(smiles_list,delta_calc)
        print ("number of structures/sec: %s"%stat_calc)
        s3_upload(csv_header)

# if file is defined through env
# by default, parallelism expands to available cores exposed to the container
    elif os.getenv('INPUT_SMILES'):
        env_SMILES=os.environ['INPUT_SMILES']
        infile=open(env_SMILES,"r")
        pool = mp.Pool()
        start_calc=datetime.datetime.now()
        smiles_list=list(infile)
        pool.map(smiles_desc,smiles_list)
        pool.close()
        end_calc=datetime.datetime.now()
        delta_calc=end_calc-start_calc
        stat_calc=calc_perf(smiles_list,delta_calc)
        print ("number of structures/sec: %s"%stat_calc)
        s3_upload(csv_header)

# if file is in the S3 import bucket, specify an AWS BATCH job def
# by default, parallelism expands to available cores exposed to the container
    elif os.getenv('INPUT_SMILES_S3'):
        env_S3_SMILES=os.environ['INPUT_SMILES_S3']
        s3_dn=subprocess.Popen("aws s3 cp %s /data" %env_S3_SMILES, shell=True)
        s3_dn.communicate()
        s3_file=env_S3_SMILES.split('/')[-1]
        infile=open("/data/%s" %s3_file,"r")
        pool = mp.Pool()
        start_calc=datetime.datetime.now()
        smiles_list=list(infile)
        pool.map(smiles_desc,smiles_list)
        pool.close()
        end_calc=datetime.datetime.now()
        delta_calc=end_calc-start_calc
        stat_calc=calc_perf(smiles_list,delta_calc)
        print ("number of structures/sec: %s"%stat_calc)
        s3_upload(csv_header)

    else:
        print ("No SMILES FOUND...")
        sys.exit(0)

