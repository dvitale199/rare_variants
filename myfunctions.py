# myfunctions

import pandas as pd
import sys
import time
import subprocess


# Created for inquiry for Lana

def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)
    # start = time.time()
    res=subprocess.run(command.split(), stdout=subprocess.PIPE)
    # end = time.time()
    # sys.stdout.write('EXEC_TIME in sec: '+ str(round(end - start, 3)) + ' : ')
    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))


def get_refGENE(gene, annovar_refgene):
    """
    give gene name and path to annovar_refGene, then pull refGene results 
    NEED annovar loaded
    """
    refgene = pd.read_csv(annovar_refgene, delim_whitespace=True, header=None)
    # Only on ordinary chromosome
    refgeneChr = refgene[refgene[2].isin(['chr'+str(i) for i in range(1,23)])].copy()

    # Get the Gene in refGene
    refgeneChrTarget = refgeneChr[refgeneChr[12]==gene].copy()
    refgeneChrTarget.to_csv(f'{gene}.refGene', header=False)
    N = refgeneChrTarget.shape[0]
    if N ==0:
        print(f'_C_: refGene does not have {gene}')
        return('')

    # Get the table of exone, intron for the longest defined one. 
    refgeneChrTarget['Length']=refgeneChrTarget[5]-refgeneChrTarget[4]
    refgeneChrTargetSorted = refgeneChrTarget.sort_values('Length', ascending=False)
    df = refgeneChrTargetSorted.iloc[0, :].copy()

    accessions = list(refgeneChrTargetSorted[1])
    selected = accessions[0]
    print(f'_C_: selected={selected}, length={df.Length/1000}K')

    if N > 1:
        dropped = accessions[1:]
        print(f'_C_: dropped={dropped}. \n_C_: ..See {gene}.refGene in detail')



    # reshape the ref table
    # exon dataframe
    exons=[(i, j) for i, j in zip(df[9].split(','),df[10].split(',')) if not (i, j)==('', '')]
    exons = pd.DataFrame(exons, columns=['start', 'stop'])
    Nrow = exons.shape[0]
    if df[3]=='-':
        exons['type'] = ['exon' + str(Nrow - i) for i in range(Nrow)]
    if df[3]=='+':
        exons['type'] = ['exon' + str(i + 1) for i in range(Nrow)]

    # intron dataframe
    introns=[(i, j) for i, j in zip(exons.stop[:Nrow-1], exons.start[1:])]
    introns = pd.DataFrame(introns, columns=['start', 'stop'])
    Nrow = introns.shape[0]
    if df[3]=='-':
        introns['type'] = ['intron' + str(Nrow - i) for i in range(Nrow)]
    if df[3]=='+':
        introns['type'] = ['intron' + str(i + 1) for i in range(Nrow)]


    # combine exon and intron dataframe
    d = pd.concat([exons, introns], axis=0, ignore_index=True).sort_values('start')
    d['gene'] = gene
    d['accession'] = df[1]
    d['chr'] = int(df[2].replace('chr', ''))
    d['length'] = d.stop.astype('int') - d.start.astype('int')
    d = d[['gene', 'accession', 'chr', 'start', 'stop', 'type', 'length']]
    d.to_csv(f'{gene}.refgene.reshape', index=False)

# get_refGENE(gene = GENE, annovar_refgene = ANNOVAR_REFGENE)


def getGeneVariants(gene, annovar_refgene, keep_id, bpmargin=1000000):
    '''
    take gene, path to annovar_refgene, path to keep IDs, and bpmargin.
    NEEDS plink2, annovar loaded.
    RETURNS 
     GENE.log(plink log)
     GENE.raw(plink export)
     GENE.refGene (refgene extract)
     GENE.refGene.reshape (re-shaped GENE.refGene)
    '''
    get_refGENE(gene = gene, annovar_refgene=annovar_refgene)
    dft = pd.read_csv(f'{gene}.refgene.reshape')
    chr = 'chr'+str(dft.chr[0])
    cut_start = min(dft.start - bpmargin)
    cut_stop = max(dft.stop + bpmargin)
    print(f'_C_: Cut {chr}:{cut_start}-{cut_stop}')

    pfile = f'/data/CARD/PD/GENOMES/august19/genotypes/pd.june2019.{chr}.sqc'
    plinkcmd = f'\
    plink\
     --pfile {pfile}\
     --geno 0.05\
     --maf 0.005\
     --keep {keep_id}\
     --chr {dft.chr[0]}\
     --from-bp {cut_start}\
     --to-bp {cut_stop}\
     --export A\
     --make-pgen\
     --out {gene}'

    shell_do(plinkcmd)

def getGeneExpression(gene, ensg, keep_id, log2expression):
	'''
	Input: gene name, ensg id, keep id, path to log2 expression of RNA
	Output: gene.log2expression
	'''
	IIDs = pd.read_csv(keep_id, sep='\t').IID
	headers = shell_do(f'head -n1 {log2expression}', 
		return_log=True).split(',')
	expvalues = shell_do(f'grep {ensg} {log2expression}', 
		return_log=True).split(',')
	df = pd.DataFrame(data={headers[0]:headers[1:], expvalues[0]:expvalues[1:]})
	dfkeep = df.loc[df.ID.isin(IIDs), :]
	dfkeep.to_csv(f'{gene}.rnaLog2')

def getMethylationBeta(keep_id, queryFile, MethylBetaFile, SampleFile):
	'''
	Input: queryFile created by "getCloseGeneSite.ipynb", MethylBetaFile, SampleFile
	Output: gene_cgID.csv for each cgID
	'''
	IIDs = pd.read_csv(keep_id, sep='\t').IID
	# Sample file to join
	sample = pd.read_csv(SampleFile, sep='\t')
	sample['IDkey'] = sample['Sentrix ID'].astype('str') + '_' + sample['Sentrix Position']
	sample['IID'] = 'PPMISI' + sample.PATNO.replace('A', '') # 3794A


	# Get header of Methylation beta file
	## Need to change the last header
	headers = shell_do(f'head -n1 {MethylBetaFile}', return_log=True).split(',')
	headers = [i.replace('\r\n', '') for i in headers ]

	queryFile = pd.read_csv(queryFile)

	for cgID, gene in zip(queryFile.Name, queryFile.UCSC_RefGene_Name):
	    betas = shell_do(f'grep -m1 {cgID} {MethylBetaFile}', return_log=True).split(',')
	    
	    if len(betas)==1:
	        print(f'{cgID} not found')
	        continue
	    
	    dfbeta = pd.DataFrame(data = {'IDkey':headers, 'beta':betas[1:]})
	    df = pd.merge(sample, dfbeta, on='IDkey', how='inner')
	    df.loc[df.IID.isin(IIDs), ['IID', 'beta']].to_csv(f'{gene}_{cgID}.methyl', index=False)


# Created for Inquiry from Alice

def pullByRSID(rsid, chrnum, genomefolder, idfile, outfile):
    '''
    INPUT: rsid, chrnum (number), genomefolder, idfile
    OUTPUT: outfile-contains the rsid,position,ref,alt
    '''
    fileToSearch = f'{genomefolder}/annotation.{chrnum}.hg38_multianno.txt'
    with open(f'{outfile}', 'a') as f:
        script = f'grep -w -m1 {rsid} {fileToSearch}'
        grepRes = shell_do(script, return_log=True).split('\n') # Can be more than 1 (multi-allelic)
        N_alt = len(grepRes)-1 # Last one is ''
        if N_alt==0:
            f.write(f'{rsid}\tNA\n')
        else:
            for i in range(N_alt):
                variantID = grepRes[i].split('\t')[21]
                bp = grepRes[i].split('\t')[20]
                f.write(f'{rsid}\t{variantID}\n')
                pfile = f'/data/CARD/PD/GENOMES/august19/genotypes/pd.june2019.chr{chrnum}.sqc'
                plinkcmd = f'plink\
                --pfile {pfile}\
                --chr {chrnum}\
                --keep {idfile}\
                --from-bp {bp}\
                --to-bp {bp}\
                --export A\
                --make-pgen\
                --out pulled_{rsid}'
                shell_do(plinkcmd)

def pullByRSIDfile(rsidfile, genomefolder, idfile, outfile):
    df = pd.read_csv(rsidfile)
    for x,y in zip(df.rsid, df.chr):
        pullByRSID(rsid = x, chrnum = y, 
                   genomefolder=genomefolder, 
                   idfile=idfile,
                   outfile=outfile)
