import cPickle
from sgRNA_learning import *


# # genome fasta file
GENOME_FASTA =  '/gstore/data/omni/crispr/piru/CRISPRai/input_files/lifted_hg38/hg38.fa'

# # paths to genome and bigWig files containing chromatin data of interest
# # see https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeOpenChromDnase    
CHROMATIN_DNASE = '/gstore/data/omni/crispr/piru/CRISPRai/input_files/lifted_hg38/wgEncodeOpenChromDnaseK562BaseOverlapSignalV2_lifted_hg38.bigWig'
CHROMATIN_FAIRE = '/gstore/data/omni/crispr/piru/CRISPRai/input_files/lifted_hg38/wgEncodeOpenChromFaireK562Sig_lifted_hg38.bigWig'
CHROMATIN_MNASE = '/gstore/data/omni/crispr/piru/CRISPRai/input_files/lifted_hg38/wgEncodeSydhNsomeK562Sig_lifted_hg38.bigWig'

### debug with manually converted files ###
LIBRARY_TABLE_TRAINING = '/Users/perampap/Documents/For publishing/sonata_trainingdata_libraryTable_hg38.txt'
SGRNA_TABLE_TRAINING = '/Users/perampap/Documents/For publishing/sonata_example/sonata_trainingdata_sgRNAInfoTable_hg38.txt'
TSS_DATA = '/Users/perampap/Documents/For publishing/sonata_example/sonata_human_tssTable_hg38.txt'
P1P2_DATA = '/Users/perampap/Documents/For publishing/sonata_example/sonata_p1p2Table_hg38.txt'  

modality = 'CRISPRa'

def predictWeissmanScore(tssTable, p1p2Table, sgrnaTable, libraryTable, modality, verbose):

    # trained model pickle file
    pickleFile = 'trained_models/' + modality + '_estimator_weissman_hg19.pkl'
   
    
    # open pickle file to continue from previously trained session/model
    try:
        with open(pickleFile) as infile:
            fitTable, estimators, scaler, reg, transformedParams_train_header = cPickle.load(infile)
    except:
        raise Exception('Trained model file not found.') 

    paramTable = getParamTable(tssTable, p1p2Table, sgrnaTable, libraryTable, fastaFile=GENOME_FASTA, CHROMATIN_DNASE, CHROMATIN_FAIRE, CHROMATIN_MNASE, verbose = verbose)
    
    transformedParams_new = getTransformedParams(paramTable, fitTable, estimators, verbose = verbose)

    print 'Predicting sgRNA scores...'
    predictedScores = pd.Series(reg.predict(scaler.transform(transformedParams_new.loc[:, transformedParams_train_header.columns].fillna(0).values)), index=transformedParams_new.index)
    
    return predictedScores


def getParamTable(tssTable, p1p2Table, sgrnaTable, libraryTable, fastaFile, CHROMATIN_DNASE, CHROMATIN_FAIRE, CHROMATIN_MNASE, verbose):

    
    try:
        genomeDict=loadGenomeAsDict(fastaFile)
    except:
        raise Exception("Genome FASTA file not found.")

    if verbose == True:
        print "Loading chromatin data..."

    try:
        bwhandleDict = {'dnase':BigWigFile(open(CHROMATIN_DNASE)), 'faire':BigWigFile(open(CHROMATIN_FAIRE)), 'mnase':BigWigFile(open(CHROMATIN_MNASE))}
    except:
        raise Exception("Could not load chromatin data.")

    # parse primary TSS and secondary TSS
    p1p2Table['primary TSS'] = p1p2Table['primary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]), int(tupString.strip('()').split(', ')[1].split('.')[0])))
    p1p2Table['secondary TSS'] = p1p2Table['secondary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]),int(tupString.strip('()').split(', ')[1].split('.')[0])))

    if verbose == True:
        print "Calculating parameters..."

    #try:
    #    paramTable = generateTypicalParamTable(libraryTable, sgrnaTable, tssTable, p1p2Table, genomeDict, bwhandleDict)
    #except:
    #    raise Exception("Error generating parameter table.")

    #paramTable = generateTypicalParamTable(libraryTable, sgrnaTable, tssTable, p1p2Table, genomeDict, bwhandleDict)

    print "library table ....."
    print libraryTable.head()
    print "sgRNA info table ....."
    print sgrnaTable.head()
    print "TSS table ....."
    print tssTable.head()
    print "p1p2 table ....."
    print p1p2Table.head()
    
    print "trying param table...."
    
    paramTable = generateTypicalParamTable(libraryTable = libraryTable, sgInfoTable = sgrnaTable, tssTable = tssTable, p1p2Table = p1p2Table, genomeDict = genomeDict, bwFileHandleDict = bwhandleDict, transcripts=False)
    
    
    return paramTable


def getTransformedParams(paramTable, fitTable, estimators, verbose):
    
    if verbose == True:
        print 'Transforming parameters...'

    try:
        transformedParams_new = transformParams(paramTable, fitTable, estimators)
    except:
        raise Exception("Error transforming parameters.")

    # reconcile differences in column headers
    colTups = []
    for (l1, l2), col in transformedParams_new.iteritems():
        colTups.append((l1,str(l2)))
    transformedParams_new.columns = pd.MultiIndex.from_tuples(colTups)

    return transformedParams_new


tssTable = TSS_DATA
p1p2Table = P1P2_DATA
sgrnaTable = SGRNA_TABLE_TRAINING
libraryTable = LIBRARY_TABLE_TRAINING

results =  predictWeissmanScore(tssTable, p1p2Table, sgrnaTable, libraryTable, verbose)
print "******DONE******"