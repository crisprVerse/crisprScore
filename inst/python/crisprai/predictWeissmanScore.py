# python script for weissman/horlbeck algorithm for sgRNA predictions
# adapted from https://github.com/mhorlbeck/CRISPRiaDesign/blob/master/CRISPRiaDesign_example_notebook.md


import cPickle
from sgRNA_learning3 import *    # contains functions to load genome and empirical data


### Data input variables ###

genome_type = "lifted_hg38"
TRAINED_DIR = 'hg19_trained/'
PREDICTION_DIR = 'hg38_predictions/'

# file output prefix
OUTPUT_PREFIX = 'CRISPRa_July19_'

# genome fasta file
GENOME_FASTA =  'input_files/lifted_hg38/hg38.fa'

# paths to genome and bigWig files containing chromatin data of interest
# see https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeOpenChromDnase    
CHROMATIN_DNASE = 'input_files/lifted_hg38/wgEncodeOpenChromDnaseK562BaseOverlapSignalV2_lifted_hg38.bigWig'
CHROMATIN_FAIRE = 'input_files/lifted_hg38/wgEncodeOpenChromFaireK562Sig_lifted_hg38.bigWig'
CHROMATIN_MNASE = 'input_files/lifted_hg38/wgEncodeSydhNsomeK562Sig_lifted_hg38.bigWig'


## file inputs for testing ##
# use sonata_hg38 input files if sonata flag is set to True
LIBRARY_TABLE_TRAINING = 'input_files/sonata_hg38/min_input/new2_trx_sonata_libraryTable.txt'
SGRNA_TABLE_TRAINING = 'input_files/sonata_hg38/min_input/new2_trx_sonata_sgrnaInfoTable.txt'
TSS_DATA = 'input_files/sonata_hg38/min_input/new2_trx_sonata_tssTable.txt'
P1P2_DATA = 'input_files/sonata_hg38/min_input/new2_trx_sonata_p1p2Table.txt'



def predictWeissmanScore(tssTable, p1p2Table, sgrnaTable, libraryTable, modality, verbose = False):

    # trained model pickle file
    PICKLE_FILE = 'trained_models/' + modality + '_estimator_weissman_hg19.pkl'
    
    # open pickle file to continue from previously trained session/model
    with open(PICKLE_FILE) as infile:
        fitTable, estimators, scaler, reg, transformedParams_train_header = cPickle.load(infile)

    paramTable = getParamTable(tssTable, p1p2Table, sgrnaTable, libraryTable, verbose = verbose)

    transformedParams_new = getTransformedParams(paramTable, fitTable, estimators, verbose = verbose)

    print 'Predicting sgRNA scores...'
    predictedScores = pd.Series(reg.predict(scaler.transform(transformedParams_new.loc[:, transformedParams_train_header.columns].fillna(0).values)), index=transformedParams_new.index)

    saveData(scoreTable = predictedScores, verbose = verbose)
    
    return predictedScores


def getParamTable(tssTable, p1p2Table, sgrnaTable, libraryTable, verbose = False):

    # load genome files
    genomeDict=loadGenomeAsDict(GENOME_FASTA)

    if verbose == True:
        print "Loading chromatin data..."

    bwhandleDict = {'dnase':BigWigFile(open(CHROMATIN_DNASE)), 'faire':BigWigFile(open(CHROMATIN_FAIRE)), 'mnase':BigWigFile(open(CHROMATIN_MNASE))}

    # df contains both primary and secondary in different columns, so need to split into seprate dfs
    p1p2Table['primary TSS'] = p1p2Table['primary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]), int(tupString.strip('()').split(', ')[1].split('.')[0])))
    p1p2Table['secondary TSS'] = p1p2Table['secondary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]),int(tupString.strip('()').split(', ')[1].split('.')[0])))

    if verbose == True:
        print "Calculating parameters..."

    paramTable = generateTypicalParamTable(libraryTable, sgrnaTable, tssTable, p1p2Table, genomeDict, bwhandleDict)

    return paramTable


def getTransformedParams(paramTable, fitTable, estimators, verbose = False):
    
    if verbose == True:
        print 'Transform and predict scores...'

    transformedParams_new = transformParams(paramTable, fitTable, estimators)

    # reconcile differences in column headers
    colTups = []
    for (l1, l2), col in transformedParams_new.iteritems():
        colTups.append((l1,str(l2)))
    transformedParams_new.columns = pd.MultiIndex.from_tuples(colTups)

    return transformedParams_new


## save results for testing ##
def saveData(scoreTable, verbose = False):

    PREDICTED_SCORES_PATH = 'debug/Aug9_predictWeissmanScore_min_Sonata.csv'

    if verbose == True:
        print 'Saving data...'

    scoreTable.to_csv(PREDICTED_SCORES_PATH, sep='\t')


def loadTestData():

    print "Loading TSS data..."
    tssTable = pd.read_csv(TSS_DATA, sep='\t', index_col=range(2))

    print "Loading P1P2 data..."
    p1p2Table = pd.read_csv(P1P2_DATA, sep='\t', header=0, index_col=range(2))

    print "Loading sgRNA info data..."
    sgInfoTable = pd.read_csv(SGRNA_TABLE_TRAINING, sep='\t', index_col=0)

    print "Loading libraryTable data..."
    libraryTable = pd.read_csv(LIBRARY_TABLE_TRAINING, sep='\t', index_col = 0)

    return tssTable, p1p2Table, sgInfoTable, libraryTable


### testing

def testfunc(passedval):
    print "The passed value is .. " + str(passedval)

def testfuncpd(df):
    print df.head()

# tssTable, p1p2Table, sgInfoTable, libraryTable = loadTestData()

# results = predictWeissmanScore(tssTable=tssTable, p1p2Table=p1p2Table, sgrnaTable=sgInfoTable, libraryTable=libraryTable, modality="CRISPRa", verbose = False)

