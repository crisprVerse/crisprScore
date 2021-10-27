# python script for weissman/horlbeck algorithm for sgRNA predictions
# adapted from https://github.com/mhorlbeck/CRISPRiaDesign/blob/master/CRISPRiaDesign_example_notebook.md


import cPickle
from sgRNA_learning import *    # contains functions to load genome and empirical data


# genome fasta file
GENOME_FASTA =  'input_files/lifted_hg38/hg38.fa'

# paths to genome and bigWig files containing chromatin data of interest
# see https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeOpenChromDnase    
CHROMATIN_DNASE = 'input_files/lifted_hg38/wgEncodeOpenChromDnaseK562BaseOverlapSignalV2_lifted_hg38.bigWig'
CHROMATIN_FAIRE = 'input_files/lifted_hg38/wgEncodeOpenChromFaireK562Sig_lifted_hg38.bigWig'
CHROMATIN_MNASE = 'input_files/lifted_hg38/wgEncodeSydhNsomeK562Sig_lifted_hg38.bigWig'


def predictWeissmanScore(tssTable, p1p2Table, sgrnaTable, libraryTable, modality, verbose = False):

    # trained model pickle file
    PICKLE_FILE = 'trained_models/' + modality + '_estimator_weissman_hg19.pkl'
    
    # open pickle file to continue from previously trained session/model
    try:
        with open(PICKLE_FILE) as infile:
            fitTable, estimators, scaler, reg, transformedParams_train_header = cPickle.load(infile)
    except:
        raise Exception('Trained model file not found.') 

    paramTable = getParamTable(tssTable, p1p2Table, sgrnaTable, libraryTable, verbose = verbose)
    
    transformedParams_new = getTransformedParams(paramTable, fitTable, estimators, verbose = verbose)

    print 'Predicting sgRNA scores...'
    predictedScores = pd.Series(reg.predict(scaler.transform(transformedParams_new.loc[:, transformedParams_train_header.columns].fillna(0).values)), index=transformedParams_new.index)
    
    return predictedScores


def getParamTable(tssTable, p1p2Table, sgrnaTable, libraryTable, verbose = False):

    try:
        genomeDict=loadGenomeAsDict(GENOME_FASTA)
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

    try:
        paramTable = generateTypicalParamTable(libraryTable, sgrnaTable, tssTable, p1p2Table, genomeDict, bwhandleDict)
    except:
        raise Exception("Error generating parameter table.")

    return paramTable


def getTransformedParams(paramTable, fitTable, estimators, verbose = False):
    
    if verbose == True:
        print 'Transforming parameters...'

    try:
        transformedParams_new = transformParams(paramTable, fitTable, estimators)
    except:
        raise Exception("Error transforming parameters.")

    # reconcil e differences in column headers
    colTups = []
    for (l1, l2), col in transformedParams_new.iteritems():
        colTups.append((l1,str(l2)))
    transformedParams_new.columns = pd.MultiIndex.from_tuples(colTups)

    return transformedParams_new