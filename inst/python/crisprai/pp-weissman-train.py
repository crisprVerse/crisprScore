# python script for weissman/horlbeck algorithm for sgRNA predictions
# adapted from https://github.com/mhorlbeck/CRISPRiaDesign/blob/master/CRISPRiaDesign_example_notebook.md


from sgRNA_learning3 import *    # contains functions to load genome and empirical data
import cPickle
import  sys


### Data input variables ###

genome_type = "hg19"
# genome_type = "lifted_hg38"

# specify i or a for CRISPRi or CRISPRa, respectively
MODALITY = 'CRISPRi'

# genome fasta file
if genome_type == 'hg19':
    GENOME_FASTA =  'input_files/' + genome_type + '/hg19.fa'
elif genome_type ==  'lifted_hg38':
    GENOME_FASTA =  'input_files/' + genome_type + '/hg38.fa'

# paths to genome and bigWig files containing chromatin data of interest
# see https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeOpenChromDnase    
CHROMATIN_DNASE = 'input_files/' + genome_type + '/wgEncodeOpenChromDnaseK562BaseOverlapSignalV2_' + genome_type + '.bigWig'
CHROMATIN_FAIRE = 'input_files/' + genome_type + '/wgEncodeOpenChromFaireK562Sig_' + genome_type + '.bigWig'
CHROMATIN_MNASE = 'input_files/' + genome_type + '/wgEncodeSydhNsomeK562Sig_' + genome_type + '.bigWig'

# paths to model training data from CRISPRi data from https://github.com/mhorlbeck/CRISPRiaDesign/blob/master/Library_design_walkthrough.md
# can also swap these to the provided CRISPRa data
if MODALITY == 'CRISPRi':
    LIBRARY_TABLE_TRAINING= 'input_files/' + genome_type + '/CRISPRi_trainingdata_libraryTable.txt'
    SGRNA_TABLE_TRAINING = 'input_files/' + genome_type + '/CRISPRi_trainingdata_sgRNAInfoTable_' + genome_type + '.txt'
    TRAINING_SETS = 'input_files/hg19/CRISPRi_trainingdata_traintestsets.txt'
    ACTIVITY_SCORES = 'input_files/hg19/CRISPRi_trainingdata_activityScores.txt'
elif MODALITY == 'CRISPRa':
    LIBRARY_TABLE_TRAINING= 'input_files/' + genome_type + '/CRISPRa_trainingdata_libraryTable.txt'
    SGRNA_TABLE_TRAINING = 'input_files/' + genome_type + '/CRISPRa_trainingdata_sgRNAInfoTable_' + genome_type + '.txt'
    TRAINING_SETS = 'input_files/hg19/CRISPRa_trainingdata_traintestsets.txt'
    ACTIVITY_SCORES = 'input_files/hg19/CRISPRa_trainingdata_activityScores.txt'

# paths to annotation files
TSS_DATA = 'input_files/' + genome_type + '/human_tssTable_' +  genome_type + '.txt'
P1P2_DATA = 'input_files/' + genome_type + '/human_p1p2Table_' + genome_type + '.txt'


### Data output variables ###

# trained model pickle file
PICKLE_FILE = 'trained_models/' + MODALITY + '_estimator_weissman_' +  genome_type + '.pkl'


### Loading data (genome, annoations) ###

# load genome files
genomeDict=loadGenomeAsDict(GENOME_FASTA)

print "loading chromatin data of interest..."
# load chromatin data of interest
bwhandleDict = {'dnase':BigWigFile(open(CHROMATIN_DNASE)), 'faire':BigWigFile(open(CHROMATIN_FAIRE)), 'mnase':BigWigFile(open(CHROMATIN_MNASE))}

print "loading training data..."
# Loading training data 
libraryTable_training = pd.read_csv(LIBRARY_TABLE_TRAINING, sep='\t', index_col = 0)
sgInfoTable_training = pd.read_csv(SGRNA_TABLE_TRAINING, sep='\t', index_col=0)

print "loading human TSS data..."
# load human TSS data
tssTable = pd.read_csv(TSS_DATA, sep='\t', index_col=range(2))

print "loading activity score  data..."
# load activity score data
activityScores = pd.read_csv(ACTIVITY_SCORES, sep='\t',index_col=0, header=None).iloc[:,0]

print "loading human p1p2 data..."
# load human primary (p1) and secondary (p2) TSS data
p1p2Table = pd.read_csv(P1P2_DATA, sep='\t', header=0, index_col=range(2)) 

print "splitting p1p2 primary and secondary into separate columns..."
# df contains both primary and secondary in different columns, so need to split into seprate dfs
p1p2Table['primary TSS'] = p1p2Table['primary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]), int(tupString.strip('()').split(', ')[1].split('.')[0])))
p1p2Table['secondary TSS'] = p1p2Table['secondary TSS'].apply(lambda tupString: (int(tupString.strip('()').split(', ')[0].split('.')[0]),int(tupString.strip('()').split(', ')[1].split('.')[0])))

print "calculating parameters..."
# calculate paramters
paramTable_trainingGuides = generateTypicalParamTable(libraryTable_training, sgInfoTable_training, tssTable, p1p2Table, genomeDict, bwhandleDict)


### Training model ###

print 'Training model...'

print "loading 5-fold cross-validation splits.."
# load in the 5-fold cross-validation splits used to generate the model
with open(TRAINING_SETS) as infile:
    geneFold_train, geneFold_test, fitTable = cPickle.load(infile)

print "fitting parameters..."
# fit parameters
transformedParams_train, estimators = fitParams(paramTable_trainingGuides.loc[activityScores.dropna().index].iloc[geneFold_train], activityScores.loc[activityScores.dropna().index].iloc[geneFold_train], fitTable)

print "transforming parameters..."
# transform parameters
transformedParams_test = transformParams(paramTable_trainingGuides.loc[activityScores.dropna().index].iloc[geneFold_test], fitTable, estimators)

print "calculating regression..."
reg = linear_model.ElasticNetCV(l1_ratio=[.5, .75, .9, .99, 1], n_jobs=16, max_iter=2000)
print "scaler..."
scaler = preprocessing.StandardScaler()
print "fitting model..."
reg.fit(scaler.fit_transform(transformedParams_train), activityScores.loc[activityScores.dropna().index].iloc[geneFold_train])
print "predictedScores..."
predictedScores = pd.Series(reg.predict(scaler.transform(transformedParams_test)), index=transformedParams_test.index)
print "testScores..."
testScores = activityScores.loc[activityScores.dropna().index].iloc[geneFold_test]

print "Printing outputs..."
# print outputs
print 'Prediction AUC-ROC:', metrics.roc_auc_score((testScores >= .75).values, np.array(predictedScores.values,dtype='float64'))
print 'Prediction R^2:', reg.score(scaler.transform(transformedParams_test), testScores)
print 'Regression parameters:', reg.l1_ratio_, reg.alpha_
coefs = pd.DataFrame(zip(*[abs(reg.coef_),reg.coef_]), index = transformedParams_test.columns, columns=['abs','true'])
print 'Number of features used:', len(coefs) - sum(coefs['abs'] < .00000000001)


## Saving trained model ###

print 'Saving model...'

# save the model
transformedParams_train_header = transformedParams_train.head()
estimatorString = cPickle.dumps((fitTable, estimators, scaler, reg, transformedParams_train_header))

with open(PICKLE_FILE,'w') as outfile:
    outfile.write(estimatorString)
