import os
import sys
import subprocess
import tempfile
import multiprocessing
import numpy as np 
import scipy as sp 
import pandas as pd
from ConfigParser import SafeConfigParser
from Bio import Seq, SeqIO
import pysam
from sklearn import linear_model, svm, ensemble, preprocessing, grid_search, metrics
from expt_config_parser import parseExptConfig, parseLibraryConfig
from bx.bbi.bigwig_file import BigWigFile


def generateSgrnaDistanceTable_p1p2Strategy(sgInfoTable, libraryTable, p1p2Table, transcripts=False):

	sgDistanceSeries = []

	if transcripts == False: # when sgRNAs weren't designed based on the p1p2 strategy

		for name, group in sgInfoTable['pam_coordinate'].groupby(libraryTable['gene']):
			
			try:

				if name in p1p2Table.index:

					tssRow = p1p2Table.loc[name]

					if len(tssRow) == 1:
						tssRow = tssRow.iloc[0]
						for sgId, pamCoord in group.iteritems():
							if (tssRow['strand'] == '+'):
								sgDistanceSeries.append((sgId, name, tssRow.name,
									pamCoord - tssRow['primary TSS'][0],
									pamCoord - tssRow['primary TSS'][1],
									pamCoord - tssRow['secondary TSS'][0],
									pamCoord - tssRow['secondary TSS'][1]))
							else:
								sgDistanceSeries.append((sgId, name, tssRow.name,
									(pamCoord - tssRow['primary TSS'][1]) * -1,
									(pamCoord - tssRow['primary TSS'][0]) * -1,
									(pamCoord - tssRow['secondary TSS'][1]) * -1,
									(pamCoord - tssRow['secondary TSS'][0]) * -1))

					else:

						for sgId, pamCoord in group.iteritems():
							closestTssRow = tssRow.loc[tssRow.apply(lambda row: abs(pamCoord - row['primary TSS'][0]), axis=1).idxmin()]

							if isinstance(closestTssRow, pd.core.frame.DataFrame):
								
								if len(closestTssRow) > 1:
									for i in range(len(closestTssRow)):

										# convert current row to pd Series
										current_closestTssRow = closestTssRow.iloc[i]

										if (current_closestTssRow['strand'] == '+'):
											sgDistanceSeries.append(((sgId + str(i)), name, current_closestTssRow.name,
												pamCoord - current_closestTssRow['primary TSS'][0],
												pamCoord - current_closestTssRow['primary TSS'][1],
												pamCoord - current_closestTssRow['secondary TSS'][0],
												pamCoord - current_closestTssRow['secondary TSS'][1]))
										else:
											sgDistanceSeries.append(((sgId + str(i)), name, current_closestTssRow.name,
												(pamCoord - current_closestTssRow['primary TSS'][1]) * -1,
												(pamCoord - current_closestTssRow['primary TSS'][0]) * -1,
												(pamCoord - current_closestTssRow['secondary TSS'][1]) * -1,
												(pamCoord - current_closestTssRow['secondary TSS'][0]) * -1))
							
							else:

								if (closestTssRow['strand'] == '+'):
									sgDistanceSeries.append((sgId, name, closestTssRow.name,
										pamCoord - closestTssRow['primary TSS'][0],
										pamCoord - closestTssRow['primary TSS'][1],
										pamCoord - closestTssRow['secondary TSS'][0],
										pamCoord - closestTssRow['secondary TSS'][1]))
								else:
									sgDistanceSeries.append((sgId, name, closestTssRow.name,
										(pamCoord - closestTssRow['primary TSS'][1]) * -1,
										(pamCoord - closestTssRow['primary TSS'][0]) * -1,
										(pamCoord - closestTssRow['secondary TSS'][1]) * -1,
										(pamCoord - closestTssRow['secondary TSS'][0]) * -1))
					
			except:
				raise Exception("Error: cannnot match to promoter table.")		

	else:
		for name, group in sgInfoTable['pam_coordinate'].groupby([libraryTable['gene'],libraryTable['transcripts']]):
			if name in p1p2Table.index:
				tssRow = p1p2Table.loc[[name]]

				if len(tssRow) == 1:
					tssRow = tssRow.iloc[0]
					for sgId, pamCoord in group.iteritems():
						if (tssRow['strand'] == '+'):
							sgDistanceSeries.append((sgId, tssRow.name[0], tssRow.name[1],
								pamCoord - tssRow['primary TSS'][0],
								pamCoord - tssRow['primary TSS'][1],
								pamCoord - tssRow['secondary TSS'][0],
								pamCoord - tssRow['secondary TSS'][1]))
						else:
							sgDistanceSeries.append((sgId, tssRow.name[0], tssRow.name[1],
								(pamCoord - tssRow['primary TSS'][1]) * -1,
								(pamCoord - tssRow['primary TSS'][0]) * -1,
								(pamCoord - tssRow['secondary TSS'][1]) * -1,
								(pamCoord - tssRow['secondary TSS'][0]) * -1))

				else:
					print name, tssRow
					raise ValueError('all gene/trans pairs should be unique')

	return pd.DataFrame(sgDistanceSeries, columns=['sgId', 'gene', 'transcript', 'primary TSS-Up', 'primary TSS-Down', 'secondary TSS-Up', 'secondary TSS-Down']).set_index(keys=['sgId'])

def generateSgrnaLengthSeries(libraryTable):
	lengthSeries =  libraryTable.apply(lambda row: len(row['sequence']),axis=1)
	lengthSeries.name = 'length'
	return lengthSeries

def generateRelativeBasesAndStrand(sgInfoTable, tssTable, libraryTable, genomeDict):
	relbases = []
	strands = []
	sgIds = []
	for gene, sgInfoGroup in sgInfoTable.groupby(libraryTable['gene']):
		tssRowGroup = tssTable.loc[gene]

		if len(set(tssRowGroup['chromosome'].values)) == 1:
			chrom = tssRowGroup['chromosome'].values[0]
		else:
			raise ValueError('multiple annotated chromosomes for ' + gene)

		if len(set(tssRowGroup['strand'].values)) == 1:
			strand = tssRowGroup['strand'].values[0]
		else:
			raise ValueError('multiple annotated strands for ' + gene)

		for sg, sgInfo in sgInfoGroup.iterrows():
			sgIds.append(sg)
			geneTup = (sgInfo['gene_name'],','.join(sgInfo['transcript_list']))
			strands.append(True if sgInfo['strand'] == strand else False)

			baseMatrix = []
			for pos in np.arange(-30,10):
				baseMatrix.append(getBaseRelativeToPam(chrom, sgInfo['pam_coordinate'],sgInfo['length'], sgInfo['strand'], pos, genomeDict))
			relbases.append(baseMatrix)

	relbases = pd.DataFrame(relbases, index = sgIds, columns = np.arange(-30,10)).loc[libraryTable.index]
	strands = pd.DataFrame(strands, index = sgIds, columns = ['same strand']).loc[libraryTable.index]

	return relbases, strands

def generateBooleanBaseTable(baseTable):
	relbases_bool = []
	for base in ['A','G','C','T']:
		relbases_bool.append(baseTable.applymap(lambda val: val == base))

	return pd.concat(relbases_bool, keys=['A','G','C','T'], axis=1)

def generateBooleanDoubleBaseTable(baseTable):
	doubleBaseTable = []
	tableCols = []
	for b1 in ['A','G','C','T']:
		for b2 in ['A','G','C','T']:
			for i in np.arange(-30,8):
				doubleBaseTable.append(pd.concat((baseTable[i] == b1, baseTable[i+1] == b2),axis=1).all(axis=1))
				tableCols.append(((b1,b2),i))
	return pd.concat(doubleBaseTable, keys=tableCols, axis=1)

def getBaseRelativeToPam(chrom, pamPos, length, strand, relPos, genomeDict):
	rc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

	pamPos = int(pamPos)

	if strand == '+':
		return rc[genomeDict[chrom][pamPos - relPos].upper()]
	elif strand == '-':
		return genomeDict[chrom][pamPos + relPos].upper()
	else:
		raise ValueError()

def getMaxLengthHomopolymer(sequence, base):
	sequence = sequence.upper()
	base = base.upper()
	
	maxBaseCount = 0
	curBaseCount = 0
	for b in sequence:
		if b == base:
			curBaseCount += 1
		else:
			maxBaseCount = max((curBaseCount, maxBaseCount))
			curBaseCount = 0
		
	return max((curBaseCount, maxBaseCount))

def getFractionBaseList(sequence, baseList):
	baseSet = [base.upper() for base in baseList]
	counter = 0.0
	for b in sequence.upper():
		if b in baseSet:
			counter += 1.0
			
	return counter / len(sequence)

def getRNAfoldingTable(libraryTable):
	tempfile_fa = tempfile.NamedTemporaryFile('w+t', delete=False)
	tempfile_rnafold = tempfile.NamedTemporaryFile('w+t', delete=False)

	for name, row in libraryTable.iterrows():
		tempfile_fa.write('>' + name + '\n' + row['sequence'] + '\n')

	tempfile_fa.close()
	tempfile_rnafold.close()

	subprocess.call('RNAfold --noPS < %s > %s' % (tempfile_fa.name, tempfile_rnafold.name), shell=True)

	mfeSeries_noScaffold = parseViennaMFE(tempfile_rnafold.name, libraryTable)
	isPaired = parseViennaPairing(tempfile_rnafold.name, libraryTable)

	tempfile_fa = tempfile.NamedTemporaryFile('w+t', delete=False)
	tempfile_rnafold = tempfile.NamedTemporaryFile('w+t', delete=False)

	with open(tempfile_fa.name,'w') as outfile:
		for name, row in libraryTable.iterrows():
			outfile.write('>' + name + '\n' + row['sequence'] + 'GTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT\n')

	tempfile_fa.close()
	tempfile_rnafold.close()

	subprocess.call('RNAfold --noPS < %s > %s' % (tempfile_fa.name, tempfile_rnafold.name), shell=True)

	mfeSeries_wScaffold = parseViennaMFE(tempfile_rnafold.name, libraryTable)

	return pd.concat((mfeSeries_noScaffold, mfeSeries_wScaffold, isPaired), keys=('no scaffold', 'with scaffold', 'is Paired'), axis=1)

def parseViennaMFE(viennaOutputFile, libraryTable):
	mfeList = []
	with open(viennaOutputFile) as infile:
		for i, line in enumerate(infile):
			if i%3 == 2:
				mfeList.append(float(line.strip().strip('.() ')))
	return pd.Series(mfeList, index=libraryTable.index, name='RNA minimum free energy')

def parseViennaPairing(viennaOutputFile, libraryTable):
	paired = []
	with open(viennaOutputFile) as infile:
		for i, line in enumerate(infile):
			if i%3 == 2:
				foldString = line.strip().split(' ')[0]
				paired.append([char != '.' for char in foldString[-18:]])
	return pd.DataFrame(paired, index=libraryTable.index, columns = range(-20,-2))

def getChromatinDataSeriesByGene(bigwigFileHandle, libraryTable, sgInfoTable, p1p2Table, sgrnaDistanceTable_p1p2, colname = '', naValue = 0, normWindow = 1000):
	bwindex = bigwigFileHandle #BigWigFile(open(bigwigFile))

	chromatinScores = []
	for (gene, transcript), sgInfoGroup in sgInfoTable.groupby([sgrnaDistanceTable_p1p2['gene'], sgrnaDistanceTable_p1p2['transcript']]):   

		tssRow = p1p2Table.loc[[(gene, transcript)]].iloc[0,:]

		chrom = tssRow['chromosome']

		normWindowArray = bwindex.get_as_array(chrom, max(0, tssRow['primary TSS'][0] - normWindow), tssRow['primary TSS'][0] + normWindow)
		normWindowArray = np.nan_to_num(normWindowArray)    # replace nans with 0.0
		if normWindowArray is not None:
			normFactor = np.nanmax(normWindowArray)
		else:
			normFactor = 1

		windowMin = max(0, min(sgInfoGroup['pam_coordinate']) - max(sgInfoGroup['length']) - 10)
		windowMax = max(sgInfoGroup['pam_coordinate']) + max(sgInfoGroup['length']) + 10
		chromatinWindow = bwindex.get_as_array(chrom, windowMin, windowMax)
		chromatinWindow = np.nan_to_num(chromatinWindow)    # replace nans with 0.0

		chromatinScores.append(sgInfoGroup.apply(lambda row: getChromatinData(row, chromatinWindow, windowMin, normFactor), axis=1))


	chromatinSeries = pd.concat(chromatinScores)

	return chromatinSeries.fillna(naValue)

def getChromatinData(sgInfoRow, chromatinWindowArray, windowMin, normFactor):
	if sgInfoRow['strand'] == '+':
		sgRange = int(sgInfoRow['pam_coordinate']) + int(sgInfoRow['length'])
	else:
		sgRange = int(sgInfoRow['pam_coordinate']) - int(sgInfoRow['length'])


	if chromatinWindowArray is not None:# and len(chromatinWindowArray) > 0:
		chromatinArray = chromatinWindowArray[min(int(sgInfoRow['pam_coordinate']), sgRange) - int(windowMin): max(int(sgInfoRow['pam_coordinate']), sgRange) - int(windowMin)]
		return np.nanmean(chromatinArray)/normFactor
	else: #often chrY when using K562 data..
		return np.nan

def generateTypicalParamTable(libraryTable, sgInfoTable, tssTable, p1p2Table, genomeDict, bwFileHandleDict, verbose, transcripts=False):
	
	printNow('\nGenerating sgRNA length series...', verbose)
	lengthSeries = generateSgrnaLengthSeries(libraryTable)
	
	printNow('\nGenerating sgRNA distance table...', verbose)
	sgrnaPositionTable_p1p2 = generateSgrnaDistanceTable_p1p2Strategy(sgInfoTable, libraryTable, p1p2Table, transcripts)

	printNow('\nGenerating relative base and strand table...', verbose)
	baseTable, strand = generateRelativeBasesAndStrand(sgInfoTable, tssTable, libraryTable, genomeDict)

	printNow('\nGenerating boolean base tables...', verbose)
	booleanBaseTable = generateBooleanBaseTable(baseTable)
	doubleBaseTable = generateBooleanDoubleBaseTable(baseTable)

	printNow('\nGenerating homopolymer table...', verbose)
	baseList = ['A','G','C','T']
	homopolymerTable = pd.concat([libraryTable.apply(lambda row: np.floor(getMaxLengthHomopolymer(row['sequence'], base)), axis=1) for base in baseList],keys=baseList,axis=1)

	printNow('\nGenerating base fractions...', verbose)
	baseFractions = pd.concat([libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['A']),axis=1),
							libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['G']),axis=1),
							libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['C']),axis=1),
							libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['T']),axis=1),
							libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['G','C']),axis=1),
							libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['G','A']),axis=1),
							libraryTable.apply(lambda row: getFractionBaseList(row['sequence'], ['C','A']),axis=1)],keys=['A','G','C','T','GC','purine','CA'],axis=1)


	printNow('\nGenerating DNase table...', verbose)
	dnaseSeries = getChromatinDataSeriesByGene(bwFileHandleDict['dnase'], libraryTable, sgInfoTable, p1p2Table, sgrnaPositionTable_p1p2)

	printNow('\nGenerating FAIRE table...', verbose)
	faireSeries = getChromatinDataSeriesByGene(bwFileHandleDict['faire'], libraryTable, sgInfoTable, p1p2Table, sgrnaPositionTable_p1p2)

	printNow('\nGenerating MNase table...', verbose)
	mnaseSeries = getChromatinDataSeriesByGene(bwFileHandleDict['mnase'], libraryTable, sgInfoTable, p1p2Table, sgrnaPositionTable_p1p2)

	printNow('\nGenerating RNA folding table...', verbose)
	rnafolding = getRNAfoldingTable(libraryTable)

	printNow('\nMerging tables...', verbose)

	return pd.concat([lengthSeries,
		   sgrnaPositionTable_p1p2.iloc[:,2:],
		   homopolymerTable,
		   baseFractions,
		   strand,
		   booleanBaseTable['A'],
		   booleanBaseTable['T'],
		   booleanBaseTable['G'],
		   booleanBaseTable['C'],
		   doubleBaseTable,
		   pd.concat([dnaseSeries,faireSeries,mnaseSeries],keys=['DNase','FAIRE','MNase'], axis=1),
		   rnafolding['no scaffold'],
		   rnafolding['with scaffold'],
		   rnafolding['is Paired']],keys=['length',
		   'distance',
		   'homopolymers',
		   'base fractions',
		   'strand',
		   'base table-A',
		   'base table-T',
		   'base table-G',
		   'base table-C',
		   'base dimers',
		   'accessibility',
		   'RNA folding-no scaffold',
		   'RNA folding-with scaffold',
		   'RNA folding-pairing, no scaffold'],axis=1)

def transformParams(paramTable, fitTable, estimators):
	transformedParams = []
	
	for i, (name, col) in enumerate(paramTable.iteritems()):
		fitRow = fitTable.iloc[i]

		# replace nans with zeroes
		col = col.fillna(0)

		if fitRow['type'] == 'binary':
			transformedParams.append(col)
		elif fitRow['type'] == 'continuous':
			col_reshape = col.values.reshape(len(col),1)
			transformedParams.append(pd.Series(estimators[i].predict(col_reshape), index=col.index, name=name))
		elif fitRow['type'] == 'binnable':
			binStats = estimators[i]
			assignedBins = applyBins(col, binStats.index.values)
			transformedParams.append(assignedBins.apply(lambda binVal: binStats.loc[binVal]))

		elif fitRow['type'] == 'binnable_onehot':
			binStats = estimators[i]
			
			assignedBins = applyBins(col, binStats.index.values)
			binGroups = col.groupby(assignedBins)
			
			oneHotFrame = pd.DataFrame(np.zeros((len(assignedBins),len(binGroups))), index = assignedBins.index, \
				columns=pd.MultiIndex.from_tuples([(name[0],', '.join([name[1],key])) for key in sorted(binGroups.groups.keys())]))

			for groupName, group in binGroups:
				oneHotFrame.loc[group.index, (name[0],', '.join([name[1],groupName]))] = 1

			transformedParams.append(oneHotFrame)
			
	return pd.concat(transformedParams, axis=1)

def applyBins(column, binStrings):
	leftLabel = ''
	rightLabel = ''
	binTups = []
	for binVal in binStrings:
		if binVal[0] == '<':
			leftLabel = binVal
		elif binVal[0] == '>':
			rightLabel = binVal
			rightBound = float(binVal[3:])
		else:
			binTups.append((float(binVal),binVal))
			
	binTups.sort()

	leftBound = binTups[0][0]
	if leftLabel == '':
		leftLabel = binTups[0][1]
	
	if rightLabel == '':
		rightLabel = binTups[-1][1]
		rightBound = binTups[-1][0]

	def binFunc(val):
		return leftLabel if val < leftBound else (rightLabel if val >= rightBound else [tup[1] for tup in binTups if val >= tup[0]][-1])

	return column.apply(binFunc)

def loadGenomeAsDict(genomeFasta, verbose):
	printNow('\nLoading genome file...', verbose)
	genomeDict = SeqIO.to_dict(SeqIO.parse(genomeFasta,'fasta'))
	return genomeDict

def printNow(outputString, verbose = False):
	if verbose == True:
		sys.stdout.write(outputString)
		sys.stdout.flush()