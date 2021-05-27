from sgrna_modeler import datasets as da
from sgrna_modeler import models as sg
from sgrna_modeler import enzymes as en

def getEnPAMGB(sequences):	
	model = sg.SklearnSgrnaModel()
	model_weights = sg.get_enpam_gb()
	model.load_model(model_weights, en.cas12a, 'enPAM_GB')
	results = model.predict_seqs(sequences)
	return results

#model = sg.SklearnSgrnaModel()
#model_weights = sg.get_enpam_gb()
#model.load_model(model_weights, en.cas12a, 'enPAM_GB')
#sequences=["TGGTTTTAAAACAGAATATACAGTCTAAAAAACC","CATGTTTTTTTGGGAACCAATCGATAATCACATT"]
#prediction = model.predict_seqs(sequences)

# Existing data:
#kim_2018_train = da.load_kim_2018_train()
#model.fit(kim_2018_train)
#kim_2018_test = da.load_kim_2018_test()
#results = model.predict(kim_2018_test)