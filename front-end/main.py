import math
import pickle
import pprint
import pandas as pd
import statistics as stats

from rdkit import Chem
from rdkit.Chem import Descriptors as props,rdmolops

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split

file = open("proctor.file", "rb")
data = pickle.load(file)

def createRandomForest(s):
	file = pd.read_csv(s)
	X = file.iloc[:,1:17]
	l = file.iloc[:,17:18]
	l = l.T
	y = [l[i][0] for i in range(0,len(X))]

	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,random_state =30)
	
	clf=RandomForestClassifier(n_estimators=500)
	clf.fit(X_train,y_train)

	return clf


def proctor(SMILE,targets,type="full"):
	try:
		featureset=getAllFeatures(SMILE,targets)
	except:
		d = {"features":{},"prediction":None,"probability":None,"score":None,"error":1}
		return d
	for key,value in featureset.items():
		featureset[key] = round(featureset[key],6)
	df = pd.DataFrame(tuple(featureset.values())).T
	
	clf = createRandomForest("./forest_input_to_py.csv")
	prob = clf.predict_proba(df)[0][1]
	pred = ("Safe","Toxic")[int(prob<0.5)]
	d = { "features": featureset,"prediction": "\nPrediction:\t"+str(pred),"probability" : "\nProbability:\t"+str(prob),"score":  "\nProctor Score:\t"+str(math.log2(prob/(1-prob)))}
	return d

def get_PC_value(self,targets):

	#target features of genes
	expr = pd.DataFrame([self.xs(i) for i in targets])

	#Median of target features of genes
	expr_median = expr.apply(stats.median,0) 
	
	#Standard Scaling
	scale 		= StandardScaler()
	pc 			= data["target_data"]["pc"]
	data_scaled = scale.fit_transform(pc)
	scaler 		= StandardScaler().fit(pc)
	pc_center 	= scaler.mean_ #pc$center
	pc_scale 	= scaler.scale_ #pc$scale
	
	#Principal Component Analysis
	pca = PCA(n_components=30)
	data_fit_transform = pca.fit_transform(data_scaled)

	#rotation matrix
	loadings = pd.DataFrame(pca.components_.T)   #$pc$rotation

	#centering and scaling
	m = pd.DataFrame(list(expr_median))
	l = m.T.apply(lambda row: row - pc_center, axis=1)
	k = l.apply(lambda row: row / pc_scale, axis=1)

	#Prediction for PCA
	pc_preds = k.dot(loadings) #pc_preds 
	pc_preds = list(pc_preds.xs(0))
	return pc_preds[:3]

def getAllFeatures(SMILE,targets):

	m 		= Chem.MolFromSmiles(SMILE)
	
	MW 		= props.MolWt(m) # molecular weight
	XlogP 	= props.MolLogP(m) # octanol-water partition coefficient log P
	HBD 	= props.NumHDonors(m) #hydrogen bond donor count
	HBA 	= props.NumHAcceptors(m) #hydrogen bond acceptor count
	PSA 	= props.TPSA(m) #polar surface area
	FC 		= rdmolops.GetFormalCharge(m)#formal charge
	RBC 	= props.NumRotatableBonds(m) #rotatable bonds count
	refr 	= props.MolMR(m) # refractivity
	alogP 	= None #
	nA 		= m.GetNumAtoms() #sum(atomcountMA(mol,addH=FALSE)) #number atoms
	AROMs 	= props.NumAromaticRings(m)
	nALERTS = len([1 for i in data["ualerts"] if m.HasSubstructMatch(Chem.MolFromSmarts(i))])

	Ro5 	= int(MW < 500 and HBD < 5 and HBA < 10 and XlogP < 5)
	Veber 	= int(RBC <= 10 and PSA <= 140)
	Ghose 	= int(PSA <140 and (-0.4 <= XlogP < 5.6) and (160 <= MW < 480) and (20 <= nA < 70))
	QED 	= Chem.QED.qed(m)
	
	t_set  	= set(targets)
	(lof, btwn, degree) = (data['target_data']['lof'],
                           data['target_data']['btwn'],
                           data['target_data']['degree'])

	lossFreq  = max([lof.xs(i)[1] / lof.xs(i)[2] for i in lof.index if i in t_set])
	maxBtwn   = max([btwn.xs(i)[0] for i in btwn.index if i in t_set])
	maxDegree = max([degree.xs(i)[0] for i in degree.index if i in t_set])

	pc = get_PC_value(data["target_data"]["expr"],targets)

	return {
        'MolecularWeight': MW,
        'XLogP': XlogP,
        'HydrogenBondDonorCount': HBD,
        'HydrogenBondAcceptorCount': HBA,
        'PolarSurfaceArea': PSA,
        'FormalCharge': FC,
        'NumRings': AROMs,
        'RotatableBondCount': RBC,
        'Refractivity': refr,
        'lossFreq': lossFreq,
        'maxBtwn': maxBtwn,
        'maxDegree': maxDegree,
        # 'Ro5': Ro5,
        # 'Ghose': Ghose,
        # 'Veber': Veber,
        'wQED': QED,
        'PC1': pc[0],
        'PC2': pc[1],
        'PC3': pc[2],
        }


SMILE   = "[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C"
targets = ["NR3C1","NR0B1","ANXA1","NOS2"] #dexamethasone

SMILE2 = "CCN(CC)C(=O)C1(CC1CN)C1=CC=CC=C1"
targets2 = ["SLC6A4","SLC6A2"] #Milnacipran s

SMILE3   = "CN(C)CCOC(C1=CC=C(Cl)C=C1)C1=CC=CC=N1"
targets3 = ["HRH1"] ##carbinoxamine s

SMILE4="CC1=C(C)C=C2N(C[C@H](O)[C@H](O)[C@H](O)CO)C3=NC(=O)NC(=O)C3=NC2=C1" #Riboflavin
targets4=["RFK","BLVRB"]


SMILE5 = "[H][C@@]1(CC[C@@]2([H])\C(CCC[C@]12C)=C\C=C1C[C@@H](O)C[C@H](O)C1)[C@H](C)\C=C\[C@H](C)C(C)(C)O"
targets5 = ["VDR"] #Paricalcitol s

SMILE6  = "CC(C)C[C@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)C1=CN=CC=N1)B(O)O" #bortezomib
targets6 = ["PSMB5","PSMB1"]# highly toxic

SMILE7 = "CCCCC1=NC2(CCCC2)C(=O)N1CC1=CC=C(C=C1)C1=CC=CC=C1C1=NNN=N1"
targets7 = ["AGTR1","JUN"] #Irbesartan s


SMILE8 = "[H][C@@]12C=C(C)CC[C@@]1([H])C(C)(C)OC1=C2C(O)=CC(CCCCC)=C1"
targets8 = ["CNR1","CNR2"] # dronabinol t

SMILE9   = "CN1CCC2=C(C1)C(C1=CC=CC=C21)C1=CC=CC=C1"
targets9 = ["HRH1"]  #phenindamine s

#abiraterone not existing in drugbank
def execute(drug,target):
	return proctor(drug,target)
