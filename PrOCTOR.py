import pickle
import pandas as pd
import statistics as stats

from rdkit import Chem
from rdkit.Chem import Descriptors as props,rdmolops
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pprint

file = open("proctor.file", "rb")
data = pickle.load(file)

def proctor(SMILE,targets,type="full"):

	featureset=getAllFeatures(SMILE,targets)
	pprint.pp(featureset)
	return 0

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
        'Ro5': Ro5,
        'Ghose': Ghose,
        'Veber': Veber,
        'wQED': QED,
        'PC1': pc[0],
        'PC2': pc[1],
        'PC3': pc[2],
        }


SMILE   = "[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C"
targets = ["NR3C1","NR0B1","ANXA1","NOS2"]

proctor(SMILE,targets)