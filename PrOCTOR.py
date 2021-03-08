####
# PrOCTOR.R
# Description: Run PrOCTOR for any given structure/target pair
# Last Updated: October 20, 2016
# Update: Fixed bug in plotting function. Added predictions for target/structure specific models and a plot all function.
# Functions:
#   1) PrOCTOR(SMILE,targets,type)
#         - applies PrOCTOR to inputted drug
#         - type = type of PrOCTOR model to use: full, target, structure (default=full)
#         - example: Dexamethasone
#             --> PrOCTOR(SMILE="[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",targets=c("NR3C1","NR0B1","ANXA1","NOS2"))
#
#   2) PrOCTOR_plot(SMILE,targetstype)
#         - applies PrOCTOR and displays a graph that compares the inputted drug to the training set
#
#   2) PrOCTOR_plotall(SMILE,targets)
#         - applies PrOCTOR and displays graphs that compares the inputted drug to the training set
#           for full, target, and structural version of PrOCTOR


import pickle
import pandas
import statistics as stats

from rdkit import Chem
from rdkit.Chem import Descriptors as props,rdmolops
import rpy2.robjects as robjects


readRDS = robjects.r['readRDS']
df = readRDS('proctor.rds')
ualerts = list(df.rx2("unwantedALERTS"))

ualerts.remove("(C(=O)O[C,H1]).(C(=O)O[C,H1]).(C(=O)O[C,H1])")

def proctor(SMILE,targets,type="full"):
	featureset=getAllFeatures(SMILE,targets)
	return 0


def getAllFeatures(SMILE,targets):
	m = Chem.MolFromSmiles(SMILE)
	MW= props.MolWt(m) # molecular weight
	XlogP=props.MolLogP(m) # octanol-water partition coefficient log P
	HBD=props.NumHDonors(m) #hydrogen bond donor count
	HBA=props.NumHAcceptors(m) #hydrogen bond acceptor count
	PSA=props.TPSA(m) #polar surface area
	FC= rdmolops.GetFormalCharge(m)#formal charge
	RBC=props.NumRotatableBonds(m) #rotatable bonds count
	refr=props.MolMR(m) # refractivity
	alogP=None #
	nA= m.GetNumAtoms() #sum(atomcountMA(mol,addH=FALSE)) #number atoms
	AROMs=props.NumAromaticRings(m)
	nALERTS = len([1 for i in ualerts if m.HasSubstructMatch(Chem.MolFromSmarts(i))])
	nALERTS=0
	Ro5=int(MW<500 and HBD<5 and HBA<10 and XlogP<5)
	Veber=int(RBC<=10 and PSA<=140)
	Ghose=int(PSA<140 and (-0.4<=XlogP<5.6) and (160<=MW<480) and (20<=nA<70))
	QED = Chem.QED.qed(m)

	print(MW,XlogP,HBD,HBA,PSA,FC,RBC,refr,alogP,nA,AROMs,nALERTS,Ro5,Veber,Ghose,QED)
	# return(c(MolecularWeight=MW,XLogP=XlogP,HydrogenBondDonorCount=HBD,HydrogenBondAcceptorCount=HBA,PolarSurfaceArea=PSA,FormalCharge=FC,NumRings=AROMs,RotatableBondCount=RBC,Refractivity=refr,lossFreq=max(data$target_data$lof$deleterious[rownames(data$target_data$lof) %in% targets]/data$target_data$lof$Total[rownames(data$target_data$lof) %in% targets]),maxBtwn=max(data$target_data$btwn[names(data$target_data$btwn) %in% targets]),maxDegree=max(data$target_data$degree[names(data$target_data$degree) %in% targets]),Ro5=Ro5,Ghose=Ghose,Veber=Veber,wQED=QED$wQED,get_PC_value(targets)))


def get_PC_value(self,targets):
	expr = pandas.DataFrame([self.xs(i) for i in targets])
	expr.apply(stats.median,1)
	# print(expr)

with open("data$target_data$expr.file", "rb") as file:
	df = pickle.load(file)

SMILE="[H][C@@]12C[C@@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C"
targets = ["NR3C1","NR0B1","ANXA1","NOS2"]
get_PC_value(df,targets)

proctor(SMILE,targets)