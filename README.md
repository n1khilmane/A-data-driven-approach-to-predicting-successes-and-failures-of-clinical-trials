We introduce a new data-driven approach that directly predicts the likelihood of toxicity in clinical trials. The model integrates properties of a compound’s targets and its structure to provide a new measure. Drug target network connectivity and expression levels, along with molecular weight, were identified as important indicators of adverse clinical events. Altogether, our method provides a data-driven broadly applicable strategy to identify drugs likely to possess manageable toxicity in clinical trials and will help drive the design of therapeutic agents with less toxicity.

## Prerequisite:

Python 3.8 installed and running.


## Installation:
```
git clone https://github.com/S-P-RAJAT/final-year-project.git
python3 -m venv venv
. venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

## Execution:


### Running the front-end:

```
cd front-end/
export FLASK_APP=ui
export FLASK_ENV=development
flask run
```
### Running the model:
```
cd code/
python3 PrOCTOR.py
python3 randomForest.py
