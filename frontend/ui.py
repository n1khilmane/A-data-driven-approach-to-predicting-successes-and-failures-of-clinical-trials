from flask import Flask,render_template,request
from main import *
import pandas as pd
app = Flask(__name__)
file = pd.read_csv("targetsList.csv")

@app.route('/', methods =["GET", "POST"])
def hello_world(name=None):
    data = {"features":{},"prediction":None,"probability":None,"score":None}
    s = ""
    method = "GET"
    drug, targets= "",[]
    if request.method == "POST":
        drug = request.form.get("drug") 
        targets = request.form.getlist("targets")
        data = execute(drug,targets)
        if drug == "":
            data["error"] = 1
        if data["error"] != 1:
            return render_template('index1.html',name=name,data = data,drug = drug,tar_string = targets,method = method,file = file)+s
        method = "POST"
    return render_template('index.html', name=name,data = data,drug = drug,tar_string = targets,method = method,file = file)+s