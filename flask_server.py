import sys 
import os 
from flask import Flask, flash, jsonify, render_template, redirect, url_for, request
from flaskwebgui import FlaskUI
import time
import tensorflow as tf
import tensorflow_addons as tfa
from keras.models import load_model
from tensorflow.keras.preprocessing import sequence
import numpy as np
import pandas as pd
import easygui


max_sequence_size = 2000
# target_function = '0005524' # ATP binding
# target_function = '0046872' # metal ion binding
# target_function = '0008270' # zinc ion binding
# target_function = '0003677' # DNA binding
# target_function = '0003676' # nucleic acid binding
target_function_index = {0: '0005524',
                         1: '0046872',
                         2: '0008270',
                         3: '0003677',
                         4: '0003676'}

acid_letters = ['_', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y']

# ========================================================================== #
def sequence_to_indices(sequence):
    try:
        indices = [acid_letters.index(c) for c in list(sequence)]
        return indices
    except Exception:
        print (sequence)
        raise Exception
# ========================================================================== #


def load_file():
    path = easygui.fileopenbox(filetypes=["*.txt", "*.fa"])
    with open(path) as f:
        
        X, X_all, infos = [], [], []
        rang = []
        pos = 0
        data = f.read().split('>')[1:]
        for element in data:
            pos += 1
            prot = element.split('\n')
            info = prot[0]
            seq = ''
            for i in range(1, len(prot)):
                seq += prot[i]

            # we're doing this to reduce input size
            if len(seq) > max_sequence_size:
                continue

            try:                                            
                X.append(seq.upper())
                infos.append(info)
                rang.append(pos)
                
            except KeyError:
                pass # For some proteins, we don't have annotations. skip these
    
    for i in range(len(X)): # need to do this in loop to ensure ordering
        x = sequence_to_indices(X[i])        
        X_all.append(x)
    
    ############################################
    #x = sequence_to_indices(X[0])        
    #X_all.append(x)
    ############################################
    X_all = np.asarray(X_all)
    X_all = sequence.pad_sequences(X_all, maxlen=max_sequence_size)
    
    return (rang, infos, X_all)



def prediction(inputt):
    global model
    
    inputt = pd.DataFrame(inputt)
    predict = model.predict(inputt)
    #val = np.argmax(predict)
    val = np.array(predict)
    
    return val



if getattr(sys, 'frozen', False): 
    template_folder = os.path.join(sys._MEIPASS, 'templates') 
    static_folder = os.path.join(sys._MEIPASS, 'static') 
    print(template_folder) 
    print(static_folder) 
    app = Flask(__name__, template_folder=template_folder, static_folder=static_folder) 
    ui = FlaskUI(app, width=900, height=650)
else: 
    app = Flask(__name__)
    #ui = FlaskUI(app, width=900, height=650, app_mode=True)
    ui = FlaskUI(app, width=900, height=650)



@app.route("/") 
def hello(): 
    return render_template("gui.html") 


@app.route("/home") 
def home():
    return render_template("home.html") 


@app.route("/loadModel")
def loadModel():
    global model
    model = load_model('model/model_221.h5')
    time.sleep(3.0)
    print('\n****************************************\n')
    print(model.summary())
    print('\n****************************************\n')
    print('Model Loaded')
    return jsonify("oh so slow")


@app.route('/submit_form', methods=['POST'])
def submit_form():
    
    seq = (request.form['sequence']).upper()
    x = [sequence_to_indices(seq)]
    inputt = sequence.pad_sequences(x, maxlen=max_sequence_size)
    
    output = prediction(inputt)
    output = output*100
    output = np.around(output, 2)
    columns = ['ATP binding', 'metal ion binding', 'zinc ion binding', 'DNA binding', 'nucleic acid binding']
    print ('\n############################################', output, '\n############################################')
    
    data = dict(seq=seq, a=output[0][0], b=output[0][1], c=output[0][2], d=output[0][3], e=output[0][4])
    
    return render_template('predict.html', data=data, columns=columns)


@app.route('/upload_predict')
def upload_predict():
        
    (rang, infos, inputt) = load_file()
    output = prediction(inputt)
    print(output)
    
    #df = pd.DataFrame(output)
    output = np.array(output)
    output = output*100
    output = np.around(output, 2)
    columns = ['ATP binding', 'metal ion binding', 'zinc ion binding', 'DNA binding', 'nucleic acid binding']
    liste = []
    compter = -1
    for row in output:
        compter += 1
        infos[compter] = infos[compter].rstrip('\n')        
        db = str(infos[compter].split('|')[0])
        go = str(infos[compter].split('|')[1])
        name = str(infos[compter].split('|')[2])
        if db == 'sp':
            db = '(SP) Swiss-Prot'
        values = dict(rang=rang[compter], a=row[0], b=row[1], c=row[2], d=row[3], e=row[4], db=db, go=go, name=name)
        print('\n/////////////////////////////////////\n',values,'\n///////////////////////////////////////\n')
        liste.append(values)
    index = rang
    return render_template('upload_predict.html', data=liste, columns=columns)
    
    #return "File uploaded !"
  

if __name__ == "__main__": 
    model = None 
    
    app.run(host='127.0.0.1', port=5000, debug=True) 
    #ui.run()

