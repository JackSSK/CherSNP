import os
import re
from sklearn import svm
import gen_tools as t
import observer as obs

import numpy as np

# For sklearn < 0.23 versions joblib is there
try:
    from sklearn.externals import joblib
except:
    import joblib


class Format_error(Exception):
    pass

# This class is for loading Classifiers from files or Trainer directly
class Classifiers:
    def __init__(self, init = None, term = None, entCDS = None, outCDS = None,
        dict = None, filenames = "None"):
        # Initialization
        self.init = 0
        self.term = 0
        self.entCDS = 0
        self.outCDS = 0
        self.dict = {}

        # What if not loading from files
        if filenames == "None":
            try:
                self.init = init
                self.term = term
                self.entCDS = entCDS
                self.outCDS = outCDS
                self.dict = dict
            except:
                raise Format_error("Incorrect Classifiers format")

        # Loading from other files
        else:
            if len(filenames) != 5:
                raise Format_error("Predict: Incorrect Classifiers load file format")
            self._load(filenames)

    def _load(self, filenames):
        self.init = joblib.load(filenames[0])
        self.term = joblib.load(filenames[1])
        self.entCDS = joblib.load(filenames[2])
        self.outCDS = joblib.load(filenames[3])
        self.dict = t.decode_json(filenames[4])

# This class is to train out SVM classifiers using sklearn package
class Trainer:
    def __init__(self, seq_file, gff_file, save = True,
        mode = "hasRNA", filenames = "None"):
        print("Trace On!")
        self.observ = obs.Observer(seq_file, gff_file, mode)
        # test to see dictionary and observations in a JSON file

        print(" I am the bone of my swords")
        # init part
        init_clf = self._fit("init")

        print(" Steel is my body, and fire is my blood")
        # term part
        term_clf = self._fit("term")

        print(" I haved created over thousands blades")
        # entCDS part
        entCDS_clf = self._fit("entCDS")

        print(" Unknown to death, nor known to life")
        # outCDS part
        outCDS_clf = self._fit("outCDS")

        print(" Have withstood pain to creat many weapons")
        self.clfs = Classifiers(init_clf, term_clf, entCDS_clf, outCDS_clf,
            self.observ.dict)

        print(" Yet, those hands will never hold anything")
        if filenames == "None":
            print("Trainer: One or more filename missing, using default")
            filenames = [
                'default/Katyusha.pkl',
                'default/Erika.pkl',
                'default/Nadeshiko.pkl',
                'default/Juliet.pkl',
                'default/Jeanne.js'
            ]

        self.filenames = filenames

        if save:
            self._save()
        print(" So as I pray, Unlimited Blade Works!")
    # This function is to build classifier
    def _fit(self, name):
        clf = svm.SVC(kernel='rbf',class_weight={0: 10})
        correct = len(self.observ.subj[name]["correct"])
        wrong = len(self.observ.subj[name]["wrong"])
        notation = [1] * correct + [0] * wrong
        weight = [2] * correct + [1] * wrong
        obser = self.observ.subj[name]["correct"] + self.observ.subj[name]["wrong"]
        try:
            clf.fit(obser,notation,sample_weight=weight)
        except:
            print("\nWARNING: no obervation\n")

        return clf

    #Save classifiers into files
    def _save(self):
        joblib.dump(self.clfs.init, self.filenames[0])
        joblib.dump(self.clfs.term, self.filenames[1])
        joblib.dump(self.clfs.entCDS, self.filenames[2])
        joblib.dump(self.clfs.outCDS, self.filenames[3])
        t.encode_json(self.observ.dict, out = self.filenames[4])
