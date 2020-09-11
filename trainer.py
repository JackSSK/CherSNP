import os
import re
from sklearn import svm
import gen_tools as t
import observer as obs
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

        elif filenames == "Default":
            filenames = [
                'default/Katyusha.pkl',
                'default/Erika.pkl',
                'default/Nadeshiko.pkl',
                'default/Mary.pkl',
                'default/Jeanne.js'
            ]
            self._load(filenames)

        else:
            if len(filenames) != 5:
                raise Format_error("Incorrect Classifiers load file format")
            self._load(filenames)

    def _load(self, filenames):
        self.init = joblib.load(filenames[0])
        self.term = joblib.load(filenames[1])
        self.entCDS = joblib.load(filenames[2])
        self.outCDS = joblib.load(filenames[3])
        self.dict = t.decode_json(filenames[4])

# This class is to train out SVM classifiers using sklearn package
class Trainer:
    def __init__(self, seq_file, gff_file, save = True, filenames = "None"):
        self.observ = obs.Observer(seq_file, gff_file)
        # test to see dictionary and observations in a JSON file
        # t.encode_json([self.observ.dict, self.observ.subj])

        # init part
        init_clf = self._fit("init")

        # term part
        term_clf = self._fit("term")

        # entCDS part
        entCDS_clf = self._fit("entCDS")

        # outCDS part
        outCDS_clf = self._fit("outCDS")

        self.clfs = Classifiers(init_clf, term_clf, entCDS_clf, outCDS_clf,
            self.observ.dict)

        if len(filenames) != 5:
            print('Warning: Invalid filenames, Will use default filenames')
            filenames = "None"
        if filenames == "None":
            filenames = [
                'default/Katyusha.pkl',
                'default/Erika.pkl',
                'default/Nadeshiko.pkl',
                'default/Mary.pkl',
                'default/Jeanne.js'
            ]

        self.filenames = filenames

        if save:
            self._save()

    # This function is to build classifier
    def _fit(self, name):
        clf = svm.SVC(kernel='rbf')
        correct = len(self.observ.subj[name]["correct"])
        wrong = len(self.observ.subj[name]["wrong"])
        notation = [1] * correct + [0] * wrong
        weight = [2] * correct + [1] * wrong
        obser = self.observ.subj[name]["correct"] + self.observ.subj[name]["wrong"]
        clf.fit(obser,notation,sample_weight=weight)
        # clf.fit(obser,notation)

        # test here to see accuracy with data just trained with
        # answer = clf.predict(obser)
        # correct = 0
        # for i in range(len(answer)):
        #     if answer[i] == notation[i]:
        #         correct += 1
        # print(name, ': ', correct/len(answer))

        return clf

    #Save classifiers into files
    def _save(self):
        joblib.dump(self.clfs.init, self.filenames[0])
        joblib.dump(self.clfs.term, self.filenames[1])
        joblib.dump(self.clfs.entCDS, self.filenames[2])
        joblib.dump(self.clfs.outCDS, self.filenames[3])
        t.encode_json(self.observ.dict, out = self.filenames[4])
