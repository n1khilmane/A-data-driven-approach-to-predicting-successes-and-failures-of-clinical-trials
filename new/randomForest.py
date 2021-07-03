import pandas as pd
file = pd.read_csv("./forest_input_to_py.csv")
X = file.iloc[:,1:17]
l = file.iloc[:,17:18]

l = l.T
y = []
for i in range(0,len(X)):
    y = y+[l[i][0]]

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,random_state = 30) # 70% training and 30% test



#Import Random Forest Model
from sklearn.ensemble import RandomForestClassifier


#Create a Gaussian Classifier
clf=RandomForestClassifier(n_estimators=500)

#Train the model using the training sets y_pred=clf.predict(X_test)
clf.fit(X_train,y_train)

y_pred=clf.predict(X_test)

#Import scikit-learn metrics module for accuracy calculation
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score,plot_confusion_matrix
import matplotlib.pyplot as plt

# Model Accuracy, how often is the classifier correct?
print('Accuracy: %.3f' % accuracy_score(y_test, y_pred))
print('Precision: %.3f'% precision_score(y_test, y_pred,pos_label="passed"))
print('Recall: %.3f' % recall_score(y_test, y_pred,pos_label="passed"))
print('F1 Score: %.3f' % f1_score(y_test, y_pred,pos_label="passed"))


titles_options = [("Confusion matrix, without normalization", None),
                  ("Normalized confusion matrix", 'true')]
for title, normalize in titles_options:
    disp = plot_confusion_matrix(clf, X_test, y_test,
                                 cmap=plt.cm.Blues,
                                 normalize=normalize)
    disp.ax_.set_title(title)

    print(title)
    print(disp.confusion_matrix)

plt.show()