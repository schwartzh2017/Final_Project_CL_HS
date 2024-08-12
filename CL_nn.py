import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import copy
import torch.optim as optim
import tqdm
from sklearn.model_selection import StratifiedKFold, train_test_split

df1 = pd.read_csv('datasets/bcsc_risk_factors_expanded1.csv')
df2 = pd.read_csv('datasets/bcsc_risk_factors_expanded2.csv')
df3 = pd.read_csv('datasets/bcsc_risk_factors_expanded3.csv')

#slight cleaning
df = pd.concat([df1, df2, df3])
df = df[df.ne(9).all(1)] #drop unknowns (9s)
df.drop(['year'], axis=1, inplace=True)

history_counts = df.breast_cancer_history.value_counts()
print(history_counts[1]/history_counts.sum())

X = df.iloc[:, 0:df.shape[1]-1]
y = df.iloc[:, df.shape[1]-1]

X = torch.tensor(X.values, dtype=torch.float32)
y = torch.tensor(y.values, dtype=torch.float32).reshape(-1, 1)

class Deep(nn.Module):
    def __init__(self):
        super().__init__()
        self.layer1 = nn.Linear(10, 10)
        self.act1 = nn.ReLU()
        self.layer2 = nn.Linear(10, 10)
        self.act2 = nn.ReLU()
        self.layer3 = nn.Linear(10, 10)
        self.act3 = nn.ReLU()
        self.output = nn.Linear(10, 1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = self.act1(self.layer1(x))
        x = self.act2(self.layer2(x))
        x = self.act3(self.layer3(x))
        x = self.sigmoid(self.output(x))
        return x

def model_train(model, X_train, y_train, X_val, y_val):
    # loss function and optimizer
    loss_fn = nn.BCELoss()  # binary cross entropy
    optimizer = optim.Adam(model.parameters(), lr=0.0001)

    n_epochs = 250   # number of epochs to run
    batch_size = 10  # size of each batch
    batch_start = torch.arange(0, len(X_train), batch_size)

    # Hold the best model
    best_acc = - np.inf   # init to negative infinity
    best_weights = None

    for epoch in range(n_epochs):
        # print(epoch)
        model.train()
        with tqdm.tqdm(batch_start, unit="batch", mininterval=0, disable=True) as bar:
            bar.set_description(f"Epoch {epoch}")
            running_loss = 0.0
            for start in bar:
                # take a batch
                X_batch = X_train[start:start+batch_size]
                y_batch = y_train[start:start+batch_size]
                # forward pass
                y_pred = model(X_batch)
                loss = loss_fn(y_pred, y_batch)
                # backward pass
                optimizer.zero_grad()
                loss.backward()
                # update weights
                optimizer.step()
                # print progress
                acc = (y_pred.round() == y_batch).float().mean()
                bar.set_postfix(
                    loss=float(loss),
                    acc=float(acc)
                )
                running_loss += loss.item()
            epoch_loss = running_loss / len(bar)
            print(f"Epoch {epoch + 1}/{n_epochs}, Loss: {epoch_loss:.4f}")
        # evaluate accuracy at end of each epoch
        model.eval()
        y_pred = model(X_val)
        acc = (y_pred.round() == y_val).float().mean()
        acc = float(acc)
        if acc > best_acc:
            best_acc = acc
            best_weights = copy.deepcopy(model.state_dict())
    # restore model and return best accuracy
    model.load_state_dict(best_weights)
    return best_acc

# train-test split: Hold out the test set for final model evaluation
# define 5-fold cross validation test harness
kfold = StratifiedKFold(n_splits=5, shuffle=True)
cv_scores = []
for train, test in kfold.split(X, y):
    # create model, train, and get accuracy
    model = Deep()
    acc = model_train(model, X[train], y[train], X[test], y[test])
    print("Accuracy (wide): %.2f" % acc)
    cv_scores.append(acc)

# evaluate the model
acc = np.mean(cv_scores)
std = np.std(cv_scores)
print("Model accuracy: %.2f%% (+/- %.2f%%)" % (acc*100, std*100))
