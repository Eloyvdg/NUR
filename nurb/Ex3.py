import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
import sys
from downhill import downhill_simplex

# Question 3: Spiral and elliptical galaxies
# Problem 3.a
# data = np.genfromtxt('galaxy_data.txt')
data=np.genfromtxt(os.path.join(sys.path[0],"galaxy_data.txt"))


features = data[:,:4]

for i in range(features.shape[1]): 
    features[:,i] = (features[:,i] - np.mean(features[:,i])) / np.std(features[:,i])
    
np.savetxt('galaxy_data_fs.txt', features)

bin_class = data[:,4]
fig, ax = plt.subplots(2,2, figsize=(10,8))
ax[0,0].hist(features[:,0], bins=20)
ax[0,0].set(ylabel='N', xlabel=r'$\kappa_{CO}$')
ax[0,1].hist(features[:,1], bins=20)
ax[0,1].set(xlabel='Color')
ax[1,0].hist(features[:,2], bins=20)
ax[1,0].set(ylabel='N', xlabel='Extended')
ax[1,1].hist(features[:,3], bins=20)
ax[1,1].set(xlabel='Emission line flux')
plt.savefig("fig3a.png")
plt.close()

# # Problem 3.b

def sigmoid(z):
    return 1 / (1 + np.exp(-z))

def cost_function(weights_input, features_input, labels_input):
    z = np.sum(weights_input * features_input, axis = 1)
    h_theta = sigmoid(z)
    cost = -1/len(h_theta) * np.sum(labels_input * np.log(h_theta) + (1 - labels_input)* np.log(1 - h_theta))
    return cost

def cost_function2(weights_input, features_input): 
    z = np.sum(weights_input * features_input, axis = 1)
    return sigmoid(z)

def precision(TP, FP): 
    return TP / (TP + FP)
    
def recall(TP, FN): 
    return TP / (TP + FN)

def f_score(P, R, beta):
    return (1+beta**2) * (P * R) / (beta**2 * P + R)
    
initial_simplex = np.array([[-1.0, 0.6],
                            [-0.7, 1.1],
                            [-0.8, 1.5]])

best_weights = np.zeros((4, 4, 2))
confusion_matrix = np.zeros((4, 4, 4))

names = [r'$\kappa_{CO}$', 'Color', 'Extended', 'Emission line flux']

fig, ax  = plt.subplots(1,1, figsize=(10,5), constrained_layout=True)
for i in range(4):
    for j in range(i+1, 4): 
        feature_columns = [i,j]
        features_min = features[:,feature_columns]
        cost_minimize = lambda weights_input: cost_function(weights_input, features_min, bin_class)
        minimum, value, best = downhill_simplex(initial_simplex, cost_minimize, 100, 1e-10)
        best_weights[i, j, :] = minimum
        
        results = cost_function2(minimum, features_min)
        results[results < 0.5] = 0
        results[results >= 0.5] = 1
        
        TP, TN, FP, FN = 0, 0, 0, 0
        for k in range(len(results)): 
            if results[k] == 1:
                if bin_class[k] == 1: 
                    TP += 1
                else: 
                    FP += 1
            else: 
                if bin_class[k] == 0:
                    TN += 1
                else: 
                    FN += 1
                    
        arr_confusion = np.array([TP, TN, FP, FN])         
        confusion_matrix[i,j,:] = arr_confusion
        
        ax.plot(np.arange(0,len(best)), best, label=f'{names[i]}+{names[j]}')
ax.set(xlabel='Number of iterations', ylabel='Cost function')
plt.legend(loc=(1.05,0))
plt.savefig("fig3b.png")
plt.close()

# Problem 3.c
fig, ax = plt.subplots(3,2,figsize=(10,15))
plot_idx = [[0,0], [0,1], [1,0], [1,1], [2,0], [2,1]]
for i, comb in enumerate(itertools.combinations(np.arange(0,4), 2)):
    
    ax[plot_idx[i][0],plot_idx[i][1]].scatter(features[:,comb[0]], features[:,comb[1]], c=bin_class)
    x1_boundary = np.linspace(np.min(features[:,comb[0]]), np.max(features[:,comb[0]]), 100)
    x2_boundary = -(best_weights[comb][0] * x1_boundary) / best_weights[comb][1]
    ax[plot_idx[i][0],plot_idx[i][1]].plot(x1_boundary, x2_boundary, c='black')
    
    ax[plot_idx[i][0],plot_idx[i][1]].set(xlabel=names[comb[0]], ylabel=names[comb[1]])
    ax[plot_idx[i][0],plot_idx[i][1]].plot([0.5,0.5],[0,1], 'k--')
    

plt.savefig("fig3c.png")
plt.close()

for i in range(4):
    for j in range(i+1, 4): 
        confusion = confusion_matrix[i,j,:]
        P = precision(confusion[0], confusion[2])
        R = recall(confusion[0], confusion[3])
        F = f_score(P, R, 1)
        print('---------------------------------------------------------------------')
        print(f'''For the combination of {names[i]} and {names[j]}, the results are: 
              TP = {confusion[0]}, TN = {confusion[1]}, FP = {confusion[2]}, FN = {confusion[3]}
              P = {P}, R = {R}, f = {F}''')
