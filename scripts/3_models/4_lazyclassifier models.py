############################## import data ##############################


############################## other models ##############################
from lazypredict.Supervised import LazyClassifier

# Fit all models
clf = LazyClassifier(predictions=True)
models, predictions = clf.fit(X_train, X_test, y_train, y_test)
models

#Plot graph with 2 y axes
fig, ax1 = plt.subplots()

#Plot bars
ax1.bar(models.index.tolist(), 
        models['Balanced Accuracy'], alpha=0.3)
ax1.bar_label(ax1.containers[0], label_type='center', rotation=90, color='black', fmt='%.2f', padding=10)
ax1.set_xlabel('Classification Method')
ax1.set_xticklabels(models.index.tolist(), rotation=90)

# Make the y-axis label and tick labels match the line color.
ax1.set_ylabel('Balanced Accuracy', color='b')
[tl.set_color('b') for tl in ax1.get_yticklabels()]


#Set up ax2 to be the second y axis with x shared
ax2 = ax1.twinx()
#Plot a line
ax2.plot(models.index.tolist(), models['Time Taken'], 'r-')

# Make the y-axis label and tick labels match the line color.
ax2.set_ylabel('Time Taken', color='r')
[tl.set_color('r') for tl in ax2.get_yticklabels()]

plt.show()
