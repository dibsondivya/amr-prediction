############################## import data ##############################


############################## glm model ##############################
import statsmodels.api as sm
linear_model=sm.GLM(endog = y_train, exog = X_train) # used Gaussian
result=linear_model.fit(xname=readfile_by_kmer_df.columns.tolist())
#print(result.summary())
results_as_html = result.summary().tables[1].as_html()
coefs_df = pd.read_html(results_as_html, header=0, index_col=0)[0]
coefs_df['kmer'] = readfile_by_kmer_df.columns.tolist()
coefs_df = coefs_df.set_index('kmer')
coefs_df

significant_coefs = coefs_df[coefs_df['P>|z|'] < 0.05] # 94 rows left
ribo_list = []
for kmer in significant_coefs.index.tolist():
    ribotype = library_kmer_2_ribotype[kmer]
    ribo_list.append(ribotype)
significant_coefs['library ribo'] = ribo_list
significant_coefs['kmer'] = significant_coefs.index
significant_coefs

# significant plot
import seaborn as sns
ax = sns.barplot(x='kmer', y='coef', data=significant_coefs, hue='library ribo', dodge=False)
plt.title(" K-mers with significant linear coefficients (P<0.05) ")
plt.xticks(rotation=90, fontsize=5)

plt.ylabel("Linear Coefficient")
plt.xlabel("k-mer")

plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

# plot dist of library ribotypes
significant_coefs['library ribo'].value_counts().plot(kind='bar')

# try plot of confidence intervals
#for lower,upper,y in zip(significant_coefs['[0.025'],significant_coefs['0.975]'],range(len(significant_coefs))):
#    plt.plot((lower,upper),(y,y),'ro-',color='orange')
#plt.yticks(range(len(significant_coefs)),list(significant_coefs.index))

# get linear plot
from statsmodels.graphics.api import abline_plot

fig, ax = plt.subplots()
yhat = result.mu
colors = ['tab:blue'] * len(yhat)
colors[231] = 'r'

ax.scatter(yhat, y_train,c=colors)
ax.scatter([],[],c='r',label='Outlier')
ax.legend(loc='best')

line_fit = sm.OLS(y_train, sm.add_constant(yhat, prepend=True)).fit()
abline_plot(model_results=line_fit, ax=ax)

ax.set_title('Linear Model Fit Plot on Training Set')
ax.set_ylabel('True Ribotypes')
ax.set_xlabel('Predicted Ribtoypes');
ax.set_yticks([1, 14, 27, 78])
ax.set_xticks([1, 14, 27, 58, 78])

## plot on test set
fig, ax = plt.subplots()
ypred = result.predict(exog=X_test)
colors = ['tab:blue'] * len(ypred)

ax.scatter(ypred, y_test,c=colors)
ax.legend(loc='best')

line_fit = sm.OLS(y_train, sm.add_constant(yhat, prepend=True)).fit()
abline_plot(model_results=line_fit, ax=ax)

ax.set_title('Linear Model Fit Plot on Test Set')
ax.set_ylabel('True Ribotypes')
ax.set_xlabel('Predicted Ribtoypes');
ax.set_yticks([1, 14, 27, 78])

# Residual Dependence Plot
fig, ax = plt.subplots()

ax.scatter(yhat, result.resid_pearson)
ax.hlines(0, 0, 1)
ax.set_xlim(0, 1)
ax.set_title('Residual Dependence Plot')
ax.set_ylabel('Pearson Residuals')
ax.set_xlabel('Fitted values')

# QQ plot
from statsmodels import graphics
graphics.gofplots.qqplot(result.resid_deviance.copy(), line='r')

## remove outlier
import numpy as np
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]
find_nearest(yhat, 60) # 58.31744936865806
itemindex = np.where(yhat == 58.31744936865806)
itemindex # 231

del_linear_model=sm.GLM(endog = np.delete(y_train, 231, 0), exog = np.delete(X_train, 231, 0)) # used Gaussian
del_result=del_linear_model.fit(xname=readfile_by_kmer_df.columns.tolist())

from statsmodels import graphics
graphics.gofplots.qqplot(del_result.resid_deviance.copy(), line='r')
