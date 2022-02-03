import numpy as np
from nimfa.models.nmf_std import Nmf_std
from sklearn.decomposition import NMF
import nimfa
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

scaler = MinMaxScaler(feature_range=(0, 1))
kernel = "logarithm"
kernels = ["logarithm"]#, "square", "exponential","log-sigma-1-mm-3D", "log-sigma-2-mm-3D",
           #"log-sigma-3-mm-3D", "log-sigma-4-mm-3D", "log-sigma-5-mm-3D"]

for kernel in kernels:
    kernel = "logarithm"
    print(kernel)
    w_sparsness = []
    h_sparsness = []
    features_of_interest = [
      kernel + "_glcm_SumEntropy",
      kernel + "_glcm_Idn",
      kernel + "_glcm_DifferenceEntropy",
      kernel + "_glcm_DifferenceVariance",
      kernel + "_glcm_DifferenceAverage",
      kernel + "_glszm_SmallAreaEmphasis",
      kernel + "_glcm_Idmn",
      kernel + "_glrlm_GrayLevelNonUniformityNormalized",
      kernel + "_glrlm_LongRunEmphasis",
      kernel + "_gldm_SmallDependenceEmphasis",
      kernel + "_gldm_LargeDependenceEmphasis",

      # Added featuress
      kernel + "_glrlm_RunPercentage",
      kernel + "_glrlm_ShortRunEmphasis",
      kernel + "_glszm_LargeAreaEmphasis",
      kernel + "_glcm_Imc1",
      kernel + "_glcm_Autocorrelation",
      "source"]

    df1 = pd.read_csv("../data/Moffitt.csv")
    df1['source'] = "Moffitt"
    df1 = df1.loc[:, features_of_interest]
    df2 = pd.read_csv("../data/HarvardRT.csv")
    df2['source'] = "HarvardRT"
    df2 = df2.loc[:, features_of_interest]
    df3 = pd.read_csv("../data/M-SPORE.csv")
    df3['source'] = "M-SPORE"
    df3 = df3.loc[:, features_of_interest]
    df4 = pd.read_csv("../data/Maastro.csv")
    df4['source'] = "Maastro"
    df4 = df4.loc[:, features_of_interest]
    df5 = pd.read_csv("../data/MUMC.csv")
    df5['source'] = "MUMC"
    df5 = df5.loc[:, features_of_interest]
    df6 = pd.read_csv("../data/Radboud.csv")
    df6['source'] = "Radboud"
    df6 = df6.loc[:, features_of_interest]

    df = pd.concat([df1,df2,df3,df4,df5,df6], ignore_index=True)
    df.loc[:, features_of_interest[:-1]] = scaler.fit_transform(df.loc[:, features_of_interest[:-1]])
    df = df.dropna()
    V = np.array(df.loc[:, features_of_interest[:-1]], dtype="float32")

    # lsnmf = nimfa.Nmf(V,
    #                   seed="random",
    #                   rank=10,
    #                   max_iter=1000,
    #                   update='euclidean',
    #                   objective='fro',
    #                  min_residuals=1e-6)
    #lsnmf_fit = lsnmf()
    #dict = lsnmf_fit.fit.estimate_rank(rank_range=[2,3,4,5,6,7,8,9,10])

    for k in range(2, 10):

        lsnmf = nimfa.Nmf(V,
                          seed="random",
                          rank=k,
                          max_iter=1000,
                          update='euclidean',
                          objective='fro',
                          min_residuals=1e-4)
        lsnmf_fit = lsnmf()

        # print("K = ", k)
        # print('Rss: %5.4f' % lsnmf_fit.fit.rss())
        # print('Accuracy: %5.4f' % lsnmf_fit.fit.evar())
        # print('K-L divergence: %5.4f' % lsnmf_fit.distance(metric='kl'))
        # print('Coph cor: %5.4f' % lsnmf_fit.fit.coph_cor())
        # print('Sparseness, W: %5.4f, H: %5.4f' % lsnmf_fit.fit.sparseness())
        w_s, h_s = lsnmf_fit.fit.sparseness()
        h_sparsness.append(h_s)
        w_sparsness.append(w_s)
        print("\t K = ", k)
        print("\tH Sparsness: ", h_s)
        print("\tW Sparsness: ", w_s)
        break

    print("\tMean H Sparsness: ",np.mean(h_sparsness))
    print("\tMean W Sparsness: ", np.mean(w_sparsness))


    break

model = NMF(n_components=3,
                init='random',
                random_state=27,
                max_iter=1000,
                tol=1e-4,
                solver="cd", # cd or mu
                beta_loss="frobenius") #frobenius, itakura-saito, kullback-leibler
W = model.fit_transform(V)
H = model.components_


H = pd.DataFrame(H)
H.columns = features_of_interest[:-1]

H_m = H.reset_index().melt(id_vars=['index'])
plt.rcParams["figure.figsize"] = (20,10)
sns.barplot(x="variable", y="value", hue="index", data=H_m)
plt.title(kernel)
plt.xticks(rotation=90)
plt.show()

tracerx = pd.read_excel("/Users/au589901/GenomeDK_local/CancerEvolution/phd/Datasets/TRACERx/Radiomics/Radiomics20200801/PyRadiomics_20200718/logarithm.xlsx")
tracerx.loc[:, features_of_interest[:-1]] = scaler.fit_transform(tracerx.loc[:, features_of_interest[:-1]])

x = tracerx.loc[:, features_of_interest[:-1]]
x = x.dropna()
x = x.astype("float32")

res = model.transform(X=x)
res = pd.DataFrame(res)
res['sampleid'] = tracerx['sampleid']

print(res)
res.to_csv("tracerx_predicted.csv")