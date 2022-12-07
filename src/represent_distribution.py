"""
distribution globale AG AC confondues
distribution AG 
distribution AB
distribution AG accessible > 40 
distribution AB accessible > 40
distribution AG accessible > 70 
distribution AB accessible > 70
distribution AG inaccessible < 40 
distribution AB inaccessible < 40
frequence théorique contact > 40
frequence théorique contact > 70

def frequence_theorique_2ways(AB,AG):
    \"""
    deux manieres de calculer:
    - par complexe:
        pour chaque complexe on recupere les fréquences calculées pour chaque AA, puis on les multiplies (400 valeurs)
        Ensuite ces 400 valeurs sont moyennées pour chaque complexe sur une matrice 20*20
        . On obtient une matrice moyenne des frequences de contact
    - sur toute l'ABDB:
        pour chaque AA on moyenne toutes les frequences, tout complexe confondu puis on les multiplies (400 valeurs)
        . On obtient une matrice des frequences moyennes de contact 
    \"""
    AG.loc[np.isnan(AG['Frequency']),'Frequency'] = 0
    AB.loc[np.isnan(AB['Frequency']),'Frequency'] = 0
    # POUR ABDB
    # recuperer les indexes des AA identiques pour AG et AB
    AAorder = pd.unique(AG['TypeAA'] )
    AAIAG = np.array([ AG.loc[AG['TypeAA'] == aa, 'Frequency' ] for aa in AAorder ])
    AAIAB = np.array([ AB.loc[AB['TypeAA'] == aa, 'Frequency' ] for aa in AAorder ])    
    # moyenne des frequences par AA -> multiplication des frequences entre iAB et iAG
    
    meanAG = np.mean(AAIAG,axis=1)
    meanAB = np.mean(AAIAB,axis=1)
    matAA = pd.DataFrame( index= AAorder, columns= AAorder)
    matAA = matAA.fillna(0)
    # matrice 20*20
    for i,aai in enumerate(AAorder):
        for j,aaj in enumerate(AAorder):
            matAA.loc[aai,aaj] = meanAB[i]*meanAG[j]
    # print(matAA)
    # PAR COMPLEXE
    # recuperer les indexes des complexes identiques pour AG et AB 
    pdb_unique = pd.unique( AG['PDBID'])
    cplxIAG = np.array([ AG.loc[AG['PDBID'] == pdbid, 'Frequency' ] for pdbid in pdb_unique ])
    cplxIAB = np.array([ AB.loc[AB['PDBID'] == pdbid, 'Frequency' ] for pdbid in pdb_unique ])
    # pour chaque indexe iAG iAB multiplier la frequence et donner une valeur iAAAG iAAAB
    multi = np.array([ np.array([ cplxIAB[i][k]*f for k in range(20) for f in e ]) for i,e in enumerate(cplxIAG)])
    meanMulti = np.mean(multi, axis=0)
    # moyenne de ces valeurs (division par le nombre de complexes) dans une matrice 20*20
    matCPLX = pd.DataFrame( index= AAorder, columns= AAorder)
    matCPLX = matCPLX.fillna(0)
    k = 0
    for aai in AAorder:
        for aaj in AAorder:
            matCPLX.loc[aai,aaj] = meanMulti[k]
            k+=1

    return (matAA,matCPLX)
"""
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def barplot(df):

    sns.barplot(x='TypeAA',y='Frequency', data=df)
    plt.show()


def barplotHue(ndf,outputpathName,sectionName):
    ax = sns.barplot(x='TypeAA',y='Frequency', hue='Section',data=ndf)
    # plt.show()
    name = 'Distribution of Amino acids for {}'.format(sectionName)
    ax.set_title(name)
    plt.xticks(rotation=45)
    plt.savefig('{}{}'.format(outputpathName,'_'.join(name.split())),dpi=1000)
    plt.close()

def matrix_representation(matrix,outputpathName,name):
    """represent a matrixF in an hitmap
    """

    sns.heatmap(matrix, cmap = 'Blues')
    plt.xlabel('Antigen')
    plt.ylabel('Antibody')
    
    plt.title(name)
    plt.savefig('{}{}'.format(outputpathName,'_'.join(name.split())),dpi=1000)
    plt.close()

def transform(data,colnames,colsumname,sectionName):

    D = pd.melt(data,id_vars='PDBID', value_vars = colnames,  value_name='NbAA',var_name = 'TypeAA', ignore_index=False )
    D['Total'] = data.loc[D.index,colsumname]
    D['Frequency'] = D['NbAA'] / D['Total']
    D['Section'] = sectionName
    D['TypeAA'] = [aat.split('_')[-1] for aat in D['TypeAA']]
    D.reset_index(drop=True, inplace=True)
    iAB = [i for i in D.index if D.loc[i,'PDBID'].split('_')[-1] == 'AB']
    AB = D.loc[ iAB ,: ]
    AG = D.drop( iAB )
    AB['PDBID'] = [ id[:-3] for id in AB['PDBID']]
    AG['PDBID'] = [ id[:-3] for id in AG['PDBID']]
    return (D,AB,AG)


def frequence_theorique(AB,AG):
    """
    deux manieres de calculer:
    - par complexe:
        pour chaque complexe on recupere les fréquences calculées pour chaque AA, puis on les multiplies (400 valeurs)
        Ensuite ces 400 valeurs sont moyennées pour chaque complexe sur une matrice 20*20
        . On obtient une matrice moyenne des frequences de contact
    - sur toute l'ABDB:
        pour chaque AA on moyenne toutes les frequences, tout complexe confondu puis on les multiplies (400 valeurs)
        . On obtient une matrice des frequences moyennes de contact 
    """
    AG.loc[np.isnan(AG['Frequency']),'Frequency'] = 0
    AB.loc[np.isnan(AB['Frequency']),'Frequency'] = 0
    # POUR ABDB
    # recuperer les indexes des AA identiques pour AG et AB
    AAorder = pd.unique(AG['TypeAA'] )
    AAIAG = np.array([ AG.loc[AG['TypeAA'] == aa, 'Frequency' ] for aa in AAorder ])
    AAIAB = np.array([ AB.loc[AB['TypeAA'] == aa, 'Frequency' ] for aa in AAorder ])    
    # moyenne des frequences par AA -> multiplication des frequences entre iAB et iAG
    
    meanAG = np.mean(AAIAG,axis=1)
    meanAB = np.mean(AAIAB,axis=1)
    matAA = pd.DataFrame( index= AAorder, columns= AAorder)
    matAA = matAA.fillna(0)
    # matrice 20*20
    for i,aai in enumerate(AAorder):
        for j,aaj in enumerate(AAorder):
            matAA.loc[aai,aaj] = meanAB[i]*meanAG[j]
    return matAA


            

def organizer(distribFile,outputpathName):
    
    distribution = pd.read_csv(distribFile, sep = '\t', engine = 'c')
    header = distribution.columns
    globalcol = header[1:21]
    S40col = header[22:42]
    S70col = header[43:63]
    s40col = header[64:84]

    GlobalD,Global_AB,Global_AG = transform(distribution,globalcol,'AA number','Global')
    S40colD,S40colD_AB,S40colD_AG = transform(distribution,S40col,'S40 number','Sup_40')
    S70colD,S70colD_AB,S70colD_AG = transform(distribution,S70col,'S70 number','Sup_70')
    s40colD,s40colD_AB,s40colD_AG = transform(distribution,s40col,'s40 number','Sub_40')

    # S40AA,S40CP = frequence_theorique(S40colD_AB,S40colD_AG)
    # matrix_representation(S40AA,outputpathName,'Sup 40 mean frequency matrix of contacts')
    # matrix_representation(S40CP,outputpathName,'Sup 40 frequency matrix of mean contacts')

    # S70AA,S70CP = frequence_theorique(S70colD_AB,S70colD_AG)
    # matrix_representation(S70AA,outputpathName,'Sup 70 mean frequency matrix of contacts')
    # matrix_representation(S70CP,outputpathName,'Sup 70 frequency matrix of mean contacts')

    s40AA = frequence_theorique(s40colD_AB,s40colD_AG)
    matrix_representation(s40AA,outputpathName,'Sub 40 mean frequency matrix of contacts')

    # GAA,GCP = frequence_theorique(Global_AB,Global_AG)
    # matrix_representation(GAA,outputpathName,'Global mean frequency matrix of contacts')
    # matrix_representation(GCP,outputpathName,'Global frequency matrix of mean contacts')
    # ndf = pd.concat([GlobalD,S40colD,S70colD,s40colD])
    # barplotHue(ndf,outputpathName,'Antibody and Antigen')

    # ndf = pd.concat([Global_AB,S40colD_AB,S70colD_AB,s40colD_AB])
    # barplotHue(ndf,outputpathName,'Antibodies')

    # ndf = pd.concat([Global_AG,S40colD_AG,S70colD_AG,s40colD_AG])
    # barplotHue(ndf,outputpathName,'Antigens')

if __name__ == '__main__':

    try:
        distribFile = sys.argv[1]
        outputpathName = sys.argv[2]

    except:
        sys.exit(__doc__)

    organizer(distribFile,outputpathName)