from __future__ import print_function
import numpy as np
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import allel
import sys
import os
import subprocess
import cPickle as pickle
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPClassifier
from sklearn.decomposition import PCA
from sklearn import preprocessing
from pyfaidx import Fasta 

def hyperparameter_tuning(training_labels,training_ped,map_file,reference_fasta):
    """ Tune parameters for the MLP
    training_labels ; str : path to the training labels file 
    training_ped ; str : path to the training ped file
    map_file ; str : path to the mapfile
    reference_fasta ; str : path to the reference fasta file
    """    
    print ("Scaling the input....\n")
    population_labels,training_samples = return_populations(training_labels)    
    #training_genotypes = get_genotypes_eigenstrat(training_genofile,nsamples)
    training_genotypes = convert_ped_to_geno(training_ped,map_file,reference_fasta)
    training_genotypes[training_genotypes == 9] = 0
    scaler = allel.stats.preprocessing.get_scaler("standard",True,2)
    transformed_training = scaler.fit_transform(training_genotypes)
    training_genotypes = transformed_training.transpose()
    
    print ("Running The pipeline PCA and Neural Net !!")
    pipe = Pipeline([
                     ('pca',PCA(n_components=6)),
                     ('nnet',MLPClassifier(solver='lbfgs',activation='logistic',alpha=0.1,hidden_layer_sizes=(6,)))])

    # Split the dataset and hold out 20% for classification
    X_train, X_test, y_train, y_test = train_test_split(
            training_genotypes,population_labels,test_size=0.2)
    
    ## Tune some parameters .. 
    tuned_parameters = [{'nnet__alpha':[1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1],'nnet__hidden_layer_sizes':[(6,),(8,),(12,),(16,),(24,)]}]
    scores = ['precision','recall']

    for score in scores:
        print("# Tuning hyper-parameters for %s" % score)
        print()
        clf = GridSearchCV(pipe,tuned_parameters,cv=10,scoring='%s_macro' % score,n_jobs=10)
        clf.fit(X_train, y_train)
        print("Best parameters set found on development set:")    
        print()
        print(clf.best_params_)
        print()
        print("Grid scores on development set:")
        print()
        means = clf.cv_results_['mean_test_score']
        stds = clf.cv_results_['std_test_score']
        for mean, std, params in zip(means, stds, clf.cv_results_['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))
        print()            
        print("Detailed classification report:")
        print()
        print("The model is trained on the full development set.")
        print("The scores are computed on the full evaluation set.")
        print()
        y_true, y_pred = y_test, clf.predict(X_test)
        print(classification_report(y_true, y_pred))
        print()            

def convert_ped_to_geno(pedfile,mapfile,ref_fasta):
    """ Convert a ped file to a 0,1 or 2 coded geno file
    returns a 2d numpy array (n_snps x n_samples)
    
    pedfile : str ; path to the pedfile
    mapfile : str ; path to the mapfile
    reference_fasta : str ; path to the reference fasta file
    """    
    genome = Fasta(ref_fasta)
    positions = []
    ref_alleles = []
    with open(mapfile,'r') as IN:
        for line in IN:
            line=line.strip('\n')
            pos = int(line.split('\t')[3])
            chrom = line.split('\t')[0]
            ref_allele = get_ref_allele(genome,chrom,pos)
            positions.append(pos)
            ref_alleles.append(ref_allele)            
    
    geno_mat = np.zeros([len(positions),10000]) ## A new genotype matrix,the number of samples are arbitrarily declared
    i = 0 
    with open(pedfile,'r') as IN:
        for line in IN:
            line=line.strip('\n')
            temp_genotypes = line.split(" ")[6:]
            span = 2
            genotypes = ["".join(temp_genotypes[k:k+span]) for k in range(0,len(temp_genotypes),span)]
            j = 0
            coded_genos = [] ## Store the coded genotypes
            
            for geno in genotypes:
                count = 0
                for allele in list(geno):
                    if allele == '0': ## Missing
                        count = 9
                        break
                    if allele != ref_alleles[j]:
                        count+=1
                j+=1
                coded_genos.append(count)                
            geno_mat[:,i] = coded_genos
            i+=1
            
    return geno_mat[:,0:i]

def get_ref_allele(genome,chrom,pos):
    """ Get the reference site from a pyfasta object
    
    genome : a pyfasta object ; the reference genome
    chrom : str ; the chromosome 
    pos : int ; the position
    
    return : str ; the reference allele 
    """
    return str(genome[chrom][pos-1:pos]) ## The fasta file will be 0-based,
                                         ## so pos will correspond to pos-1 in the fasta

def get_genotypes_eigenstrat(genofile,nsamples):
    """Read an eigenstrat genotype file (snps x samples)
    and store in a 2d numpy array

    Values are : 0 = homref, 1 = het and 2 = hom alt

    genofile : str; path to the genotype file
    nsamples : int; the number of samples in your file 

    return : genotype_array ; a snps x samples numpy float array
    """
    i=0
    nsamples = int(nsamples)
    ## Store an arbitrarily large array in the future
    genotype_array = np.zeros((29961,nsamples),dtype=np.float)
    
    with open(genofile,'r') as IN:
        for line in IN:
            line=line.strip('\n')
            genotype_array[i,:] = list(line)
            i+=1    
    print ("Read {0} snps".format(i))
    return genotype_array[0:i,0:nsamples]

def return_populations(labels_file):
    """ Return a list of training lables ordered by the samples
    in the genotype file

    labels_file : str ; the path to the training labels file
    returns : tuple ; (list,list) ; the population_labels and the
                                    training sample names 
    """
    population_labels = []
    training_samples = []
    with open(labels_file,'r') as L:
        for line in L:
            line=line.strip('\n')
            contents=line.split('\t')
            population_labels.append(contents[1])
            training_samples.append(contents[0])
            
    return (population_labels,training_samples)

def train_model(training_pedfile,map_file,training_labels,reference_fasta,output_model,l2_alpha,hidden_layers):
    """
    Train the model 
    training_genotfile : str ; the path to the training genotype file 
    training_labels : str ; the path to the training class labels
    output_model : str ; path to the save the training model
    """
    pipe = Pipeline([
        ('scaling',preprocessing.StandardScaler()),
        ('pca',PCA(n_components=6)),
        ('nnet',MLPClassifier(solver='lbfgs',activation='logistic',alpha=l2_alpha,hidden_layer_sizes=(hidden_layers,),random_state=0))])
            
    population_labels,training_samples = return_populations(training_labels)    
    #training_genotypes = get_genotypes_eigenstrat(training_genofile,nsamples)
    training_genotypes = convert_ped_to_geno(training_pedfile,map_file,reference_fasta)
    training_genotypes[training_genotypes == 9] = 0
    transposed_genotypes = training_genotypes.transpose()
    print ("Fitting the model...")
    pipe.fit(transposed_genotypes,population_labels)
    print ("Saving the model...")
    with open(output_model,"wb") as out_file:
        pickle.dump(pipe,out_file)

def predict_new_samples(testing_pedfile,map_file,reference_fasta,output_file,training_model):
    """ Predict new samples
    testing_pedfile : str ; path to the pedfile 
    map_file : str ; path to the mapfile
    reference_fasta : str ; path to the reference genome file
    output_file : str ; path to the output file containing the probabilities
    training_model : str ; path to the pickled python object (pca+MLP training model) 
    """    
    testing_genotypes = convert_ped_to_geno(testing_pedfile,map_file,reference_fasta)
    ## Get genotyping rate
    genotyping_rate = float((testing_genotypes==9).sum(axis=0))/float(np.shape(testing_genotypes)[0])
    testing_genotypes[testing_genotypes == 9] = 0
    transposed_genotypes = testing_genotypes.transpose()
    print ("Opening the model...")
    with open(training_model,"rb") as input_model:
        model = pickle.load(input_model)
        print ("Outputing the probabilities...")
        with open(output_file,'w') as OUT:
            i = 0
            for e in model.predict_proba(transposed_genotypes):
                output_line = '\t'.join([str(prob) for prob in e])
                geno_rate = str(1.00 - genotyping_rate)
                output_line = output_line+'\t'+geno_rate
                OUT.write(output_line+'\n')
                i+=1
