import argparse
import sys
from helper_functions import hyperparameter_tuning,train_model,predict_new_samples,convert_ped_to_geno

def main(args):
    """  The main function
    """
    if args.tuning == True: ## For creating your own model
        hyperparameter_tuning(args.training_labels,args.training_ped,args.mapfile,args.reference_fasta)
        print ("Exiting ... \nRun the program again specifying the tuned parameters for training\n")
        sys.exit()
        
    if args.train == True: ## For training your own model
        train_model(args.training_ped,args.mapfile,args.training_labels,args.reference_fasta,args.output_model,args.l2_alpha,args.hidden_layers)
        
    else: ## Predict new samples
        predict_new_samples(args.testing_ped,args.mapfile,args.reference_fasta,args.output_file,args.training_model)
        
   
if __name__ == '__main__':
    parser = argparse.ArgumentParser('model_ancestry.py')
    parser.add_argument('--training-ped',dest='training_ped',help='the input ped file',required=False)
    parser.add_argument('--testing-ped',dest='testing_ped',help='the input ped file',required=False)
    parser.add_argument('--mapfile',dest='mapfile',help='the ethnicity markers in plink map format',required=True,default='data/filtered.37MB.master.training.map')
    parser.add_argument('--training-labels',dest='training_labels',help='the training labels file, a tsv file with sample label',required=False)
    parser.add_argument('--output-model',dest='output_model',help='the output pca file, is a pickleded python object',required=False)
    parser.add_argument('--train',dest='train',help='To Specify whether to train a model',action='store_true',default=False)
    parser.add_argument('--input-model',dest='training_model',help='the input pickled training model',required=False)
    parser.add_argument('--alpha',dest='l2_alpha',help='the alpha value obtained from model tuning',default=1e-05,type=float)
    parser.add_argument('--hidden-layers',dest='hidden_layers',help='the number of hidden layers obtained from model tuning',default=6,type=int)
    parser.add_argument('--output-prob-file',dest='output_file',help='the output file with the predict probabilities')
    parser.add_argument('--tuning',dest='tuning',help='To Specify whether to tune model parameters',action='store_true',default=False)
    parser.add_argument('--reference',dest='reference_fasta',help='The reference fasta',default='/scratch/HS_Build37/BWA_INDEX_hs37d5/hs37d5.fa',required=False)
    args = parser.parse_args()
    if args.train == False and args.training_model == None:
        parser.error("\n\nPlease Enter the path to the training model using --input-model \n")
    if args.train == True and (args.training_ped == None or args.training_labels == None):
        parser.error("\n\nPlease Enter the path to both the training pedfile and the training labels using the appropriate options \n")
    if args.train == False and (args.output_file == None or args.testing_ped == None):
        parser.error("\n\nPlease Enter the path to both the testing pedfile and the output probability file using the appropriate options \n")
    print ("Running program with the arguments : {0}".format(args))
    main(args)
