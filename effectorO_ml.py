from argparse import ArgumentParser, FileType
from os import makedirs
from os.path import exists
from time import sleep
from numpy import array, round
import pandas as pd
from joblib import load

from gen_src.fasta_content import Fasta_Content
from ml_src.get_average_features import get_average_features

## take in:
##    1) secreted proteins fasta file
##    2) (optional) model file path (will default to best Random Forest)

## RUN LIKE THIS:
##    Using conda: python3 effectorO.py -i {INPUT_FASTA_PATH}
##    Using requirements.txt: python3.6 effectorO.py {INPUT_FASTA_PATH}

## output:
##    1) csv conatining (by column):
##			 - protein_id from FASTA file
##			 - sequence from FASTA file
##			 - prediction - 0=non_effector | 1=effector
##			 - probability that a sequence is an effector
##			 - meaning - verbose string version of prediction
##    2) fasta file of predicted effectors


def import_ml_model(model_path:str):
	if model_path == None: model_path = "trained_models/RF_88_best.sav"
	with open(model_path, 'rb') as fmodel:
		return load(fmodel)


def main():
	parser = ArgumentParser(prog="A Python script to run EffectorO.")
	parser.add_argument("-i", "--input_fasta", type=FileType('r'), help="Input FASTA file", required=True)
	parser.add_argument("-m", "--model_path", type=str, help="Path to trained model to be used [default: provided Random Forest model]", required=False)
	args = parser.parse_args()

	print("Parsing FASTA file...")
	FASTA_CONTENT_TO_PREDICT = Fasta_Content()
	FASTA_CONTENT_TO_PREDICT.parse_fasta_file(args.input_fasta)
	print("Calculating sequence features...")
	seq_features = array([get_average_features(seq) for seq in FASTA_CONTENT_TO_PREDICT.get_sequences()])

	print("Importing ML model...")
	model_path = str(args.model_path) if (args.model_path) else None
	if model_path != None:
		model_name =  model_path[model_path.rfind('/')+1:model_path.rfind('.')]
	MODEL = import_ml_model(model_path)

	print("Sequences to run secreted oomycete-trained " +
			 	f"{'Random Forest' if (model_path == None) else (model_name)} " +
				f"effector classifier on: {len(seq_features)}")
	if model_path == None:
		print("\n**NOTES**:\n\n" +
					"Positive training dataset: " +
					"~100 experimentally validated oomycete avirulence effectors\n" +
					"Negative training dataset: " +
					"~100 secreted orthologous oomycete genes\n" +
					"\n**END OF NOTES**\n")

	print("Predicting effectors amongst FASTA sequences...")
	predictions = MODEL.predict(seq_features)
	probabilities = MODEL.predict_proba(seq_features)
	PREDICTION_MAP = {'0': "predicted_non-effector", '1': "predicted_effector"}
	meanings = array([PREDICTION_MAP[pred] for pred in predictions])

	result_df = pd.DataFrame({"protein_id": FASTA_CONTENT_TO_PREDICT.get_ids(),
														"sequence": FASTA_CONTENT_TO_PREDICT.get_sequences(),
														"prediction": predictions,
														"probability": probabilities[:,1],
														"meaning": meanings})
	result_df['probability'] = round(result_df['probability'], 2)

	print(f"\nCounts of predicted classes:\n{result_df['meaning'].value_counts().to_string(header=False)}\n")

	print("Creating a directory for EffectorO outputs...")
	OUTPUT_DIR = "effectoro_results/"
	if exists(OUTPUT_DIR):
		print("Warning: directory 'effectoro_results' already exists. Replacing its contents in 3 seconds...")
		sleep(3)
	makedirs(OUTPUT_DIR, exist_ok=True)

	print("Constructing a CSV containing information on sequences from the FASTA file...")
	csv_filename = f"{OUTPUT_DIR}/{FASTA_CONTENT_TO_PREDICT.get_fasta_filename()}.effector_classification_table.csv"
	result_df.to_csv(csv_filename, index=False)

	print("Constructing a FASTA file with just predicted effectors...")
	out_filename = f"{OUTPUT_DIR}/{FASTA_CONTENT_TO_PREDICT.get_fasta_filename()}.predicted_effectors.fasta"
	with open(out_filename, 'w') as outfile:
		predicted_effectors = result_df[result_df['meaning'].astype(str) == "predicted_effector"]
		seqs_to_write = ('>' + predicted_effectors["protein_id"].astype(str) + ' ' +
									 	 predicted_effectors["meaning"].astype(str) + ' ' +
									 	 "probability=" + predicted_effectors["probability"].astype(str) + '\n' +
										 predicted_effectors["sequence"].astype(str) + '\n').to_list()
		outfile.writelines(seqs_to_write)

	print("\n\nOutput fasta of predicted effectors available in:\n" +
			 	out_filename + "\n\n" +
			 	"Detailed CSV with fasta IDs | sequences | predictions | probabilities in:\n" +
			  csv_filename + "\n\n")


if __name__ == "__main__":
	main()
