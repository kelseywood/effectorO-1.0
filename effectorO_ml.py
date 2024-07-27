import sys
import os
from argparse import ArgumentParser
from time import sleep
from numpy import ndarray, array, round
import pandas as pd
from joblib import load

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)

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


def _analyze_fasta(fasta_input:str):
	print("Parsing FASTA file...")
	fasta_content = Fasta_Content()
	fasta_content.parse_fasta_file(fasta_input)
	print("Calculating sequence features...")
	seq_features = array([get_average_features(seq) for seq in fasta_content.get_sequences()])
	return (fasta_content, seq_features)

def __import_ml_model(model_path:str):
	with open(model_path, 'rb') as fmodel:
		return load(fmodel)

def _analyze_model_path(model_path_input, seq_features:ndarray):
	print("Importing ML model...")
	model_path = str(model_path_input) if (model_path_input) else os.path.join(SCRIPT_DIR, "trained_models", "RF_88_best.sav")
	model_name = os.path.splitext(os.path.basename(model_path))[0]
	model = __import_ml_model(model_path)
	print("Sequences to run secreted oomycete-trained " +
			 	f"{'Random Forest' if (model_name == 'RF_88_best') else (model_name)} " +
				f"effector classifier on: {len(seq_features)}")
	if model_path_input == None:
		print("\n**NOTES**:\n\n" +
					"Positive training dataset: " +
					"~100 experimentally validated oomycete avirulence effectors\n" +
					"Negative training dataset: " +
					"~100 secreted orthologous oomycete genes\n" +
					"\n**END OF NOTES**\n")
	return model

def _predict_effectors(model, seq_features:ndarray, fasta_content:Fasta_Content):
	print("Predicting effectors amongst FASTA sequences...")
	predictions = model.predict(seq_features)
	probabilities = model.predict_proba(seq_features)
	PREDICTION_MAP = {'0': "predicted_non-effector", '1': "predicted_effector"}
	meanings = array([PREDICTION_MAP[pred] for pred in predictions])
	result_df = pd.DataFrame({"protein_id": fasta_content.get_ids(),
														"sequence": fasta_content.get_sequences(),
														"prediction": predictions,
														"probability": probabilities[:,1],
														"meaning": meanings})
	result_df['probability'] = round(result_df['probability'], 2)
	print(f"\nCounts of predicted classes:\n{result_df['meaning'].value_counts().to_string(header=False)}\n")
	return result_df

def _create_output_dir(output_dir:str):
	print("Creating a directory of EffectorO outputs...")
	if not output_dir: output_dir = "effectoro_results"
	if os.path.exists(output_dir):
		print("Warning: directory 'effectoro_results' already exists. Placing contents into the existing directory...")
	os.makedirs(output_dir, exist_ok=True)
	return output_dir

def _create_csv(predictions_df:pd.DataFrame, filepath:str):
	print("Constructing a CSV containing information on sequences from the FASTA file...")
	predictions_df.to_csv(filepath, index=False)

def _create_fasta(predictions_df:pd.DataFrame, filepath:str):
	print("Constructing a FASTA file with just predicted effectors...")
	with open(filepath, 'w') as outfile:
		predicted_effectors = predictions_df[predictions_df['meaning'].astype(str) == "predicted_effector"]
		seqs_to_write = [
				f">{row['protein_id']} {row['meaning']} probability={row['probability']}\n{row['sequence']}\n"
				for _, row in predicted_effectors.iterrows()
		]
		outfile.writelines(seqs_to_write)

def run_effectoro_ml(input_fasta:str, model_path, output_dir:str):
	fasta_content_to_predict, seq_features = _analyze_fasta(input_fasta)
	model = _analyze_model_path(model_path, seq_features)
	predictions_df = _predict_effectors(model, seq_features, fasta_content_to_predict)
	output_dir = _create_output_dir(output_dir)
	csv_filepath = os.path.join(output_dir, fasta_content_to_predict.get_fasta_filename() + ".effector_classification_table.csv")
	_create_csv(predictions_df, csv_filepath)
	fasta_filepath = os.path.join(output_dir, fasta_content_to_predict.get_fasta_filename() + ".predicted_effectors.fasta")
	_create_fasta(predictions_df, fasta_filepath)

	print(f"\nOutput fasta of predicted effectors available in:\n{fasta_filepath}")
	print(f"\nDetailed CSV with fasta IDs | sequences | predictions | probabilities in:\n{csv_filepath}\n")


if __name__ == "__main__":
	parser = ArgumentParser(prog="A Python script to run EffectorO.")
	parser.add_argument("-i", "--input_fasta", type=str, help="Input FASTA file", required=True)
	parser.add_argument("-m", "--model_path", type=str, help="Path to trained model to be used [default: provided Random Forest model]", required=False)
	parser.add_argument("-o", "--output_dir", type=str, help="Output directory")
	parser.add_argument("-s", "--suppress_prints", action="store_true", help="Suppress all print statements")
	args = parser.parse_args()

	# suppress print statements if indicated
	if args.suppress_prints:
		sys.stdout = open(os.devnull, 'w')
		sys.stderr = open(os.devnull, 'w')

	# run the effectoro-ML script
	run_effectoro_ml(args.input_fasta, args.model_path, args.output_dir)
